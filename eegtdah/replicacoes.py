"""Replicações pré-registradas: referência (R0), alegações espectrais (R1 a R7), vazamento, alegações espaciais e aquisição."""
import itertools
import json
import math
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import eigh

from . import estatistica as est
from . import sinais as s
from .config import BORDA8, CANAIS, FIGURAS, FS, POST, RAIZ, RESULTADOS, indices, isolar_threads

FB = {'fb1': (0.6, 11.6), 'fb2': (11.6, 21.6), 'fb3': (21.6, 31.6), 'fb4': (31.6, 41.6)}
BANDAS = {'larga': (1, 40), 'teta': (4, 8), 'alfa': (8, 13), 'beta': (13, 25), **FB}
BANDAS_E2 = ((0.5, 4), (4, 8), (8, 13), (13, 30), (30, 45))
REFS = ('nativa', 'car', 'rest', 'csd')
PACOTES = (('Fp1', 'Fp2'), ('F7', 'F3', 'Fz', 'F4', 'F8'), ('T7', 'C3', 'Cz', 'C4', 'T8'),
           ('P7', 'P3', 'Pz', 'P4', 'P8'), ('O1', 'O2'))
IU = np.triu_indices(19, 1)
PAR = {(CANAIS[i], CANAIS[j]): k for k, (i, j) in enumerate(zip(*IU))}
CZ, IPOST = CANAIS.index('Cz'), indices(POST)
INTERNOS = [i for i, c in enumerate(CANAIS) if c not in BORDA8]
N_REP, N_PERM, N_BOOT, MIN_A = 5, 100, 1000, (30, 8)
PARES_M3 = (('adhdata_tarefa', 'hbn_fechado'), ('tdbrain_fechado', 'hbn_fechado'), ('tdbrain_aberto', 'hbn_aberto'))
E3 = (('E3', 'teta_frontal_homologo', (('Fp1', 'Fp2'), ('F3', 'F4'), ('F7', 'F8')), 4, 8, 1),
      ('E3', 'alfa_frontal_longo', (('F7', 'F8'), ('Fp1', 'F8'), ('Fp2', 'F7'), ('F3', 'F8'), ('F4', 'F7')), 8, 13, -1))
CURTOS = (('F3', 'Fz'), ('F4', 'Fz'), ('Fp1', 'F3'), ('Fp2', 'F4'), ('F7', 'F3'), ('F8', 'F4'), ('F3', 'C3'),
          ('F4', 'C4'), ('Fz', 'Cz'))
LONGOS = (('F7', 'F8'), ('Fp1', 'F8'), ('Fp2', 'F7'), ('F7', 'F4'), ('F8', 'F3'))
E4 = tuple(('E4', f'alfa_{n}_{k}', p, lo, lo + 1, sg) for n, lo, sg in (('baixo', 8, 1), ('alto', 10, -1))
           for k, p in (('curtos', CURTOS), ('longos', LONGOS), ('todos', CURTOS + LONGOS)))
NOME = {'adhdata_tarefa': 'adhdata (tarefa)', 'hbn_fechado': 'HBN OF', 'hbn_aberto': 'HBN OA',
        'tdbrain_fechado': 'TDBRAIN OF', 'tdbrain_aberto': 'TDBRAIN OA'}
_CFG = {}


def _usar(cfg):
    _CFG.update(cfg)


def _dir(*p):
    return Path(_CFG['saida_dir']).joinpath(*p)


def _meta(cfg, bloco):
    m = pd.read_csv(Path(cfg['pacote_dir']) / bloco / 'metadados.csv')
    m = m[m['grupo'].isin(['TDAH', 'controle'])]
    if bloco.startswith('tdbrain'):
        m = m[m['idade'] >= 18]
    m = m.reset_index(drop=True)
    m['y'] = (m['grupo'] == 'TDAH').astype(int)
    m['sexo_f'] = np.where(m['sexo'] == 'F', 1.0, np.where(m['sexo'] == 'M', 0.0, np.nan))
    m['log_dur'] = np.log(m['duracao_s'].astype(float))
    return m


def _taghi_sub(m, bloco):
    if not bloco.startswith('hbn'):
        return set(m['codigo'])
    ctrl = m[m.grupo == 'controle']
    return set(ctrl['codigo']) | set(m[m.grupo == 'TDAH'].sample(len(ctrl), random_state=0)['codigo'])


def _chis(W, hi=40.0):
    P, f = s.espectro(W)
    Pm = np.median(P, 0)
    return np.array([m.aperiodic_params_[-1] if (m := s.fooof(f, q, hi)) else np.nan for q in Pm]), Pm, f


def _alfa_post(f, p, hi=40.0):
    from scipy.signal import savgol_filter
    m = s.fooof(f, p, hi)
    if m is None:
        return np.full(3, np.nan)
    ff, r = s.razao_pa(m)
    pa = savgol_filter(r, 5, 2)
    j = np.flatnonzero((ff >= 7) & (ff <= 13))
    e = np.clip(pa[j] - 1, 0, None)
    k = j[np.argmax(pa[j])]
    return np.array([np.nan if k in (j[0], j[-1]) else ff[k], (ff[j] * e).sum() / e.sum() if e.sum() > 0 else np.nan,
                     pa[j].max() - 1])


def _covs(W, bandas, n=15):
    from mne.filter import filter_data
    if len(W) < n:
        return {b: np.full((W.shape[1],) * 2, np.nan) for b in bandas}
    return {b: np.mean([w @ w.T / w.shape[1] for w in filter_data(W[:n].astype(np.float64), FS, lo, min(hi, FS / 2 - 1),
                                                                  verbose='ERROR')], 0) for b, (lo, hi) in bandas.items()}


def _coerencia(W):
    f = np.fft.rfftfreq(W.shape[-1], 1 / FS)
    sel = (f >= 1) & (f <= 45)
    X = np.fft.rfft(W * np.hanning(W.shape[-1]), axis=-1)[:, :, sel]
    S = np.einsum('eif,ejf->ijf', X, X.conj()) / len(W)
    P = np.real(np.einsum('iif->if', S))
    return (np.abs(S[IU[0], IU[1]]) ** 2 / (P[IU[0]] * P[IU[1]] + 1e-30)).astype(np.float32), f[sel]


def _entropia_perm(W, m=3):
    code = (np.argsort(np.lib.stride_tricks.sliding_window_view(W, m, axis=-1), -1) * m ** np.arange(m)).sum(-1)
    n_ep, c, L = code.shape
    cnt = np.bincount((code.reshape(-1, L) + (np.arange(n_ep * c) * m ** m)[:, None]).ravel(),
                      minlength=n_ep * c * m ** m).reshape(n_ep, c, m ** m) / L
    with np.errstate(divide='ignore', invalid='ignore'):
        return -np.nansum(np.where(cnt > 0, cnt * np.log(cnt), 0), -1) / np.log(math.factorial(m))


def _katz(x):
    d = np.abs(np.diff(x))
    L, dist = d.sum(), np.abs(x - x[0]).max()
    n = L / d.mean()
    return float(np.log10(n) / (np.log10(n) + np.log10(dist / L))) if dist > 0 and L > 0 else np.nan


def _higuchi(x, kmax=64):
    N, ks, lk = len(x), np.arange(1, kmax + 1), []
    for k in ks:
        lk.append(np.mean([np.abs(np.diff(x[m::k])).sum() * (N - 1) / ((len(x[m::k]) - 1) * k) / k
                           for m in range(k) if len(x[m::k]) >= 2]))
    return float(-np.polyfit(np.log(ks), np.log(lk), 1)[0])


def _apen(x, m=2, rf=0.22):
    r = rf * np.std(x)

    def phi(mm):
        X = np.lib.stride_tricks.sliding_window_view(x, mm)
        C = np.concatenate([(np.abs(X[i:i + 256, None] - X[None]).max(-1) <= r).mean(1) for i in range(0, len(X), 256)])
        return np.log(C).mean()
    return float(phi(m) - phi(m + 1))


def _lle(x, m=2, passos=20):
    X = np.lib.stride_tricks.sliding_window_view(x, m)
    n = len(X) - passos
    if n < 50:
        return np.nan
    D = np.abs(X[:n, None] - X[None, :n]).max(-1)
    for i in range(n):
        D[i, max(0, i - 10):i + 11] = np.inf
    viz = D.argmin(1)
    div = []
    for j in range(1, passos):
        d = np.linalg.norm(X[np.arange(n) + j] - X[viz + j], axis=1)
        div.append(np.log(d[d > 0]).mean() if (d > 0).any() else np.nan)
    div, ok = np.array(div), np.isfinite(div)
    return float(np.polyfit(np.arange(1, passos)[ok], div[ok], 1)[0] * FS) if ok.sum() > 5 else np.nan


def _taghi(xc, segs):
    from mne.filter import filter_data
    n, nj = int(30 * FS), int(10 * FS)
    E = [xc[:, i:i + n] for a, b in segs for i in range(int(a), int(b) - n + 1, n)]
    if not E:
        return {}
    E = np.stack(E)
    cov, nl = np.zeros((len(E), 4, 19, 19)), np.full((len(E), 4, 3, 6, 4), np.nan)
    pac = [indices(p) for p in PACOTES]
    for b, (lo, hi) in enumerate(FB.values()):
        Ef = filter_data(E, FS, lo, min(hi, FS / 2 - 1), verbose='ERROR')
        for e in range(len(E)):
            cov[e, b] = Ef[e] @ Ef[e].T / n
            for w in range(3):
                med = np.array([[_lle(q), _apen(q), _katz(q), _higuchi(q)] for q in Ef[e, :, w * nj:(w + 1) * nj]])
                nl[e, b, w] = np.vstack([np.nanmean(med[ix], 0) for ix in pac] + [np.nanmean(med, 0)])
    chi = [np.nanmedian(_chis(s.epocas(E[e], [(0, n)]), 30.0)[0]) for e in range(len(E))]
    return {'cov30': cov.astype(np.float32), 'nl': nl.astype(np.float32), 'chi_ep': np.array(chi)}


def _uma(args):
    pacote, saida, bloco, cod, taghi = args
    alvo = Path(saida) / bloco / f'{cod}.npz'
    if alvo.exists():
        return None
    try:
        from mne.filter import filter_data
        x, segs = s.carregar(pacote, bloco, cod)
        W = {r: s.epocas(s.referenciar(x, r), segs) for r in REFS}
        Wc, out = W['car'], {}
        for r in REFS:
            out[f'chi_{r}'], Pm, f = _chis(W[r])
            if r in ('car', 'nativa'):
                out[f'psd_{r}'] = Pm[:, f <= 45].astype(np.float32)
            if r == 'car':
                out['chi30'] = np.array([m.aperiodic_params_[-1] if (m := s.fooof(f, q, 30.0)) else np.nan for q in Pm])
                out['alfa_post'], out['f_psd'] = _alfa_post(f, np.median(Pm[IPOST], 0)), f[f <= 45]
        h = len(Wc) // 2
        out['chi_metades'] = np.full((2, 19), np.nan)
        out['alfa_metades'] = np.full((2, 3), np.nan)
        for k, Wh in enumerate((Wc[:h], Wc[h:])):
            if len(Wh) >= 3:
                out['chi_metades'][k], Pm, f = _chis(Wh)
                out['alfa_metades'][k] = _alfa_post(f, np.median(Pm[IPOST], 0))
        out.update({f'cov_{b}': c for b, c in _covs(Wc, BANDAS).items()})
        for r in ('rest', 'csd'):
            out[f'cov_{r}'] = _covs(W[r], {'l': (1, 40)})['l']
        out['coh_car'], out['f_coh'] = _coerencia(Wc)
        out['coh_csd'] = _coerencia(W['csd'])[0]
        Wb = filter_data(Wc.astype(np.float64), FS, 1, 40, verbose='ERROR')
        out['cov_ep'] = (np.einsum('eit,ejt->eij', Wb, Wb) / Wb.shape[-1]).astype(np.float32)
        Pw, fw = np.abs(np.fft.rfft(Wc * np.hanning(Wc.shape[-1]), axis=-1)) ** 2, np.fft.rfftfreq(Wc.shape[-1], 1 / FS)
        out['e2'] = np.stack([np.stack([np.log(Pw[:, :, (fw >= lo) & (fw < hi)].mean(-1) + 1e-30), _entropia_perm(
            filter_data(Wc.astype(np.float64), FS, lo, min(hi, FS / 2 - 1), verbose='ERROR'))], -1)
            for lo, hi in BANDAS_E2], 2).astype(np.float32)
        v, (P, f) = Wc.var(-1).mean(0), s.espectro(Wc)
        p = np.median(P, 0).mean(0)
        out['A1'] = np.array([np.log(v.sum()), np.log(v[indices(BORDA8)].mean() / v[INTERNOS].mean()),
                              np.log10(p[(f >= 25) & (f <= 32)].mean() / p[(f >= 8) & (f <= 13)].mean())])
        if taghi:
            out.update(_taghi(s.car(x), segs))
        alvo.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(alvo.with_suffix('.tmp.npz'), **out)
        alvo.with_suffix('.tmp.npz').replace(alvo)
        return None
    except Exception as exc:
        return f'{bloco}: {type(exc).__name__}'


def extrair(cfg, blocos, jobs=8, limite=None, taghi=True):
    isolar_threads()
    tarefas = []
    for b in blocos:
        m = _meta(cfg, b)
        sub = _taghi_sub(m, b)
        cods = list(m['codigo']) if not limite else [c for g in (1, 0) for c in m.loc[m.y == g, 'codigo'][:limite]]
        tarefas += [(cfg['pacote_dir'], Path(cfg['saida_dir']) / 'replicacoes', b, c, taghi and c in sub) for c in cods]
    with Pool(jobs) as pool:
        erros = [e for e in pool.imap_unordered(_uma, tarefas, chunksize=2) if e]
    print(f'replicações: {len(tarefas)} registros, {len(erros)} erros', pd.Series(erros).value_counts().to_dict(), flush=True)


@lru_cache(None)
def _tab(bloco):
    m = _meta(_CFG, bloco)
    return m[[_dir('replicacoes', bloco, f'{c}.npz').exists() for c in m['codigo']]].reset_index(drop=True)


@lru_cache(None)
def _campo(bloco, chave):
    L = []
    for c in _tab(bloco)['codigo']:
        with np.load(_dir('replicacoes', bloco, f'{c}.npz')) as z:
            L.append(z[chave] if chave in z.files else None)
    try:
        return np.stack(L)
    except (ValueError, TypeError):
        return L


def _disponiveis():
    return [b for b in NOME if _dir('replicacoes', b).exists() and len(_tab(b)) >= 20 and _tab(b).y.nunique() == 2]


def _auc(y, sc):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, sc))


def _ic(y, s1, s2=None, q=(2.5, 97.5)):
    rng, v = np.random.default_rng(0), []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) == 2:
            v.append(_auc(y[i], s1[i]) - (_auc(y[i], s2[i]) if s2 is not None else 0))
    return tuple(np.percentile(v, q))


def _g_cond(y, v, Z, n=500):
    def g(yy, vv, ZZ):
        X = np.column_stack([np.ones(len(yy)), yy, ZZ])
        b = np.linalg.lstsq(X, vv, rcond=None)[0]
        df = len(yy) - X.shape[1]
        sd = np.sqrt(((vv - X @ b) ** 2).sum() / df)
        return (1 - 3 / (4 * df - 1)) * b[1] / sd if sd > 0 else np.nan
    ok = np.isfinite(v)
    y, v, Z = y[ok], v[ok], Z[ok]
    rng = np.random.default_rng(0)
    bs = [g(y[i], v[i], Z[i]) for i in (rng.integers(0, len(y), len(y)) for _ in range(n)) if len(set(y[i])) == 2]
    return float(g(y, v, Z)), *np.percentile(bs, [2.5, 97.5]), int(ok.sum())


def _resid(X, Z, tr, y):
    if Z.shape[1] == 0:
        return X
    c = tr[y[tr] == 0]
    B = np.linalg.lstsq(np.column_stack([np.ones(len(c)), Z[c]]), X[c], rcond=None)[0]
    return X - np.column_stack([np.ones(len(X)), Z]) @ B


def _padr(X, tr):
    return (X - X[tr].mean(0)) / (X[tr].std(0) + 1e-9)


def _modelo(clf):
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.svm import SVC
    return {'LR': lambda: LogisticRegression(C=1.0, max_iter=3000, class_weight='balanced'),
            'LDA': lambda: LinearDiscriminantAnalysis(solver='lsqr', shrinkage='auto'),
            'KNN3': lambda: KNeighborsClassifier(3), 'SVM': lambda: SVC(kernel='linear', C=1.0),
            'QSVM': lambda: SVC(kernel='poly', degree=2, C=1.0),
            'LR1': lambda: LogisticRegression(C=1.0, max_iter=2000),
            'RBF': lambda: SVC(kernel='rbf', C=1.0, gamma='scale', class_weight='balanced')}[clf]()


def _esc(m, X):
    return m.decision_function(X) if hasattr(m, 'decision_function') else m.predict_proba(X)[:, 1]


def _cv(Xfn, y, Z, parts, clf='LR', resid=True):
    dec = np.full((max(1, len(parts) // 10), len(y)), np.nan)
    for k, (tr, te) in enumerate(parts):
        X = Xfn(tr, y)
        X = _padr(_resid(X, Z, tr, y) if resid else X, tr)
        dec[k // 10, te] = _esc(_modelo(clf).fit(X[tr], y[tr]), X[te])
    a = [_auc(y, r) for r in dec]
    return float(np.mean(a)), dec.mean(0), a


def _csp_W(C, y, m):
    _, V = eigh(C[y == 1].mean(0), C[y == 1].mean(0) + C[y == 0].mean(0))
    return np.concatenate([V[:, :m], V[:, -m:]], 1)


def _csp_f(C, W):
    v = np.einsum('ncd,ck,dk->nk', C, W, W)
    return np.log(v / v.sum(1, keepdims=True))


def _tang(C):
    L = est.logm_spd(est.regularizar(C))
    return lambda tr, y: est.vetorizar(L, tr)


def _fbcsp(cov, vazado=None, m=2):
    Cs = {b: est.regularizar(cov[b]) for b in FB}
    if vazado is not None:
        X = np.concatenate([_csp_f(Cs[b], _csp_W(Cs[b], vazado, m)) for b in FB], 1)
        return lambda tr, y: X
    return lambda tr, y: np.concatenate([_csp_f(Cs[b], _csp_W(Cs[b][tr], y[tr], m)) for b in FB], 1)


def _d(bloco, fonte):
    t = _tab(bloco)
    d = {'y': t.y.values, 'Z': est.C0(t)}
    if fonte in ('bandas', 'ref'):
        ch = BANDAS if fonte == 'bandas' else ('larga', 'rest', 'csd')
        d['cov'] = {b: _campo(bloco, f'cov_{b}').astype(np.float64) for b in ch}
        ok = np.all([np.isfinite(v).all((1, 2)) for v in d['cov'].values()], 0)
        d = {'y': d['y'][ok], 'Z': d['Z'][ok], 'cov': {k: v[ok] for k, v in d['cov'].items()}, 'n_excl': int((~ok).sum())}
    elif fonte == 'coh':
        d.update(coh={r: _campo(bloco, f'coh_{r}') for r in ('car', 'csd')}, f=_campo(bloco, 'f_coh')[0])
    elif fonte == 'ep':
        p = Path(_CFG['pacote_dir']) / bloco / 'emg_70_95.csv'
        e = pd.read_csv(p).set_index('codigo')['emg_70_95'] if p.exists() else pd.Series(dtype=float)
        d.update(cov_ep=_campo(bloco, 'cov_ep'), emg=e.reindex(t['codigo']).values.astype(float))
    else:
        d['e2'] = _campo(bloco, 'e2')
    return d


def _yp(d, p):
    return d['y'] if p is None else np.random.default_rng(1000 + p).permutation(d['y'])


def _parts(y, p):
    return est.particoes(y, N_REP if p is None else 1, 0 if p is None else p, k=10)


def _lin(id_, b, estim, **kw):
    return {'id': id_, 'bloco': b, 'coluna': NOME.get(b, b), 'estimador': estim, **kw}


def _e1(b, p=None):
    d = _d(b, 'bandas')
    y, Z = _yp(d, p), d['Z']
    parts, fn = _parts(y, p), _fbcsp(d['cov'])
    if p is not None:
        return {c: _cv(fn, y, Z, parts, c)[0] for c in ('KNN3', 'LR')}
    aZ = _cv(lambda tr, yy: Z, y, Z, parts, resid=False)[0] if Z.shape[1] else 0.5
    out = []
    for c in ('KNN3', 'LR'):
        a, sc, reps = _cv(fn, y, Z, parts, c)
        aXZ = _cv(lambda tr, yy: np.column_stack([fn(tr, yy), Z]), y, Z, parts, c, False)[0]
        lo, hi = _ic(y, sc)
        out.append(_lin('E1', b, f'fbcsp_m2_{c}', auc=a, ic_lo=lo, ic_hi=hi, auc_dp_rep=np.std(reps),
                        auc_base=aZ, auc_mais_base=aXZ, ganho=aXZ - aZ, n=len(y), n_tdah=int(y.sum())))
    return out


def _e7(b, p=None):
    d = _d(b, 'bandas')
    y, Z = _yp(d, p), d['Z']
    parts = _parts(y, p)
    (a_c, s_c, _), (a_v, s_v, _) = (_cv(_fbcsp(d['cov'], v), y, Z, parts) for v in (None, y))
    if p is not None:
        return {'correto': a_c, 'vazado': a_v}
    lo, hi = _ic(y, s_v, s_c)
    return [_lin('E7', b, 'fbcsp_correto', auc=a_c, n=len(y)),
            _lin('E7', b, 'fbcsp_vazado', auc=a_v, n=len(y), delta_vazado_menos_correto=a_v - a_c, delta_lo=lo, delta_hi=hi)]


def _media_riemann(Cs, it=30, tol=1e-7):
    M = Cs.mean(0)
    for _ in range(it):
        w, V = np.linalg.eigh(M)
        Mih, Mh = (V * w ** -0.5) @ V.T, (V * w ** 0.5) @ V.T
        L = est.logm_spd(np.einsum('ij,njk,kl->nil', Mih, Cs, Mih)).mean(0)
        w2, V2 = np.linalg.eigh(L)
        M = Mh @ ((V2 * np.exp(w2)) @ V2.T) @ Mh
        if np.linalg.norm(L) < tol:
            break
    return M


def _dist_riemann(Cs, M):
    w, V = np.linalg.eigh(M)
    Mih = (V * w ** -0.5) @ V.T
    return np.linalg.norm(est.logm_spd(np.einsum('ij,njk,kl->nil', Mih, Cs, Mih)), axis=(1, 2))


def _m2(b, p=None):
    d = _d(b, 'bandas')
    y, Z = _yp(d, p), d['Z']
    parts, C = _parts(y, p), est.regularizar(d['cov']['larga'])
    dec = np.full((len(parts) // 10, len(y)), np.nan)
    for k, (tr, te) in enumerate(parts):
        dec[k // 10, te] = (_dist_riemann(C[te], _media_riemann(C[tr][y[tr] == 0]))
                            - _dist_riemann(C[te], _media_riemann(C[tr][y[tr] == 1])))
    a_m, s_m = float(np.mean([_auc(y, r) for r in dec])), dec.mean(0)
    a_f, s_f, _ = _cv(_fbcsp(d['cov']), y, Z, parts, 'LDA', resid=False)
    if p is not None:
        return {'mdrm': a_m, 'fbcsp_lda': a_f}
    lo, hi = _ic(y, s_m, s_f, (5, 95))
    ml, mh = _ic(y, s_m)
    return [_lin('M2', b, 'mdrm', auc=a_m, ic_lo=ml, ic_hi=mh, n=len(y)),
            _lin('M2', b, 'fbcsp_m2_LDA', auc=a_f, n=len(y), delta_mdrm_menos_fbcsp=a_m - a_f, delta_ic90_lo=lo,
                 delta_ic90_hi=hi)]


def _m5(b, p=None):
    d = _d(b, 'ref')
    y, Z = _yp(d, p), d['Z']
    parts = _parts(y, p)
    r = {k: _cv(_tang(d['cov'][c]), y, Z, parts)[:2] for k, c in (('car', 'larga'), ('rest', 'rest'), ('csd', 'csd'))}
    if p is not None:
        return {k: v[0] for k, v in r.items()}
    out = []
    for k in r:
        kw = {}
        for o in (o for o in r if o != k):
            kw[f'delta_vs_{o}'] = r[k][0] - r[o][0]
            kw[f'delta_vs_{o}_lo'], kw[f'delta_vs_{o}_hi'] = _ic(y, r[k][1], r[o][1])
        out.append(_lin('M5', b, f'tangente_{k}', auc=r[k][0], n=len(y), n_excluidos=d['n_excl'], **kw))
    return out


def _m3_fora(Cs, ys, Cd, yd, modo):
    if modo == 'csp':
        W = _csp_W(Cs, ys, 3)
        Xs, Xd = _csp_f(Cs, W), _csp_f(Cd, W)
    else:
        Ls, Ld = est.logm_spd(Cs), est.logm_spd(Cd)
        Xs = est.vetorizar(Ls, np.arange(len(Ls)))
        Xd = est.vetorizar(np.concatenate([Ld, Ls]), np.arange(len(Ld)) if modo == 'ts_recentrado' else
                           np.arange(len(Ld), len(Ld) + len(Ls)))[:len(Ld)]
    mu, sd = Xs.mean(0), Xs.std(0) + 1e-9
    m = _modelo('LR').fit((Xs - mu) / sd, ys)
    Xd = (Xd - Xd.mean(0)) / (Xd.std(0) + 1e-9) if modo == 'ts_recentrado' else (Xd - mu) / sd
    return _auc(yd, m.decision_function(Xd))


def _m3(par, p=None):
    out = [] if p is None else {}
    for src, dst in (par, par[::-1]):
        ds, dd = _d(src, 'bandas'), _d(dst, 'bandas')
        Cs, Cd = est.regularizar(ds['cov']['larga']), est.regularizar(dd['cov']['larga'])
        ys = _yp(ds, p)
        if p is not None:
            out.update({f'{src}->{dst}:{m}': _m3_fora(Cs, ys, Cd, dd['y'], m) for m in ('csp', 'ts', 'ts_recentrado')})
            continue
        pt = _parts(ys, None)
        dentro = {'csp': _cv(lambda tr, y: _csp_f(Cs, _csp_W(Cs[tr], y[tr], 3)), ys, ds['Z'], pt, resid=False)[0],
                  'ts': _cv(_tang(ds['cov']['larga']), ys, ds['Z'], pt, resid=False)[0]}
        for m in ('csp', 'ts', 'ts_recentrado'):
            fora, dn = _m3_fora(Cs, ys, Cd, dd['y'], m), dentro['csp' if m == 'csp' else 'ts']
            out.append({'id': 'M3', 'bloco': f'{src}->{dst}', 'coluna': NOME[dst], 'origem': src, 'destino': dst,
                        'estimador': m, 'auc_dentro': dn, 'auc_fora': fora, 'queda': dn - fora, 'n_origem': len(ys),
                        'n_destino': len(dd['y'])})
    return out


def _filtro_principal(C, y):
    Cr = est.regularizar(C)
    lam, V = eigh(Cr[y == 1].mean(0), Cr[y == 1].mean(0) + Cr[y == 0].mean(0))
    j = int(np.argmax(np.abs(lam - 0.5)))
    return V[:, j] / np.linalg.norm(V[:, j]), float(lam[j])


def _e6(b, p=None):
    d, out = _d(b, 'ep'), []
    for igualar in (False, True):
        sel = np.flatnonzero([len(c) >= 10 for c in d['cov_ep']]) if igualar else np.arange(len(d['y']))
        y, ce = d['y'][sel], [d['cov_ep'][i].astype(np.float64) for i in sel]
        nm = min(len(c) for c in ce) if igualar else None

        def treino(tr):
            C = [ce[i][:nm] for i in tr]
            return np.concatenate(C), np.concatenate([[y[i]] * len(c) for i, c in zip(tr, C)]), np.concatenate(
                [c[:, [0, 1], [0, 1]].sum(1) for c in C])
        coss, lams = [], []
        for tr, _ in _parts(y, None):
            C, yy, fp = treino(tr)
            w, lam = _filtro_principal(C, yy)
            mask = np.arange(len(C)) != int(np.argmax(fp))
            coss.append(abs(float(w @ _filtro_principal(C[mask], yy[mask])[0])))
            lams.append(lam)
        w, _ = _filtro_principal(*treino(np.arange(len(y)))[:2])
        proj = np.array([np.log(np.einsum('eij,i,j->e', est.regularizar(c[:nm]), w, w)).mean() for c in ce])
        vfp = np.array([np.log(c[:nm][:, [0, 1], [0, 1]].sum(1)).mean() for c in ce])
        emg = d['emg'][sel]
        ok = np.isfinite(emg)
        coss = np.array(coss)
        out.append(_lin('E6', b, 'csp_principal' + ('_epocas_igualadas' if igualar else ''),
                        frac_cos_abaixo_0_90=float((coss < 0.9).mean()), cos_mediano=float(np.median(coss)),
                        lambda_mediano=float(np.median(lams)), r_variancia_frontal=float(np.corrcoef(proj, vfp)[0, 1]),
                        r_emg=float(np.corrcoef(proj[ok], emg[ok])[0, 1]) if ok.sum() > 10 else np.nan, n=len(y),
                        n_epocas_media=float(np.mean([len(c) for c in ce]))))
    return out


def _e2(b, p=None):
    from sklearn.metrics import accuracy_score
    from sklearn.model_selection import GroupKFold, KFold
    d = _d(b, 'e2')
    yp = _yp(d, p)
    X = np.nan_to_num(np.concatenate([e[:15].reshape(len(e[:15]), -1) for e in d['e2']]).astype(np.float64))
    g = np.concatenate([[i] * len(e[:15]) for i, e in enumerate(d['e2'])])
    yy, seed, r = yp[g], (p or 0), {}
    for esq in ('epoca', 'participante'):
        if esq == 'epoca':
            dobras = KFold(10, shuffle=True, random_state=seed).split(X)
        else:
            mapa = {c: i for i, c in enumerate(np.random.default_rng(seed).permutation(np.unique(g)))}
            dobras = GroupKFold(10).split(X, yy, np.array([mapa[c] for c in g]))
        pred, dec = np.full(len(yy), -1), np.full(len(yy), np.nan)
        for tr, te in dobras:
            Xs = _padr(X, tr)
            m = _modelo('RBF').fit(Xs[tr], yy[tr])
            pred[te], dec[te] = m.predict(Xs[te]), m.decision_function(Xs[te])
        r[f'acc_epoca_{esq}'] = float(accuracy_score(yy, pred))
        r[f'auc_participante_{esq}'] = _auc(yp, np.array([dec[g == c].mean() for c in range(len(yp))]))
    return r if p is not None else [_lin('E2', b, 'psd_entropia_190_svmrbf', n=len(yp), n_epocas=len(yy), **r)]


def _z_coh(v):
    return np.arctanh(np.sqrt(np.clip(v.astype(np.float64), 0, 1 - 1e-6)))


def _e3e4(b, p=None):
    d, out = _d(b, 'coh'), []
    for ref in ('car', 'csd'):
        for id_, nome, pares, lo, hi, sg in E3 + E4:
            idx = [PAR.get(q, PAR.get(q[::-1])) for q in pares]
            v = _z_coh(d['coh'][ref][:, idx][:, :, (d['f'] >= lo) & (d['f'] < hi)].mean((1, 2)))
            g, glo, ghi, n = _g_cond(d['y'], v, d['Z'])
            out.append(_lin(id_, b, f'{nome}_{ref}', medida=nome, ref=ref, g=g, g_lo=glo, g_hi=ghi, sinal_previsto=sg,
                            n=n, passa=bool(np.sign(g) == sg and (glo > 0 if sg > 0 else ghi < 0))))
    return out


def _e5(b, p=None):
    from sklearn.decomposition import PCA
    d = _d(b, 'coh')
    y, Z, f = _yp(d, p), d['Z'], d['f']
    X = _z_coh(np.stack([d['coh']['car'][:, :, (f >= lo) & (f < lo + 2)].mean(-1) for lo in range(1, 33, 2)], -1)
               .reshape(len(y), -1))

    def fn(tr, yy):
        Xr = _padr(_resid(X, Z, tr, yy), tr)
        return PCA(min(30, len(tr) - 1), random_state=0).fit(Xr[tr]).transform(Xr)
    parts = _parts(y, p)
    a, sc, reps = _cv(fn, y, Z, parts, 'LDA', resid=False)
    if p is not None:
        return {'lda': a}
    aZ = _cv(lambda tr, yy: Z, y, Z, parts, resid=False)[0] if Z.shape[1] else 0.5
    aXZ = _cv(lambda tr, yy: np.column_stack([fn(tr, yy), Z]), y, Z, parts, resid=False)[0]
    lo, hi = _ic(y, sc)
    return [_lin('E5', b, 'coerencia_171x16_pca30_lda', auc=a, ic_lo=lo, ic_hi=hi, auc_base=aZ, auc_mais_base=aXZ,
                 ganho=aXZ - aZ, n=len(y), auc_dp_rep=float(np.std(reps)))]


TESTES = {'e1': _e1, 'e7': _e7, 'm2': _m2, 'm5': _m5, 'm3': _m3, 'e6': _e6, 'e3e4': _e3e4, 'e5': _e5, 'e2': _e2}
SEM_NULO, LIMITE_PERM = {'e3e4', 'e6'}, {'e2': 30}


def _nome(u):
    return u if isinstance(u, str) else '__'.join(u)


def _init(cfg):
    isolar_threads()
    _usar(cfg)


def _tarefa(args):
    t, u, p = args
    try:
        if p is None:
            return t, u, p, TESTES[t](u)
        alvo = _dir('agregados', 'perm', f'{t}_{_nome(u)}_{p:03d}.json')
        if not alvo.exists():
            alvo.parent.mkdir(parents=True, exist_ok=True)
            alvo.write_text(json.dumps(TESTES[t](u, p)))
        return t, u, p, None
    except Exception as exc:
        return t, u, p, f'erro: {exc!r}'[:300]


def _nulos(t, u):
    L = [json.loads(a.read_text()) for a in sorted(_dir('agregados', 'perm').glob(f'{t}_{_nome(u)}_*.json'))]
    return {k: np.array([x[k] for x in L if k in x]) for k in (L[0] if L else {})}


def _f(v, nd=2):
    return f'{v:.{nd}f}'.replace('.', ',') if np.isfinite(v) else 'nd'


def _p95(v):
    return float(np.percentile(v, 95)) if len(v) else np.nan


def _consolidar(R):
    st = []

    def add(id_, b, estado, texto, base, ok):
        st.append({'id': id_, 'coluna': NOME.get(b, b), 'estado': estado, 'texto': texto, 'base': base,
                   'criterio_passou': bool(ok)})
    for b in NOME:
        def q(t, est_):
            r = R.get(t, {}).get(b)
            return None if r is None else r[r.estimador == est_].iloc[0]
        k, nl = q('e1', 'fbcsp_m2_KNN3'), _nulos('e1', b)
        if k is not None and nl:
            n95 = _p95(nl['KNN3'])
            ok = k.auc >= 0.7 and k.ic_lo > 0.5 and n95 <= 0.55 and k.ganho >= 0.05
            add('E1', b, 'replica' if ok else 'não replica', f'KNN-3 {_f(k.auc)} [{_f(k.ic_lo)}; {_f(k.ic_hi)}]; ganho '
                f'{_f(k.ganho)}', f'nulo p95 {_f(n95)}; base {_f(k.auc_base)}', ok)
        v, c, nl = q('e7', 'fbcsp_vazado'), q('e7', 'fbcsp_correto'), _nulos('e7', b)
        if v is not None and nl:
            ok = v.delta_vazado_menos_correto >= 0.05
            add('E7', b, 'replica' if ok else 'não replica', f'vazado {_f(v.auc)} × correto {_f(c.auc)}; Δ '
                f'{_f(v.delta_vazado_menos_correto)}', f'nulo p95 vazado {_f(_p95(nl["vazado"]))}; correto '
                f'{_f(_p95(nl["correto"]))}', ok)
        m, f_, nl = q('m2', 'mdrm'), q('m2', 'fbcsp_m2_LDA'), _nulos('m2', b)
        if m is not None and nl:
            n95 = _p95(nl['mdrm'])
            eq = abs(f_.delta_mdrm_menos_fbcsp) <= 0.03 and f_.delta_ic90_lo >= -0.03 and f_.delta_ic90_hi <= 0.03
            aval = not (n95 > 0.55 or max(m.auc, f_.auc) <= n95)
            add('M2', b, ('replica' if eq else 'não replica') if aval else 'não avaliável', f'MDRM {_f(m.auc)} × FBCSP+LDA '
                f'{_f(f_.auc)}; Δ {_f(f_.delta_mdrm_menos_fbcsp)}', f'nulo p95 {_f(n95)}; TOST ±0,03', aval and eq)
        a, nl = {r: q('m5', f'tangente_{r}') for r in ('car', 'rest', 'csd')}, _nulos('m5', b)
        if a['car'] is not None and nl:
            n95 = {r: _p95(nl[r]) for r in a}
            dm = max(abs(a['car'].delta_vs_rest), abs(a['car'].delta_vs_csd), abs(a['rest'].delta_vs_csd))
            aval = not (max(n95.values()) > 0.55 or all(a[r].auc <= n95[r] for r in a))
            add('M5', b, ('replica' if dm <= 0.01 else 'não replica') if aval else 'não avaliável',
                '; '.join(f'{r.upper()} {_f(a[r].auc)}' for r in a) + f'; |Δ| máx {_f(dm, 3)}', 'critério |Δ| ≤ 0,01',
                aval and dm <= 0.01)
        a, g = q('e6', 'csp_principal'), q('e6', 'csp_principal_epocas_igualadas')
        if a is not None:
            rr = np.nanmax([abs(a.r_variancia_frontal), abs(a.r_emg) if np.isfinite(a.r_emg) else 0])
            ok = a.frac_cos_abaixo_0_90 >= 0.2 or rr >= 0.5
            add('E6', b, 'replica' if ok else 'não replica', f'cos < 0,90 em {_f(100 * a.frac_cos_abaixo_0_90, 0)}%; r frontal '
                f'{_f(a.r_variancia_frontal)}; r EMG {_f(a.r_emg)}', f'épocas igualadas: {_f(100 * g.frac_cos_abaixo_0_90, 0)}%',
                ok)
        m, nl = q('e2', 'psd_entropia_190_svmrbf'), _nulos('e2', b)
        if m is not None and nl:
            n95 = _p95(nl['auc_participante_participante'])
            num, gen = m.acc_epoca_epoca >= 0.95, m.auc_participante_participante >= 0.7 and n95 <= 0.55
            vaz = m.acc_epoca_epoca - m.acc_epoca_participante >= 0.15
            add('E2', b, 'replica' if num and gen else 'em parte' if num else 'não replica', f'acc época '
                f'{_f(m.acc_epoca_epoca)} × por participante {_f(m.acc_epoca_participante)}; AUC '
                f'{_f(m.auc_participante_participante)}', f'nulo p95 {_f(n95)}; vazamento '
                f'{"confirmado" if vaz else "não confirmado"}', num and gen)
        r = R.get('e3e4', {}).get(b)
        if r is not None:
            for id_, ms in (('E3', ('teta_frontal_homologo', 'alfa_frontal_longo')), ('E4', ('alfa_baixo_todos',
                                                                                          'alfa_alto_todos'))):
                qq = {ref: r[(r.id == id_) & (r.ref == ref) & r.medida.isin(ms)].set_index('medida') for ref in ('car', 'csd')}
                n_ok = int(qq['car'].passa.sum())
                add(id_, b, 'replica' if n_ok == 2 else 'em parte' if n_ok else 'não replica',
                    '; '.join(f'g {_f(qq["car"].loc[x].g)} [{_f(qq["car"].loc[x].g_lo)}; {_f(qq["car"].loc[x].g_hi)}]'
                              for x in ms), 'CSD: ' + '; '.join(f'g {_f(qq["csd"].loc[x].g)}' for x in ms), n_ok == 2)
        m, nl = q('e5', 'coerencia_171x16_pca30_lda'), _nulos('e5', b)
        if m is not None and nl:
            n95 = _p95(nl['lda'])
            ok = m.auc >= 0.7 and n95 <= 0.55 and m.ganho >= 0.05
            add('E5', b, 'replica' if ok else 'não replica', f'AUC {_f(m.auc)} [{_f(m.ic_lo)}; {_f(m.ic_hi)}]; ganho '
                f'{_f(m.ganho)}', f'nulo p95 {_f(n95)}; base {_f(m.auc_base)}', ok)
    for par, r in R.get('m3', {}).items():
        nl = _nulos('m3', par)
        for (src, dst), q_ in r.groupby(['origem', 'destino']):
            q_ = q_.set_index('estimador')
            ok = q_.loc['ts'].auc_fora - q_.loc['csp'].auc_fora >= 0.05 and q_.loc['ts'].queda < q_.loc['csp'].queda
            add('M3', dst, 'replica' if ok else 'não replica', f'de {src}: TS {_f(q_.loc["ts"].auc_fora)} × CSP '
                f'{_f(q_.loc["csp"].auc_fora)} fora', 'nulos p95 ' + ', '.join(
                    f'{m_} {_f(_p95(nl.get(f"{src}->{dst}:{m_}", [])))}' for m_ in ('csp', 'ts', 'ts_recentrado')), ok)
    return pd.DataFrame(st)


def _espaciais(jobs, n_perm):
    blocos = _disponiveis()
    unid = {t: [p for p in PARES_M3 if set(p) <= set(blocos)] if t == 'm3' else blocos for t in TESTES}
    reais = [(t, u, None) for t in TESTES for u in unid[t]]
    nulos = [(t, u, p) for t in TESTES if t not in SEM_NULO for u in unid[t] for p in range(min(n_perm, LIMITE_PERM.get(t, n_perm)))]
    R = {}
    with Pool(jobs, initializer=_init, initargs=(dict(_CFG),)) as pool:
        for t, u, _, r in pool.imap_unordered(_tarefa, reais + nulos, chunksize=1):
            if isinstance(r, str):
                print('ERRO', t, _nome(u), r, flush=True)
            elif r is not None:
                R.setdefault(t, {})[u] = pd.DataFrame(r)
    todos = [x for d in R.values() for x in d.values()]
    return (pd.concat(todos, ignore_index=True) if todos else pd.DataFrame()), _consolidar(R)


def _vaz_dados(b):
    cov, nl, chi, cods = _campo(b, 'cov30'), _campo(b, 'nl'), _campo(b, 'chi_ep'), _tab(b)
    ii = [i for i in range(len(cods)) if cov[i] is not None]
    n = [len(cov[i]) for i in ii]
    return {'cov': np.concatenate([cov[i] for i in ii]).astype(np.float64),
            'nl': np.concatenate([nl[i].reshape(len(nl[i]), -1) for i in ii]).astype(np.float64),
            'chi': np.concatenate([chi[i] for i in ii])[:, None], 'y': np.repeat(cods.y.values[ii], n),
            'g': np.repeat(np.arange(len(ii)), n)}


def _top_k(X, y, k=10):
    sc = [abs(_auc(y[o], X[o, j]) - 0.5) if (o := np.isfinite(X[:, j])).sum() > 10 and len(set(y[o])) == 2 else 0
          for j in range(X.shape[1])]
    return np.argsort(sc)[::-1][:k]


def _protocolo(args):
    from sklearn.model_selection import GroupKFold, KFold
    b, pipe, V, sem, nulo = args
    d = _VAZ[b]
    V1, V2, V3, V4 = V
    y, g = d['y'].copy(), d['g']
    if nulo:
        cri = np.unique(g)
        perm = dict(zip(cri, np.random.default_rng(10_000 + sem).permutation([y[g == c][0] for c in cri])))
        y = np.array([perm[c] for c in g])
    if V1:
        dobras = list(KFold(5, shuffle=True, random_state=sem).split(y))
    else:
        mapa = {c: i for i, c in enumerate(np.random.default_rng(sem).permutation(np.unique(g)))}
        dobras = list(GroupKFold(5).split(y, y, np.array([mapa[c] for c in g])))

    def feats(ref):
        return d['chi'] if pipe == 'chi' else np.concatenate([_csp_vaz(d['cov'], y, ref), d['nl']], 1)

    def norm(X, tr):
        return np.nan_to_num((X - np.nanmean(X[tr], 0)) / (np.nanstd(X[tr], 0) + 1e-9))
    if V2:
        Xa = norm(feats(np.arange(len(y))), np.arange(len(y)))
        sa = _top_k(Xa, y) if pipe != 'chi' else np.arange(Xa.shape[1])
    clfs = ('LR1', 'KNN3', 'SVM', 'QSVM') if V3 else ('LR1',)
    esc = {k: np.full(len(y), np.nan) for k in clfs}
    for tr, te in dobras:
        if V2:
            X, sel = Xa, sa
        else:
            X = norm(feats(tr), tr)
            sel = _top_k(X[tr], y[tr]) if pipe != 'chi' else np.arange(X.shape[1])
        for k in clfs:
            esc[k][te] = _esc(_modelo(k).fit(X[tr][:, sel], y[tr]), X[te][:, sel])
    cri = np.unique(g)
    aucs = {k: _auc(y, sc) if V4 else _auc([y[g == c][0] for c in cri], [np.nanmean(sc[g == c]) for c in cri])
            for k, sc in esc.items() if np.isfinite(sc).all()}
    melhor = max(aucs, key=aucs.get) if V3 else 'LR1'
    return {'bloco': b, 'pipeline': pipe, 'V1_divisao_epoca': V1, 'V2_selecao_toda': V2, 'V3_clf_no_teste': V3,
            'V4_metrica_epoca': V4, 'semente': sem, 'nulo': nulo, 'auc': aucs[melhor], 'classificador': melhor}


def _csp_vaz(cov, y, ref, m=2, reg=1e-6):
    Ws = []
    for b in range(cov.shape[1]):
        c, yr = cov[ref, b], y[ref]
        S1, S = c[yr == 1].mean(0), c[yr == 1].mean(0) + c[yr == 0].mean(0)
        V = eigh(S1, S + reg * np.trace(S) / len(S) * np.eye(len(S)))[1]
        Ws.append(np.concatenate([V[:, :m], V[:, -m:]], 1))
    return np.concatenate([np.log(np.clip((v := np.einsum('ecd,ck,dk->ek', cov[:, b], W, W)) / v.sum(1, keepdims=True),
                                          1e-12, None)) for b, W in enumerate(Ws)], 1)


_VAZ = {}


def _init_vaz(cfg, blocos):
    _init(cfg)
    for b in blocos:
        _VAZ[b] = _vaz_dados(b)


def _vazamento(jobs, sementes=5):
    blocos = [b for b in _disponiveis() if any(c is not None for c in _campo(b, 'cov30'))]
    if not blocos:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    tarefas = [(b, p, V, sm, n) for b in blocos for p in ('taghi', 'chi') for V in itertools.product((0, 1), repeat=4)
               for sm in range(sementes) for n in (False, True)]
    with Pool(jobs, initializer=_init_vaz, initargs=(dict(_CFG), blocos)) as pool:
        df = pd.DataFrame(pool.map(_protocolo, tarefas, chunksize=4))
    ef = []
    for (b, p, nulo), d in df.groupby(['bloco', 'pipeline', 'nulo']):
        m = d.groupby(['V1_divisao_epoca', 'V2_selecao_toda', 'V3_clf_no_teste', 'V4_metrica_epoca']).auc.mean()
        X, lin = np.array(list(m.index)) * 2 - 1, {'bloco': b, 'pipeline': p, 'nulo': nulo,
                                                    'auc_tudo_correto': m.loc[(0, 0, 0, 0)], 'auc_tudo_vazando': m.loc[(1, 1, 1, 1)]}
        for r in (1, 2):
            for c in itertools.combinations(range(4), r):
                lin[':'.join(f'V{i + 1}' for i in c)] = float((X[:, list(c)].prod(1) * m.values).mean() * 2)
        ef.append(lin)
    esc = {'tudo vazando': ('taghi', (1, 1, 1, 1)), 'só divisão por época': ('taghi', (1, 0, 0, 0)),
           'só seleção antes da divisão': ('taghi', (0, 1, 0, 0)), 'só classificador escolhido no teste': ('taghi', (0, 0, 1, 0)),
           'só métrica por época': ('taghi', (0, 0, 0, 1)), 'tudo correto': ('taghi', (0, 0, 0, 0)),
           'nosso χ, tudo correto': ('chi', (0, 0, 0, 0))}
    cols = ['V1_divisao_epoca', 'V2_selecao_toda', 'V3_clf_no_teste', 'V4_metrica_epoca']
    escada = [{'bloco': b, 'protocolo': k, 'nulo': nulo, 'auc_media': d.auc.mean()} for k, (p, V) in esc.items()
              for (b, nulo), d in df[(df.pipeline == p) & (df[cols].values == V).all(1)].groupby(['bloco', 'nulo'])]
    return df, pd.DataFrame(ef), pd.DataFrame(escada)


def _referencias():
    out = []
    for b in _disponiveis():
        t = _tab(b)
        y, Z = t.y.values, est.C0(t)
        for r in REFS:
            g, lo, hi = est.g_boot(y, np.nanmedian(_campo(b, f'chi_{r}'), 1), Z, B=500, seed=42)
            out.append({'bloco': b, 'tipo': 'chi', 'ref': r, 'U40': 1, 'chave': 'global', 'g': g, 'g_lo': lo, 'g_hi': hi})
    return pd.DataFrame(out)


def _icc(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 20:
        return np.nan
    X = np.column_stack([a[ok], b[ok]])
    n, gm = len(X), X.mean()
    msr = 2 * ((X.mean(1) - gm) ** 2).sum() / (n - 1)
    mse = ((X - X.mean(1, keepdims=True) - X.mean(0, keepdims=True) + gm) ** 2).sum() / (n - 1)
    return float((msr - mse) / (msr + mse))


def _banda(f, P, lo, hi):
    return P[..., (f >= lo) & (f < hi)].mean(-1)


def _espectrais():
    rep, esp, icc, alfa = {}, {}, [], {}
    for b in _disponiveis():
        t = _tab(b)
        y, Z, idade, ctrl = t.y.values, est.C0(t), t.idade.values.astype(float), t.y.values == 0
        chi, chi30 = np.nanmedian(_campo(b, 'chi_car'), 1), np.nanmedian(_campo(b, 'chi30'), 1)
        ap, f = _campo(b, 'alfa_post'), _campo(b, 'f_psd')[0]
        iaf, cg, alt = ap.T
        Pc, Pn = _campo(b, 'psd_car').astype(float), _campo(b, 'psd_nativa').astype(float)
        cm, am = _campo(b, 'chi_metades'), _campo(b, 'alfa_metades')
        icc += [{'bloco': b, 'medida': 'chi_global_splithalf', 'icc': _icc(np.nanmedian(cm[:, 0], 1), np.nanmedian(cm[:, 1], 1))},
                {'bloco': b, 'medida': 'IAF_splithalf', 'icc': _icc(am[:, 0, 0], am[:, 1, 0])},
                {'bloco': b, 'medida': 'altura_alfa_splithalf', 'icc': _icc(am[:, 0, 2], am[:, 1, 2])}]
        alfa[b] = dict(zip(t.codigo, np.log10(_banda(f, Pc[:, IPOST], 8, 13).mean(1))))
        if np.isfinite(idade).all():
            for nome, v in (('chi', chi), ('IAF', iaf)):
                ok = ctrl & np.isfinite(v)
                rho, se = stats.spearmanr(idade[ok], v[ok])[0], 1.06 / np.sqrt(max(ok.sum() - 3, 1))
                rep[f'R1R2|{b}|{nome}'] = {'rho_idade_ctrl': float(rho), 'ic': [float(np.tanh(np.arctanh(rho) + k * 1.96 * se))
                                                                                 for k in (-1, 1)], 'n': int(ok.sum())}
                r, tq = stats.linregress(idade[ok], v[ok]), stats.t.ppf(.975, ok.sum() - 2)
                rep[f'R2|{b}|{nome}'] = {'inclinacao': float(r.slope), 'ic': [float(r.slope - tq * r.stderr),
                                                                               float(r.slope + tq * r.stderr)],
                                         'n': int(ok.sum()), 'faixa_idade': [float(idade[ok].min()), float(idade[ok].max())]}
            lo_q, hi_q = np.percentile(idade[ctrl], [20, 80])
            for nome, sel in (('novos', ctrl & (idade <= lo_q)), ('velhos', ctrl & (idade >= hi_q))):
                esp[f'{b}|{nome}'] = np.log10(Pc[sel][:, IPOST].mean(1)).mean(0)
                rep[f'R1|{b}|{nome}'] = {'idade_media': float(idade[sel].mean()), 'n': int(sel.sum())}
            esp['f'] = f
        for refn, P in (('car', Pc), ('nativa', Pn)):
            tbr = np.log10(_banda(f, P[:, CZ], 4, 8) / _banda(f, P[:, CZ], 13, 30))
            ok = np.isfinite(tbr) & np.isfinite(chi) & np.isfinite(cg)
            r4 = {'n': int(ok.sum())}
            for nome, extra in (('so_covariaveis', []), ('mais_chi', [chi]), ('mais_chi_iaf', [chi, cg])):
                g, glo, ghi = est.g_boot(y[ok], tbr[ok], np.column_stack([Z[ok]] + [e[ok] for e in extra]), B=500, seed=0)
                r4[nome] = {'g': g, 'ic': [glo, ghi]}
            oka = ok & np.isfinite(alt)
            r4.update(rho_tbr_chi=stats.spearmanr(tbr[ok], chi[ok])[0], rho_tbr_iaf=stats.spearmanr(tbr[ok], cg[ok])[0],
                      rho_tbr_altura_alfa=stats.spearmanr(tbr[oka], alt[oka])[0])
            rep[f'R4R5|{b}|{refn}'] = r4
        g, lo, hi = est.g_boot(y, chi30, Z, B=500, seed=0)
        rep[f'R3|{b}'] = {'g_chi_fixed_1_30': g, 'ic': [lo, hi], 'replica_panda': bool(lo < 0 < hi and abs(g) < 0.2)}
        dif = (chi - chi30)[np.isfinite(chi - chi30)]
        rng = np.random.default_rng(0)
        bm = [np.median(rng.choice(dif, len(dif))) for _ in range(1000)]
        rep[f'R7|{b}'] = {'mediana_chi_U40_menos_U30': float(np.median(dif)), 'ic': list(np.percentile(bm, [2.5, 97.5])),
                          'chi_U30_mediana': float(np.nanmedian(chi30)), 'chi_U40_mediana': float(np.nanmedian(chi)),
                          'wilcoxon_p': float(stats.wilcoxon(dif)[1]), 'n': len(dif)}
    for banco in ('hbn', 'tdbrain'):
        bf, ba = f'{banco}_fechado', f'{banco}_aberto'
        if bf in alfa and ba in alfa:
            comum = sorted(set(alfa[bf]) & set(alfa[ba]))
            dif = np.array([alfa[bf][c] - alfa[ba][c] for c in comum])
            dif = dif[np.isfinite(dif)]
            rep[f'R6|{banco}'] = {'n_pares': len(dif), 'log10_OF_menos_OA_mediana': float(np.median(dif)),
                                  'frac_OF_maior': float((dif > 0).mean()), 'wilcoxon_p': float(stats.wilcoxon(dif)[1])}
            for nome, k in (('IAF', 0), ('altura_alfa', 2)):
                va = {c: v[k] for c, v in zip(_tab(bf).codigo, _campo(bf, 'alfa_post'))}
                vb = {c: v[k] for c, v in zip(_tab(ba).codigo, _campo(ba, 'alfa_post'))}
                icc.append({'bloco': banco, 'medida': f'{nome}_OAxOF', 'icc': _icc(np.array([vb[c] for c in comum]),
                                                                                   np.array([va[c] for c in comum]))})
    return rep, esp, pd.DataFrame(icc)


def _vereditos(rep):
    V = []

    def add(rid, b, valor, ok):
        V.append({'id': rid, 'banco': NOME.get(b, b), 'nosso_valor': valor,
                  'veredito': ok if isinstance(ok, str) else 'replica' if ok else 'não replica'})
    for b in NOME:
        if f'R1R2|{b}|chi' in rep:
            c, i = rep[f'R1R2|{b}|chi'], rep[f'R1R2|{b}|IAF']
            add('R1', b, f'ρχ {c["rho_idade_ctrl"]:+.2f}; ρIAF {i["rho_idade_ctrl"]:+.2f}', c['ic'][1] < 0 and i['ic'][0] > 0)
            for nome, autor in (('chi', -0.0097), ('IAF', 0.216)):
                r = rep[f'R2|{b}|{nome}']
                cri = r['ic'][1] < 0 if autor < 0 else r['ic'][0] > 0
                cont = r['ic'][0] <= autor <= r['ic'][1] or r['ic'][0] <= autor / 2 <= r['ic'][1]
                v = 'fora da faixa (adultos)' if r['faixa_idade'][0] >= 18 else 'em parte' if cri and not cont else cri
                add(f'R2_{nome}', b, f'{r["inclinacao"]:+.4f} [{r["ic"][0]:+.4f}, {r["ic"][1]:+.4f}]', v)
        if f'R3|{b}' in rep:
            r = rep[f'R3|{b}']
            add('R3', b, f'g {r["g_chi_fixed_1_30"]:+.2f} [{r["ic"][0]:+.2f}, {r["ic"][1]:+.2f}]', r['replica_panda'])
            r = rep[f'R4R5|{b}|car']
            add('R4', b, f'g TBR corrigida por χ {r["mais_chi"]["g"]:+.2f}', r['mais_chi']['ic'][0] <= 0 <= r['mais_chi']['ic'][1])
            add('R5', b, f'ρ(TBR, χ) {r["rho_tbr_chi"]:+.2f}', abs(r['rho_tbr_chi']) >= 0.3)
            r = rep[f'R7|{b}']
            add('R7', b, f'χ(U40) − χ(U30) {r["mediana_chi_U40_menos_U30"]:+.3f}', r['ic'][0] > 0 or r['ic'][1] < 0)
    for banco in ('hbn', 'tdbrain'):
        if f'R6|{banco}' in rep:
            r = rep[f'R6|{banco}']
            add('R6', banco, f'{100 * r["frac_OF_maior"]:.1f}% com alfa OF > OA', r['frac_OF_maior'] >= 0.9)
    return pd.DataFrame(V)


def _lote(t):
    return pd.Series((t['duracao_s'].values < 120.5).astype(float), index=t['codigo'])


def _subestudo(t):
    try:
        mp = pd.read_csv(_CFG['mapa_ids_csv'])
        mp = mp[mp.banco.str.startswith('tdbrain')].drop_duplicates('codigo').set_index('codigo').id_original
        x = pd.read_excel(_CFG['tdbrain_xlsx']).drop_duplicates('TDBRAIN_ID').set_index('TDBRAIN_ID')['Dataset']
        return pd.Series(x.reindex(mp.reindex(t['codigo']).values).astype(str).str.contains('NF', case=False).astype(float)
                         .values, index=t['codigo'])
    except Exception:
        return pd.Series(np.nan, index=t['codigo'])


def _aninhado(fe, Z, y, parts):
    dec = {k: np.full((len(parts) // 10, len(y)), np.nan) for k in ('esp', 'esp_Z', 'Z')}
    for k, (tr, te) in enumerate(parts):
        S = fe(tr, y)
        for nome, X in (('esp', S), ('esp_Z', np.column_stack([S, Z])), ('Z', Z)):
            if X.shape[1]:
                X = _padr(X, tr)
                dec[nome][k // 10, te] = _modelo('LR').fit(X[tr], y[tr]).decision_function(X[te])
    return {k: float(np.mean([_auc(y, d) for d in v])) for k, v in dec.items() if np.isfinite(v).all()}


def _braco_A():
    out = []
    for b in _disponiveis():
        d, t = _d(b, 'bandas'), _tab(b)
        ok = np.all([np.isfinite(_campo(b, f'cov_{k}').astype(float)).all((1, 2)) for k in BANDAS], 0)
        t = t[ok].reset_index(drop=True)
        y, Z, cov = d['y'], d['Z'], d['cov']['larga']
        todos = [i for i in range(19) if not (b.startswith('hbn') and i == CZ)]

        def tang(C):
            return _tang(C[:, todos][:, :, todos])
        A1 = _campo(b, 'A1')[ok]
        if b.endswith(('fechado', 'tarefa')) and np.isfinite(A1).all(1).sum() > MIN_A[0]:
            A1 = np.where(np.isfinite(A1), A1, np.nanmedian(A1, 0))
            r = _aninhado(lambda tr, yy: A1, Z, y, est.particoes(y, 5, k=10))
            out.append({'teste': 'A1', 'bloco': b, 'representacao': 'assinatura_nao_cerebral', 'canais': '—', 'n': len(y),
                        **r, 'ganho_sobre_C0': r['esp_Z'] - r['Z']})
        if not b.startswith('tdbrain'):
            continue
        L = pd.DataFrame({'lote_120s': _lote(t), 'subestudo_adhd_nf': _subestudo(t)}).reindex(t['codigo'])
        for alvo in ('lote_120s', 'subestudo_adhd_nf'):
            for grp in (1, 0):
                sel = (y == grp) & L[alvo].notna().values
                yy = L[alvo].values[sel].astype(int)
                if min((yy == 0).sum(), (yy == 1).sum()) >= MIN_A[1]:
                    r = _aninhado(tang(cov[sel]), np.zeros((sel.sum(), 0)), yy, est.particoes(yy, 5, k=10))
                    out.append({'teste': 'A5a', 'bloco': b, 'representacao': f'tangente_prediz_{alvo}',
                                'canais': f'grupo={grp}', 'n': int(sel.sum()), **r})
        rng, idade, sexo = np.random.default_rng(0), t.idade.values, t.sexo_f.values
        faixa, esc = np.floor(idade / 5), []
        for (fx, sx), _ in pd.DataFrame({'f': faixa, 's': sexo})[y == 1].groupby(['f', 's']):
            cand = np.flatnonzero((y == 0) & (faixa == fx) & (sexo == sx))
            esc += list(rng.choice(cand, min(len(cand), int(((y == 1) & (faixa == fx) & (sexo == sx)).sum())), replace=False))
        sel = np.isin(np.arange(len(y)), esc) | ((y == 1) & (idade <= 55.6))
        if min(y[sel].sum(), (1 - y[sel]).sum()) >= MIN_A[1]:
            r = _aninhado(tang(cov[sel]), Z[sel], y[sel], est.particoes(y[sel], 5, k=10))
            out.append({'teste': 'A5c_pareado', 'bloco': b, 'representacao': 'tangente', 'canais': 'pareado_idade_sexo',
                        'n': int(sel.sum()), 'n_tdah': int(y[sel].sum()), **r, 'ganho_sobre_C0': r['esp_Z'] - r['Z']})
    return pd.DataFrame(out)


def rodar(cfg, jobs=8, n_perm=None, publicar=True, vazamento=True, espaciais=True):
    isolar_threads()
    _usar(cfg)
    _tab.cache_clear()
    _campo.cache_clear()
    destinos = [_dir('agregados')] + ([RESULTADOS / 'replicacoes'] if publicar else [])
    saidas = {'referencias_chi.csv': _referencias(), 'braco_A.csv': _braco_A()}
    rep, esp, icc = _espectrais()
    saidas.update({'icc.csv': icc, 'vereditos.csv': _vereditos(rep), 'espectros.csv': pd.DataFrame(esp)})
    if vazamento:
        saidas['vazamento_todos.csv'], saidas['efeitos.csv'], saidas['escada.csv'] = _vazamento(jobs)
    if espaciais:
        saidas['espaciais.csv'], saidas['status_em.csv'] = _espaciais(jobs, N_PERM if n_perm is None else n_perm)
    for d in destinos:
        d.mkdir(parents=True, exist_ok=True)
        (d / 'replicacoes_espectrais.json').write_text(json.dumps(rep, indent=1, ensure_ascii=False, default=float),
                                                       encoding='utf-8')
        for nome, df in saidas.items():
            df.to_csv(d / nome, index=False)
    print('replicações agregadas em', ', '.join(str(d) for d in destinos), flush=True)


def _ler(cfg, nome):
    for d in (Path(cfg['saida_dir']) / 'agregados', RESULTADOS / 'replicacoes'):
        if (d / nome).exists():
            return d / nome


def figuras(cfg, destino=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman', 'DejaVu Serif'], 'font.size': 8.5,
                         'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 300, 'pdf.fonttype': 42})
    destino = Path(destino or FIGURAS)
    destino.mkdir(parents=True, exist_ok=True)
    cor = {'nativa': '#D55E00', 'car': '#0072B2', 'rest': '#009E73', 'csd': '#E69F00'}
    from matplotlib.ticker import FuncFormatter, ScalarFormatter
    virgula = FuncFormatter(lambda v, _: f'{v:g}'.replace('.', ',').replace('-', '−'))

    def br(v, nd=2):
        return f'{v:.{nd}f}'.replace('.', ',').replace('-', '−')

    def salvar(fig, nome):
        for a in fig.axes:
            for eixo in (a.xaxis, a.yaxis):
                if isinstance(eixo.get_major_formatter(), ScalarFormatter):
                    eixo.set_major_formatter(virgula)
        fig.savefig(destino / f'{nome}.png', bbox_inches='tight')
        fig.savefig(destino / f'{nome}.pdf', bbox_inches='tight')
        plt.close(fig)
    if p := _ler(cfg, 'referencias_chi.csv'):
        d = pd.read_csv(p)
        blocos = [b for b in NOME if b in set(d.bloco)]
        fig, ax = plt.subplots(figsize=(6.77, 3.3))
        for i, b in enumerate(blocos):
            for j, (r, mk) in enumerate(zip(REFS, 'os^D')):
                g, lo, hi = d[(d.bloco == b) & (d.ref == r)].iloc[0][['g', 'g_lo', 'g_hi']]
                y = i + (j - 1.5) * 0.17
                ax.plot([lo, hi], [y, y], color=cor[r], lw=1.2)
                ax.plot(g, y, mk, color=cor[r], ms=5.5, mfc=cor[r] if r == 'nativa' else 'white',
                        label=(r if r == 'nativa' else r.upper()) if i == 0 else None)
        ax.axvline(0, color='k', lw=0.8)
        ax.axhspan(2.5, 4.5, color='#7F7F7F', alpha=0.08, lw=0)
        ax.set_yticks(range(len(blocos)), [NOME[b] for b in blocos])
        ax.set_ylim(len(blocos) - 0.5, -0.5)
        ax.set_xlabel('g de Hedges do χ (TDAH − controle; IC 95 %)')
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=4, frameon=False)
        salvar(fig, 'R0_referencia')
    pe, pj = _ler(cfg, 'espectros.csv'), _ler(cfg, 'replicacoes_espectrais.json')
    if pe and pj and pe.stat().st_size > 5:
        e, rep = pd.read_csv(pe), json.loads(pj.read_text(encoding='utf-8'))
        tro = np.genfromtxt(RAIZ / 'referencias' / 'literatura' / 'troendle2022_fig2.csv', delimiter=';', names=True)
        ft, f = 2 + 0.25 * np.arange(len(tro)), e['f'].values
        sel = (f >= 2) & (f <= 40)
        blocos = [b for b in ('hbn_fechado', 'hbn_aberto', 'tdbrain_fechado', 'tdbrain_aberto') if f'{b}|novos' in e]
        fig, ax = plt.subplots(2, max(len(blocos), 1), figsize=(6.77, 4.3), sharex=True, squeeze=False,
                               gridspec_kw={'height_ratios': [1.25, 1]})
        for j, b in enumerate(blocos):
            for g, c, rot in (('novos', '#E69F00', '20% mais novos'), ('velhos', '#0072B2', '20% mais velhos')):
                ax[0, j].plot(f[sel], e[f'{b}|{g}'][sel], color=c, lw=1.3, label=rot)
            rc, ri = rep[f'R1R2|{b}|chi'], rep[f'R1R2|{b}|IAF']
            ok = rc['ic'][1] < 0 and ri['ic'][0] > 0
            ax[0, j].set_title(f'{NOME[b]}\n' + ('replica' if ok else 'adultos' if b.startswith('tdbrain') else 'não replica'),
                               fontsize=8)
            dn = (e[f'{b}|velhos'] - e[f'{b}|novos']).values[sel]
            dt = np.interp(f[sel], ft, tro['totalPsdYOlder'] - tro['totalPsdYounger'])
            ax[1, j].axhline(0, color='0.75', lw=0.6)
            ax[1, j].plot(f[sel], dn, color='#3C3C3C', lw=1.3, label='nosso: velhos − novos')
            ax[1, j].plot(f[sel], dt, color='#3C3C3C', lw=1.0, ls=(0, (4, 2)), label='Tröndle et al.: velhos − novos')
            ax[1, j].set_title(f'forma: r = {br(np.corrcoef(dn, dt)[0, 1])}', loc='left', fontsize=8)
            ax[1, j].set_xlabel('Hz')
        ax[0, 0].set_ylabel('log₁₀ potência')
        ax[1, 0].set_ylabel('Δ log₁₀ potência')
        h = [x for a in (ax[0, 0], ax[1, 0]) for x in zip(*a.get_legend_handles_labels())]
        fig.legend([x[0] for x in h], [x[1] for x in h], loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.02))
        fig.tight_layout(rect=(0, 0.05, 1, 1))
        salvar(fig, 'R1_espectro_idade_troendle2022')
    pe, pf = _ler(cfg, 'escada.csv'), _ler(cfg, 'efeitos.csv')
    if pe and pf and pe.stat().st_size > 5:
        e, ef = pd.read_csv(pe), pd.read_csv(pf)
        prot = list(dict.fromkeys(e.protocolo))
        blocos = [b for b in ('adhdata_tarefa', 'hbn_fechado', 'tdbrain_fechado', 'tdbrain_aberto') if b in set(e.bloco)]
        fig = plt.figure(figsize=(6.77, 4.9))
        gs = fig.add_gridspec(2, max(len(blocos), 1), height_ratios=[1, 0.6], hspace=0.38, wspace=0.1, left=0.215)
        for k, b in enumerate(blocos):
            a = fig.add_subplot(gs[0, k])
            d = e[e.bloco == b].pivot_table(index='protocolo', columns='nulo', values='auc_media')
            a.axvline(0.5, color='#9AA3AA', lw=0.7)
            a.axvline(0.83, color='k', lw=0.8, ls='--')
            for y, pr in enumerate(prot):
                real, nulo = d.loc[pr, False], d.loc[pr, True]
                a.plot([nulo, real], [y, y], color='#C9CFD4', lw=2.2)
                a.scatter(nulo, y, s=26, facecolor='white', edgecolor='#8A96A0')
                a.scatter(real, y, s=30, color='#137E83')
            a.set_xlim(0.3, 1.03)
            a.set_ylim(len(prot) - 0.5, -0.6)
            a.set_yticks(range(len(prot)), prot if k == 0 else [], fontsize=7.4)
            a.set_title(NOME[b], fontsize=7.8)
        t = ef[ef.pipeline == 'taghi'].groupby('nulo')[['V1', 'V2', 'V3', 'V4']].mean()
        ab, y = fig.add_subplot(gs[1, :]), np.arange(4)
        ab.barh(y - 0.19, t.loc[True].values, height=0.36, color='#C4452C', label='rótulo permutado: AUC criada sem sinal')
        ab.barh(y + 0.19, t.loc[False].values, height=0.36, color='#137E83', alpha=0.75, label='rótulo real')
        ab.axvline(0, color='k', lw=0.7)
        ab.set_yticks(y, ['divisão por época', 'seleção antes da divisão', 'classificador escolhido no teste',
                          'métrica por época'], fontsize=7.2)
        ab.invert_yaxis()
        ab.set_xlabel('mudança na AUC ao ligar o vazamento (média dos blocos)')
        ab.legend(loc='lower right', fontsize=6.8, frameon=False)
        salvar(fig, 'V1_vazamento_artigo')
