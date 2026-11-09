"""Ensaios do artigo: fenótipos, impressão digital, pareamento, complementaridade, integração e tabuleiro."""
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from . import estatistica as e
from .caracteristicas import tabela
from .config import BANCO, BLOCOS, BORDA8, C18, RESULTADOS, isolar_threads

ART = ['emg_t', 'razao_emg', 'log_var', 'borda_int', 'ocular']
FEATS = ['chi', 'chi_post', 'offset', 'iaf', 'alfa_pw', 'tbr']
VARS_PROC = ['rej_emg', 'rej_olho', 'ica_musc', 'ica_olho', 'internos', 'dur60', 'nativa', 'U30']
VARS_COV = ['base', 'rej_emg', 'rej_olho', 'ica_musc', 'ica_olho', 'internos', 'dur60']
MODELOS = ['C0', 'chi', 'per', 'tbr', 'esp', 'art', 'chi+esp', 'chi+per+esp']
APRESENTACOES = [('desatento', 'F90.0'), ('hiperativo', 'F90.1'), ('combinado', 'F90.2')]
N_MIN = 5


def suprimir(cont, chaves=('banco', 'variavel'), col='n'):
    out, suprimidas = [], []
    for k, g in cont.groupby(list(chaves), sort=False):
        pos = g[g[col] >= 1].sort_values(col, kind='stable')
        n = int(((pos[col] >= 1) & (pos[col] < N_MIN)).sum())
        if not n:
            out.append(g)
            continue
        while n < len(pos) and (pos[col].iloc[:n].sum() < N_MIN or n < 2):
            n += 1
        if n >= len(pos) or pos[col].iloc[:n].sum() < N_MIN:
            suprimidas.append(k)
            continue
        junta = pos.iloc[:n]
        out += [g.drop(junta.index), junta.iloc[:1].assign(nivel='outros: ' + ' + '.join(junta['nivel']), **{col: junta[col].sum()})]
    return pd.concat(out, ignore_index=True), suprimidas


def sem_pequenos(df, col='n'):
    return df[~((df[col] >= 1) & (df[col] < N_MIN))].reset_index(drop=True)


def subtipo(cods):
    for c in cods:
        for p, n in (('F90.0', 'desatento'), ('F90.1', 'hiperativo'), ('F90.2', 'combinado')):
            if c.startswith(p):
                return n
    return 'outro_nao_especificado' if any(c.startswith('F90') for c in cods) else 'sem_F90'


def fenotipos(cfg, destino=RESULTADOS):
    mapa, pac, dua, cont = pd.read_csv(cfg['mapa_ids_csv']), Path(cfg['pacote_dir']), Path(cfg['saida_dir']), []

    def contar(b, nome, s):
        cont.extend({'banco': b, 'variavel': nome, 'nivel': str(k), 'n': int(v)} for k, v in s.value_counts(dropna=False).items())

    def ids(banco, bloco):
        mp = mapa[mapa['banco'].str.startswith(banco)][['id_original', 'codigo']].drop_duplicates('codigo')
        return pd.read_csv(pac / bloco / 'metadados.csv').merge(mp, on='codigo', how='left')

    def eid(df, col, unico=True):
        df['EID'] = df[col].astype(str).str.split(',').str[0].str.strip().str.replace('sub-', '', regex=False)
        return df.drop_duplicates('EID') if unico else df
    m = eid(ids('hbn', 'hbn_fechado'), 'id_original', False)
    f = eid(pd.read_csv(cfg['hbn_fenotipo_csv'], low_memory=False), 'Identifiers')
    dxc = [c for c in f.columns if c.startswith('Diagnosis_ClinicianConsensus,DX_') and c.endswith('_Code')]
    dxcat = [c for c in f.columns if c.startswith('Diagnosis_ClinicianConsensus,DX_') and c.endswith('_Cat')]
    cods = f[dxc].apply(lambda r: [str(v) for v in r if pd.notna(v)], axis=1)
    f['subtipo'] = cods.map(subtipo)
    f['n_comorb'] = cods.map(lambda cs: sum(c.strip() not in ('', 'nan') and not c.startswith('F90')
                                            and not c.upper().startswith('NO') for c in cs))
    f['amp_pos'] = pd.to_numeric(f.get('DrugScreen,AMP500'), errors='coerce')
    cats, codes = f[dxcat].astype(str).agg(' | '.join, axis=1).str.lower(), f[dxc].astype(str).agg(' | '.join, axis=1)
    f['comorb_ansiedade'], f['comorb_tea'] = cats.str.contains('anxiety').astype(int), codes.str.contains('F84').astype(int)
    f['comorb_aprend'], f['comorb_humor'] = codes.str.contains('F81').astype(int), cats.str.contains('depress|bipolar|mood').astype(int)
    o = eid(pd.read_csv(cfg['hbn_sintomas_csv'], low_memory=False), 'Identifiers')
    o = o[['EID', 'SWAN,SWAN_Total', 'SDQ,SDQ_Hyperactivity', 'CGAS,CGAS_Score']].set_axis(['EID', 'swan', 'sdq_hiper', 'cgas'], axis=1)
    comorb = ['comorb_ansiedade', 'comorb_tea', 'comorb_aprend', 'comorb_humor']
    h = m.merge(f[['EID', 'subtipo', 'n_comorb', 'amp_pos'] + comorb], on='EID', how='left').merge(o, on='EID', how='left')
    h.drop(columns=['id_original', 'EID']).to_csv(dua / 'fenotipos_hbn.csv', index=False)
    for g in ('TDAH', 'controle'):
        s = h[h['grupo'] == g]
        for nome, v in [('subtipo', s['subtipo']), ('n_comorb', s['n_comorb'].clip(upper=3)), ('amp_pos', s['amp_pos']),
                        ('sexo', s['sexo'])] + [(c, s[c]) for c in comorb] + [('swan_presente', s['swan'].notna())]:
            contar(f'hbn_{g}', nome, v)
    idade = h.groupby(['grupo', 'subtipo'])['idade'].agg(['mean', 'std', 'median', 'count']).reset_index()
    idade[idade['count'] >= N_MIN].to_csv(destino / 'e0_idade_por_apresentacao.csv', index=False)
    t = ids('tdbrain', 'tdbrain_fechado')
    x = pd.read_excel(cfg['tdbrain_xlsx']).drop_duplicates('TDBRAIN_ID')
    cols = ['ADHD_pre_Att_leading', 'ADHD_pre_Hyp_leading', 'coffee', 'cigarette', 'drugs', 'sleep', 'education', 'Dataset',
            'indication']
    t = t.merge(x[['TDBRAIN_ID', 'sessTime'] + cols], left_on='id_original', right_on='TDBRAIN_ID', how='left')
    t['subestudo_nf'] = t['Dataset'].astype(str).str.contains('NF', case=False).astype(int)
    t['sesstime_ausente'] = t['sessTime'].isna().astype(int)
    t.drop(columns=['TDBRAIN_ID', 'id_original']).to_csv(dua / 'fenotipos_tdbrain.csv', index=False)
    for g in ('TDAH', 'controle'):
        s = t[(t['grupo'] == g) & (t['idade'] >= 18)]
        for c in cols[:6] + ['Dataset', 'subestudo_nf']:
            contar(f'tdbrain_{g}', c, s[c])
    tab, suprimidas = suprimir(pd.DataFrame(cont))
    tab.to_csv(destino / 'e0_contagens.csv', index=False)
    return suprimidas


def dados(cfg, bloco):
    banco = BANCO[bloco][0]
    d, covs = tabela(cfg, bloco)
    d = d.rename(columns={'base__chi': 'chi', 'base__alfa_pa': 'alfa', 'base__iaf_cg': 'iaf', 'aq__tbr_fc': 'tbr_fc',
                          **{f'aq__{a}': a for a in ART}})
    fen = Path(cfg['saida_dir']) / f'fenotipos_{banco}.csv'
    if fen.exists():
        d = d.merge(pd.read_csv(fen).drop(columns=['grupo', 'idade', 'sexo', 'duracao_s'], errors='ignore'), on='codigo', how='left')
        for c in ('swan', 'sdq_hiper', 'cgas', 'amp_pos', 'n_comorb'):
            if c in d:
                d[c] = pd.to_numeric(d[c], errors='coerce')
    C = covs['base']
    ok = (np.isfinite(d['chi']) & np.isfinite(C).all((1, 2))).values
    return d[ok].reset_index(drop=True), C[ok]


def familias(d, L):
    return {'C0': ([e.C0(d)], None), 'chi': ([d[['chi']].values], None), 'per': ([d[['alfa', 'iaf']].values], None),
            'tbr': ([d[['tbr_fc']].values], None), 'esp': ([], L), 'art': ([d[ART].values], None),
            'chi+esp': ([d[['chi']].values], L), 'chi+per+esp': ([d[['chi', 'alfa', 'iaf']].values], L)}


def rho_boot(a, b, n=300):
    rng = np.random.default_rng(0)
    bs = [stats.spearmanr(a[i], b[i])[0] for i in (rng.choice(len(a), len(a)) for _ in range(n))]
    return np.nanpercentile(bs, [2.5, 97.5])


def p1(D):
    linhas = []
    for bloco, (d, C) in D.items():
        ctrl = d['y'].values == 0
        V = e.vetorizar(e.logs_tangente(C), np.flatnonzero(ctrl))
        Vc = V - V[ctrl].mean(0)
        feats = {'chi': d['chi'].values, 'alfa': d['alfa'].values, 'iaf': d['iaf'].values, 'tbr': d['tbr_fc'].values,
                 'esp_dist': np.linalg.norm(V, axis=1), 'esp_pc1': Vc @ np.linalg.svd(Vc[ctrl], full_matrices=False)[2][0]}
        confs = {c: d[c].values for c in ['idade', 'sexo_f', 'log_dur'] + ART}
        for fn, fv in feats.items():
            for cn, cv in confs.items():
                ok = ctrl & np.isfinite(fv) & np.isfinite(cv)
                if ok.sum() < 20 or np.nanstd(cv[ok]) < 1e-9:
                    continue
                r, p = stats.spearmanr(fv[ok], cv[ok])
                lo, hi = rho_boot(fv[ok], cv[ok])
                linhas.append({'banco': BANCO[bloco][0], 'cond': BANCO[bloco][1], 'feature': fn, 'confundidor': cn, 'rho': r,
                               'p': p, 'ic_lo': lo, 'ic_hi': hi, 'n_controles': int(ok.sum())})
    return pd.DataFrame(linhas)


def p1_estado(D, banco):
    (dF, CF), (dA, CA) = D[f'{banco}_fechado'], D[f'{banco}_aberto']
    m = dF[['codigo']].reset_index().merge(dA[['codigo']].reset_index(), on='codigo', suffixes=('_f', '_a'))
    i, j = m['index_f'].values, m['index_a'].values
    out = []
    for fn in ('chi', 'alfa', 'iaf', 'tbr_fc'):
        a, b = dF[fn].values[i], dA[fn].values[j]
        ok = np.isfinite(a) & np.isfinite(b)
        out.append({'banco': banco, 'feature': fn, 'rho_OF_OA': stats.spearmanr(a[ok], b[ok])[0],
                    'delta_med_sobre_dp': float(np.median(a[ok] - b[ok]) / np.std(a[ok])), 'n': int(ok.sum())})
    Lf, La = e.logs_tangente(CF[i]), e.logs_tangente(CA[j])
    entre = np.linalg.norm(Lf - La[np.random.default_rng(0).permutation(len(j))], axis=(1, 2))
    out.append({'banco': banco, 'feature': 'esp_tangente', 'rho_OF_OA': np.nan, 'n': len(i),
                'delta_med_sobre_dp': float(np.median(np.linalg.norm(Lf - La, axis=(1, 2))) / np.median(entre)),
                'nota': 'razão distância dentro da pessoa (OF×OA) / entre pessoas'})
    return out


_D = {}


def _init(D):
    _D.update({k: (d, C, e.logs_tangente(C)) for k, (d, C) in D.items()})


def designs_hbn(d):
    tdah, ctrl, idade = np.flatnonzero(d['y'] == 1), np.flatnonzero(d['y'] == 0), d['idade'].values
    faixas = [(5, 8), (8, 11), (11, 22)]
    nmin = int(min(((idade[ctrl] >= lo) & (idade[ctrl] < hi)).sum() for lo, hi in faixas) * 0.9)
    out = [('global', tdah, ctrl, None)]
    out += [(f'idade_{lo}-{hi}', tdah[(idade[tdah] >= lo - 0.5) & (idade[tdah] < hi + 0.5)],
             ctrl[(idade[ctrl] >= lo) & (idade[ctrl] < hi)], nmin) for lo, hi in faixas]
    sx = d['sexo_f'].values
    nf = int((sx[ctrl] == 1).sum() * 0.9)
    out += [(n, tdah[sx[tdah] == s], ctrl[sx[ctrl] == s], nf) for s, n in ((0.0, 'so_meninos'), (1.0, 'so_meninas'))]
    if 'subtipo' in d:
        st = d['subtipo'].values
        out += [(f'subtipo_{s}', tdah[st[tdah] == s], ctrl, 75) for s in ('desatento', 'combinado', 'hiperativo')]
        out.append(('desatento_x_combinado', tdah[st[tdah] == 'desatento'], tdah[st[tdah] == 'combinado'], 300))
    if 'n_comorb' in d:
        nc = d['n_comorb'].values
        out += [('tdah_sem_comorbidade', tdah[nc[tdah] == 0], ctrl, 150), ('tdah_com_comorbidade', tdah[nc[tdah] >= 1], ctrl, 150)]
    return out


def _sorteio(args):
    nome, bloco, a, b, seed, n_max, perm = args
    d, _, L = _D[bloco]
    pares = e.parear(d, a, b, seed, n_max=n_max)
    if len(pares) < 15:
        return None
    idx = np.array([i for p in pares for i in p])
    dd, y = d.iloc[idx].reset_index(drop=True), np.tile([1, 0], len(pares))
    if perm:
        y = np.random.default_rng(10_000 + seed).permutation(y)
    parts, fam = e.particoes(y, rep=3, seed=seed), familias(dd, L[idx])
    return {'design': nome, 'banco': 'hbn', 'cond': BANCO[bloco][1], 'seed': seed, 'n_pares': len(pares), 'permutado': perm,
            'idade_media': float(dd['idade'].mean()), 'frac_f': float(np.nanmean(dd['sexo_f'])),
            **{m: e.auc_modelo(y, fam[m][0], parts, fam[m][1]) for m in MODELOS}}


def p2(D, jobs, sorteios=20):
    so_oa = ('global', 'subtipo_desatento', 'subtipo_combinado', 'desatento_x_combinado', 'so_meninos', 'so_meninas')
    tarefas = [(nome, bloco, a, b, s, nmax, perm) for bloco in ('hbn_fechado', 'hbn_aberto')
               for nome, a, b, nmax in designs_hbn(D[bloco][0]) if bloco == 'hbn_fechado' or nome in so_oa
               for perm, seeds in ((False, range(sorteios)), (True, range(100, 120))) for s in seeds]
    with Pool(jobs, initializer=_init, initargs=(D,)) as pool:
        res = pd.DataFrame([r for r in pool.imap_unordered(_sorteio, tarefas, chunksize=1) if r])
    linhas = []
    for (c, des), g in res.groupby(['cond', 'design']):
        real, p95 = g[~g['permutado']], np.nanpercentile(g.loc[g['permutado'], MODELOS].values, 95)
        linhas += [{'cond': c, 'design': des, 'modelo': m, 'n_pares': int(real['n_pares'].median()), 'auc': real[m].mean(),
                    'lo': np.percentile(real[m], 2.5), 'hi': np.percentile(real[m], 97.5), 'nulo_p95': p95,
                    'frac_acima': (real[m] > p95).mean()} for m in MODELOS]
    return res, pd.DataFrame(linhas)


FILTROS = {'lote_121s': lambda d: d['duracao_s'] >= 120.5, 'suporte_18_56': lambda d: d['idade'] <= 56}


def _p3(args):
    bloco, nome, filtro, canais, modelos, perm = args
    d, C, _ = _D[bloco]
    if filtro:
        m = FILTROS[filtro](d).values
        d, C = d[m].reset_index(drop=True), C[m]
    y = np.random.default_rng(perm).permutation(d['y'].values) if perm else d['y'].values
    L = e.logs_tangente(C, [C18.index(c) for c in canais] if canais else None)
    parts, c0, chi, art = e.particoes(y), e.C0(d), d[['chi']].values, d[ART].values
    combos = {'C0': ([c0], None), 'C0+chi': ([c0, chi], None), 'C0+esp': ([c0], L), 'C0+chi+esp': ([c0, chi], L),
              'C0+art': ([c0, art], None), 'C0+art+esp': ([c0, art], L), 'C0+art+chi': ([c0, art, chi], None),
              'esp': ([], L), 'art': ([art], None)}
    return {'banco': BANCO[bloco][0], 'cond': BANCO[bloco][1], 'analise': nome, 'n_tdah': int(d['y'].sum()),
            'n_ctrl': int((1 - d['y']).sum()), 'perm': perm, **{m: e.auc_modelo(y, combos[m][0], parts, combos[m][1]) for m in modelos}}


def p3(D, jobs, n_perm=30):
    todos = ['C0', 'C0+chi', 'C0+esp', 'C0+chi+esp', 'C0+art', 'C0+art+esp', 'C0+art+chi', 'esp', 'art']
    internos = [c for c in C18 if c not in BORDA8]
    t = []
    for bloco in D:
        t.append((bloco, 'completo', None, None, todos, 0))
        t += [(bloco, 'completo', None, None, ['C0+esp', 'esp', 'C0+chi', 'C0+chi+esp'], p) for p in range(1, n_perm + 1)]
        if not bloco.startswith('hbn') or bloco == 'hbn_fechado':
            t += [(bloco, 'so_internos', None, internos, ['esp', 'C0+esp'], 0),
                  (bloco, 'so_borda', None, list(BORDA8), ['esp', 'C0+esp'], 0)]
        if bloco.startswith('tdbrain'):
            t += [(bloco, f, f, None, todos, 0) for f in FILTROS]
    with Pool(jobs, initializer=_init, initargs=(D,)) as pool:
        return pd.DataFrame(pool.map(_p3, t, chunksize=1))


def p4(D):
    d, C = D['hbn_fechado']
    t = d['y'].values == 1
    d = d.assign(esp_dist=np.linalg.norm(e.vetorizar(e.logs_tangente(C), np.flatnonzero(~t)), axis=1))
    contrastes = {'amp_positivo_vs_negativo': (d['amp_pos'] == 1, d['amp_pos'] == 0),
                  'com_comorbidade_vs_sem': (d['n_comorb'] >= 1, d['n_comorb'] == 0),
                  'ansiedade_vs_sem': (d['comorb_ansiedade'] == 1, d['n_comorb'] == 0),
                  'tea_vs_sem': (d['comorb_tea'] == 1, d['n_comorb'] == 0),
                  'hiperativo_vs_desatento': (d['subtipo'] == 'hiperativo', d['subtipo'] == 'desatento'),
                  'combinado_vs_desatento': (d['subtipo'] == 'combinado', d['subtipo'] == 'desatento')}
    rng, g_linhas, s_linhas = np.random.default_rng(0), [], []
    for nome, (a, b) in contrastes.items():
        m = t & (a.values | b.values)
        dd, y = d[m].reset_index(drop=True), a.values[m].astype(int)
        Z = e.C0(dd)
        for fn in ('chi', 'alfa', 'iaf', 'tbr_fc', 'esp_dist', 'emg_t'):
            v = dd[fn].values.astype(float)
            gb = [e.hedges(y[bi], v[bi], Z[bi]) for bi in e.reamostras(y, 500, rng)]
            g_linhas.append({'contraste': nome, 'feature': fn, 'n_a': int(y.sum()), 'n_b': int((1 - y).sum()),
                             'g': e.hedges(y, v, Z), 'ic_lo': np.nanpercentile(gb, 2.5), 'ic_hi': np.nanpercentile(gb, 97.5)})
    Z = e.C0(d)
    for esc in ('swan', 'sdq_hiper'):
        for fn in ('chi', 'alfa', 'iaf', 'tbr_fc', 'esp_dist'):
            ok = (np.isfinite(d[fn]) & np.isfinite(d[esc])).values
            X = np.column_stack([np.ones(ok.sum()), Z[ok]])
            r1, r2 = (d[c].values[ok] - X @ np.linalg.lstsq(X, d[c].values[ok], rcond=None)[0] for c in (fn, esc))
            r, p = stats.spearmanr(r1, r2)
            s_linhas.append({'escala': esc, 'feature': fn, 'rho_parcial': r, 'p': p, 'n': int(ok.sum())})
    return pd.DataFrame(g_linhas), pd.DataFrame(s_linhas)


def _fit(Xtr, ytr, Xte, C):
    from sklearn.linear_model import LogisticRegression
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    clf = LogisticRegression(C=C, max_iter=2000).fit((Xtr - mu) / sd, ytr)
    return clf.decision_function((Xte - mu) / sd), clf.coef_[0]


def _resid(V, Z, tr, y):
    if Z.shape[1] == 0:
        return V
    ctrl = tr[y[tr] == 0]
    B = np.linalg.lstsq(np.column_stack([np.ones(len(ctrl)), Z[ctrl]]), V[ctrl], rcond=None)[0]
    return V - np.column_stack([np.ones(len(V)), Z]) @ B


def integrar(d, L, y, rep=5, seed=0):
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    Z, parts = e.C0(d), e.particoes(y, rep=rep, seed=seed)
    dec = {k: np.full((rep, len(y)), np.nan) for k in ('C0', 'concat', 'stack', 'chi', 'per', 'esp')}
    pesos = []
    for k, (tr, te) in enumerate(parts):
        r = k // 5
        if Z.shape[1]:
            dec['C0'][r, te] = _fit(Z[tr], y[tr], Z[te], 1.0)[0]
        X = np.column_stack([Z, d[['chi', 'alfa', 'iaf']].values, e.vetorizar(L, tr)])
        X = np.where(np.isfinite(X), X, np.nanmedian(X[tr], axis=0))
        dec['concat'][r, te] = _fit(X[tr], y[tr], X[te], 0.05)[0]
        s_tr, s_te = {}, {}
        for nome, V in (('chi', d[['chi']].values), ('per', d[['alfa', 'iaf']].values), ('esp', e.vetorizar(L, tr))):
            R = _resid(np.where(np.isfinite(V), V, np.nanmedian(V[tr], axis=0)), Z, tr, y)
            c, s_tr[nome] = (0.05 if R.shape[1] > 20 else 1.0), np.zeros(len(tr))
            for itr, ite in StratifiedKFold(5, shuffle=True, random_state=0).split(tr, y[tr]):
                s_tr[nome][ite] = _fit(R[tr[itr]], y[tr[itr]], R[tr[ite]], c)[0]
            s_te[nome] = dec[nome][r, te] = _fit(R[tr], y[tr], R[te], c)[0]
        Mtr, Mte = (np.column_stack([Z[i]] + [s[f] for f in ('chi', 'per', 'esp')]) for i, s in ((tr, s_tr), (te, s_te)))
        dec['stack'][r, te], coef = _fit(Mtr, y[tr], Mte, 1.0)
        pesos.append(coef[-3:])
    res = {k: float(np.mean([roc_auc_score(y, dd) for dd in v])) for k, v in dec.items() if np.isfinite(v).all()}
    return {**res, **dict(zip(('peso_chi', 'peso_per', 'peso_esp'), np.mean(pesos, 0)))}


def _e4(args):
    bloco, seed, pareado = args
    d, _, L = _D[bloco]
    if pareado:
        pares = e.parear(d, np.flatnonzero(d['y'] == 1), np.flatnonzero(d['y'] == 0), seed)
        idx = np.array([i for p in pares for i in p])
        d, L, y, rep = d.iloc[idx].reset_index(drop=True), L[idx], np.tile([1, 0], len(pares)), 3
    else:
        y, rep = d['y'].values, 5
    return {'banco': BANCO[bloco][0], 'cond': BANCO[bloco][1], 'pareado': pareado, 'seed': seed, 'n': len(y),
            **integrar(d, L, y, rep, seed)}


def e4(D, jobs, sorteios=20):
    t = [(b, 0, False) for b in ('hbn_fechado', 'hbn_aberto', 'tdbrain_fechado', 'adhdata_tarefa')]
    with Pool(jobs, initializer=_init, initargs=(D,)) as pool:
        return pd.DataFrame(pool.map(_e4, t + [('hbn_fechado', s, True) for s in range(sorteios)], chunksize=1))


def dados_tabuleiro(cfg, bloco):
    d, covs = tabela(cfg, bloco)
    d = d.rename(columns={'aq__log_var': 'log_var'})
    covs = {k: C for k, C in covs.items() if np.isfinite(C).any()}
    valido = np.ones(len(d), bool)
    for C in covs.values():
        tr = np.trace(C, axis1=1, axis2=2)
        valido &= np.isfinite(C).all((1, 2)) & np.isfinite(tr) & (np.abs(tr) > 1e-10)
    return d, covs, valido


def auc_1d(y, v, Z):
    return e.auc_modelo(y, [Z, np.where(np.isfinite(v), v, np.nanmedian(v))[:, None]], e.particoes(y))


def distancia(C, y, valido, loo=False):
    iv, out = np.flatnonzero(valido), np.full(len(y), np.nan)
    ctrl, L = y[iv] == 0, e.logs_tangente(C[iv])
    out[iv] = np.linalg.norm(e.vetorizar(L, np.flatnonzero(ctrl)), axis=1)
    if loo:
        iu = np.triu_indices(L.shape[-1])
        X = L[:, iu[0], iu[1]] * np.where(iu[0] == iu[1], 1.0, np.sqrt(2))
        nc, soma = ctrl.sum(), X[ctrl].sum(0)
        for i in np.flatnonzero(ctrl):
            out[iv[i]] = np.linalg.norm(X[i] - (soma - X[i]) / (nc - 1))
    return out


def tabuleiro(cfg, bloco):
    d, covs, valido = dados_tabuleiro(cfg, bloco)
    y, Z, ctrl = d['y'].values, e.C0(d), d['y'].values == 0
    esp = distancia(covs['base'], y, valido)
    conf = {c: d[k].values for c, k in (('idade', 'idade'), ('sexo_f', 'sexo_f'), ('log_dur', 'log_dur'), ('emg_t', 'base__emg_t'),
                                         ('ocular', 'base__ocular'), ('log_var', 'log_var'))}
    longo, linhas = [], []
    for fn, fv in {**{f: d[f'base__{f}'].values for f in FEATS}, 'espacial': esp}.items():
        for cn, cv in conf.items():
            ok = ctrl & np.isfinite(fv) & np.isfinite(cv)
            if ok.sum() >= 20 and np.nanstd(cv[ok]) >= 1e-9:
                longo.append({'bloco': bloco, 'feature': fn, 'coluna': cn, 'tipo': 'rho_controles',
                              'valor': stats.spearmanr(fv[ok], cv[ok])[0]})
    c0 = e.auc_modelo(y, [Z], e.particoes(y)) if Z.shape[1] else np.nan
    for var in ['base'] + VARS_PROC:
        for f in FEATS:
            v = d[f'{var}__{f}'].values.astype(float)
            if not np.isfinite(v).any():
                continue
            g, lo, hi = e.g_boot(y, v, Z)
            linhas.append({'bloco': bloco, 'variante': var, 'feature': f, 'g': g, 'g_lo': lo, 'g_hi': hi,
                           'auc_c0_mais_f': auc_1d(y, v, Z), 'auc_c0': c0, 'falha': float(np.isnan(v).mean())})
        for f in ('emg_t', 'ocular'):
            linhas.append({'bloco': bloco, 'variante': var, 'feature': f'artefato_{f}', 'auc_c0': c0,
                           'g': np.nanmedian(d[f'{var}__{f}'].values - d[f'base__{f}'].values)})
    for var in VARS_COV:
        if var in covs and valido.sum() >= 30:
            L, yy = e.logs_tangente(covs[var][valido]), y[valido]
            parts = e.particoes(yy)
            linhas.append({'bloco': bloco, 'variante': var, 'feature': 'espacial', 'auc_c0': c0, 'falha': float((~valido).mean()),
                           'auc_c0_mais_f': e.auc_modelo(yy, [Z[valido]], parts, L), 'auc_so_f': e.auc_modelo(yy, [], parts, L)})
    return longo, linhas


def tabuleiro_todos(cfg, blocos=BLOCOS, destino=RESULTADOS):
    with Pool(len(blocos)) as pool:
        out = pool.starmap(tabuleiro, [(cfg, b) for b in blocos])
    g = pd.DataFrame([l for _, b in out for l in b])
    g.to_csv(destino / 'e6_g_auc_variantes.csv', index=False)
    base, rows = g[g['variante'] == 'base'].set_index(['bloco', 'feature']), []
    for _, r in g[(g['variante'] != 'base') & ~g['feature'].str.startswith('artefato')].iterrows():
        b = base.loc[(r['bloco'], r['feature'])]
        esp = r['feature'] == 'espacial'
        rows.append({'bloco': r['bloco'], 'feature': r['feature'], 'coluna': r['variante'], 'tipo': 'delta_auc' if esp else 'delta_g',
                     'valor': r['auc_so_f'] - b['auc_so_f'] if esp else r['g'] - b['g']})
    for (bl, f), b in base[~base.index.get_level_values(1).str.startswith('artefato')].iterrows():
        rows.append({'bloco': bl, 'feature': f, 'coluna': 'sinal', 'tipo': 'auc', 'valor': b['auc_so_f']} if f == 'espacial' else
                    {'bloco': bl, 'feature': f, 'coluna': 'sinal', 'tipo': 'g', 'valor': b['g'], 'exclui0': bool(b['g_lo'] > 0 or b['g_hi'] < 0)})
    pd.concat([pd.DataFrame([l for a, _ in out for l in a]), pd.DataFrame(rows)], ignore_index=True).to_csv(
        destino / 'e6_tabuleiro_longo.csv', index=False)


def subtipos(cfg, destino=RESULTADOS):
    fen = pd.read_csv(Path(cfg['saida_dir']) / 'fenotipos_hbn.csv')[['codigo', 'subtipo']]
    longo, resumo = [], []
    for bloco in ('hbn_fechado', 'hbn_aberto'):
        d, covs, valido = dados_tabuleiro(cfg, bloco)
        d = d.merge(fen, on='codigo', how='left')
        for st, cid in APRESENTACOES:
            sel = ((d['grupo'] == 'controle') | (d['subtipo'] == st)).values
            ds = d[sel].reset_index(drop=True)
            y, Z = ds['y'].values, e.C0(ds)
            if y.sum() < N_MIN:
                continue
            base = {'bloco': bloco, 'cid': cid, 'apresentacao': st, 'n_tdah': int(y.sum()), 'n_controle': int((y == 0).sum())}
            for f in FEATS:
                g, lo, hi = e.g_boot(y, ds[f'base__{f}'].values.astype(float), Z)
                longo.append({'bloco': bloco, 'feature': f, 'coluna': f'sinal_{cid}', 'tipo': 'g', 'valor': g,
                              'exclui0': bool(lo > 0 or hi < 0)})
                resumo.append({**base, 'feature': f, 'g': g, 'g_lo': lo, 'g_hi': hi, 'auc_espacial': np.nan})
            ok = valido[sel]
            yy = y[ok]
            a = e.auc_modelo(yy, [], e.particoes(yy), e.logs_tangente(covs['base'][sel][ok]))
            longo.append({'bloco': bloco, 'feature': 'espacial', 'coluna': f'sinal_{cid}', 'tipo': 'auc', 'valor': a})
            resumo.append({**base, 'n_tdah': int(yy.sum()), 'n_controle': int((yy == 0).sum()), 'feature': 'espacial',
                           'g': np.nan, 'g_lo': np.nan, 'g_hi': np.nan, 'auc_espacial': a})
    pd.DataFrame(longo).to_csv(destino / 'e8_subtipos_longo.csv', index=False)
    pd.DataFrame(resumo).to_csv(destino / 'e8_subtipos_resumo.csv', index=False)


def amostra_potencia(cfg, bloco):
    d, covs, valido = dados_tabuleiro(cfg, bloco)
    d['esp_dist'] = distancia(covs['base'], d['y'].values, valido)
    d = d[(np.isfinite(d['log_var']) & np.isfinite(d['base__emg_t']) & np.isfinite(d['base__ocular'])).values]
    d = d.reset_index(drop=True)
    Z0, pot = e.C0(d), d['log_var'].values
    return d, {'C0': Z0, 'C0+pot': np.column_stack([Z0, pot]),
               'C0+pot+art': np.column_stack([Z0, pot, d['base__emg_t'].values, d['base__ocular'].values])}


def ajuste_potencia(cfg, blocos=BLOCOS, destino=RESULTADOS):
    linhas = []
    for bloco in blocos:
        d, Zs = amostra_potencia(cfg, bloco)
        y, Z0, pot = d['y'].values, Zs['C0'], d['log_var'].values
        n = {'n_tdah': int(y.sum()), 'n_controle': int((y == 0).sum())}
        g, lo, hi = e.g_boot(y, pot, Z0)
        linhas.append({'bloco': bloco, 'feature': 'potencia_total', 'ajuste': 'C0', 'g': g, 'g_lo': lo, 'g_hi': hi,
                       'auc_c0_mais_f': auc_1d(y, pot, Z0), 'auc_c0': e.auc_modelo(y, [Z0], e.particoes(y)), **n})
        for f in FEATS:
            for nome, Z in Zs.items():
                g, lo, hi = e.g_boot(y, d[f'base__{f}'].values.astype(float), Z)
                linhas.append({'bloco': bloco, 'feature': f, 'ajuste': nome, 'g': g, 'g_lo': lo, 'g_hi': hi,
                               'auc_c0_mais_f': np.nan, 'auc_c0': np.nan, **n})
    pd.DataFrame(linhas).to_csv(destino / 'e9_ajuste_potencia.csv', index=False)


def _g_dois(y, v, Zs):
    def ajuste(Z):
        X = np.column_stack([np.ones(len(v)), y.astype(float), Z])
        b = np.linalg.lstsq(X, v, rcond=None)[0]
        r = v - X @ b
        return b[1], np.sqrt(r @ r / (len(v) - X.shape[1]))
    k, s0 = 1 - 3 / (4 * len(v) - 9), ajuste(Zs['C0'])[1]
    return {n: (b1 / s * k, b1 / s0 * k) for n, (b1, s) in ((n, ajuste(Z)) for n, Z in Zs.items())}


def denominador_fixo(cfg, blocos=BLOCOS, destino=RESULTADOS):
    linhas = []
    for bloco in blocos:
        d, Zs = amostra_potencia(cfg, bloco)
        for f in ('chi', 'offset', 'alfa_pw', 'tbr', 'chi_post', 'esp_dist'):
            v = d[f if f == 'esp_dist' else f'base__{f}'].values.astype(float)
            m = np.isfinite(v)
            y, v, Zm = d['y'].values[m], v[m], {k: Z[m] for k, Z in Zs.items()}
            pt = _g_dois(y, v, Zm)
            B = [_g_dois(y[b], v[b], {k: Z[b] for k, Z in Zm.items()}) for b in e.reamostras(y, 300, 0)]
            for k in Zs:
                A = np.array([r[k] for r in B])
                (la, ha), (lb, hb) = np.nanpercentile(A[:, 0], [2.5, 97.5]), np.nanpercentile(A[:, 1], [2.5, 97.5])
                linhas.append({'bloco': bloco, 'feature': f, 'ajuste': k, 'g_a': pt[k][0], 'g_a_lo': la, 'g_a_hi': ha,
                               'g_b': pt[k][1], 'g_b_lo': lb, 'g_b_hi': hb, 'dif_a_menos_b': pt[k][0] - pt[k][1],
                               'n_tdah': int(y.sum()), 'n_controle': int((y == 0).sum())})
    pd.DataFrame(linhas).to_csv(destino / 'e14_denominador_fixo.csv', index=False)


def fluxo(cfg, blocos=BLOCOS, destino=RESULTADOS):
    L = []

    def add(bloco, etapa, grupo, n):
        L.append({'bloco': bloco, 'etapa': etapa, 'grupo': grupo, 'n': int(n)})
    for bloco in blocos:
        meta = pd.read_csv(Path(cfg['pacote_dir']) / bloco / 'metadados.csv')
        add(bloco, '1_metadados_todos', 'todos', len(meta))
        for g, c in meta['grupo'].value_counts().items():
            add(bloco, '1_metadados_por_grupo', str(g), c)
        with np.load(Path(cfg['saida_dir']) / 'caracteristicas' / f'{bloco}.npz') as z:
            extraidos = set(z['codigo'])
        for g in ('TDAH', 'controle'):
            mg = meta[meta['grupo'] == g]
            for etapa in ('2_extraidos_e1', '2_extraidos_e5'):
                add(bloco, etapa, g, len(set(mg['codigo']) & extraidos))
            if bloco.startswith('tdbrain'):
                add(bloco, '2b_idade_ausente', g, mg['idade'].isna().sum())
        d, covs, valido = dados_tabuleiro(cfg, bloco)
        okart = (np.isfinite(d['log_var']) & np.isfinite(d['base__emg_t']) & np.isfinite(d['base__ocular'])).values
        for g, yv in (('TDAH', 1), ('controle', 0)):
            s = d['y'].values == yv
            add(bloco, '3_apos_filtro_idade_tdbrain' if bloco.startswith('tdbrain') else '3_analisados_e6', g, s.sum())
            add(bloco, '4_cov_valida_todas_variantes', g, (s & valido).sum())
            add(bloco, '5_amostra_e9_e14', g, (s & okart).sum())
            add(bloco, '6_alfa_pw_zero_base', g, (s & (d['base__alfa_pw'].values == 0)).sum())
            for f in ('iaf', 'alfa_pw', 'chi', 'chi_post', 'offset', 'tbr'):
                add(bloco, f'6_{f}_ausente_base', g, (s & ~np.isfinite(d[f'base__{f}'].values)).sum())
        for v, C in covs.items():
            tr = np.trace(C, axis1=1, axis2=2)
            add(bloco, f'7_cov_invalida_{v}', 'TDAH+controle',
                (~(np.isfinite(C).all((1, 2)) & np.isfinite(tr) & (np.abs(tr) > 1e-10))).sum())
        add(bloco, '7_menos_de_15_epocas_e1', 'TDAH+controle', (d['base__n_ep'] < 15).sum())
    sem_pequenos(pd.DataFrame(L)).to_csv(destino / 'e16_fluxo.csv', index=False)


def lote_tdbrain(cfg, destino=RESULTADOS):
    fen = pd.read_csv(Path(cfg['saida_dir']) / 'fenotipos_tdbrain.csv')[['codigo', 'subestudo_nf', 'sesstime_ausente']]
    L = []
    for bloco in ('tdbrain_fechado', 'tdbrain_aberto'):
        m = pd.read_csv(Path(cfg['pacote_dir']) / bloco / 'metadados.csv')
        d = m[m['grupo'].isin(['TDAH', 'controle'])].merge(fen, on='codigo', how='left')
        d['lote_120s'] = (d['duracao_s'].astype(float) < 120.5).astype(int)
        for faixa, dd in (('todas_idades', d), ('adultos_18+', d[d['idade'] >= 18])):
            for g, x in dd.groupby('grupo'):
                lote, nf, curto = x['lote_120s'], x['subestudo_nf'], x['duracao_s'].astype(float) < 120.5
                med = [('n', len(x)), ('sem_linha_no_lote', nf.isna().sum()), ('lote_120s', (lote == 1).sum()),
                       ('lote_121s', (lote == 0).sum()), ('subestudo_adhd_nf', (nf == 1).sum()),
                       ('sesstime_ausente', (x['sesstime_ausente'] == 1).sum())]
                for s, rot in ((1, 'dentro'), (0, 'fora')):
                    med += [(f'lote_120s_{rot}_adhd_nf', ((nf == s) & (lote == 1)).sum()), (f'n_{rot}_adhd_nf', (nf == s).sum())]
                med += [(f'lote{a}_e_dur{b}', ((lote == int(a == 120)) & (curto == (b == '<120,5'))).sum())
                        for a in (120, 121) for b in ('<120,5', '>=120,5')]
                L += [{'bloco': bloco, 'faixa': faixa, 'grupo': g, 'medida': k, 'n': int(v)} for k, v in med]
    sem_pequenos(pd.DataFrame(L)).to_csv(destino / 'e17_lote_tdbrain.csv', index=False)


RHO = ['idade', 'sexo_f', 'log_dur', 'emg_t', 'ocular', 'log_var']
VARS_TAB = ['rej_emg', 'rej_olho', 'ica_olho', 'ica_musc', 'internos', 'dur60', 'nativa', 'U30']
NAO_SE_APLICA = {('tbr', 'U30'), ('espacial', 'dur60'), ('espacial', 'U30'), ('espacial', 'nativa')}


def tabela_numerica(res=RESULTADOS):
    lg, gv, linhas = pd.read_csv(Path(res) / 'e6_tabuleiro_longo.csv'), pd.read_csv(Path(res) / 'e6_g_auc_variantes.csv'), []
    ordem = ['adhdata_tarefa', 'hbn_fechado', 'tdbrain_fechado', 'hbn_aberto', 'tdbrain_aberto']
    for bl in [b for b in ordem if b in set(lg['bloco'])] + [b for b in pd.unique(lg['bloco']) if b not in ordem]:
        for f in ('chi', 'chi_post', 'offset', 'alfa_pw', 'iaf', 'tbr', 'espacial'):
            s = lg[(lg['bloco'] == bl) & (lg['feature'] == f)]
            r = {'bloco': bl, 'estado': BANCO.get(bl, (bl, bl))[1], 'feature': f, 'metrica_proc': 'delta_auc' if f == 'espacial' else 'delta_g'}
            for c in RHO:
                v = s.loc[(s['tipo'] == 'rho_controles') & (s['coluna'] == c), 'valor']
                r[f'rho_{c}'] = v.iloc[0] if len(v) else np.nan
            for var in VARS_TAB:
                v = s.loc[s['tipo'].isin(['delta_g', 'delta_auc']) & (s['coluna'] == var), 'valor']
                r[f'd_{var}'] = np.nan if (f, var) in NAO_SE_APLICA or not len(v) else v.iloc[0]
            b = gv[(gv['bloco'] == bl) & (gv['feature'] == f) & (gv['variante'] == 'base')]
            g = b[['g', 'g_lo', 'g_hi']].iloc[0].tolist() if len(b) and f != 'espacial' else [np.nan] * 3
            r.update(zip(('g', 'g_lo', 'g_hi'), g), auc_base=b['auc_so_f'].iloc[0] if len(b) and f == 'espacial' else np.nan)
            linhas.append(r)
    return pd.DataFrame(linhas)


def escore(t):
    grupos = {'amostra': ['rho_idade', 'rho_sexo_f', 'rho_log_dur'], 'artefato_aquisicao': ['rho_emg_t', 'rho_ocular', 'rho_log_var'],
              'processamento': [f'd_{v}' for v in VARS_TAB]}
    L = []
    padrao = {'OF': ['adhdata_tarefa', 'hbn_fechado', 'tdbrain_fechado'], 'OA': ['hbn_aberto', 'tdbrain_aberto']}
    outros = {b: [b] for b in t['bloco'].unique() if b not in sum(padrao.values(), [])}
    for est, blocos in [(k, v) for k, v in padrao.items() if t['bloco'].isin(v).any()] + list(outros.items()):
        s = t[t['bloco'].isin(blocos)]
        for f in t['feature'].unique():
            tot = np.zeros(3, int)
            for nome, cols in grupos.items():
                lim = (0.2, 0.4) if nome != 'processamento' else (0.03, 0.07) if f == 'espacial' else (0.10, 0.20)
                v = np.abs(s.loc[s['feature'] == f, cols].values.ravel().astype(float))
                v = v[np.isfinite(v)]
                c = np.array([len(v), ((v >= lim[0]) & (v < lim[1])).sum(), (v >= lim[1]).sum()])
                tot += c
                L.append({'estado': est, 'feature': f, 'grupo_linhas': nome, 'casas': c[0], 'amarelo': c[1], 'vermelho': c[2],
                          'escore': (c[1] + 2 * c[2]) / (2 * c[0]) if c[0] else np.nan})
            L.append({'estado': est, 'feature': f, 'grupo_linhas': 'total', 'casas': tot[0], 'amarelo': tot[1],
                      'vermelho': tot[2], 'escore': (tot[1] + 2 * tot[2]) / (2 * tot[0])})
    return pd.DataFrame(L)


def tabela_e_escore(res=RESULTADOS, destino=RESULTADOS):
    t = tabela_numerica(res)
    t.round(2).to_csv(Path(destino) / 'e18_tabuleiro_numerico.csv', index=False)
    escore(t).round(2).to_csv(Path(destino) / 'e20_escore_cuidado.csv', index=False)


def espacial_loo(cfg, blocos=BLOCOS, destino=RESULTADOS):
    L = []
    for bloco in blocos:
        d, covs, valido = dados_tabuleiro(cfg, bloco)
        y = d['y'].values
        antes, depois = distancia(covs['base'], y, valido), distancia(covs['base'], y, valido, loo=True)
        ctrl = (y == 0) & valido
        for c, col in (('idade', 'idade'), ('sexo_f', 'sexo_f'), ('log_dur', 'log_dur'), ('emg_t', 'base__emg_t'),
                       ('ocular', 'base__ocular'), ('log_var', 'log_var')):
            cv = d[col].values.astype(float)
            ok = ctrl & np.isfinite(cv)
            if ok.sum() >= 20 and np.nanstd(cv[ok]) >= 1e-9:
                ra, rd = stats.spearmanr(antes[ok], cv[ok])[0], stats.spearmanr(depois[ok], cv[ok])[0]
                L.append({'bloco': bloco, 'coluna': c, 'rho_antes': ra, 'rho_loo': rd, 'mudanca': rd - ra,
                          'n_controles': int(ok.sum())})
    pd.DataFrame(L).to_csv(destino / 'e21_espacial_loo.csv', index=False)


def desequilibrio(cfg, blocos=BLOCOS, destino=RESULTADOS):
    linhas = []
    for bloco in blocos:
        d, _, _ = dados_tabuleiro(cfg, bloco)
        y, Z = d['y'].values, e.C0(d)
        for k, col in (('idade', 'idade'), ('sexo_f', 'sexo_f'), ('log_dur', 'log_dur'), ('emg_t', 'base__emg_t'),
                       ('ocular', 'base__ocular'), ('log_var', 'log_var')):
            v = d[col].values.astype(float)
            g = e.g_boot(y, v, Z) if k in ('emg_t', 'ocular', 'log_var') else (np.nan,) * 3
            linhas.append({'bloco': bloco, 'variavel': k, 'smd': e.smd(y, v), 'g_c0': g[0], 'g_lo': g[1], 'g_hi': g[2]})
    pd.DataFrame(linhas).to_csv(destino / 'e10_desequilibrio.csv', index=False)


PARTES = ('fenotipos', 'e2', 'e4', 'tabuleiro', 'subtipos', 'potencia', 'desequilibrio', 'denominador', 'fluxo', 'lote',
          'loo', 'tabela')


def rodar(cfg, partes=PARTES, jobs=8, sorteios=20, destino=RESULTADOS, blocos=BLOCOS):
    isolar_threads()
    destino = Path(destino)
    destino.mkdir(parents=True, exist_ok=True)
    if 'fenotipos' in partes:
        fenotipos(cfg, destino)
    if {'e2', 'e4'} & set(partes):
        D = {b: dados(cfg, b) for b in BLOCOS}
        (destino / 'e2_amostras.json').write_text(json.dumps({
            '_'.join(BANCO[b]): {'n_tdah': int(d['y'].sum()), 'n_ctrl': int((1 - d['y']).sum())} for b, (d, _) in D.items()},
            indent=1), encoding='utf-8')
    if 'e2' in partes:
        p1(D).to_csv(destino / 'e2_p1_impressao_digital.csv', index=False)
        pd.DataFrame(p1_estado(D, 'hbn') + p1_estado(D, 'tdbrain')).to_csv(destino / 'e2_p1_estado.csv', index=False)
        g, s = p4(D)
        g.to_csv(destino / 'e2_p4_dentro_tdah.csv', index=False)
        s.to_csv(destino / 'e2_p4_sintomas.csv', index=False)
        p3(D, jobs).to_csv(destino / 'e2_p3_complementaridade.csv', index=False)
        pares, resumo = p2(D, jobs, sorteios)
        pares.to_csv(destino / 'e2_p2_pareados.csv', index=False)
        resumo.to_csv(destino / 'e2_p2_resumo.csv', index=False)
    if 'e4' in partes:
        e4(D, jobs, sorteios).to_csv(destino / 'e4_integracao.csv', index=False)
    if 'tabuleiro' in partes:
        tabuleiro_todos(cfg, blocos, destino)
    if 'subtipos' in partes:
        subtipos(cfg, destino)
    if 'potencia' in partes:
        ajuste_potencia(cfg, blocos, destino)
    if 'desequilibrio' in partes:
        desequilibrio(cfg, blocos, destino)
    if 'denominador' in partes:
        denominador_fixo(cfg, blocos, destino)
    if 'fluxo' in partes:
        fluxo(cfg, blocos, destino)
    if 'lote' in partes:
        lote_tdbrain(cfg, destino)
    if 'loo' in partes:
        espacial_loo(cfg, blocos, destino)
    if 'tabela' in partes:
        tabela_e_escore(destino, destino)
