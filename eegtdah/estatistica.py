"""Covariáveis, efeito ajustado, validação por participante, espaço tangente e pareamento."""
import numpy as np

REG = 0.1


def C0(d):
    cols = []
    if np.isfinite(d['idade']).all():
        a = ((d['idade'] - d['idade'].mean()) / d['idade'].std()).values
        s = np.nan_to_num(d['sexo_f'].values, nan=0.5)
        cols += [a, a ** 2, s, a * s]
    du = d['log_dur'].values
    if du.std() > 1e-9:
        cols.append((du - du.mean()) / du.std())
    return np.column_stack(cols) if cols else np.zeros((len(d), 0))


def hedges(y, v, Z):
    ok = np.isfinite(v)
    y, v, Z = y[ok], v[ok], Z[ok]
    X = np.column_stack([np.ones(len(v)), y.astype(float), Z])
    b, *_ = np.linalg.lstsq(X, v, rcond=None)
    res = v - X @ b
    sd = np.sqrt(res @ res / (len(v) - X.shape[1]))
    return float(b[1] / sd * (1 - 3 / (4 * len(v) - 9)))


def reamostras(y, B, seed=0):
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    return [np.concatenate([rng.choice(np.flatnonzero(y == k), (y == k).sum()) for k in (0, 1)]) for _ in range(B)]


def g_boot(y, v, Z, B=300, seed=0, n_min=20):
    ok = np.isfinite(v)
    if ok.sum() < n_min:
        return np.nan, np.nan, np.nan
    y, v, Z = y[ok], v[ok], Z[ok]
    gb = [hedges(y[b], v[b], Z[b]) for b in reamostras(y, B, seed)]
    return (hedges(y, v, Z), *np.nanpercentile(gb, [2.5, 97.5]))


def smd(y, v):
    ok = np.isfinite(v)
    a, b = v[ok & (y == 1)], v[ok & (y == 0)]
    if len(a) < 5 or len(b) < 5:
        return np.nan
    s = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / s if s > 0 else np.nan


def particoes(y, rep=5, seed=0, k=5):
    from sklearn.model_selection import RepeatedStratifiedKFold
    return list(RepeatedStratifiedKFold(n_splits=k, n_repeats=rep, random_state=seed).split(np.zeros(len(y)), y))


def regularizar(C):
    n = C.shape[-1]
    return (1 - REG) * C / np.trace(C, axis1=-2, axis2=-1)[..., None, None] + REG * np.eye(n) / n


def logm_spd(C):
    w, V = np.linalg.eigh(C)
    return np.einsum('...ij,...j,...kj->...ik', V, np.log(np.clip(w, 1e-12, None)), V)


def logs_tangente(C, idx=None):
    if idx is not None:
        C = C[:, idx][:, :, idx]
    return logm_spd(regularizar(C))


def vetorizar(L, tr):
    L = L - L[tr].mean(0)
    iu = np.triu_indices(L.shape[-1])
    return L[:, iu[0], iu[1]] * np.where(iu[0] == iu[1], 1.0, np.sqrt(2))


def auc_modelo(y, blocos, parts, L=None, C_lr=None, k=5):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    dec = np.full((len(parts) // k, len(y)), np.nan)
    for i, (tr, te) in enumerate(parts):
        X = [b for b in blocos if b.shape[1]] + ([vetorizar(L, tr)] if L is not None else [])
        X = np.column_stack(X)
        X = np.where(np.isfinite(X), X, np.nanmedian(X[tr], axis=0))
        X = (X - X[tr].mean(0)) / (X[tr].std(0) + 1e-9)
        c = C_lr if C_lr is not None else (0.05 if X.shape[1] > 20 else 1.0)
        dec[i // k, te] = LogisticRegression(C=c, max_iter=2000).fit(X[tr], y[tr]).decision_function(X[te])
    return float(np.mean([roc_auc_score(y, dd) for dd in dec]))


def parear(d, pool_a, pool_b, seed, tol=0.5, n_max=None):
    rng = np.random.default_rng(seed)
    sx, idade = d['sexo_f'].values, d['idade'].values
    livres, pares = set(pool_a.tolist()), []
    for c in rng.permutation(pool_b):
        cand = [t for t in livres if sx[t] == sx[c] and abs(idade[t] - idade[c]) <= tol]
        if not cand:
            continue
        dif = np.array([abs(idade[t] - idade[c]) for t in cand])
        menores = [cand[j] for j in np.flatnonzero(dif <= dif.min() + 0.05)]
        t = menores[rng.integers(len(menores))]
        livres.discard(t)
        pares.append((t, c))
        if n_max and len(pares) >= n_max:
            break
    return pares
