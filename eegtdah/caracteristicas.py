"""Características por participante na base e nas variantes de processamento (o tabuleiro)."""
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from . import sinais as s
from .config import BORDA, BORDA8, C18, CANAIS, FC5, POST, indices, isolar_threads

VARIANTES = ('base', 'rej_emg', 'rej_olho', 'ica_musc', 'ica_olho', 'internos', 'dur60', 'nativa', 'U30', 'rest', 'csd')
FEATS = ('chi', 'offset', 'chi_post', 'iaf', 'alfa_pw', 'tbr', 'emg_t', 'ocular', 'n_ep')
COVS = ('base', 'rej_emg', 'rej_olho', 'ica_musc', 'ica_olho', 'internos', 'dur60')
AQ = ('emg_t', 'razao_emg', 'log_var', 'borda_int', 'ocular', 'tbr_fc')


def aquisicao(x, segs):
    W = s.epocas(s.car(x[indices(C18)]), segs)
    P, f = s.espectro(W)
    Pm = np.median(P, 0)
    cov = s.covariancia(W, C18)

    def bp(nomes, lo, hi):
        return s.potencia(f, Pm, indices(nomes, C18), lo, hi)
    v = np.var(np.concatenate(list(W), axis=-1), axis=-1)
    internos = [c for c in C18 if c not in BORDA8]
    return {'emg_t': np.log10(bp(('T7', 'T8'), 30, 45)), 'razao_emg': np.log10(bp(C18, 25, 32) / bp(C18, 8, 13)),
            'log_var': np.log10(np.trace(cov)) if np.isfinite(cov).all() else np.nan,
            'borda_int': np.log10(v[indices(BORDA8, C18)].mean() / v[indices(internos, C18)].mean()),
            'ocular': np.log10(bp(('Fp1', 'Fp2'), 1, 4) / bp(('O1', 'O2'), 1, 4)),
            'tbr_fc': np.log10(bp(FC5, 4, 8) / bp(FC5, 13, 30))}


def uma(args):
    pacote, bloco, cod, variantes = args
    try:
        x, segs = s.carregar(pacote, bloco, cod)
        res, covs, extra = {}, {}, {f'aq__{k}': v for k, v in aquisicao(x, segs).items()}
        Wb = s.epocas(s.car(x), segs)
        res['base'], (Pb, fb) = s.medir(Wb)
        covs['base'] = s.covariancia(Wb)
        extra['base__alfa_pa'], extra['base__iaf_cg'] = s.alfa_pa(fb, np.median(np.median(Pb, 0)[indices(POST)], 0))
        emg = s.banda(fb, Pb[:, indices(('T7', 'T8'))], 30, 45).mean(1)
        olho = s.banda(fb, Pb[:, indices(('Fp1', 'Fp2'))], 1, 4).mean(1)
        for nome, e in (('rej_emg', emg), ('rej_olho', olho)):
            if nome in variantes:
                W = Wb[e <= np.percentile(e, 90)]
                res[nome], _ = s.medir(W)
                covs[nome] = s.covariancia(W)
        if {'ica_musc', 'ica_olho'} & set(variantes):
            try:
                for nome, (y, n) in zip(('ica_musc', 'ica_olho'), s.ica_limpar(x).values()):
                    W = s.epocas(s.car(y), segs)
                    res[nome], _ = s.medir(W)
                    covs[nome] = s.covariancia(W)
                    extra[f'n_{nome}'] = n
            except Exception:
                pass
        if 'internos' in variantes:
            ci = tuple(c for c in CANAIS if c not in BORDA)
            Wi = s.epocas(s.car(x[indices(ci)]), segs)
            res['internos'], _ = s.medir(Wi, ci)
            covs['internos'] = s.covariancia(Wi, ci)
        if 'dur60' in variantes:
            res['dur60'], _ = s.medir(Wb[:15])
            covs['dur60'] = s.covariancia(Wb[:15])
        if 'U30' in variantes:
            res['U30'], _ = s.medir(Wb, hi=30.0)
        for ref in ('nativa', 'rest', 'csd'):
            if ref in variantes:
                res[ref], _ = s.medir(s.epocas(s.referenciar(x, ref), segs))
        return cod, res, covs, extra, None
    except Exception as exc:
        return cod, None, None, None, repr(exc)


def extrair(cfg, blocos, jobs=8, variantes=VARIANTES, limite=None, refazer=False):
    isolar_threads()
    saida = Path(cfg['saida_dir']) / 'caracteristicas'
    saida.mkdir(parents=True, exist_ok=True)
    for bloco in blocos:
        alvo = saida / f'{bloco}.npz'
        if alvo.exists() and not refazer:
            print(bloco, 'já existe', flush=True)
            continue
        meta = pd.read_csv(Path(cfg['pacote_dir']) / bloco / 'metadados.csv')
        cods = list(meta.loc[meta['grupo'].isin(['TDAH', 'controle']), 'codigo'])[:limite]
        t0 = time.time()
        with Pool(jobs) as pool:
            res = pool.map(uma, [(cfg['pacote_dir'], bloco, c, tuple(variantes)) for c in cods], chunksize=2)
        ok = [r for r in res if r[4] is None]
        F = {f'{v}__{k}': np.array([r[1].get(v, {}).get(k, np.nan) for r in ok], float) for v in VARIANTES for k in FEATS}
        E = {k: np.array([r[3].get(k, np.nan) for r in ok], float) for k in ok[0][3]} if ok else {}
        E.update({k: np.array([r[3].get(k, np.nan) for r in ok], float) for k in ('n_ica_musc', 'n_ica_olho')})
        C = {f'cov_{v}': np.stack([r[2].get(v, np.full((12, 12) if v == 'internos' else (18, 18), np.nan)) for r in ok])
             for v in COVS}
        np.savez_compressed(alvo, codigo=np.array([r[0] for r in ok]), **F, **E, **C)
        print(f'{bloco}: {len(ok)} ok, {len(res) - len(ok)} erros, {time.time() - t0:.0f} s', flush=True)


def tabela(cfg, bloco):
    meta = pd.read_csv(Path(cfg['pacote_dir']) / bloco / 'metadados.csv')
    meta = meta[meta['grupo'].isin(['TDAH', 'controle'])]
    with np.load(Path(cfg['saida_dir']) / 'caracteristicas' / f'{bloco}.npz') as z:
        d = pd.DataFrame({k: z[k] for k in z.files if not k.startswith('cov_')})
        pos = {c: i for i, c in enumerate(z['codigo'])}
        d = d.merge(meta, on='codigo')
        if bloco.startswith('tdbrain'):
            d = d[d['idade'] >= 18]
        d = d.reset_index(drop=True)
        ii = np.array([pos[c] for c in d['codigo']], int)
        covs = {k[4:]: z[k][ii] for k in z.files if k.startswith('cov_')}
    d['y'] = (d['grupo'] == 'TDAH').astype(int)
    d['sexo_f'] = np.where(d['sexo'] == 'F', 1.0, np.where(d['sexo'] == 'M', 0.0, np.nan))
    d['log_dur'] = np.log(d['duracao_s'].astype(float))
    return d, covs
