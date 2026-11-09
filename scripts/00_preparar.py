import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from eegtdah import config, preparar, sinais


def validar(cfg, bloco, n=5):
    mapa = pd.read_csv(cfg['mapa_ids_csv'])
    meta = pd.read_csv(Path(cfg['pacote_dir']) / bloco / 'metadados.csv')
    banco, cond = bloco.split('_')
    linhas = []
    for cod in meta['codigo'][:n]:
        sid = mapa.loc[(mapa['codigo'] == cod) & mapa['banco'].str.startswith(banco), 'id_original'].iloc[0]
        ref, segs_ref = sinais.carregar(cfg['pacote_dir'], bloco, cod)
        if banco == 'tdbrain':
            x, canais = preparar.tdbrain(cfg['tdbrain_dir'], sid, cond)
            segs = [(0, x.shape[1])]
        else:
            x, canais, segs = preparar.hbn(cfg['hbn_dir'], sid)[cond]
        x = x[[canais.index(c) for c in config.CANAIS]]
        m = min(x.shape[1], ref.shape[1])
        r = [np.corrcoef(x[i, :m], ref[i, :m])[0, 1] if ref[i].std() > 0 else np.nan for i in range(len(x))]
        erro = [np.abs(x[i, :m] - ref[i, :m]).max() / ref[i].std() if ref[i].std() > 0 else np.nan for i in range(len(x))]
        linhas.append({'bloco': bloco, 'dif_amostras': x.shape[1] - ref.shape[1], 'segs_iguais': np.array_equal(
            np.asarray(segs), segs_ref), 'r_min': np.nanmin(r), 'r_mediana': np.nanmedian(r),
            'n_canais_r_menor_099': int(np.sum(np.array(r) < 0.99)), 'n_canais_planos': int(np.isnan(r).sum()),
            'erro_max_sobre_dp': np.nanmax(erro)})
    return linhas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bancos', default='adhdata,tdbrain,hbn')
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--baixar-hbn', nargs='*', metavar='participants.tsv')
    ap.add_argument('--validar', nargs='*', metavar='BLOCO')
    ap.add_argument('--meu-banco', action='store_true', help='monta o pacote dos seus dados (chaves meu_banco_* do config.toml)')
    ap.add_argument('--n', type=int, default=5)
    a = ap.parse_args()
    cfg = config.ler()
    if a.meu_banco:
        print('bloco %s: %d registros' % preparar.montar_generico(cfg))
    elif a.baixar_hbn:
        preparar.baixar_hbn(cfg, a.baixar_hbn, a.jobs)
    elif a.validar is not None:
        blocos = a.validar or ['tdbrain_fechado', 'tdbrain_aberto', 'hbn_fechado', 'hbn_aberto']
        t = pd.DataFrame([l for b in blocos for l in validar(cfg, b, a.n)])
        print(t.groupby('bloco').agg(n=('r_min', 'size'), r_mediana=('r_mediana', 'median'), r_min=('r_min', 'min'),
                                      canais_r_menor_099=('n_canais_r_menor_099', 'sum'),
                                      canais_planos=('n_canais_planos', 'sum'), dif_amostras_max=('dif_amostras', 'max'),
                                      dif_amostras_min=('dif_amostras', 'min'), segs_iguais=('segs_iguais', 'all'),
                                      erro_max_sobre_dp=('erro_max_sobre_dp', 'max')).round(5).to_string())
    else:
        preparar.montar(cfg, tuple(a.bancos.split(',')), a.jobs)


if __name__ == '__main__':
    main()
