import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from eegtdah import config, ensaios, replicacoes

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Ensaios agregados do artigo e, opcionalmente, as replicações pré-registradas.')
    ap.add_argument('--partes', nargs='*', default=list(ensaios.PARTES))
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--sorteios', type=int, default=20)
    ap.add_argument('--destino', default=str(config.RESULTADOS))
    ap.add_argument('--replicacoes', action='store_true')
    ap.add_argument('--blocos', nargs='*', default=list(config.BLOCOS))
    ap.add_argument('--limite', type=int)
    a = ap.parse_args()
    cfg = config.ler()
    if a.replicacoes:
        replicacoes.extrair(cfg, list(config.BLOCOS), a.jobs, a.limite)
        replicacoes.rodar(cfg, a.jobs)
        replicacoes.figuras(cfg)
    else:
        ensaios.rodar(cfg, a.partes, a.jobs, a.sorteios, a.destino, a.blocos)
