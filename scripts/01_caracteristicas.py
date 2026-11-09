import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from eegtdah import caracteristicas, config

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Características por participante em cada variante (fica em saida_dir, sob DUA).')
    ap.add_argument('--blocos', nargs='*', default=list(config.BLOCOS))
    ap.add_argument('--variantes', nargs='*', default=list(caracteristicas.VARIANTES))
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--limite', type=int)
    ap.add_argument('--refazer', action='store_true')
    a = ap.parse_args()
    caracteristicas.extrair(config.ler(), a.blocos, a.jobs, a.variantes, a.limite, a.refazer)
