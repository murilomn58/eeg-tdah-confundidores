import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from eegtdah import figuras

ap = argparse.ArgumentParser(description='Gera as figuras do artigo a partir de resultados_agregados/.')
ap.add_argument('--so', nargs='*', choices=list(figuras.FIGS), help='só estas figuras')
ap.add_argument('--res', help='pasta dos resultados agregados (padrão: resultados_agregados/)')
ap.add_argument('--saida', help='pasta de saída (padrão: figuras/)')
ap.add_argument('--braco', help='braco_A.csv das replicações (padrão: <saida_dir>/agregados/braco_A.csv)')
a = ap.parse_args()
figuras.gerar(a.res, a.saida, a.so, a.braco)
