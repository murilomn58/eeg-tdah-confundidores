import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from eegtdah import documento

ap = argparse.ArgumentParser(description='Monta o artigo em .docx (curto ou longo).')
ap.add_argument('versao', nargs='?', choices=list(documento.VERSOES), default='curto')
ap.add_argument('--saida')
ap.add_argument('--figuras', default=str(documento.FIGURAS))
a = ap.parse_args()
saida, n, pendentes = documento.montar(a.versao, a.saida, a.figuras)
print(f'{saida}: {n} referências citadas, {len(pendentes)} pendentes ou não verificadas')
