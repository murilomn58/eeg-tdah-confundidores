"""Caminhos (config.toml) e constantes fixas do pipeline."""
import os
import tomllib
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
RESULTADOS = RAIZ / 'resultados_agregados'
FIGURAS = RAIZ / 'figuras'
TEXTO = RAIZ / 'texto'
FS = 128.0
CANAIS = ('Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'T7', 'C3', 'Cz', 'C4', 'T8',
          'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'O2')
C18 = tuple(c for c in CANAIS if c != 'Cz')
BORDA = ('Fp1', 'Fp2', 'F7', 'F8', 'T7', 'T8')
BORDA8 = BORDA + ('P7', 'P8')
POST = ('P3', 'Pz', 'P4', 'O1', 'O2')
FC = ('Fz', 'F3', 'F4', 'C3', 'C4', 'Cz')
FC5 = ('Fz', 'F3', 'F4', 'C3', 'C4')
BLOCOS = ('adhdata_tarefa', 'hbn_fechado', 'hbn_aberto', 'tdbrain_fechado', 'tdbrain_aberto')
BANCO = {'adhdata_tarefa': ('adhdata', 'tarefa'), 'hbn_fechado': ('hbn', 'OF'), 'hbn_aberto': ('hbn', 'OA'),
         'tdbrain_fechado': ('tdbrain', 'OF'), 'tdbrain_aberto': ('tdbrain', 'OA')}
FOOOF = dict(aperiodic_mode='fixed', peak_width_limits=(1.0, 8.0), max_n_peaks=6, min_peak_height=0.0,
             peak_threshold=2.0, verbose=False)
CHAVES = ('pacote_dir', 'saida_dir', 'adhdata_csv', 'tdbrain_dir', 'tdbrain_xlsx', 'hbn_dir', 'hbn_fenotipo_csv',
          'hbn_sintomas_csv', 'mapa_ids_csv', 'meu_banco_csv', 'meu_banco_mapa')
VALORES = ('meu_banco_bloco', 'meu_banco_rede')


def ler(caminho=None):
    p = Path(caminho or os.environ.get('EEGTDAH_CONFIG') or RAIZ / 'config.toml')
    if not p.exists():
        raise FileNotFoundError(f'{p} não existe: copie config.example.toml para config.toml e ajuste os caminhos')
    c = tomllib.loads(p.read_text(encoding='utf-8'))
    return {**{k: Path(c[k]) if c.get(k) else None for k in CHAVES}, **{k: c.get(k) or None for k in VALORES}}


def indices(nomes, canais=CANAIS):
    return [list(canais).index(c) for c in nomes if c in canais]


def isolar_threads():
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(v, '1')
