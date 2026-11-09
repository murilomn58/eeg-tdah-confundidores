import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from eegtdah import estatistica, preparar, sinais


def test_fooof_recupera_expoente_e_offset():
    f = np.arange(0.25, 63.25, 0.25)
    for chi, off in ((1.0, 1.5), (1.6, 2.2), (2.2, 0.8)):
        m = sinais.fooof(f, 10 ** off * f ** -chi * (1 + 0.8 * np.exp(-0.5 * ((f - 10) / 1.2) ** 2)))
        assert abs(m.aperiodic_params_[1] - chi) < 0.02 and abs(m.aperiodic_params_[0] - off) < 0.02


def test_tangente_invariante_a_ganho_global():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(4, 18, 200))
    C = np.einsum('nct,ndt->ncd', A, A) / 200
    assert np.allclose(estatistica.logs_tangente(C), estatistica.logs_tangente(37.5 * C), atol=1e-10)


def test_quantizacao_int16_ida_e_volta():
    x = np.random.default_rng(1).normal(scale=[[5.0], [80.0], [0.3]], size=(3, 5000))
    q, escala = preparar.quantizar(x)
    assert q.dtype == np.int16 and np.abs(q).max() <= 32000
    assert np.all(np.abs(q * escala[:, None] - x).max(axis=1) <= escala / 2 + 1e-12)


if __name__ == '__main__':
    for nome, fn in list(globals().items()):
        if nome.startswith('test_'):
            fn()
            print('ok', nome)
