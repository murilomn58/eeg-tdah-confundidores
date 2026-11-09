import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from eegtdah import config, preparar, sinais

FS0, DUR = 256, 60
ROTULOS = ['EEG Fp1-REF', 'EEG Fp2-REF', 'EEG F7-REF', 'EEG F3-REF', 'EEG Fz-REF', 'EEG F4-REF', 'EEG F8-REF', 'EEG T3-REF',
           'EEG C3-REF', 'EEG Cz-REF', 'EEG C4-REF', 'EEG T4-REF', 'EEG T5-REF', 'EEG P3-REF', 'EEG Pz-REF', 'EEG P4-REF',
           'EEG T6-REF', 'EEG O1-REF', 'EEG O2-REF', 'ECG']


def escrever_edf(caminho, uv, rotulos, fs):
    ns, nrec, pmin, pmax = len(rotulos), uv.shape[1] // fs, -500.0, 500.0
    campo = lambda v, n: f'{v}'[:n].ljust(n).encode('ascii')
    cab = b''.join([campo(0, 8), campo('X', 80), campo('X', 80), campo('01.01.01', 8), campo('00.00.00', 8),
                    campo(256 * (ns + 1), 8), campo('', 44), campo(nrec, 8), campo(1, 8), campo(ns, 4)])
    for valores, n in ((rotulos, 16), (['AgAgCl'] * ns, 80), (['uV'] * ns, 8), ([pmin] * ns, 8), ([pmax] * ns, 8),
                       ([-32768] * ns, 8), ([32767] * ns, 8), ([''] * ns, 80), ([fs] * ns, 8), ([''] * ns, 32)):
        cab += b''.join(campo(v, n) for v in valores)
    dig = np.round((uv - pmin) / (pmax - pmin) * 65535 - 32768).clip(-32768, 32767).astype('<i2')
    corpo = dig[:, :nrec * fs].reshape(ns, nrec, fs).transpose(1, 0, 2).tobytes()
    Path(caminho).write_bytes(cab + corpo)


def test_leitor_generico_edf(tmp_path):
    t = np.arange(FS0 * DUR) / FS0
    rng = np.random.default_rng(0)
    alfa = 20 * np.sin(2 * np.pi * 10 * t)
    uv = np.stack([alfa * (1 + 0.05 * k) + 15 * np.sin(2 * np.pi * 50 * t) + rng.normal(0, 3, t.size)
                   for k in range(len(ROTULOS))])
    escrever_edf(tmp_path / 'p1.edf', uv, ROTULOS, FS0)
    escrever_edf(tmp_path / 'p2.edf', uv[::-1], ROTULOS[::-1], FS0)
    x, canais, rede = preparar.generico(tmp_path / 'p1.edf', rede=50)
    assert canais == list(config.CANAIS) and x.shape == (19, DUR * 128) and rede == 50.0
    from mne.time_frequency import psd_array_welch
    p, f = psd_array_welch(x, 128.0, n_fft=512, verbose=False)
    assert p[:, np.argmin(abs(f - 50))].mean() < 0.01 * p[:, np.argmin(abs(f - 10))].mean()
    pd.DataFrame({'arquivo': ['p1.edf', 'p2.edf'], 'grupo': ['TDAH', 'controle'], 'idade': [9.0, 10.0],
                  'sexo': ['M', 'F']}).to_csv(tmp_path / 'lista.csv', index=False)
    cfg = {'meu_banco_csv': tmp_path / 'lista.csv', 'meu_banco_bloco': 'meu_banco', 'meu_banco_rede': None,
           'meu_banco_mapa': None, 'pacote_dir': tmp_path / 'pacote', 'mapa_ids_csv': tmp_path / 'mapa.csv'}
    assert preparar.montar_generico(cfg) == ('meu_banco', 2)
    meta = pd.read_csv(tmp_path / 'pacote' / 'meu_banco' / 'metadados.csv')
    assert list(meta.columns) == ['codigo', 'grupo', 'idade', 'sexo', 'duracao_s'] and meta['duracao_s'].tolist() == [60.0, 60.0]
    a, segs = sinais.carregar(tmp_path / 'pacote', 'meu_banco', meta['codigo'][0])
    b, _ = sinais.carregar(tmp_path / 'pacote', 'meu_banco', meta['codigo'][1])
    assert np.corrcoef(a[0], b[0])[0, 1] > 0.999 and sinais.epocas(a, segs).shape == (15, 19, 512)
