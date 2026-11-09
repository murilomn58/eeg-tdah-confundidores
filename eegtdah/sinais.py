"""Sinal do pacote: épocas, espectro, FOOOF, referências, ICA e covariância."""
import json
import warnings
from pathlib import Path

import numpy as np

from .config import CANAIS, FOOOF, FS, indices

warnings.filterwarnings('ignore')
REST = np.array(json.loads(Path(__file__).with_name('rest19.json').read_text(encoding='utf-8'))['R'])


def carregar(pacote_dir, bloco, cod):
    with np.load(Path(pacote_dir) / bloco / f'{cod}.npz') as z:
        return z['x'].astype(np.float64) * z['escala'][:, None], z['segs']


def car(x):
    return x - x.mean(-2, keepdims=True)


def epocas(x, segs, dur=4.0):
    n = int(dur * FS)
    w = [x[:, i:i + n] for a, b in segs for i in range(int(a), int(b) - n + 1, n)]
    if len(w) < 5:
        raise ValueError(f'só {len(w)} épocas de {dur:g} s')
    return np.stack(w)


def espectro(W, multitaper=False):
    from mne.time_frequency import psd_array_multitaper, psd_array_welch
    n = W.shape[-1]
    if multitaper:
        return psd_array_multitaper(W, FS, fmin=0.0, fmax=63.0, bandwidth=2.0, adaptive=False, low_bias=True,
                                    normalization='full', verbose=False)
    return psd_array_welch(W, sfreq=FS, fmin=0.0, fmax=63.0, n_fft=n, n_per_seg=n, n_overlap=0, window='hann',
                           verbose=False)


def banda(f, P, lo, hi):
    s = (f >= lo) & (f <= hi)
    return np.log10(P[..., s].mean(axis=-1))


def potencia(f, Pm, idx, lo, hi):
    s = (f >= lo) & (f <= hi)
    return Pm[idx][:, s].mean()


def fooof(f, p, hi=40.0, modo='fixed'):
    from fooof import FOOOF as Modelo
    if not np.all(np.isfinite(p)) or np.any(p <= 0):
        return None
    m = Modelo(**{**FOOOF, 'aperiodic_mode': modo})
    try:
        m.fit(f, p, [1.0, hi])
    except Exception:
        return None
    return m if m.has_model else None


def razao_pa(m):
    from fooof.sim.gen import gen_aperiodic
    return m.freqs, 10 ** (m.power_spectrum - gen_aperiodic(m.freqs, m.aperiodic_params_))


def alfa_pa(f, p, hi=40.0, lo_a=7, hi_a=13):
    from scipy.signal import savgol_filter
    m = fooof(f, p, hi)
    if m is None:
        return np.nan, np.nan
    ff, r = razao_pa(m)
    pa = savgol_filter(r, 5, 2)
    j = (ff >= lo_a) & (ff <= hi_a)
    e = np.clip(pa[j] - 1, 0, None)
    return float(pa[j].max() - 1), float((ff[j] * e).sum() / e.sum()) if e.sum() > 0 else np.nan


def referenciar(x, ref):
    if ref == 'nativa':
        return x
    if ref == 'car':
        return car(x)
    if ref == 'rest':
        return REST @ x
    if ref == 'csd':
        import mne
        info = mne.create_info(list(CANAIS), FS, 'eeg')
        info.set_montage(mne.channels.make_standard_montage('standard_1020'))
        raw = mne.io.RawArray(x * 1e-6, info, verbose='ERROR')
        raw.set_eeg_reference('average', verbose='ERROR')
        return mne.preprocessing.compute_current_source_density(raw, verbose='ERROR').get_data() * 1e6
    raise ValueError(ref)


def ica_limpar(x, p=0.8):
    import mne
    from mne_icalabel import label_components
    info = mne.create_info(list(CANAIS), FS, 'eeg')
    info.set_montage(mne.channels.make_standard_montage('standard_1020'))
    plano = np.std(x, axis=1) < 1e-9
    raw = mne.io.RawArray(x * 1e-6, info, verbose='ERROR')
    raw.info['bads'] = [CANAIS[i] for i in np.flatnonzero(plano)]
    raw.set_eeg_reference('average', verbose='ERROR')
    ajuste = raw.copy().filter(1.0, FS / 2 - 1, verbose='ERROR')
    ica = mne.preprocessing.ICA(n_components=len(CANAIS) - 1 - int(plano.sum()), method='infomax',
                                fit_params=dict(extended=True), random_state=0, max_iter='auto', verbose='ERROR')
    ica.fit(ajuste, verbose='ERROR')
    rot = label_components(ajuste, ica, method='iclabel')
    fora = {k: [i for i, (l, q) in enumerate(zip(rot['labels'], rot['y_pred_proba'])) if l == k and q > p]
            for k in ('muscle artifact', 'eye blink')}
    return {k: (ica.apply(raw.copy(), exclude=v, verbose='ERROR').get_data() * 1e6, len(v)) for k, v in fora.items()}


def covariancia(W, canais=CANAIS, n_ep=15, lo=1.0, hi=40.0):
    from mne.filter import filter_data
    sem_cz = [i for i, c in enumerate(canais) if c != 'Cz']
    if len(W) < n_ep:
        return np.full((len(sem_cz), len(sem_cz)), np.nan)
    Wc = car(W[:n_ep][:, sem_cz]).astype(np.float64)
    Wf = filter_data(Wc, FS, lo, min(hi, FS / 2 - 1), verbose='ERROR')
    return np.mean([w @ w.T / w.shape[1] for w in Wf], axis=0)


def medir(W, canais=CANAIS, hi=40.0):
    from .config import FC, POST
    P, f = espectro(W)
    Pm = np.median(P, 0)
    ms = [fooof(f, Pm[j], hi) for j in range(len(canais))]
    chis = np.array([m.aperiodic_params_[-1] if m else np.nan for m in ms])
    offs = np.array([m.aperiodic_params_[0] if m else np.nan for m in ms])
    ip, ifc = indices(POST, canais), indices(FC, canais)
    it, ifp, io = indices(('T7', 'T8'), canais), indices(('Fp1', 'Fp2'), canais), indices(('O1', 'O2'), canais)
    out = {'chi': np.nanmedian(chis), 'offset': np.nanmedian(offs), 'chi_post': np.nanmedian(chis[ip]) if ip else np.nan,
           'iaf': np.nan, 'alfa_pw': np.nan, 'n_ep': len(W)}
    if ip:
        m = fooof(f, np.median(Pm[ip], 0), hi)
        if m is not None and len(m.peak_params_):
            pk = [q for q in m.peak_params_ if 7.0 <= q[0] <= 14.0]
            out['iaf'], out['alfa_pw'] = max(pk, key=lambda q: q[1])[:2] if pk else (np.nan, 0.0)
    out['tbr'] = np.log10(potencia(f, Pm, ifc, 4, 8) / potencia(f, Pm, ifc, 13, 30)) if ifc else np.nan
    out['emg_t'] = np.log10(potencia(f, Pm, it, 30, 45)) if it else np.nan
    out['ocular'] = np.log10(potencia(f, Pm, ifp, 1, 4) / potencia(f, Pm, io, 1, 4)) if ifp and io else np.nan
    return out, (P, f)
