"""Do registro bruto ao pacote: 19 canais, passa-altas 0,5 Hz, notch na rede, 128 Hz, int16 com escala por canal."""
import json
import re
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .config import CANAIS, FS

REDES = (50.0, 60.0)
ALVOS = ('instructed_toOpenEyes', 'instructed_toCloseEyes')
CORTES = ('boundary', 'break cnt')
MAPA_EGI = json.loads(Path(__file__).with_name('mapa_egi_1020.json').read_text(encoding='utf-8'))
ANTIGOS = {'T3': 'T7', 'T4': 'T8', 'T5': 'P7', 'T6': 'P8'}
NEMAR = {'R1': 'on005505', 'R2': 'on005506', 'R3': 'on005507', 'R4': 'on005508', 'R5': 'on005509', 'R6': 'on005510',
         'R7': 'on005511', 'R8': 'on005512', 'R9': 'on005514', 'R10': 'on005515', 'R11': 'on005516', 'NC': 'nm000103'}


def detectar_rede(raw):
    from mne.time_frequency import psd_array_welch
    sf = raw.info['sfreq']
    viaveis = [f for f in REDES if f < sf / 2 * 0.9]
    d = raw.get_data(picks='eeg')
    if not viaveis or d[:, :int(sf * 60)].std(axis=1).max() * 1e6 < 1e-3:
        return None
    d = d[::max(1, len(d) // 20)][:20]
    ini = int(10 * sf) if int(10 * sf) < d.shape[1] else 0
    d = d[:, ini:ini + int(60 * sf)]
    psd, f = psd_array_welch(d, sfreq=sf, fmin=1.0, fmax=min(sf / 2 * 0.95, 100.0), n_fft=min(int(sf * 4), d.shape[1]),
                             verbose=False)
    psd, r = psd.mean(0), {}
    for f0 in viaveis:
        dentro = (f >= f0 - 1) & (f <= f0 + 1)
        fora = ((f >= f0 - 5) & (f <= f0 - 2)) | ((f >= f0 + 2) & (f <= f0 + 5))
        if dentro.any() and fora.any() and psd[dentro].max() > 0 and np.median(psd[fora]) > 0:
            r[f0] = 10 * np.log10(psd[dentro].max() / np.median(psd[fora]))
    o = sorted(r.items(), key=lambda p: p[1], reverse=True)
    if not o or o[0][1] < 3.0 or (len(o) > 1 and o[0][1] - o[1][1] < 3.0):
        return None
    return o[0][0]


def filtrar(raw, rede=None):
    rede = detectar_rede(raw) if rede is None else rede
    out = raw.copy().filter(l_freq=0.5, h_freq=None, verbose=False)
    harm = [rede * k for k in (1, 2, 3) if rede and rede * k < raw.info['sfreq'] / 2 * 0.9]
    if harm:
        out.notch_filter(freqs=harm, notch_widths=1.0, verbose=False)
    return out, rede


def raw_de(uv, canais, fs):
    import mne
    return mne.io.RawArray(np.asarray(uv, float) * 1e-6, mne.create_info(list(canais), float(fs), 'eeg'), verbose='ERROR')


def quantizar(x):
    escala = np.maximum(np.abs(x).max(axis=1), 1e-9) / 32000.0
    return np.round(x / escala[:, None]).astype(np.int16), escala


def gravar(destino, bloco, cod, x, canais, segs):
    x = np.asarray(x, float)[[list(canais).index(c) for c in CANAIS]]
    q, escala = quantizar(x)
    (Path(destino) / bloco).mkdir(parents=True, exist_ok=True)
    np.savez_compressed(Path(destino) / bloco / f'{cod}.npz', x=q, escala=escala, segs=np.asarray(segs, int))


def adhdata(csv):
    df = pd.read_csv(csv)
    for sid in pd.unique(df['ID']):
        s = df[df['ID'] == sid]
        filt, _ = filtrar(raw_de(s[list(CANAIS)].T.values, CANAIS, FS))
        yield sid, ('TDAH' if str(s['Class'].iloc[0]).upper() == 'ADHD' else 'controle'), filt.get_data() * 1e6


def tdbrain_arquivo(raiz, sid, cond):
    ses = sorted((p.name for p in (Path(raiz) / sid).glob('ses-*')), key=lambda s: int(s.split('-')[1]))[0]
    return Path(raiz) / sid / ses / 'eeg' / f'{sid}_{ses}_task-{ {"aberto": "restEO", "fechado": "restEC"}[cond]}_eeg.bdf'


def tdbrain(raiz, sid, cond):
    import mne
    raw = mne.io.read_raw_bdf(tdbrain_arquivo(raiz, sid, cond), preload=True, verbose='ERROR').pick(list(CANAIS))
    filt, _ = filtrar(raw, rede=50.0)
    filt = filt.resample(FS, verbose='ERROR')
    return filt.get_data() * 1e6, list(filt.ch_names)


def blocos_eventos(ev, n):
    amostra = lambda o: int(round(float(o) * FS))
    cortes = sorted({amostra(o) for o, v in zip(ev['onset'], ev['value']) if v in CORTES and 0 < amostra(o) < n})
    marcos = sorted((amostra(o), v) for o, v in zip(ev['onset'], ev['value']) if v in ALVOS and 0 <= amostra(o) < n)
    out = []
    for (ini, valor), (fim, _) in zip(marcos, marcos[1:]):
        bordas = [ini] + [c for c in cortes if ini < c < fim] + [fim]
        out += [(a, b, valor) for a, b in zip(bordas, bordas[1:]) if fim > ini]
    return out


def hbn(pasta, sid):
    z = np.load(Path(pasta) / sid / f'{sid}_rest19.npz')
    canais = [str(c) for c in z['canais']]
    filt, _ = filtrar(raw_de(z['digital'] * z['scale'][:, None] + z['offset'][:, None], canais, z['fs']))
    x = filt.resample(FS, verbose='ERROR').get_data() * 1e6
    ev = pd.read_csv(Path(pasta) / sid / f'{sid}_task-RestingState_events.tsv', sep='\t')
    blocos, out = blocos_eventos(ev, x.shape[1]), {}
    for alvo, cond in zip(ALVOS, ('aberto', 'fechado')):
        segs = [(a, b) for a, b, v in blocos if v == alvo]
        if segs:
            pos = np.cumsum([0] + [b - a for a, b in segs])
            out[cond] = (np.concatenate([x[:, a:b] for a, b in segs], axis=1), canais, list(zip(pos[:-1], pos[1:])))
    return out


def _tdbrain(args):
    raiz, destino, sid, cod, cond = args
    x, canais = tdbrain(raiz, sid, cond)
    gravar(destino, f'tdbrain_{cond}', cod, x, canais, [(0, x.shape[1])])
    return x.shape[1] / FS


def _hbn(args):
    pasta, destino, sid, cod = args
    dur = {}
    for cond, (x, canais, segs) in hbn(pasta, sid).items():
        gravar(destino, f'hbn_{cond}', cod, x, canais, segs)
        dur[cond] = x.shape[1] / FS
    return dur


def montar(cfg, bancos=('adhdata', 'tdbrain', 'hbn'), jobs=4):
    destino, mapa = Path(cfg['pacote_dir']), []
    if 'adhdata' in bancos:
        meta = []
        for k, (sid, grupo, x) in enumerate(adhdata(cfg['adhdata_csv'])):
            cod = f'A{k + 1:04d}'
            gravar(destino, 'adhdata_tarefa', cod, x, CANAIS, [(0, x.shape[1])])
            meta.append({'codigo': cod, 'grupo': grupo, 'idade': np.nan, 'sexo': '', 'duracao_s': x.shape[1] / FS})
            mapa.append({'banco': 'adhdata', 'id_original': sid, 'codigo': cod})
        pd.DataFrame(meta).to_csv(destino / 'adhdata_tarefa' / 'metadados.csv', index=False)
    if 'tdbrain' in bancos:
        x = pd.read_excel(cfg['tdbrain_xlsx'])
        x['age'] = pd.to_numeric(x['age'].astype(str).str.replace(',', '.'), errors='coerce')
        p = x.sort_values('sessID').drop_duplicates('TDBRAIN_ID')
        p = p[p['indication'].isin(['ADHD', 'HEALTHY'])].set_index('TDBRAIN_ID')
        raiz = Path(cfg['tdbrain_dir'])
        ids = [d.name for d in sorted(raiz.glob('sub-*')) if d.name in p.index
               and all(tdbrain_arquivo(raiz, d.name, c).exists() for c in ('aberto', 'fechado'))]
        cods = {s: f'T{k + 1:04d}' for k, s in enumerate(ids)}
        for cond in ('fechado', 'aberto'):
            with ProcessPoolExecutor(jobs) as ex:
                durs = list(ex.map(_tdbrain, [(raiz, destino, s, cods[s], cond) for s in ids]))
            pd.DataFrame([{'codigo': cods[s], 'grupo': 'TDAH' if p.at[s, 'indication'] == 'ADHD' else 'controle',
                           'idade': float(p.at[s, 'age']), 'sexo': {0: 'F', 1: 'M'}.get(p.at[s, 'gender'], ''),
                           'duracao_s': d} for s, d in zip(ids, durs)]).to_csv(
                destino / f'tdbrain_{cond}' / 'metadados.csv', index=False)
        mapa += [{'banco': 'tdbrain', 'id_original': s, 'codigo': c} for s, c in cods.items()]
    if 'hbn' in bancos:
        pasta = Path(cfg['hbn_dir'])
        inv = pd.read_csv(pasta / 'inventario_v7.csv')
        inv = inv[[(pasta / s / f'{s}_rest19.npz').exists() for s in inv['participant_id']]]
        cods = {s: f'H{k + 1:04d}' for k, s in enumerate(inv['participant_id'])}
        with ProcessPoolExecutor(jobs) as ex:
            durs = dict(zip(cods, ex.map(_hbn, [(pasta, destino, s, c) for s, c in cods.items()])))
        for cond in ('fechado', 'aberto'):
            pd.DataFrame([{'codigo': cods[r.participant_id], 'grupo': r.grupo, 'idade': r.age, 'sexo': r.sex,
                           'duracao_s': durs[r.participant_id][cond]} for r in inv.itertuples()
                          if cond in durs[r.participant_id]]).to_csv(destino / f'hbn_{cond}' / 'metadados.csv', index=False)
        mapa += [{'banco': 'hbn', 'id_original': s, 'codigo': c} for s, c in cods.items()]
    if mapa and cfg.get('mapa_ids_csv'):
        pd.DataFrame(mapa).to_csv(cfg['mapa_ids_csv'], index=False)


def inventario_hbn(cfg, participants_tsv):
    fen = pd.read_csv(cfg['hbn_sintomas_csv'], low_memory=False)
    fen = fen[fen['Identifiers'].str.match(r'[A-Z]{4}[A-Z0-9]{8}', na=False)].copy()
    fen['EID'] = fen['Identifiers'].str.split(',').str[0]
    dx = fen[[f'Diagnosis_ClinicianConsensus,DX_0{i}' for i in range(1, 5)]].fillna('').astype(str)
    fen['TDAH'] = dx.agg(' | '.join, axis=1).str.contains('ADHD', case=False)
    fen['controle'] = dx.iloc[:, 0].str.strip().eq('No Diagnosis Given')
    part = pd.concat([pd.read_csv(f, sep='\t') for f in participants_tsv])
    part['EID'] = part['participant_id'].str.replace('sub-', '', regex=False)
    part = part[part['RestingState'].astype(str).str.contains('available')]
    part['age'] = pd.to_numeric(part['age'], errors='coerce')
    m = part.merge(fen[['EID', 'TDAH', 'controle']], on='EID')
    m = m[(m.TDAH | m.controle) & m.age.between(5.0, 18.0, inclusive='left') & m.release_number.isin(list(NEMAR))]
    m['grupo'] = np.where(m.TDAH, 'TDAH', 'controle')
    m = m[['participant_id', 'release_number', 'grupo', 'age', 'sex']].drop_duplicates('participant_id')
    m['oid'] = m['release_number'].map(NEMAR)
    return m


@lru_cache
def _indice(oid):
    import requests
    idx = requests.get(f'https://zarr.nemar.org/{oid}/zarr/index.json', timeout=120).json()['stores']
    return {s['path'].split('/')[0]: s['zarr'] for s in idx if 'task-RestingState' in s['path']}


def _baixar_hbn(args):
    import requests
    from numcodecs import Blosc
    pasta, sid, oid = args
    alvo = Path(pasta) / sid / f'{sid}_rest19.npz'
    if alvo.exists():
        return 'ja_existe'
    zarr = _indice(oid)[sid]
    base = f'https://nemar.s3.us-east-2.amazonaws.com/{oid}/zarr/{zarr}/eeg_250hz'
    grupo, arr = requests.get(f'{base}/zarr.json').json()['attributes'], requests.get(f'{base}/0/zarr.json').json()
    n_ch, n = arr['shape']
    shard, inner = arr['chunk_grid']['configuration']['chunk_shape'][1], arr['codecs'][0]['configuration']['chunk_shape'][1]
    rot = {('E129' if c['label'] == 'Cz' else c['label']): int(c['row_index']) for c in grupo['channels']}
    linhas, partes = [rot[MAPA_EGI[c]] for c in CANAIS], []
    for k in range(int(np.ceil(n / shard))):
        raw = requests.get(f'{base}/0/c/0/{k}', timeout=120).content
        ii = np.frombuffer(raw[-(shard // inner) * 16 - 4:-4], dtype='<u8').reshape(-1, 2)
        partes += [np.zeros((n_ch, inner), np.int16) if o == 2 ** 64 - 1 else
                   np.frombuffer(Blosc().decode(raw[o:o + b]), dtype='<i2').reshape(n_ch, inner)[:] for o, b in ii]
    dig = np.concatenate(partes, axis=1)[linhas, :n]
    ev = requests.get(f'https://data.nemar.org/{oid}/latest/{sid}/eeg/{sid}_task-RestingState_events.tsv', timeout=120)
    alvo.parent.mkdir(parents=True, exist_ok=True)
    (alvo.parent / f'{sid}_task-RestingState_events.tsv').write_bytes(ev.content)
    np.savez_compressed(alvo, digital=dig, scale=np.asarray(arr['attributes']['scale'])[linhas],
                        offset=np.asarray(arr['attributes']['offset'])[linhas], canais=np.array(CANAIS),
                        fs=float(grupo['rate']))
    return 'ok'


def baixar_hbn(cfg, participants_tsv, jobs=6):
    pasta = Path(cfg['hbn_dir'])
    inv = inventario_hbn(cfg, participants_tsv)
    pasta.mkdir(parents=True, exist_ok=True)
    inv.to_csv(pasta / 'inventario_v7.csv', index=False)
    with ThreadPoolExecutor(jobs) as ex:
        return list(ex.map(_baixar_hbn, [(pasta, s, o) for s, o in zip(inv['participant_id'], inv['oid'])]))


def nome_1020(nome, mapa=None):
    if mapa and nome in mapa:
        return mapa[nome]
    n = re.sub(r'^EEG[\s_-]*|[\s_-]+(REF|LE|AVG|AV|A1|A2|M1|M2)$', '', nome.strip(), flags=re.I).upper()
    n = ANTIGOS.get(n, n)
    return next((c for c in CANAIS if c.upper() == n), None)


def generico(arquivo, rede=None, mapa=None):
    import mne
    raw = mne.io.read_raw(arquivo, preload=True, verbose='ERROR')
    achados = {}
    for c in raw.ch_names:
        alvo = nome_1020(c, mapa)
        if alvo and alvo not in achados:
            achados[alvo] = c
    faltam = [c for c in CANAIS if c not in achados]
    if faltam:
        raise ValueError(f'{Path(arquivo).name}: canais 10-20 ausentes: {", ".join(faltam)}')
    raw = raw.pick([achados[c] for c in CANAIS]).rename_channels({achados[c]: c for c in CANAIS})
    raw.set_channel_types({c: 'eeg' for c in CANAIS})
    filt, rede = filtrar(raw, None if rede is None else float(rede))
    if filt.info['sfreq'] != FS:
        filt = filt.resample(FS, verbose='ERROR')
    return filt.get_data() * 1e6, list(filt.ch_names), rede


def montar_generico(cfg):
    lista = pd.read_csv(cfg['meu_banco_csv'])
    bloco, base = cfg.get('meu_banco_bloco') or 'meu_banco', Path(cfg['meu_banco_csv']).parent
    mapa = json.loads(Path(cfg['meu_banco_mapa']).read_text(encoding='utf-8')) if cfg.get('meu_banco_mapa') else None
    meta, ids = [], []
    for k, r in enumerate(lista.itertuples()):
        cod = f'M{k + 1:04d}'
        arq = Path(r.arquivo) if Path(r.arquivo).is_absolute() else base / r.arquivo
        x, canais, _ = generico(arq, cfg.get('meu_banco_rede'), mapa)
        gravar(cfg['pacote_dir'], bloco, cod, x, canais, [(0, x.shape[1])])
        meta.append({'codigo': cod, 'grupo': r.grupo, 'idade': getattr(r, 'idade', np.nan), 'sexo': getattr(r, 'sexo', ''),
                     'duracao_s': x.shape[1] / FS})
        ids.append({'banco': bloco, 'id_original': arq.name, 'codigo': cod})
    pd.DataFrame(meta).to_csv(Path(cfg['pacote_dir']) / bloco / 'metadados.csv', index=False)
    if cfg.get('mapa_ids_csv'):
        pd.DataFrame(ids).to_csv(cfg['mapa_ids_csv'], index=False)
    return bloco, len(meta)
