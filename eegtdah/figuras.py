"""Figuras do artigo geradas só a partir dos resultados agregados."""
import warnings
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from .config import FIGURAS, RESULTADOS

SERIF = {'font.family': 'serif', 'font.serif': ['Times New Roman', 'DejaVu Serif'], 'savefig.dpi': 300}
SEM_BORDA = {'axes.spines.top': False, 'axes.spines.right': False}
OK = {'laranja': '#D55E00', 'azul': '#0072B2', 'verde': '#009E73', 'ouro': '#E69F00', 'cinza': '#7F7F7F'}
VIRG = matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('-', '−').replace('.', ','))


def estilo(**kw):
    return plt.rc_context({**SERIF, **kw})


def salvar(fig, sai, nome, tight=True, dpi=300):
    sai.mkdir(parents=True, exist_ok=True)
    kw = {'bbox_inches': 'tight'} if tight else {}
    fig.savefig(sai / f'{nome}.png', dpi=dpi, **kw)
    fig.savefig(sai / f'{nome}.pdf', **kw)
    plt.close(fig)


def caixa_arred(ax, x, y, w, h, texto, cor, fs=8.6, alpha=0.14, bold=False):
    box = 'round,pad=0.012,rounding_size=0.02'
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=box, fc=cor, ec=cor, alpha=alpha, lw=0))
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=box, fc='none', ec=cor, lw=1.1))
    ax.text(x + w / 2, y + h / 2, texto, ha='center', va='center', fontsize=fs, fontweight='bold' if bold else 'normal',
            wrap=True)


def seta_ax(ax, x0, y0, x1, y1, cor='#555555'):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle='-|>', color=cor, lw=1.0))


FEAT3 = {'chi': 'Expoente χ', 'alfa': 'Altura do alfa', 'iaf': 'IAF', 'tbr': 'TBR clássica',
         'esp_dist': 'Espacial: distância\nao centro dos controles', 'esp_pc1': 'Espacial: 1º componente'}
CONF3 = {'idade': 'Idade', 'sexo_f': 'Sexo\n(feminino)', 'log_dur': 'Duração\ndo registro',
         'emg_t': 'Músculo\n(T7/T8,\n30–45 Hz)', 'ocular': 'Olho\n(delta\nFp/O)', 'borda_int': 'Variância\nborda/\ninternos',
         'log_var': 'Potência\ntotal (ganho,\ncontato)'}


def mapa_confundidores(res, sai, **_):
    d = pd.read_csv(res / 'e2_p1_impressao_digital.csv')
    paineis = [('hbn', 'OF', 'HBN, crianças, olhos fechados'), ('tdbrain', 'OF', 'TDBRAIN, adultos, olhos fechados'),
               ('adhdata', 'tarefa', 'adhdata, crianças, tarefa')]
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 5.0), gridspec_kw={'wspace': 0.08})
    cmap = plt.get_cmap('RdBu_r')
    for k, (ax, (b, c, tit)) in enumerate(zip(axs, paineis)):
        s = d[(d['banco'] == b) & (d['cond'] == c)].set_index(['feature', 'confundidor'])
        M = np.full((len(FEAT3), len(CONF3)), np.nan)
        sig = np.zeros_like(M, dtype=bool)
        for i, f in enumerate(FEAT3):
            for j, cf in enumerate(CONF3):
                if (f, cf) in s.index:
                    r = s.loc[(f, cf)]
                    M[i, j], sig[i, j] = r['rho'], (r['ic_lo'] > 0) or (r['ic_hi'] < 0)
        ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-0.8, vmax=0.8, aspect='auto')
        for (i, j), v in np.ndenumerate(M):
            if np.isnan(v):
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, color='#EEEEEE', zorder=0))
                ax.text(j, i, 'n/d', ha='center', va='center', fontsize=7, color='#999999')
                continue
            ax.text(j, i, f'{v:+.2f}'.replace('.', ','), ha='center', va='center', fontsize=7.5,
                    color='white' if abs(v) > 0.5 else 'black', fontweight='bold' if sig[i, j] else 'normal',
                    alpha=1 if sig[i, j] else 0.55)
        ax.axhline(3.5, color='black', lw=0.8)
        ax.set_xticks(range(len(CONF3)))
        ax.set_xticklabels(list(CONF3.values()), fontsize=7.6, rotation=0, linespacing=0.95)
        ax.set_yticks(range(len(FEAT3)))
        ax.set_yticklabels(list(FEAT3.values()) if k == 0 else [], fontsize=8.5)
        ax.set_title(f'{"ABC"[k]}  {tit}\n(controles, n = {int(d[(d["banco"] == b) & (d["cond"] == c)]["n_controles"].max())})',
                     loc='left')
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-0.8, 0.8)), ax=axs, fraction=0.015, pad=0.01)
    cb.set_label('ρ de Spearman (só controles)')
    salvar(fig, sai, 'F3_mapa_confundidores')


def pareados(res, sai, **_):
    d = pd.read_csv(res / 'e2_p2_pareados.csv')
    d = d[d['cond'] == 'OF']
    ordem = [('global', 'Amostra pareada completa'), None, ('idade_5-8', 'Idade 5–8 anos'),
             ('idade_8-11', 'Idade 8–11 anos'), ('idade_11-22', 'Idade 11–22 anos'), None,
             ('so_meninos', 'Só meninos'), ('so_meninas', 'Só meninas'), None,
             ('subtipo_desatento', 'Desatento × controle'), ('subtipo_combinado', 'Combinado × controle'),
             ('subtipo_hiperativo', 'Hiperativo × controle'), ('desatento_x_combinado', 'Desatento × combinado'), None,
             ('tdah_sem_comorbidade', 'TDAH sem comorbidade × controle'),
             ('tdah_com_comorbidade', 'TDAH com comorbidade × controle')]
    modelos = [('chi', 'Expoente χ', OK['laranja'], 'o'), ('per', 'Periódico (alfa + IAF)', OK['azul'], 's'),
               ('tbr', 'TBR', OK['ouro'], 'D'), ('esp', 'Espacial (tangente)', OK['verde'], '^'),
               ('chi+per+esp', 'χ + periódico + espacial', 'black', 'P')]
    cols = [m for m, *_ in modelos]
    fig, ax = plt.subplots(figsize=(8.6, 8.2))
    y, ticks, labs = 0, [], []
    for item in ordem:
        if item is None:
            y += 0.6
            continue
        des, rot = item
        real, nulo = d[(d['design'] == des) & (~d['permutado'])], d[(d['design'] == des) & (d['permutado'])]
        if real.empty:
            continue
        if len(nulo):
            p05, p95 = np.nanpercentile(nulo[cols].values, [5, 95])
            ax.add_patch(plt.Rectangle((p05, y - 0.42), p95 - p05, 0.84, color='#E6E6E6', zorder=0, lw=0))
        for k, (m, _, cor, mk) in enumerate(modelos):
            v = real[m].values
            mu = np.mean(v)
            ax.errorbar(mu, y + (k - 2) * 0.15, xerr=[[mu - np.percentile(v, 2.5)], [np.percentile(v, 97.5) - mu]],
                        fmt=mk, color=cor, ms=5.5, lw=1, capsize=0, mfc=cor if m != 'chi+per+esp' else 'white')
        ticks.append(y)
        labs.append(f'{rot}  (n = {int(real["n_pares"].median())} pares)')
        y += 1
    ax.axvline(0.5, color='black', lw=0.8)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labs, fontsize=8.5)
    ax.invert_yaxis()
    ax.set(xlim=(0.30, 0.80), xlabel='AUC por participante (média de 40 pareamentos; barra: 2,5–97,5 %)')
    h = [plt.Line2D([], [], marker=mk, color=cor, ls='', mfc=cor if m != 'chi+per+esp' else 'white', label=nome)
         for m, nome, cor, mk in modelos]
    h.append(plt.Rectangle((0, 0), 1, 1, color='#E6E6E6', label='faixa do nulo (5–95 %, rótulo permutado)'))
    ax.legend(handles=h, loc='lower right', fontsize=8, frameon=False)
    salvar(fig, sai, 'F4_subconjuntos_pareados')


def aquisicao_tdbrain(res, sai, braco=None, **_):
    if braco is None or not Path(braco).exists():
        warnings.warn('braco_A.csv ausente: F5 não gerada')
        return
    p3 = pd.read_csv(res / 'e2_p3_complementaridade.csv')
    p3 = p3[p3['perm'] == 0]
    braco = pd.read_csv(braco)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={'width_ratios': [1.6, 1], 'wspace': 0.55})

    def primeiro(df, col='esp'):
        return df[col].iloc[0] if len(df) else np.nan
    for cond, cor, mk, off in (('OF', OK['verde'], 'o', -0.13), ('OA', OK['azul'], 's', 0.13)):
        bl = f'tdbrain_{"fechado" if cond == "OF" else "aberto"}'
        loc = lambda an, col='esp': primeiro(p3[(p3['banco'] == 'tdbrain') & (p3['cond'] == cond) & (p3['analise'] == an)], col)
        srv = lambda t, c=None: primeiro(braco[(braco['teste'] == t) & (braco['bloco'] == bl) & ((braco['canais'] == c) if c else True)])
        linhas = [('Amostra completa', loc('completo')), ('Só o lote de 121 s', loc('lote_121s')),
                  ('Suporte comum de idade (18–56)', loc('suporte_18_56')),
                  ('Pareado 1:1 por idade e sexo (pré-registrada)', srv('A5c_pareado', 'pareado_idade_sexo')),
                  ('Só eletrodos internos', loc('so_internos')), ('Só eletrodos de borda', loc('so_borda')),
                  ('Assinatura não cerebral sozinha (pré-registrada)', srv('A1')), ('Base C0 (idade, sexo, duração)', loc('completo', 'C0'))]
        for i, (_, v) in enumerate(linhas):
            a1.plot(v, i + off, mk, color=cor, ms=7,
                    label=('olhos fechados' if cond == 'OF' else 'olhos abertos') if i == 0 else None)
        a1.set_yticks(range(len(linhas)))
        a1.set_yticklabels([r for r, _ in linhas], fontsize=8.5)
    a1.axhline(5.5, color='#BBBBBB', lw=0.8, ls='--')
    a1.axvline(0.5, color='black', lw=0.8)
    a1.invert_yaxis()
    a1.set(xlim=(0.45, 0.95), xlabel='AUC do espaço tangente (TDAH × controle)')
    a1.set_title('A  O ganho espacial no TDBRAIN resiste aos controles…', loc='left')
    a1.legend(frameon=False, fontsize=8, loc='lower right')
    rot, vals, cores = [], [], []
    for bl, nm in (('tdbrain_fechado', 'OF'), ('tdbrain_aberto', 'OA')):
        for g, gn in (('grupo=0', 'controles'), ('grupo=1', 'TDAH')):
            r = braco[(braco['teste'] == 'A5a') & (braco['bloco'] == bl) & (braco['representacao'] == 'tangente_prediz_lote_120s')
                   & (braco['canais'] == g)]
            if len(r):
                rot.append(f'{gn}, {nm}')
                vals.append(r['esp'].iloc[0])
                cores.append(OK['verde'] if nm == 'OF' else OK['azul'])
    a2.barh(range(len(vals)), vals, color=cores, height=0.6)
    for i, v in enumerate(vals):
        a2.text(v + 0.01, i, f'{v:.2f}'.replace('.', ','), va='center', fontsize=8.5)
    a2.set_yticks(range(len(vals)))
    a2.set_yticklabels(rot, fontsize=8.5)
    a2.invert_yaxis()
    a2.axvline(0.5, color='black', lw=0.8)
    a2.set(xlim=(0.4, 1.0), xlabel='AUC: tangente prevê o lote de 120 s')
    a2.set_title('B  …mas carrega o lote de gravação', loc='left')
    salvar(fig, sai, 'F5_espacial_aquisicao_tdbrain')


def ponte_artigo2(res, sai, **_):
    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis('off')
    cols = [(0.03, 'Espectral aperiódica\n(expoente χ, offset)', OK['laranja'],
             'Idade (maturação)\nMúsculo e potência total\nTarefa × repouso\nFaixa do ajuste; referência\n'
             'Medicação (literatura)*'),
            (0.36, 'Espectral periódica\n(alfa, IAF, TBR)', OK['azul'],
             'Estado (OF × OA)\nIdade (IAF; TBR)\nSexo (altura do alfa)\nMistura com o χ (TBR)\nComorbidade (dilui)'),
            (0.69, 'Espacial e riemanniana\n(covariância, tangente)', OK['verde'],
             'Aquisição (lote, ganho, contato)\nMúsculo e olho\nMontagem e Cz\nComorbidade (TEA)\nSem relação com a idade')]
    for x, tit, cor, conf in cols:
        caixa_arred(ax, x, 0.80, 0.28, 0.13, tit, cor, fs=9.5, bold=True, alpha=0.22)
        caixa_arred(ax, x, 0.50, 0.28, 0.24, 'Confundidores\n\n' + conf, cor, fs=8.4)
        seta_ax(ax, x + 0.14, 0.80, x + 0.14, 0.745)
        seta_ax(ax, x + 0.14, 0.50, 0.5 if x == 0.36 else (0.40 if x < 0.36 else 0.60), 0.385)
    ax.text(0.03, 0.975, '* não testável com os dados disponíveis (painel toxicológico do HBN não detecta metilfenidato)',
            fontsize=7.5, color='#555555')
    caixa_arred(ax, 0.20, 0.25, 0.60, 0.13,
                'Covariáveis explícitas e estratos: idade, sexo, idade × sexo, duração, condição, índices de artefato,\n'
                'banco/lote; pareamento ou estratificação no treino; validação por participante, sem vazamento',
                OK['cinza'], fs=8.5)
    seta_ax(ax, 0.5, 0.25, 0.5, 0.185)
    caixa_arred(ax, 0.20, 0.04, 0.60, 0.14,
                'Artigo 2: modelo conjunto auditável\nlogit P(TDAH) = β₀ + β₁p_aperiódico + β₂p_periódico + β₃p_espacial '
                '+ Σγ·covariáveis + interações\n(regressão logística bayesiana; peso de cada família visível ao médico)',
                'black', fs=8.6, alpha=0.06)
    salvar(fig, sai, 'F7_ponte_artigo2')


FEATS_T = [('chi', 'χ\nglobal'), ('chi_post', 'χ\nposterior'), ('offset', 'Offset\naperiódico'),
           ('alfa_pw', 'Potência\ndo alfa'), ('iaf', 'IAF'), ('tbr', 'TBR\nclássica'), ('espacial', 'Espacial\n(tangente)')]
GRUPOS_T = [('Amostra (ρ nos controles)', [('idade', 'Idade'), ('sexo_f', 'Sexo'), ('log_dur', 'Duração')]),
            ('Artefato e aquisição (ρ nos controles)',
             [('emg_t', 'Músculo'), ('ocular', 'Olho'), ('log_var', 'Potência total\n(ganho, contato)')]),
            ('Processamento (mudança no efeito)',
             [('rej_emg', 'Rejeitar épocas: músculo'), ('rej_olho', 'Rejeitar épocas: olho'), ('ica_olho', 'ICA ocular'),
              ('internos', 'Sem canais de borda'), ('dur60', 'Só 60 s'), ('nativa', 'Referência nativa'),
              ('U30', 'Ajuste 1–30 Hz')]),
            ('Diagnóstico (base)', [('sinal', 'TDAH × controle')])]
VERDE, AMAR, VERM, CINZA, AZUL, AZUL_CL = '#4DAF7C', '#F2C14E', '#D1495B', '#E3E3E3', '#1F4E79', '#BFD3E6'
NAO_SE_APLICA = {('tbr', 'U30'), ('espacial', 'dur60'), ('espacial', 'U30')}
LIMIARES = {'rho_controles': (0.2, 0.4), 'delta_g': (0.10, 0.20), 'delta_auc': (0.03, 0.07)}


def cor_casa(tipo, v, exclui0=None):
    if v is None or not np.isfinite(v):
        return CINZA
    if tipo in LIMIARES:
        a, (l1, l2) = abs(v), LIMIARES[tipo]
        return VERDE if a < l1 else AMAR if a < l2 else VERM
    if tipo == 'g':
        return AZUL if exclui0 else AZUL_CL
    return (AZUL if v >= 0.60 else AZUL_CL) if tipo == 'auc' else CINZA


def numero(v):
    if v is None or not np.isfinite(v):
        return ''
    s = f'{abs(v):.2f}'.lstrip('0')
    return ('−' if v < 0 and s != '.00' else '') + s.replace('.', ',')


def quadrado(ax, xy, w, h, c):
    if c == 'white':
        ax.add_patch(plt.Rectangle(xy, w, h, facecolor='white', edgecolor='#BBBBBB', lw=0.4, hatch='////'))
    else:
        ax.add_patch(plt.Rectangle(xy, w, h, color=c, lw=0))


def tabuleiro(res, sai, cond='OF', numeros=False, subtipos=False, nome=None, **_):
    d = pd.read_csv(res / 'e6_tabuleiro_longo.csv')
    grupos = GRUPOS_T
    if subtipos:
        d = pd.concat([d, pd.read_csv(res / 'e8_subtipos_longo.csv')], ignore_index=True)
        grupos = GRUPOS_T[:-1] + [('Diagnóstico (base)', [('sinal', 'TDAH (todos) × controle'),
                                                          ('sinal_F90.0', 'F90.0 desatento × controle'),
                                                          ('sinal_F90.1', 'F90.1 hiperativo × controle'),
                                                          ('sinal_F90.2', 'F90.2 combinado × controle')])]
    suf = 'fechado' if cond == 'OF' else 'aberto'
    blocos = [('adhdata_tarefa', 'a'), (f'hbn_{suf}', 'H'), (f'tdbrain_{suf}', 'T')]
    fig, ax = plt.subplots(figsize=(7.0, 8.0) if numeros else (6.7, 7.4))
    W, H, gap, X0 = (1.55, 0.66, 0.06, -3.3) if numeros else (1.0, 0.62, 0.06, -3.9)
    larg = (W - 2 * gap) / 3
    linhas, y, ys_grupo = [], 0.0, []
    for g, cs in grupos:
        y -= 0.5
        ys_grupo.append((g, y + 0.08))
        for c, r in cs:
            linhas.append((c, r, y))
            y -= H
    idx = d.set_index(['bloco', 'feature', 'coluna'])
    for j, (f, rot) in enumerate(FEATS_T):
        x = j * W
        for c, _, yy in linhas:
            for k, (bl, _) in enumerate(blocos):
                txt, cc = '', CINZA
                if (f, c) in NAO_SE_APLICA:
                    cc = 'white'
                elif (bl, f, c) in idx.index:
                    s = idx.loc[[(bl, f, c)]].iloc[0]
                    ex = bool(s['exclui0']) if 'exclui0' in s and pd.notna(s['exclui0']) else None
                    cc, txt = cor_casa(s['tipo'], s['valor'], ex), numero(s['valor'])
                x_q, w_q = x + gap + k * larg, larg - 0.025
                quadrado(ax, (x_q, yy - H + gap), w_q, H - 2 * gap, cc)
                if numeros and txt:
                    ax.text(x_q + w_q / 2, yy - H / 2, txt, ha='center', va='center', fontsize=5.8,
                            color='white' if cc in (VERM, AZUL) else '#1A1A1A')
        ax.text(x + W / 2, 0.36, rot, ha='center', va='bottom', fontsize=7.8, linespacing=0.95)
        for k, (_, ab) in enumerate(blocos):
            ax.text(x + gap + k * larg + larg / 2, 0.05, ab, ha='center', va='bottom', fontsize=6, color='#555555')
    for c, r, yy in linhas:
        ax.text(-0.12, yy - H / 2, r, ha='right', va='center', fontsize=7.8, linespacing=0.95)
    for g, yy in ys_grupo:
        ax.text(X0, yy, g, ha='left', va='bottom', fontsize=7.8, fontweight='bold')
        ax.plot([X0, len(FEATS_T) * W], [yy - 0.02, yy - 0.02], color='#888888', lw=0.5)
    itens = [(VERDE, 'pouco sensível'), (AMAR, 'moderado'), (VERM, 'muito sensível'), (CINZA, 'não disponível'),
             ('white', 'não se aplica'), (AZUL, 'sinal (IC exclui 0 / AUC ≥ 0,60)'), (AZUL_CL, 'sem sinal')]
    esc = 1.27 if numeros else 1.0
    pos = [(X0, y - 0.5), (X0 + 2.0 * esc, y - 0.5), (X0 + 3.75 * esc, y - 0.5), (X0 + 5.75 * esc, y - 0.5),
           (X0 + 7.85 * esc, y - 0.5), (X0, y - 0.95), (X0 + 3.75 * esc, y - 0.95)]
    for (c, tx), (xl, yl) in zip(itens, pos):
        quadrado(ax, (xl, yl), 0.28, 0.24, c)
        ax.text(xl + 0.38, yl + 0.12, tx, va='center', fontsize=7.4)
    ax.text(X0, y - 1.4, f'a = adhdata (tarefa); H = HBN; T = TDBRAIN ({"olhos fechados" if cond == "OF" else "olhos abertos"}).',
            fontsize=7.4)
    y_fim = y - 1.6
    if numeros:
        ax.text(X0, y - 1.8, 'Números: ρ de Spearman nos controles (amostra, artefato e aquisição); Δg ou ΔAUC em relação '
                'à base (processamento);\ng de Hedges ou AUC (diagnóstico). Sem o zero à esquerda: −,45 = −0,45.',
                fontsize=7.4, va='center', linespacing=1.3)
        y_fim = y - 2.15
    if subtipos:
        r = pd.read_csv(res / 'e8_subtipos_resumo.csv')
        r = r[(r['bloco'] == f'hbn_{suf}') & (r['feature'] == 'chi')].set_index('cid')
        ns = '; '.join(f"{c}: n = {int(r.loc[c, 'n_tdah'])}" for c in ('F90.0', 'F90.1', 'F90.2'))
        ax.text(X0, y_fim - 0.3, f"Apresentações (CID-10-CM, só HBN) — {ns}; controles: n = "
                f"{int(r['n_controle'].iloc[0])}.\nO adhdata e o TDBRAIN não codificam a apresentação.",
                fontsize=7.4, va='center', linespacing=1.3)
        y_fim -= 0.65
    ax.set(xlim=(X0 - 0.05, len(FEATS_T) * W + 0.05), ylim=(y_fim, 1.2))
    ax.axis('off')
    salvar(fig, sai, nome or f'F_tabuleiro_{cond}' + ('_num' if numeros else '') + ('_cid' if subtipos else ''))


def tabuleiros(res, sai, **_):
    with estilo(**{'font.size': 8}):
        for cond in ('OF', 'OA'):
            for numeros, subtipos in ((False, False), (True, False), (True, True)):
                tabuleiro(res, sai, cond, numeros, subtipos)
        tabuleiro(res, sai, 'OF', True, False, nome='F1_tabuleiro')


BLOCOS_M = [('adhdata_tarefa', 'adhdata\n(tarefa)'), ('hbn_fechado', 'HBN\n(OF)'), ('hbn_aberto', 'HBN\n(OA)'),
            ('tdbrain_fechado', 'TDBRAIN\n(OF)'), ('tdbrain_aberto', 'TDBRAIN\n(OA)')]
COR_BL = ['#7A7A7A', '#1F4E79', '#6FA3D2', '#8C2D19', '#D98C73']


def duas_metades(res, sai, v2=False, **_):
    des, aj = pd.read_csv(res / 'e10_desequilibrio.csv'), pd.read_csv(res / 'e9_ajuste_potencia.csv')
    if v2:
        e14 = pd.read_csv(res / 'e14_denominador_fixo.csv')
        e14 = e14[e14['feature'] == 'esp_dist'].rename(columns={'g_a': 'g', 'g_a_lo': 'g_lo', 'g_a_hi': 'g_hi'})
        aj = pd.concat([aj, e14[['bloco', 'feature', 'ajuste', 'g', 'g_lo', 'g_hi']]], ignore_index=True)
    vars_ = [('idade', 'Idade'), ('sexo_f', 'Sexo (F)'), ('log_dur', 'Duração'), ('emg_t', 'Músculo'),
             ('ocular', 'Olho'), ('log_var', 'Potência total')]
    feats = [('chi', 'Expoente χ'), ('offset', 'Offset'), ('alfa_pw', 'Potência do alfa')] + [('esp_dist', 'Distância espacial')] * v2
    ajs = [('C0', 'idade, sexo, duração (C0)', '#5B6770'), ('C0+pot', 'C0 + potência total', '#2F7FB8'),
           ('C0+pot+art', 'C0 + potência + músculo + olho', '#C8553D')]
    fig = plt.figure(figsize=(7.0, 6.4))
    gs = fig.add_gridspec(2, len(feats), height_ratios=[1.0, 1.25], hspace=0.55, wspace=0.12)
    ax = fig.add_subplot(gs[0, :])
    nb = len(BLOCOS_M)
    larg = 0.8 / nb
    for i, (bl, _) in enumerate(BLOCOS_M):
        s = des[des['bloco'] == bl].set_index('variavel')
        for j, (v, _) in enumerate(vars_):
            val = s.loc[v, 'smd'] if v in s.index else np.nan
            x = j - 0.4 + larg * (i + 0.5)
            if np.isfinite(val):
                ax.bar(x, val, width=larg * 0.9, color=COR_BL[i])
            else:
                ax.text(x, 0.02, 'n/d', ha='center', va='bottom', fontsize=6, color='#888888', rotation=90)
    ax.axhspan(-0.1, 0.1, color='#EEEEEE', zorder=0)
    ax.axhline(0, color='#444444', lw=0.6)
    ax.set_xticks(range(len(vars_)))
    ax.set_xticklabels([r for _, r in vars_])
    ax.set(ylabel='Diferença padronizada\n(TDAH − controle)', ylim=(-0.95, 0.75))
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in COR_BL], [r.replace('\n', ' ') for _, r in BLOCOS_M],
              ncol=5, fontsize=7.2, frameon=False, loc='lower left', bbox_to_anchor=(0.0, -0.36))
    ax.yaxis.set_major_formatter(VIRG)
    ax.set_title('A   O fator difere entre os grupos?', loc='left', fontsize=9.5, fontweight='bold')
    for k, (f, rotf) in enumerate(feats):
        ax = fig.add_subplot(gs[1, k])
        for i, (bl, _) in enumerate(BLOCOS_M):
            for m, (a, _, cor) in enumerate(ajs):
                s = aj[(aj['bloco'] == bl) & (aj['feature'] == f) & (aj['ajuste'] == a)]
                if not len(s):
                    continue
                y = i + (m - 1) * 0.22
                g, lo, hi = s[['g', 'g_lo', 'g_hi']].iloc[0]
                ax.plot([lo, hi], [y, y], color=cor, lw=1.1)
                ax.plot(g, y, 'o', color=cor, ms=3.6, mfc=cor if (lo > 0 or hi < 0) else 'white', mew=1.0)
        ax.axvline(0, color='#444444', lw=0.6)
        ax.set_yticks(range(nb))
        ax.set_yticklabels([r for _, r in BLOCOS_M] if k == 0 else [])
        ax.invert_yaxis()
        ax.set(xlim=(-1.3, 0.85), xlabel='g de Hedges' + ('\n' if v2 else ' ') + '(TDAH − controle)')
        if v2:
            ax.set_xticks([-1, -0.5, 0, 0.5])
        ax.xaxis.set_major_formatter(VIRG)
        ax.set_title(rotf, fontsize=9)
        if k == 0:
            ax.text(-0.02, 1.13, 'B   O efeito muda quando o fator é ajustado?', transform=ax.transAxes,
                    fontsize=9.5, fontweight='bold')
    fig.legend([plt.Line2D([], [], color=c, marker='o', ms=4, lw=1.1) for _, _, c in ajs], [r for _, r, _ in ajs],
               ncol=3, frameon=False, fontsize=7.4, loc='lower center', bbox_to_anchor=(0.5, -0.04 if v2 else -0.01))
    fig.text(0.5, -0.075 if v2 else -0.045, 'Ponto cheio: IC 95 % exclui zero; vazado: inclui. Faixa cinza em A: |d| < 0,1.',
             ha='center', fontsize=7.2, color='#444444')
    salvar(fig, sai, 'F6_duas_metades_v2' if v2 else 'F2_duas_metades')


def limpar_dados(res, sai, **_):
    lg = pd.read_csv(res / 'e6_tabuleiro_longo.csv')
    d = lg[(lg['tipo'] == 'delta_g') & lg['feature'].isin(['chi', 'chi_post', 'offset', 'iaf', 'alfa_pw', 'tbr'])]
    d = d[~((d['feature'] == 'tbr') & (d['coluna'] == 'U30'))]
    vs = [('rej_emg', 'Rejeição de épocas\ncom músculo', 'l'), ('rej_olho', 'Rejeição de épocas\ncom olho', 'l'),
          ('ica_olho', 'ICA ocular', 'l'), ('ica_musc', 'ICA muscular', 'l'), ('internos', '13 canais internos', 'l'),
          ('dur60', 'Só 60 s', 'l'), ('nativa', 'Referência nativa', 'p'), ('U30', 'Faixa do ajuste\n1 a 30 Hz', 'p')]
    cor = {'l': '#3B8C5A', 'p': '#E07B39'}
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    rng = np.random.default_rng(0)
    for i, (v, _, tipo) in enumerate(vs):
        x = d.loc[d['coluna'] == v, 'valor'].abs().values
        ax.scatter(x, i + rng.uniform(-0.18, 0.18, len(x)), s=9, color=cor[tipo], alpha=0.7, lw=0, zorder=3)
        ax.plot([np.median(x)] * 2, [i - 0.32, i + 0.32], color='black', lw=1.4, zorder=4)
    for xv in (0.10, 0.20):
        ax.axvline(xv, color='#888888', lw=0.7, ls='--', zorder=1)
    ax.set_yticks(range(len(vs)))
    ax.set_yticklabels([r for _, r, _ in vs], fontsize=7.5)
    ax.invert_yaxis()
    ax.set(xlim=(-0.01, 0.6), xlabel='|Δg| em relação à base')
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:.1f}'.replace('.', ',')))
    ax.legend([plt.Line2D([], [], marker='o', ls='', color=cor[k], ms=4) for k in 'lp'] + [plt.Line2D([], [], color='black', lw=1.4)],
              ['limpeza', 'processamento', 'mediana'], fontsize=7, frameon=False, loc='lower right')
    salvar(fig, sai, 'F8_limpar_processar_dados')


def potencia_esquema(res, sai, **_):
    f = np.logspace(0, np.log10(40), 300)
    ref = 10 ** (1.0 - 1.2 * np.log10(f))
    ganho, extra = ref * 3.2, ref * (1 + 6 * np.exp(-(f - 1) / 2.5))
    inc, intc = np.polyfit(np.log10(f), np.log10(extra), 1)
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    for y, c, kw in ((ref, '#8a8a8a', {}), (ganho, '#d9822b', {}), (extra, '#3f5f86', {}),
                     (10 ** (intc + inc * np.log10(f)), '#3f5f86', {'lw': 1, 'ls': '--'})):
        ax.loglog(f, y, color=c, **{'lw': 1.6, **kw})
    i, j, k = (np.argmin(abs(f - x)) for x in (30, 9, 2.2))
    ax.annotate('', xy=(f[i], ganho[i] * 0.93), xytext=(f[i], ref[i] * 1.07), arrowprops=dict(arrowstyle='-|>', color='#d9822b', lw=1.1))
    ax.text(f[j], ganho[j] * 1.45, 'ganho puro: a reta sobe inteira\n(offset muda, χ igual)', color='#d9822b', fontsize=7.5,
            ha='left', va='bottom', rotation=-17, rotation_mode='anchor')
    ax.text(1.05, extra[0] * 1.15, 'potência extra nas baixas frequências\n(olho, movimento, idade): χ muda', color='#3f5f86',
            fontsize=7.5, va='bottom')
    ax.text(f[k], ref[k] * 0.5, 'referência', color='#8a8a8a', fontsize=7.5, ha='left', va='top', rotation=-17, rotation_mode='anchor')
    ax.set(xlabel='Frequência (Hz)', ylabel='Potência (escala log)', xticks=[1, 2, 5, 10, 20, 40], yticks=[],
           ylim=(ref[-1] * 0.4, extra[0] * 4))
    ax.set_xticklabels(['1', '2', '5', '10', '20', '40'])
    ax.minorticks_off()
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    salvar(fig, sai, 'F7_potencia_esquema', tight=False)


CIN, AZC, VERDEC, LAR, ROXO, AMARC = '#f2f2f2', '#dde7f0', '#e4efe4', '#fbe9d7', '#e8e1f0', '#fff6d6'
TINTA = '#555555'


class Tela:
    def __init__(self, w, h):
        self.fig, self.ax = plt.subplots(figsize=(w, h))
        self.ax.set(xlim=(0, 100), ylim=(0, 100))
        self.ax.axis('off')
        self.fig.subplots_adjust(0, 0, 1, 1)

    def caixa(self, x, y, w, h, texto, cor=CIN, fs=6.8, titulo=None, ref=None, ha='center'):
        self.ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3,rounding_size=1.0', fc=cor, ec=TINTA, lw=0.7))
        tx = x + w / 2 if ha == 'center' else x + 1.5
        if titulo:
            self.ax.text(tx, y + h - 1.6, titulo, ha=ha, va='top', fontsize=fs + 0.5, fontweight='bold')
            self.ax.text(tx, y + (h - 3.2) / 2 - 0.2, texto, ha=ha, va='center', fontsize=fs, linespacing=1.25)
        else:
            self.ax.text(tx, y + h / 2, texto, ha=ha, va='center', fontsize=fs, linespacing=1.25)
        if ref:
            self.ax.text(x + w - 1.0, y + 0.9, ref, ha='right', va='bottom', fontsize=5.8, color='#333333', style='italic')

    def seta(self, x, y0, y1=None, x1=None):
        y1 = y0 - 2.3 if y1 is None else y1
        self.ax.add_patch(FancyArrowPatch((x, y0), (x if x1 is None else x1, y1), arrowstyle='-|>', mutation_scale=7,
                                          lw=0.7, color=TINTA, shrinkA=0, shrinkB=0))

    def linha(self, x0, y0, x1, y1):
        self.ax.plot([x0, x1], [y0, y1], color=TINTA, lw=0.7)

    def rotulo(self, y, texto):
        self.ax.text(1.0, y, texto, ha='left', va='center', fontsize=7.5, fontweight='bold', rotation=90, color='#333333')

    def nota(self, x, texto):
        self.ax.text(x, 12.5 if x == 50 else 13, texto, ha='center', va='center', fontsize=6.2, color=TINTA)

    def fechar(self, sai, nome, h):
        self.ax.set_ylim(10, 100)
        self.fig.set_size_inches(6.3, h * 0.9)
        salvar(self.fig, sai, nome, tight=False)


def fluxograma(res, sai, **_):
    with estilo(**{'font.size': 7.5}):
        t = Tela(6.3, 5.6)
        X0, XW = 9, 89
        t.rotulo(94, 'Registro')
        w = (XW - 3) / 3
        for k, (tt, s, r) in enumerate([
                ('adhdata', '19 canais, referência nas orelhas,\n128 Hz; tarefa de atenção visual', '[23]'),
                ('HBN', 'rede de 129 canais, referência Cz,\n500 Hz; repouso OF e OA', '[24]'),
                ('TDBRAIN', '26 canais, referência nas mastoides,\n500 Hz; repouso OF e OA', '[25]')]):
            x = X0 + k * (w + 1.5)
            t.caixa(x, 89, w, 10, s, AZC, fs=6.4, titulo=tt, ref=r)
            t.seta(x + w / 2, 88.7)
        t.rotulo(80.5, 'Pré-proc.')
        t.caixa(X0, 75, XW, 11,
                'Montagem nos 19 canais do sistema 10–20 (MNE-Python [32])  →  passa-altas 0,5 Hz  →  notch na rede e harmônicos\n'
                '→  reamostragem para 128 Hz  →  referência média comum (CAR); ordem das etapas conforme o PREP [36]\n'
                '→  épocas de 4 s sem sobreposição.  Variantes: referência nativa ou REST [38]; só 13 canais internos; só os primeiros 60 s;\n'
                'rejeição das épocas acima do percentil 90 em músculo ou em olho; ICA infomax estendido [40] + ICLabel [29, 30], olho removido',
                CIN, fs=6.3)
        t.seta(X0 + XW * 0.28, 74.7)
        t.seta(X0 + XW * 0.76, 74.7)
        t.rotulo(44, 'Características')
        we = XW * 0.55
        wd, xd0 = XW - we - 2, X0 + we + 2
        ramos = [(X0, we, [(62, 9.5, 'Welch por canal e época: janela de Hann, 0,25 Hz de resolução;\n'
                            'espectro do canal = mediana entre as épocas', 'Espectro', '[39]', 6.4),
                           (47, 12.5, 'modo fixed, 1–40 Hz (variante: 1–30 Hz), canal a canal, eq. (1):\n'
                            'expoente χ e offset = mediana dos 19 canais; χ posterior = mediana em P3, Pz, P4, O1, O2\n'
                            'pico alfa: novo ajuste no espectro mediano posterior; maior pico em 7–14 Hz → IAF e potência\n'
                            'acima da reta aperiódica', 'FOOOF 1.1.1 (specparam)', '[4]', 6.3),
                           (32, 12, 'TBR = log10 [P(4–8 Hz) / P(13–30 Hz)] em Fz, F3, F4, C3, C4, Cz, eq. (2)  [2, 41]\n'
                            'índice muscular = log10 P(30–45 Hz) em T7 e T8  [15, 16]\n'
                            'índice ocular = log10 [P(1–4 Hz) em Fp1/Fp2 ÷ P(1–4 Hz) em O1/O2]  [17]\n'
                            'potência total = log do traço da covariância sem normalização',
                            'TBR e índices de aquisição', None, 6.3)], VERDEC),
                 (xd0, wd, [(62, 9.5, 'banda 1–40 Hz, 15 épocas fixas de 4 s,\nmatriz 19 × 19 por participante',
                             'Covariância', None, 6.4),
                            (47, 12.5, 'C′ = 0,9 C/tr(C) + 0,1 I/19\nL = log(C′) (log-euclidiano [27])\nvetor = triângulo '
                             'superior de\nL − L̄ dos controles de treino', 'Espaço tangente', '[7]', 6.3),
                            (32, 12, 'distância = norma do vetor\n(para as correlações do tabuleiro)\nAUC = regressão '
                             'logística L2,\nvalidação cruzada 5 × 5 [35]', 'Resumo espacial', None, 6.3)], LAR)]
        for x, w, caixas, cor in ramos:
            for k, (y, h, txt, tit, ref, fs) in enumerate(caixas):
                t.caixa(x, y, w, h, txt, cor, fs=fs, titulo=tit, ref=ref)
                t.seta(x + w / 2, y - 0.3) if k < 2 else t.seta(x + w / 2, 31.7, 28.2)
        t.rotulo(23.5, 'Saída')
        t.caixa(X0, 19, XW, 9,
                'Uma linha por participante e por versão do processamento (base + 8 variantes): χ, χ posterior, offset, IAF, potência do alfa, TBR,\n'
                'distância espacial, índices muscular e ocular, potência total; mais idade, sexo, duração e grupo. Entrada das medidas da Figura 2.',
                ROXO, fs=6.4)
        t.nota(X0 + XW / 2, 'OF: olhos fechados; OA: olhos abertos. Os números entre colchetes remetem à lista de referências; '
               'as equações, à Seção 2.2.')
        t.fechar(sai, 'F0a_pipeline', 5.6)
        arvore(sai)


def arvore(sai):
    t = Tela(6.3, 6.0)
    t.caixa(2, 86, 96, 10, 'Para cada característica (χ, χ posterior, offset, IAF, potência do alfa, TBR, distância espacial)\n'
            'e cada bloco (adhdata tarefa; HBN OF e OA; TDBRAIN OF e OA)', AZC, fs=6.8, titulo='Unidade de análise')
    xs, wr = [2, 35, 68], 30
    t.linha(50, 85.7, 50, 83)
    t.linha(xs[0] + wr / 2, 83, xs[2] + wr / 2, 83)
    for x in xs:
        t.seta(x + wr / 2, 83, 80.3)
    colunas = [
        ('A. O fator acompanha a característica?', 'ρ de Spearman da característica\ncom cada fator, só nos controles',
         'Seis fatores', 'Amostra: idade · sexo · duração\nArtefato e aquisição: índice muscular ·\níndice ocular · potência total\n'
         '(adhdata: sem idade e sexo)', '|ρ| < 0,2 verde · 0,2–0,4 amarelo · ≥ 0,4 vermelho\n→ linhas Amostra e Artefato da Figura 3'),
        ('B. O fator difere entre os grupos?', 'TDAH − controle para cada fator', 'Duas medidas',
         'idade, sexo, duração: diferença\npadronizada d, eq. (3)\nmúsculo, olho, potência total: g de Hedges\n'
         'ajustado por C0, eqs. (4) e (5), IC bootstrap', '|d| < 0,1: faixa cinza (sem desequilíbrio)\n→ painel A da Figura 4')]
    for x, (tit, txt, tit2, txt2, leitura) in zip(xs, colunas):
        t.caixa(x, 70, wr, 10, txt, LAR, fs=6.5, titulo=tit)
        t.seta(x + wr / 2, 69.7)
        t.caixa(x, 53, wr, 14, txt2, AMARC, fs=6.4, titulo=tit2)
        t.seta(x + wr / 2, 52.7)
        t.caixa(x, 40, wr, 10, leitura, VERDEC, fs=6.3, titulo='Leitura')
    x2 = xs[2]
    t.caixa(x2, 70, wr, 10, 'g de Hedges ajustado, eqs. (4) e (5),\nrecalculado em cada variação', LAR, fs=6.5,
            titulo='C. O efeito de grupo muda?')
    t.linha(x2 + wr / 2, 69.7, x2 + wr / 2, 67.5)
    t.linha(x2 + 7, 67.5, x2 + wr - 7, 67.5)
    t.seta(x2 + 7, 67.5, 65.3)
    t.seta(x2 + wr - 7, 67.5, 65.3)
    wc = wr / 2 - 0.75
    xc2 = x2 + wc + 1.5
    t.caixa(x2, 46, wc, 19, 'base: CAR, 1–40 Hz\n  limpeza: rejeição músculo ·\n  rejeição olho · ICA ocular\n  canais: 13 internos\n'
            '  tempo: 60 s\n  referência: nativa\n  faixa: 1–30 Hz', AMARC, fs=6.0, titulo='C1. Processamento', ha='left')
    t.caixa(xc2, 46, wc, 19, 'C0: idade, idade²,\n  sexo, idade × sexo,\n  log da duração\nC0 + potência total\n'
            'C0 + potência total\n  + músculo + olho', AMARC, fs=6.0, titulo='C2. Covariáveis', ha='left')
    t.seta(x2 + wc / 2, 45.7)
    t.seta(xc2 + wc / 2, 45.7)
    t.caixa(x2, 33, wc, 10.5, 'Δg = g(variante) − g(base)\n|Δg| < 0,10 · 0,10–0,20 · ≥ 0,20\n(espacial: ΔAUC 0,03 · 0,07)\n'
            '→ linhas Processamento da Fig. 3', VERDEC, fs=5.9, titulo='Leitura')
    t.caixa(xc2, 33, wc, 10.5, 'g com cada conjunto; ponto cheio\nse o IC 95 % exclui zero\n→ painel B da Figura 4',
            VERDEC, fs=5.9, titulo='Leitura')
    for x, y in ((xs[0] + wr / 2, 39.7), (xs[1] + wr / 2, 39.7), (x2 + wc / 2, 32.7), (xc2 + wc / 2, 32.7)):
        t.linha(x, y, x, 30)
    t.linha(xs[0] + wr / 2, 30, xc2 + wc / 2, 30)
    t.seta(50, 30, 27.3)
    t.caixa(2, 17, 96, 10, 'Efeito de grupo na base (g com C0; espacial: AUC) na última linha da Figura 3.  '
            'Só é tratado como resultado o que se repete entre bancos, estados ou variantes;\n'
            'o que aparece num banco só é tratado como linha de ação. Sem correção para comparações múltiplas: o tabuleiro é um mapa descritivo.',
            ROXO, fs=6.3, titulo='Base diagnóstica e critério de leitura')
    t.nota(50, 'Base: referência média comum e ajuste em 1–40 Hz. ICA: componentes independentes rotulados pelo ICLabel. '
           'IC: intervalo de confiança por bootstrap (300 reamostragens).')
    t.fechar(sai, 'F0b_arvore', 6.0)


AZE, LARE, VERDEE, CZE = '#4a6a8a', '#d9822b', '#5a8a5a', '#777777'


def esquemas(res, sai, **_):
    with estilo(**{'font.size': 8}):
        rng = np.random.default_rng(3)
        fig, axs = plt.subplots(1, 6, figsize=(12, 3.2))
        t = np.linspace(0, 4, 2000)
        bruto = (np.cumsum(rng.normal(0, 1, t.size)) / 30 + 0.5 * np.sin(2 * np.pi * 10 * t)
                 + 0.25 * np.sin(2 * np.pi * 60 * t) + 0.3 * rng.normal(0, 1, t.size))
        filt = 0.5 * np.sin(2 * np.pi * 10 * t) + 0.25 * rng.normal(0, 1, t.size)
        etapas = ['EEG bruto\n(19 canais)', 'Filtro 0,5 Hz\n+ notch', 'Reamostragem\n128 Hz', 'Referência\nmédia (CAR)',
                  'Épocas de 4 s', 'Espectro (Welch)\n+ reta aperiódica']
        f = np.linspace(1, 40, 200)
        for k, ax in enumerate(axs):
            ax.set_title(etapas[k], fontsize=8.5)
            ax.axis('off')
            if k == 5:
                ax.loglog(f, 10 ** (1.2 - 1.3 * np.log10(f)) * (1 + 1.6 * np.exp(-((f - 10) / 1.5) ** 2)), color=AZE, lw=1.5)
                ax.loglog(f, 10 ** (1.2 - 1.3 * np.log10(f)), color=LARE, lw=1.2, ls='--')
                ax.text(10, 1.6, 'alfa', color=AZE, fontsize=8, ha='center')
                ax.text(22, 0.03, 'χ', color=LARE, fontsize=10)
                continue
            passo = 4 if k == 2 else 1
            for c in range(5):
                ax.plot(t[::passo], ((bruto if k == 0 else filt) * (0.8 + 0.1 * c) + 3 * c)[::passo],
                        lw=0.6 if k == 2 else 0.5, color=AZE)
            for x in ((1, 2, 3) if k == 4 else ()):
                ax.axvline(x, color=LARE, lw=1, ls='--')
            ax.annotate('', xy=(1.12, 0.5), xytext=(1.0, 0.5), xycoords='axes fraction',
                        arrowprops=dict(arrowstyle='-|>', color=CZE, lw=1))
        fig.tight_layout()
        salvar(fig, sai, 'F2_faixa_esq', tight=False)
        fig, axs = plt.subplots(1, 3, figsize=(10, 3.6))
        casos = [('HBN: fator igual nos dois grupos', 'músculo ≈ nos dois grupos\n→ ajuste não muda o efeito', False, False),
                 ('TDBRAIN: fator diferente entre os grupos', 'potência total menor no TDAH\n→ efeito do offset some ao ajustar',
                  True, False),
                 ('adhdata: fator faz parte do quadro', 'mais movimento no TDAH\n→ ajustar remove efeito verdadeiro', True, True)]
        for ax, (tit, leg, difere, parte) in zip(axs, casos):
            ax.set(xlim=(0, 10), ylim=(0, 10))
            ax.axis('off')
            ax.set_title(tit, fontsize=9)
            for x, y, txt, cor in [(2, 8, 'grupo\n(TDAH × controle)', AZE), (8, 8, 'característica\n(ex.: expoente χ)', VERDEE),
                                   (5, 3, 'fator\n(ex.: músculo)', LARE)]:
                ax.add_patch(FancyBboxPatch((x - 1.6, y - 0.9), 3.2, 1.8, boxstyle='round,pad=0.2', fc='white', ec=cor, lw=1.5))
                ax.text(x, y, txt, ha='center', va='center', fontsize=8)
            ax.add_patch(FancyArrowPatch((3.6, 8), (6.4, 8), arrowstyle='-|>', mutation_scale=12, color=CZE, lw=1.2))
            ax.text(5, 8.5, 'efeito de grupo', ha='center', fontsize=7.5, color=CZE)
            ax.add_patch(FancyArrowPatch((5.8, 3.9), (7.6, 7.0), arrowstyle='-|>', mutation_scale=12, color=LARE, lw=1.5))
            ax.text(7.4, 5.2, 'acompanha', fontsize=7.5, color=LARE)
            if difere:
                ax.add_patch(FancyArrowPatch((2.4, 7.0), (4.4, 3.9), arrowstyle='-|>' if parte else '-', mutation_scale=12,
                                             color=LARE, lw=1.5, ls='-' if parte else '--'))
                ax.text(1.6, 5.2, 'é parte do\nquadro' if parte else 'difere entre\nos grupos', fontsize=7.5, color=LARE,
                        ha='center')
            else:
                ax.text(2.8, 5.2, 'igual nos\ngrupos', fontsize=7.5, color=CZE, ha='center')
            ax.text(5, 0.6, leg, ha='center', fontsize=7.5)
        fig.tight_layout()
        salvar(fig, sai, 'F5_duas_metades_esq', tight=False)
        fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
        base = 10 ** (1.0 - 1.2 * np.log10(f))
        for ax in axs:
            ax.loglog(f, base, color=CZE, lw=1.5, label='referência')
            ax.set(xlabel='frequência (Hz)', ylabel='potência')
        b5 = base[np.argmin(abs(f - 5))]
        axs[0].loglog(f, base * 3, color=LARE, lw=1.8, label='ganho ×3')
        axs[0].set_title('Ganho puro: a reta sobe inteira\n→ muda o offset, o expoente χ não muda', fontsize=9)
        axs[0].annotate('', xy=(5, b5 * 3), xytext=(5, b5), arrowprops=dict(arrowstyle='-|>', color=LARE))
        axs[0].text(5.6, b5 * 1.6, 'offset', color=LARE)
        axs[1].loglog(f, base * (1 + 4 * np.exp(-(f - 1) / 3)), color=AZE, lw=1.8, label='potência extra em baixas frequências')
        axs[1].loglog(f, 10 ** (1.35 - 1.55 * np.log10(f)), color=AZE, lw=1, ls='--', label='nova reta (χ maior)')
        axs[1].set_title('Potência extra nas baixas frequências (olho, movimento, idade)\n→ a reta inclina, o expoente χ muda',
                         fontsize=9)
        for ax in axs:
            ax.legend(fontsize=7, loc='lower left')
        fig.tight_layout()
        salvar(fig, sai, 'F7_potencia_esq', tight=False)
        fig, ax = plt.subplots(figsize=(8, 3.6))
        itens = [('rejeição de épocas com músculo', 0.03, VERDEE), ('rejeição de épocas com olho', 0.04, VERDEE),
                 ('ICA ocular', 0.05, VERDEE), ('só 60 s', 0.04, VERDEE), ('faixa do ajuste 1–30 Hz', 0.17, LARE),
                 ('referência nativa', 0.33, LARE)]
        for i, (nome, dg, cor) in enumerate(itens):
            ax.barh(i, dg, color=cor, height=0.6)
            ax.text(dg + 0.01, i, nome, va='center', fontsize=8)
        for x, txt in ((0.10, '|Δg| = 0,10'), (0.20, '0,20')):
            ax.axvline(x, color=CZE, ls='--', lw=1)
            ax.text(x, 6.0, txt, ha='center', fontsize=7.5, color=CZE)
        ax.set_yticks([])
        ax.set(xlim=(0, 0.6), xlabel='mudança típica no efeito de grupo, |Δg| (ilustrativo)')
        ax.set_title('Limpar muda pouco (verde); escolher referência e faixa muda muito (laranja)', fontsize=9)
        for s in ('top', 'right', 'left'):
            ax.spines[s].set_visible(False)
        fig.tight_layout()
        salvar(fig, sai, 'F8_limpar_esq', tight=False)


E3 = {'font.size': 10, **SEM_BORDA, 'axes.titlesize': 10.5, 'axes.titleweight': 'bold'}
FIGS = {'mapa_confundidores': (mapa_confundidores, E3), 'pareados': (pareados, E3),
        'aquisicao_tdbrain': (aquisicao_tdbrain, E3), 'ponte_artigo2': (ponte_artigo2, E3),
        'tabuleiro': (tabuleiros, {}), 'duas_metades': (duas_metades, {'font.size': 8.5, **SEM_BORDA}),
        'duas_metades_v2': (lambda r, s, **k: duas_metades(r, s, v2=True), {'font.size': 8.5, **SEM_BORDA}),
        'limpar_dados': (limpar_dados, {'font.size': 8, **SEM_BORDA}), 'potencia_esquema': (potencia_esquema, {'font.size': 8}),
        'fluxograma': (fluxograma, {}), 'esquemas': (esquemas, {})}


def gerar(res=None, saida=None, so=None, braco=None):
    res, sai = Path(res or RESULTADOS), Path(saida or FIGURAS)
    if braco is None:
        try:
            from .config import ler
            braco = Path(ler()['saida_dir']) / 'agregados' / 'braco_A.csv'
        except Exception:
            braco = None
    for nome in so or FIGS:
        fn, rc = FIGS[nome]
        with estilo(**rc):
            fn(res, sai, braco=braco)
        print('ok', nome, flush=True)
