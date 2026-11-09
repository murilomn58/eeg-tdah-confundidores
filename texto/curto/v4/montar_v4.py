import copy
import json
import re
import sys
from pathlib import Path

import docx
from docx.oxml.ns import qn
from docx.shared import Cm, Pt
from docx.text.paragraph import Paragraph
from docx.table import Table
from PIL import Image, ImageChops

RAIZ = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
AQUI = Path(__file__).resolve().parent
FIGS = Path(sys.argv[2]) if len(sys.argv) > 2 else AQUI.parents[2] / 'figuras'
ORIG = RAIZ / 'v3_Artigo2_curto_confundidores_2026-11-08.docx'
SAIDA = RAIZ / 'v4_Artigo2_curto_confundidores_2026-11-08.docx'
PREENCHER = json.loads((AQUI / 'preencher.json').read_text(encoding='utf-8')) if (AQUI / 'preencher.json').exists() else {}

d = docx.Document(ORIG)
P = list(d.paragraphs)
CORPO = P[22]


def blocos(texto):
    out, atual = [], None
    for linha in texto.splitlines():
        if linha.startswith('@'):
            if atual:
                out.append(atual)
            atual = [linha.strip(), []]
        elif atual is not None:
            atual[1].append(linha)
    if atual:
        out.append(atual)
    return [(c, '\n'.join(l).strip()) for c, l in out]


def preencher(t):
    for k, v in PREENCHER.items():
        t = t.replace('{' + k + '}', v)
    return t


def segmentos(md):
    toks = re.split(r'(\*\*|\*|_\{|\^\{|\})', md)
    b = i = False
    va = []
    out = []
    for tk in toks:
        if tk == '**':
            b = not b
        elif tk == '*':
            i = not i
        elif tk == '_{':
            va.append('subscript')
        elif tk == '^{':
            va.append('superscript')
        elif tk == '}' and va:
            va.pop()
        elif tk:
            out.append((tk, b, i, va[-1] if va else None))
    return out


def rpr_base(p):
    for r in p.runs:
        if r.text.strip() and not r.bold:
            return copy.deepcopy(r._r.rPr) if r._r.rPr is not None else None
    r = p.runs[0] if p.runs else None
    return copy.deepcopy(r._r.rPr) if r is not None and r._r.rPr is not None else None


def escrever(p, md, rpr=None):
    rpr = rpr if rpr is not None else rpr_base(p)
    for ch in list(p._p):
        if not ch.tag.endswith('}pPr'):
            p._p.remove(ch)
    for txt, b, i, va in segmentos(preencher(md)):
        r = p.add_run(txt)
        if rpr is not None:
            r._r.insert(0, copy.deepcopy(rpr))
        r.bold = True if b else None
        r.italic = True if i else None
        if va == 'subscript':
            r.font.subscript = True
        elif va == 'superscript':
            r.font.superscript = True
        if not b:
            r.bold = False if r.bold is None and rpr is not None and rpr.find(qn('w:b')) is not None else r.bold
        if not i:
            r.italic = False if rpr is not None and rpr.find(qn('w:i')) is not None else None


def novo_paragrafo(depois_de, modelo=CORPO):
    el = copy.deepcopy(modelo._p)
    for ch in list(el):
        if not ch.tag.endswith('}pPr'):
            el.remove(ch)
    depois_de.addnext(el)
    return Paragraph(el, modelo._parent)


def trocar(p, a, b):
    txt = p.text
    assert txt.count(a) == 1, (a[:50], txt[:60])
    ini = txt.find(a)
    pos, k0, k1 = 0, None, None
    for k, r in enumerate(p.runs):
        n = len(r.text)
        if k0 is None and pos + n > ini:
            k0 = k
        if pos + n >= ini + len(a):
            k1 = k
            break
        pos += n
    rs = p.runs
    junto = ''.join(r.text for r in rs[k0:k1 + 1])
    for r in rs[k0 + 1:k1 + 1]:
        r._r.getparent().remove(r._r)
    rs[k0].text = junto.replace(a, b)
    assert b in p.text


def autocrop(arq, pad=24):
    im = Image.open(FIGS / arq).convert('RGB')
    diff = ImageChops.difference(im, Image.new('RGB', im.size, (255, 255, 255))).convert('L').point(lambda v: 255 if v > 18 else 0)
    bb = diff.getbbox()
    if bb:
        im = im.crop((max(0, bb[0] - pad), max(0, bb[1] - pad), min(im.width, bb[2] + pad), min(im.height, bb[3] + pad)))
    dst = AQUI / ('crop_' + arq)
    im.save(dst)
    return dst


def trocar_imagem(p, arq, largura_cm, cortar=True):
    for r in list(p.runs):
        r._r.getparent().remove(r._r)
    p.add_run().add_picture(str(autocrop(arq) if cortar else FIGS / arq), width=Cm(largura_cm))


ancora = {}
for cab, corpo in blocos((AQUI / 'v4_texto.md').read_text(encoding='utf-8')):
    m = re.match(r'@(p|after p|sub p|del p)(\d+)', cab)
    tipo, n = m.group(1), int(m.group(2))
    if tipo == 'p':
        escrever(P[n], corpo)
    elif tipo == 'after p':
        base = ancora.get(n, P[n]._p)
        q = novo_paragrafo(base)
        escrever(q, corpo, rpr_base(CORPO))
        ancora[n] = q._p
    elif tipo == 'sub p':
        a, b = [x.strip() for x in corpo.split('\n=>\n')]
        trocar(P[n], a, b)
    elif tipo == 'del p':
        P[n]._p.getparent().remove(P[n]._p)

trocar_imagem(P[15], 'F1_pipeline_v4_ia.png', 14.5)
trocar_imagem(P[34], 'F3_arvore_v4_ia.png', 16.0)
trocar_imagem(P[64], 'F5_duas_metades_v4c.png', 16.0)
trocar_imagem(P[70], 'F6_duas_metades_v2.png', 15.0, cortar=False)

t1, t2 = d.tables[0], d.tables[1]
for row in t1.rows:
    if row.cells[0].text.strip() == 'adhdata':
        c = row.cells[-1]
        c.paragraphs[0].runs[-1].text = c.paragraphs[0].runs[-1].text.rstrip('.') + '; C0 só com a duração'
for row in t2.rows:
    nome = row.cells[0].text.strip()
    if nome == 'Referência':
        escrever(row.cells[2].paragraphs[0], 'CAR (base); referência nativa de cada banco (variante)')
        escrever(row.cells[3].paragraphs[0], 'Yao et al. (2019)')
    if nome == 'Espaço tangente':
        escrever(row.cells[2].paragraphs[0], '18 canais, sem o Cz; 15 primeiras épocas de 4 s; traço normalizado; regularização de 0,1; métrica log-euclidiana; 171 variáveis')
    if nome == 'Classificação':
        escrever(row.cells[2].paragraphs[0], 'Regressão logística L2 (C = 0,05); validação cruzada 5 × 5 por participante, com as mesmas dobras em todas as variantes')


def sem_bordas(tbl):
    tblPr = tbl._tbl.tblPr
    for b in tblPr.findall(qn('w:tblBorders')):
        tblPr.remove(b)
    bordas = tblPr.makeelement(qn('w:tblBorders'), {})
    for lado in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
        e = bordas.makeelement(qn(f'w:{lado}'), {qn('w:val'): 'nil'})
        bordas.append(e)
    tblPr.append(bordas)


def celula_com(cel, modelo, md=None, imagem=None, largura=7.6):
    alvo = cel.paragraphs[0]
    novo = copy.deepcopy(modelo._p)
    alvo._p.addnext(novo)
    alvo._p.getparent().remove(alvo._p)
    q = Paragraph(novo, cel)
    if imagem:
        for r in list(q.runs):
            r._r.getparent().remove(r._r)
        q.add_run().add_picture(str(FIGS / imagem), width=Cm(largura))
    else:
        escrever(q, md)


leg, img, fon = P[74], P[75], P[76]
tab = d.add_table(rows=3, cols=2)
sem_bordas(tab)
for tr in tab._tbl.findall(qn('w:tr')):
    trPr = tr.get_or_add_trPr()
    trPr.append(trPr.makeelement(qn('w:cantSplit'), {}))
ancora.get(77, P[77]._p).addnext(tab._tbl)
celula_com(tab.cell(0, 0), leg, '**Figura 7 – **Por que a potência total afeta o offset e o expoente de formas diferentes')
celula_com(tab.cell(0, 1), leg, '**Figura 8 – **Quanto cada variante de processamento muda o efeito de grupo')
celula_com(tab.cell(1, 0), img, imagem='F7_potencia_esquema.png')
celula_com(tab.cell(1, 1), img, imagem='F8_limpar_processar_dados.png')
celula_com(tab.cell(2, 0), fon, '**Fonte: **Elaborada pelos autores; esquema sem dados. Um ganho que multiplica o espectro inteiro desloca a reta aperiódica sem mudar a inclinação: só o offset da eq. (1) muda. Potência adicional nas baixas frequências inclina a reta e muda o expoente χ.')
celula_com(tab.cell(2, 1), fon, '**Fonte: **Elaborada pelos autores. Cada ponto é uma casa (bloco × característica espectral) e mostra |Δg| em relação à base; traço preto: mediana; linhas tracejadas: 0,10 e 0,20. Verde: limpeza; laranja: escolhas de processamento.')
for lin in (0, 1):
    for col in (0, 1):
        for par in tab.cell(lin, col).paragraphs:
            par.paragraph_format.keep_with_next = True
for p in (P[74], P[75], P[76], P[80], P[81], P[82]):
    p._p.getparent().remove(p._p)
vazio = novo_paragrafo(tab._tbl)
escrever(vazio, '')

ESCORE = [
    ('Característica', 'Amostra', 'Artefato e aquisição', 'Processamento', 'Total OF', 'Total OA'),
    ('Offset', '0,43', '0,56', '0,08', '0,25', '0,21'),
    ('Expoente χ', '0,36', '0,39', '0,10', '0,21', '0,21'),
    ('Distância espacial', '0,00', '0,44', '0,10', '0,18', '0,23'),
    ('Expoente posterior', '0,29', '0,33', '0,06', '0,16', '0,20'),
    ('TBR', '0,29', '0,28', '0,07', '0,16', '0,17'),
    ('Potência do alfa', '0,14', '0,17', '0,06', '0,10', '0,05'),
    ('IAF', '0,14', '0,06', '0,04', '0,06', '0,20'),
]
leg_t2 = next(p for p in d.paragraphs if p.text.startswith('Tabela 2 –'))
fon_t2 = next(p for p in d.paragraphs if p.text.startswith('Fonte: elaborada pelos autores; parâmetros conferidos'))
ref = next(p for p in d.paragraphs if p.text.startswith('A Tabela 3 condensa o tabuleiro'))
leg3 = novo_paragrafo(ref._p, leg_t2)
escrever(leg3, '**Tabela 3 – **Escore de cuidado por característica (0 a 1): grupos de linhas do tabuleiro de olhos fechados (três bancos) e total de olhos abertos (HBN e TDBRAIN).', rpr_base(leg_t2))
novo_tbl = copy.deepcopy(t2._tbl)
leg3._p.addnext(novo_tbl)
t3 = Table(novo_tbl, leg3._parent)
for tr in list(novo_tbl.findall(qn('w:tr')))[len(ESCORE):]:
    novo_tbl.remove(tr)
grid = novo_tbl.find(qn('w:tblGrid'))
cols = grid.findall(qn('w:gridCol'))
while len(cols) < 6:
    grid.append(copy.deepcopy(cols[-1]))
    cols = grid.findall(qn('w:gridCol'))
for tr in novo_tbl.findall(qn('w:tr')):
    tcs = tr.findall(qn('w:tc'))
    while len(tcs) < 6:
        tr.append(copy.deepcopy(tcs[-1]))
        tcs = tr.findall(qn('w:tc'))
larg = [4.0, 2.2, 3.0, 2.6, 2.2, 2.2]
for col, w in zip(cols, larg):
    col.set(qn('w:w'), str(int(Cm(w).twips)))
for i, linha in enumerate(ESCORE):
    for j, txt in enumerate(linha):
        cel = t3.cell(i, j)
        cel.width = Cm(larg[j])
        par = cel.paragraphs[0]
        for extra in cel.paragraphs[1:]:
            extra._p.getparent().remove(extra._p)
        md = f'**{txt}**' if i == 0 else txt
        escrever(par, md)
fon3 = Paragraph(copy.deepcopy(fon_t2._p), fon_t2._parent)
novo_tbl.addnext(fon3._p)
escrever(fon3, '**Fonte: **elaborada pelos autores a partir do tabuleiro (Figura 4) e da sua versão de olhos abertos; casa amarela = 1, vermelha = 2, soma dividida pelo máximo possível em cada grupo de linhas.', rpr_base(fon_t2))

d.save(SAIDA)
pend = [k for k in ('DOI', 'N_ARQ', 'N_LOC', 'VALIDACAO') if k not in PREENCHER]
print('ok', SAIDA.name, 'pendentes:', pend)
