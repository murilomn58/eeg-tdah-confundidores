"""Monta o artigo (versão curta ou longa) em .docx a partir de texto/<versao>/*.md, figuras e referências verificadas."""
import json
import re
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH as AL
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

from .config import FIGURAS, RAIZ, TEXTO

AVISO = RGBColor(0x99, 0x33, 0x00)
VERSOES = {
    'curto': dict(metodos='02_metodos.md', entrelinha=1.15, margens=(2.5, 2.5, 2, 2), ref_pt=10.5, ref_depois=4,
                  aviso='VERSÃO CURTA (08/11/2026): rascunho para revisão dos autores; números conferidos com os arquivos '
                        'de resultados; não submeter sem revisão'),
    'longo': dict(metodos='02_metodologia.md', entrelinha=1.5, margens=(3, 3, 2, 2), ref_pt=None, ref_depois=8,
                  aviso='RASCUNHO v3 (08/11/2026): texto para revisão dos autores; números conferidos com os arquivos de '
                        'resultados; não submeter sem revisão')}
DECLARACOES = ('**Uso de ferramentas de inteligência artificial.** Um modelo de linguagem (Claude, Anthropic) foi usado como '
               'apoio na triagem e na extração da revisão de escopo, com protocolo escrito e verificação por trecho literal, na '
               'escrita do código de análise, na conferência das referências e na redação. Os autores revisaram o conteúdo e '
               'assumem a responsabilidade por ele. **Ética.** Só foram usados bancos públicos e desidentificados, sob os termos '
               'de uso de cada um (HBN e TDBRAIN); cada banco tem a aprovação ética descrita pelos seus autores. **Dados e '
               'código.** Os bancos estão disponíveis nas fontes originais; os pré-registros, os hashes e o código de análise '
               'acompanham o projeto. [PENDENTE: financiamento e conflitos de interesse, a confirmar com a orientação.]')


class Citador:
    def __init__(self):
        self.ordem = []

    def trocar(self, texto):
        def sub(m):
            ns = sorted({self.ordem.index(k) + 1 if k in self.ordem else self.ordem.append(k) or len(self.ordem)
                         for k in (c.strip().lstrip('@') for c in m.group(1).split(';'))})
            partes, i = [], 0
            while i < len(ns):
                j = i
                while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1:
                    j += 1
                partes.append(f'{ns[i]}–{ns[j]}' if j - i >= 2 else ', '.join(map(str, ns[i:j + 1])))
                i = j + 1
            return '[' + ', '.join(partes) + ']'
        return re.sub(r'\[(@[^\]]+)\]', sub, texto)


def runs(p, texto, tam=None, italico=False):
    for parte in filter(None, re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*)', texto)):
        negrito, it = parte.startswith('**'), parte.startswith('*') and not parte.startswith('**')
        r = p.add_run(parte[2:-2] if negrito else parte[1:-1] if it else parte)
        r.bold, r.italic = negrito or None, it or italico or None
        if tam:
            r.font.size = Pt(tam)


def par(doc, texto='', alinhar=AL.JUSTIFY, recuo=None, antes=None, depois=None, linha=None, manter=False):
    p = doc.add_paragraph()
    p.alignment, f = alinhar, p.paragraph_format
    if recuo is not None:
        f.first_line_indent = Cm(recuo)
    for k, v in (('space_before', antes), ('space_after', depois)):
        if v is not None:
            setattr(f, k, Pt(v))
    if linha:
        f.line_spacing = linha
    f.keep_with_next = manter or None
    if texto:
        runs(p, texto)
    return p


def titulo(doc, texto, nivel):
    r = par(doc, alinhar=AL.LEFT, antes=18 if nivel == 1 else 12, depois=6, manter=True).add_run(
        texto.upper() if nivel == 1 else texto)
    r.bold, r.font.size = True, Pt(12)


def legenda(doc, rotulo, texto, cit):
    p = par(doc, alinhar=AL.CENTER, antes=10, depois=4, linha=1.0, manter=True)
    r = p.add_run(f'{rotulo} – ')
    r.bold, r.font.size = True, Pt(10)
    runs(p, cit.trocar(texto), 10)


def fonte(doc, texto, cit):
    p = par(doc, linha=1.0, depois=12)
    r = p.add_run('Fonte: ')
    r.bold, r.font.size = True, Pt(10)
    runs(p, cit.trocar(texto), 10)


def tabela(doc, t, cit):
    rot, _, resto = t['legenda'].partition(' – ')
    legenda(doc, rot, resto, cit)
    tab = doc.add_table(rows=len(t['linhas']), cols=len(t['linhas'][0]))
    tab.alignment = WD_TABLE_ALIGNMENT.CENTER
    b = OxmlElement('w:tblBorders')
    for lado, val in (('top', 'single'), ('bottom', 'single'), ('insideH', 'single'), ('left', 'nil'), ('right', 'nil'),
                      ('insideV', 'nil')):
        el = OxmlElement(f'w:{lado}')
        for k, v in (('w:val', val), ('w:sz', '6' if val == 'single' else '0'), ('w:color', '808080')):
            el.set(qn(k), v)
        b.append(el)
    tab._tbl.tblPr.append(b)
    for i, lin in enumerate(t['linhas']):
        for j, val in enumerate(lin):
            c = tab.cell(i, j)
            c.text = ''
            c.paragraphs[0].paragraph_format.line_spacing = 1.0
            r = c.paragraphs[0].add_run(cit.trocar(val))
            r.font.size, r.bold = Pt(9), i == 0 or None
            if i == 0:
                s = OxmlElement('w:shd')
                for k, v in (('w:val', 'clear'), ('w:color', 'auto'), ('w:fill', 'E8EEF4')):
                    s.set(qn(k), v)
                c._tc.get_or_add_tcPr().append(s)
    fonte(doc, t['fonte'], cit)


def ler_tabelas(pasta):
    out, cur = {}, None
    for linha in (pasta / 'tabelas.md').read_text(encoding='utf-8').splitlines():
        if linha.startswith('=== '):
            cur = out.setdefault(linha[4:].strip(), {'linhas': []})
        elif linha.strip() and cur is not None:
            chave = next((k for k in ('LEGENDA:', 'FONTE:') if linha.startswith(k)), None)
            if chave:
                cur[chave[:-1].lower()] = linha[len(chave):].strip()
            else:
                cur['linhas'].append([c.strip() for c in linha.split('|')])
    return out


def corpo(doc, md, cit, tabelas, figs, figdir, cont):
    for b in (x.strip() for x in md.split('\n\n')):
        if not b:
            continue
        if b.startswith('#'):
            nivel = len(b) - len(b.lstrip('#'))
            titulo(doc, b[nivel + 1:], nivel)
        elif b.startswith('{{FIG:'):
            cont['fig'] += 1
            f = figs[b[6:-2]]
            legenda(doc, f'Figura {cont["fig"]}', f['titulo'], cit)
            p = par(doc, alinhar=AL.CENTER, manter=True)
            caminho = Path(figdir) / f['arquivo']
            if caminho.exists():
                p.add_run().add_picture(str(caminho), width=Cm(f['largura_cm']))
            else:
                p.add_run(f'[FIGURA AUSENTE: {f["arquivo"]}]').bold = True
            fonte(doc, f['fonte'], cit)
        elif b.startswith('{{TAB:'):
            tabela(doc, tabelas[b[6:-2]], cit)
        elif b.startswith('- ') or re.match(r'^\d+\. ', b):
            for item in re.split(r'\n(?=- |\d+\. )', b):
                num = re.match(r'^(\d+)\. ', item)
                p = par(doc, depois=4)
                p.paragraph_format.left_indent, p.paragraph_format.first_line_indent = Cm(1.25), Cm(-0.5)
                runs(p, (f'{num.group(1)}. ' if num else '• ') + cit.trocar(re.sub(r'^(- |\d+\. )', '', item).replace('\n', ' ')))
        else:
            par(doc, cit.trocar(b.replace('\n', ' ')), recuo=1.25, depois=6)


def montar(versao='curto', saida=None, figdir=FIGURAS):
    v, pasta = VERSOES[versao], TEXTO / versao
    refs = json.loads((RAIZ / 'referencias' / 'referencias_verificadas.json').read_text(encoding='utf-8'))
    figs = json.loads((pasta / 'figuras.json').read_text(encoding='utf-8'))
    cit, doc = Citador(), Document()
    st = doc.styles['Normal']
    st.font.name, st.font.size = 'Times New Roman', Pt(12)
    st.element.rPr.rFonts.set(qn('w:eastAsia'), 'Times New Roman')
    st.paragraph_format.line_spacing, st.paragraph_format.space_after, st.paragraph_format.space_before = v['entrelinha'], Pt(0), Pt(0)
    s = doc.sections[0]
    s.page_height, s.page_width = Cm(29.7), Cm(21.0)
    s.top_margin, s.left_margin, s.bottom_margin, s.right_margin = (Cm(m) for m in v['margens'])
    r = par(doc, alinhar=AL.CENTER).add_run(v['aviso'])
    r.italic, r.font.size, r.font.color.rgb = True, Pt(9), AVISO
    cab = {k.strip(): val.strip() for k, _, val in (l.partition(':') for l in (pasta / '00_resumo.md').read_text(
        encoding='utf-8').splitlines()) if k.strip().isupper()}
    r = par(doc, alinhar=AL.CENTER, antes=12, depois=12).add_run(cab['TITULO'])
    r.bold, r.font.size = True, Pt(14)
    for linha in cab['AUTORES'].split('; ') + [cab['AFILIACAO']]:
        par(doc, alinhar=AL.CENTER, linha=1.0).add_run(linha).font.size = Pt(11)
    for rot, chave, antes in (('Resumo: ', 'RESUMO', 18), ('Palavras-chave: ', 'PALAVRAS-CHAVE', 8)):
        p = par(doc, antes=antes, linha=1.0)
        p.add_run(rot).bold = True
        runs(p, cit.trocar(cab[chave]))
    cont, tabelas = {'fig': 0}, ler_tabelas(pasta)
    for nome in ('01_introducao.md', v['metodos'], '03_resultados_discussao.md', '04_conclusao.md'):
        corpo(doc, (pasta / nome).read_text(encoding='utf-8'), cit, tabelas, figs, figdir, cont)
    titulo(doc, 'Declarações', 2)
    par(doc, DECLARACOES, recuo=1.25, depois=6)
    titulo(doc, 'Referências', 1)
    pendentes = []
    for i, k in enumerate(cit.ordem, 1):
        ref = refs.get(k) or {}
        p = par(doc, alinhar=AL.LEFT, linha=1.0, depois=v['ref_depois'])
        p.paragraph_format.left_indent, p.paragraph_format.first_line_indent = Cm(0.9), Cm(-0.9)
        r = p.add_run(f'[{i}] ')
        if v['ref_pt']:
            r.font.size = Pt(v['ref_pt'])
        if ref.get('abnt'):
            runs(p, ref['abnt'], v['ref_pt'])
        nota = None if ref.get('abnt') and ref.get('verificada') else (
            '  [referência não verificada: ' + (ref.get('nota') or 'conferir') + ']' if ref.get('abnt') else f'[REFERÊNCIA PENDENTE: {k}]')
        if nota:
            r = p.add_run(nota)
            r.italic, r.bold, r.font.color.rgb = bool(ref.get('abnt')) or None, not ref.get('abnt') or None, AVISO
            pendentes.append(k)
    saida = Path(saida or RAIZ / f'artigo1_{versao}.docx')
    doc.save(saida)
    return saida, len(cit.ordem), pendentes
