import re
import subprocess
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
TEXTO = {'.py', '.md', '.txt', '.csv', '.json', '.toml', '.cff', '.ipynb', '.gitignore', '.gitattributes', '.tsv', ''}
PADROES = {
    'usuario': re.compile('Narc[ií]sius', re.I), 'caminho_windows': re.compile(r'\b[A-Za-z]:[\\/]'),
    'ip_interno': re.compile(r'161\.24\.'), 'nome_servidor': re.compile('lesson[i]a', re.I), 'rede_interna': re.compile('CCG[E]'),
    'maquina': re.compile(r'\bservido[r]', re.I), 'codigo_participante': re.compile(r'\b[HTAN]\d{4}\b'),
    'id_bids': re.compile(r'sub-[A-Z0-9]'), 'id_ndar': re.compile('NDA[R]'),
    'chave': re.compile(r'\b(api_ke[y]|toke[n]|senh[a])', re.I)}
BINARIO = ('usuario', 'ip_interno', 'nome_servidor', 'rede_interna')


def arquivos():
    out = subprocess.run(['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'], cwd=RAIZ,
                         capture_output=True).stdout.decode('utf-8').split('\0')
    return [a for a in out if a and (RAIZ / a).is_file()]


def main():
    achados, lista = [], arquivos()
    for arq in lista:
        p = RAIZ / arq
        texto = p.suffix.lower() in TEXTO
        conteudo = p.read_text(encoding='utf-8', errors='replace') if texto else p.read_bytes().decode('latin-1')
        for nome in (PADROES if texto else BINARIO):
            n = len(PADROES[nome].findall(conteudo))
            if n:
                achados.append(f'{arq}  {nome} ({n}x)')
    print('\n'.join(f'PROBLEMA  {a}' for a in achados))
    print(f'{len(lista)} arquivos versionáveis; {len(achados)} problemas')
    sys.exit(1 if achados else 0)


if __name__ == '__main__':
    main()
