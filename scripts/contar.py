import subprocess
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
arquivos = subprocess.run(['git', 'ls-files', '--cached', '--others', '--exclude-standard', '*.py'], cwd=RAIZ,
                          capture_output=True, text=True, encoding='utf-8').stdout.split()
linhas = {a: sum(1 for l in (RAIZ / a).read_text(encoding='utf-8').splitlines() if l.strip()) for a in sorted(arquivos)}
for a, n in linhas.items():
    print(f'{n:6d}  {a}')
print(f'{sum(linhas.values()):6d}  total em {len(linhas)} arquivos .py (linhas não brancas)')
