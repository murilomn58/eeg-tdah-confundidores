# eeg-tdah-confundidores

Este código mede, para cada característica de EEG usada na busca de marcadores do TDAH (expoente e offset aperiódicos, potência e frequência do pico alfa, razão teta/beta e uma distância no espaço tangente da covariância), quanto ela depende de fatores alheios ao diagnóstico: idade, sexo, duração do registro, músculo, olho, potência total e escolhas de processamento. O resultado é um **tabuleiro de variáveis** por banco e um **escore de cuidado** por característica. Ele reproduz o artigo em três bancos públicos e pode ser aplicado a qualquer conjunto de registros com os 19 canais do sistema 10-20.

## Instalação

Python 3.11:

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt      # no Windows: .venv\Scripts\pip
python -m pytest tests                         # 4 testes com sinais sintéticos
```

Copie `config.example.toml` para `config.toml` e ajuste os caminhos. A pasta `saida_dir` deve ficar fora do repositório: é nela que ficam os arquivos por participante.

## Reproduzir o artigo

| Banco | Onde obter |
|---|---|
| adhdata (crianças, tarefa) | IEEE DataPort, DOI [10.21227/rzfh-zn36](https://doi.org/10.21227/rzfh-zn36) |
| Healthy Brain Network (crianças e adolescentes, repouso) | EEG no NEMAR/OpenNeuro (ds005505 a ds005516); fenótipo sob o acordo de uso do Child Mind Institute |
| TDBRAIN (adultos, repouso) | Brainclinics, sob acordo de uso de dados |

Detalhes em [`docs/dados.md`](docs/dados.md). Com os caminhos no `config.toml`, a ordem é:

```bash
python scripts/00_preparar.py                  # registros brutos -> pacote de 19 canais a 128 Hz
python scripts/01_caracteristicas.py --jobs 8  # características por participante, base e variantes
python scripts/02_ensaios.py --jobs 8          # tabelas agregadas -> resultados_agregados/
python scripts/02_ensaios.py --replicacoes     # replicações de alegações publicadas (demorado)
python scripts/03_figuras.py                   # figuras -> figuras/
python scripts/04_documento.py curto           # artigo em .docx
```

Para baixar o HBN do NEMAR, use `python scripts/00_preparar.py --baixar-hbn participants_*.tsv` (precisa também de `requests` e `numcodecs`).

## Usar com os seus dados

**Entrada mínima.** Uma pasta com os registros (qualquer formato que o MNE-Python leia: EDF, BDF, BrainVision, EEGLAB, FIF e outros) e um CSV com uma linha por registro:

```
arquivo,grupo,idade,sexo
p001.edf,TDAH,9.5,M
p002.edf,controle,10.1,F
```

`grupo` é `TDAH` ou `controle`; `idade` e `sexo` (`F` ou `M`) podem ficar vazios, e então as covariáveis usam só a duração. Os 19 canais são encontrados pelo nome (Fp1, Fp2, F7, F3, Fz, F4, F8, T7, C3, Cz, C4, T8, P7, P3, Pz, P4, P8, O1, O2; prefixos como `EEG` e sufixos como `-REF` são ignorados, e T3, T4, T5, T6 viram T7, T8, P7, P8). Se os nomes forem outros, informe um JSON `{"nome no arquivo": "nome 10-20"}`.

**Configuração.** No `config.toml`:

```toml
pacote_dir = "/caminho/para/pacote"
saida_dir = "/caminho/para/saida"
meu_banco_csv = "/caminho/para/meus_registros/lista.csv"
meu_banco_bloco = "meu_banco"
meu_banco_rede = 60          # frequência da rede elétrica; sem esta linha ela é detectada no espectro
meu_banco_mapa = ""          # opcional: JSON com os nomes dos canais
mapa_ids_csv = "/caminho/para/saida/mapa_ids.csv"
```

**Execução e o que cada passo gera.**

```bash
python scripts/00_preparar.py --meu-banco
python scripts/01_caracteristicas.py --blocos meu_banco --jobs 8
python scripts/02_ensaios.py --blocos meu_banco --partes tabuleiro potencia desequilibrio denominador loo tabela --destino meus_resultados
```

1. `00_preparar.py --meu-banco`: lê cada registro, aplica passa-altas de 0,5 Hz e *notch* na rede, reamostra para 128 Hz e grava `pacote_dir/meu_banco/<codigo>.npz` e `metadados.csv` (formato em [`docs/pacote.md`](docs/pacote.md)). Quem já tem os dados nesse formato pode pular este passo.
2. `01_caracteristicas.py`: grava `saida_dir/caracteristicas/meu_banco.npz` com as características de cada participante na versão base (referência média comum) e nas variantes de processamento. `--variantes base rej_emg internos dur60 nativa U30` deixa de fora a ICA, que é a parte mais lenta.
3. `02_ensaios.py`: grava as tabelas agregadas, entre elas `e6_tabuleiro_longo.csv` e `e6_g_auc_variantes.csv` (o tabuleiro), `e9_ajuste_potencia.csv` e `e10_desequilibrio.csv` (as duas metades de um confundidor), `e18_tabuleiro_numerico.csv` e `e20_escore_cuidado.csv` (o escore).

**Como ler o tabuleiro.** Cada característica tem três tipos de linha. *Amostra e artefato*: correlação de Spearman (ρ) com idade, sexo, duração, músculo, olho e potência total, só nos controles; |ρ| < 0,2 é pouco sensível, 0,2 a 0,4 moderado, 0,4 ou mais muito sensível. *Processamento*: mudança no g de Hedges do efeito de grupo quando a versão base é trocada por uma variante (rejeição de épocas, ICA, só canais internos, só 60 s, referência nativa, ajuste de 1 a 30 Hz); |Δg| < 0,10, 0,10 a 0,20, 0,20 ou mais (na família espacial, mudança na AUC: 0,03 e 0,07). *Sinal*: o g de Hedges ajustado por idade, idade², sexo, idade × sexo e duração. Um fator só confunde o efeito se acompanhar a característica **e** diferir entre os grupos; a segunda metade está em `e10_desequilibrio.csv`.

**Como ler o escore.** Em `e20_escore_cuidado.csv`, cada casa moderada vale 1 e cada casa muito sensível vale 2; a soma é dividida pelo máximo possível. O escore vai de 0 (nada preocupa) a 1 (tudo preocupa) e é dado por grupo de linhas (amostra, artefato e aquisição, processamento) e no total.

## Tabelas agregadas e caderno

`resultados_agregados/` tem as tabelas do artigo, sem nenhum dado por participante, e `figuras/` as figuras geradas a partir delas. Sem baixar nenhum banco, dá para:

```bash
python scripts/03_figuras.py                   # refaz as figuras a partir de resultados_agregados/
jupyter notebook notebooks/visao_geral.ipynb   # tabuleiro, duas metades e escore de cuidado
```

Três ilustrações da versão curta do artigo (pipeline, árvore do método e duas metades) foram feitas com um modelo generativo de imagens a partir dos esquemas de `figuras.py` e não são reproduzíveis por código; `scripts/figuras_ilustradas/corrigir_fig5.py` redesenha nelas a faixa com os valores de g.

A validação do código contra os resultados do artigo está resumida em [`docs/validacao.md`](docs/validacao.md), e as impressões digitais dos arquivos de código em [`HASHES.md`](HASHES.md). Antes de publicar uma cópia com dados próprios, `python scripts/verificar_publicavel.py` procura identificadores e caminhos locais nos arquivos versionáveis.

## Como citar

Use os metadados de [`CITATION.cff`](CITATION.cff). DOI: a definir.

## Licença

MIT, ver [`LICENSE`](LICENSE).
