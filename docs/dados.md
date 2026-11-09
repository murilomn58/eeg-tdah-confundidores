# Dados: como obter e o que não é distribuído

## adhdata

IEEE DataPort, Nasrabadi et al. (2020), DOI [10.21227/rzfh-zn36](https://doi.org/10.21227/rzfh-zn36). Um CSV com os 19 canais em µV, a coluna `Class` (`ADHD` ou `Control`) e a coluna `ID`; 128 Hz; referência nas orelhas. Crianças de 7 a 12 anos em tarefa de atenção visual, sem idade e sexo por participante. Chave no `config.toml`: `adhdata_csv`.

## Healthy Brain Network (HBN)

- EEG: NEMAR, cópias dos conjuntos OpenNeuro ds005505 a ds005516, licença CC BY-SA 4.0. `python scripts/00_preparar.py --baixar-hbn participants_*.tsv` baixa só o RestingState e só os 19 canais 10-20 (a versão de 250 Hz do NEMAR), mais o `events.tsv`. Os blocos de olhos abertos e fechados são recortados pelos eventos `instructed_toOpenEyes` e `instructed_toCloseEyes`.
- Fenótipo: diagnóstico de consenso clínico e escalas (SWAN, SDQ, CGAS), sob o acordo de uso do Child Mind Institute. Chaves: `hbn_fenotipo_csv` (diagnósticos) e `hbn_sintomas_csv` (escalas).
- Entram crianças e adolescentes de 5 a 18 anos com TDAH no consenso ou sem diagnóstico.

## TDBRAIN

Brainclinics (van Dijk et al., 2022), mediante acordo de uso de dados. BIDS com `.bdf` de repouso de olhos fechados e abertos, 500 Hz, rede de 50 Hz, e a planilha `TDBRAIN_participants_V3.1.xlsx`. Entram a primeira sessão de cada pessoa e as indicações `ADHD` e `HEALTHY`; as análises usam os adultos. O lote de gravação (registros de 120 s ou de 121 s) é definido pela duração. Chaves: `tdbrain_dir`, `tdbrain_xlsx`.

## O que não é distribuído

Nada que identifique ou descreva um participante entra no repositório: o pacote pré-processado, as características e os fenótipos por participante (todos em `saida_dir`) e a tabela que liga o identificador original ao código do pacote (`mapa_ids_csv`, gerada por cada usuário). O `.gitignore` bloqueia esses arquivos.

## Células pequenas nas tabelas agregadas

Nenhuma tabela publicada tem contagem entre 1 e 4.

- `e0_contagens.csv`: dentro de cada banco e variável, as categorias com n entre 1 e 4 são somadas às próximas menores até formar uma linha `outros` com pelo menos 5 participantes e pelo menos duas categorias, de modo que nenhuma categoria pequena possa ser recuperada por diferença com o total do grupo; se isso exigir juntar todas as categorias, a variável é suprimida naquele banco.
- `e8_subtipos_*.csv`: apresentações com menos de 5 participantes são omitidas.
- `e16_fluxo.csv` e `e17_lote_tdbrain.csv`: linhas com n entre 1 e 4 são omitidas. No fluxo, essas linhas são falhas técnicas do sinal (registro sem pico alfa ou com covariância inválida, por exemplo), e algumas podem ser deduzidas pela diferença entre etapas consecutivas.
