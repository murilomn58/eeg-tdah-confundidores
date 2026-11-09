# Formato de entrada (o pacote)

Todas as análises partem de um pacote com os 19 canais do sistema 10-20 a 128 Hz. `scripts/00_preparar.py` gera o pacote a partir dos registros brutos; quem já tiver os dados nesse formato pode começar direto em `scripts/01_caracteristicas.py`.

## Pastas

```
<pacote_dir>/
  <bloco>/            um bloco por banco e condição (por exemplo meu_banco, ou hbn_fechado e hbn_aberto)
    <codigo>.npz      um arquivo por participante
    metadados.csv     uma linha por participante
```

## Cada `<codigo>.npz`

| Chave | Tipo | Conteúdo |
|---|---|---|
| `x` | int16, (19, n) | sinal quantizado, canais na ordem Fp1, Fp2, F7, F3, Fz, F4, F8, T7, C3, Cz, C4, T8, P7, P3, Pz, P4, P8, O1, O2 |
| `escala` | float, (19,) | passo de quantização de cada canal, em µV (`max(|x|) / 32000`) |
| `segs` | int, (k, 2) | trechos contínuos `[início, fim)`, em amostras a 128 Hz |

O sinal em µV é `x * escala[:, None]`. As épocas de 4 s nunca atravessam a fronteira entre trechos de `segs`; um registro contínuo tem um trecho só, `[(0, n)]`. O sinal fica na referência original do registro; as análises trocam a referência depois (média comum na versão base).

## `metadados.csv`

| Coluna | Conteúdo |
|---|---|
| `codigo` | nome do arquivo `.npz`, sem extensão |
| `grupo` | `TDAH` ou `controle` |
| `idade` | anos (pode ficar vazio) |
| `sexo` | `F` ou `M` (pode ficar vazio) |
| `duracao_s` | duração do registro em segundos |

## Pré-processamento aplicado por `00_preparar.py`

1. Seleção dos 19 canais pelo nome (no HBN, pelo mapa da rede EGI em `eegtdah/mapa_egi_1020.json`).
2. Passa-altas de 0,5 Hz no sinal contínuo.
3. *Notch* de 1 Hz de largura na rede elétrica e em até três harmônicos abaixo de 90 % de Nyquist. A rede é informada (`meu_banco_rede`) ou detectada no espectro: um pico em 50 ou 60 Hz precisa ficar 3 dB acima da vizinhança e vencer a outra candidata por 3 dB.
4. Reamostragem para 128 Hz.
5. Quantização em int16 com escala por canal.
