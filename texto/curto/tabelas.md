=== bancos
LEGENDA: Tabela 1 – Bancos de dados utilizados.
Banco | TDAH / controle | Idade | Condição | Referência nativa | Observação
adhdata | 61 / 60 | 7–12 anos | Tarefa de atenção visual | Orelhas | Sem idade e sexo por participante; metilfenidato no grupo TDAH
HBN | 1.354 / 179 | 5–21 anos | Repouso, olhos fechados e abertos | Cz (rede de 129 canais, reduzida a 19) | Controles da mesma coorte de encaminhamento, sem diagnóstico
TDBRAIN | 95 / 175 | 18 anos ou mais | Repouso, olhos fechados e abertos | Mastoides | Controles saudáveis; lote de gravação associado ao diagnóstico
FONTE: Nasrabadi et al. (2020), Alexander et al. (2017) e van Dijk et al. (2022); contagens após o pipeline.

=== implementacao
LEGENDA: Tabela 2 – Etapas do pipeline, ferramentas, parâmetros e referências.
Etapa | Ferramenta (versão) | Parâmetros | Referência
Leitura e montagem | MNE-Python 1.12.1 | 19 canais do sistema 10–20; no HBN, eletrodos da rede de 129 mapeados para o 10–20 | Gramfort et al. (2013)
Filtros e reamostragem | MNE-Python 1.12.1 | Passa-altas de 0,5 Hz; notch na rede e harmônicos; 128 Hz | Bigdely-Shamlo et al. (2015)
Referência | Implementação própria e MNE-Python | CAR (base); REST, CSD e nativa (variantes) | Yao (2001)
Espectro | MNE-Python (Welch) | Épocas de 4 s sem sobreposição, janela de Hann, resolução de 0,25 Hz; mediana entre épocas | Welch (1967)
Parametrização | FOOOF 1.1.1 (hoje specparam) | Modo fixed, 1–40 Hz; largura de pico de 1 a 8 Hz; até 6 picos; limiar de 2,0 | Donoghue et al. (2020)
ICA (variante) | MNE-Python (infomax estendido) + mne-icalabel 0.10.0 | 18 componentes, 1–63 Hz; remoção com probabilidade > 0,8 | Lee, Girolami e Sejnowski (1999); Pion-Tonachini et al. (2019)
Espaço tangente | NumPy 1.26.4 | 15 épocas de 4 s; traço normalizado; regularização de 0,1; métrica log-euclidiana | Barachant et al. (2012); Arsigny et al. (2007)
Classificação | scikit-learn 1.9.0 | Regressão logística L2 (C = 1; 0,05 com mais de 20 variáveis); validação cruzada 5 × 5 por participante | Pedregosa et al. (2011)
Estatística | NumPy 1.26.4 e SciPy 1.17.1 | g de Hedges ajustado; bootstrap estratificado (300); Spearman | Hedges (1981); Virtanen et al. (2020)
FONTE: elaborada pelos autores; parâmetros conferidos no código.
