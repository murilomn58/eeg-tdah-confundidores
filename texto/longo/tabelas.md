=== bancos
LEGENDA: Tabela 1 – Bancos de dados utilizados. Os três têm os 19 canais do sistema 10–20 e diferem na idade, na condição de registro e na referência nativa, que o pipeline harmoniza.
Banco | TDAH / controle | Idade | Condição | Referência nativa | Informação clínica disponível
adhdata | 61 / 60 | 7–12 anos | Tarefa de atenção visual | Orelhas | Sem idade e sexo por participante; metilfenidato por até 6 meses no grupo TDAH
HBN | 1.354 / 179 | 5–21 anos | Repouso, olhos fechados e abertos | Cz (rede de 129 canais) | Consenso clínico (DSM-5) com códigos CID-10-CM; controles da mesma coorte sem diagnóstico; apresentações; comorbidades; SDQ; exame toxicológico parcial
TDBRAIN | 95 / 175 (adultos) | 18–83 anos | Repouso, olhos fechados e abertos | Mastoides | Indicação clínica; subestudo de origem; medicação não registrada
FONTE: Nasrabadi et al. (2020), Alexander et al. (2017) e van Dijk et al. (2022); contagens após o pipeline.

=== familias
LEGENDA: Tabela 2 – Famílias de características avaliadas e confundidores candidatos apontados pela literatura.
Família | Características | O que mede | Confundidores candidatos
Espectral aperiódica | Expoente χ (specparam fixed, 1–40 Hz), média dos canais | Inclinação do fundo 1/f; proposto como índice do balanço excitação/inibição | Idade, medicação, tarefa × repouso, faixa do ajuste, referência, músculo
Espectral periódica | Altura do alfa e IAF (centro de gravidade) acima do fundo aperiódico | Oscilação alfa, maturação e estado de alerta | Estado (olhos fechados × abertos), idade, sexo, sonolência
Razão de bandas | TBR clássica frontocentral (4–8 / 13–30 Hz) | Mistura de atividade periódica e aperiódica | Idade, expoente, referência
Espacial e riemanniana | Espaço tangente log-euclidiano da covariância 1–40 Hz (18 canais); CSP, FBCSP, coerência | Organização espacial do sinal; conectividade funcional | Aquisição (amplificador, touca, impedância), montagem, músculo, olho, comorbidade
FONTE: elaborada pelos autores.

=== sintese
LEGENDA: Tabela 3 – Síntese por família: confundidores observados nos controles, o que sobra após o controle e a implicação para um classificador.
Família | Confundidor dominante nos dados | Após o controle | Implicação para o modelo
Expoente aperiódico | Idade (ρ de −0,24 a −0,45); potência total (ρ de 0,31 a 0,51); músculo e olho no TDBRAIN; tarefa (adhdata); referência e faixa do ajuste | ΔAUC = 0,022 sobre a base no HBN; no limite do nulo com pareamento (0,55); ausente em adultos, com sinal que muda com a referência | Entrar com idade, sexo e idade × sexo; fixar e relatar referência e faixa do ajuste; tratar a medicação como não mensurada
Offset aperiódico | Potência total (ρ de 0,62 a 0,76); idade (ρ ≈ −0,5); músculo no HBN; olho no TDBRAIN | IC exclui zero nos cinco blocos sem ajuste; no TDBRAIN de olhos fechados, some com o ajuste pela potência total | Tratar como medida de potência de banda larga; exigir covariáveis de aquisição antes de interpretar
Periódico (alfa e IAF) | Estado OF × OA; idade (IAF); sexo (potência do alfa); detecção do pico (IAF) | ΔAUC = 0,028 sobre base e expoente (reanálise com Freedman–Lane); AUC ≈ 0,58 pareado de olhos abertos; potência do alfa estável em todas as variantes no HBN (g de −0,22 a −0,32) | Modelar por condição de registro; interação com sexo; tratar como exploratório até confirmação
TBR | Idade (ρ de −0,43 a −0,51); mistura com o expoente (ρ até 0,61); referência; músculo e olho no TDBRAIN | Some nas crianças após correção pelo expoente; resta g = +0,32 só no TDBRAIN (corrigida pelo expoente); sem correção, g = +0,28 com IC incluindo zero, que sobe para +0,47 com a referência nativa | Não usar como marcador isolado; substituir pelas componentes separadas
Espacial e riemanniana | Aquisição (potência total, lote); músculo; olho; comorbidade (TEA, g = +0,27); a distância ao centro não acompanha a idade | Acaso no HBN; AUC 0,85 no TDBRAIN, que cai para 0,76 a 0,81 com ICA ocular ou sem as bordas, e a covariância prevê o lote (até 0,90) | Exigir aquisição comum aos grupos ou covariáveis de aquisição; modelar comorbidade; regularização própria
FONTE: elaborada pelos autores a partir dos resultados das Seções 3.2 a 3.8 e da Figura 6. Os ρ das Seções 3.1 a 3.6 usam as definições do pipeline principal (IAF por centro de gravidade, TBR sem o Cz) e não são idênticos aos da Figura 6, que usa o pico do FOOOF.
