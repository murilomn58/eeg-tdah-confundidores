@p0
*VERSÃO CURTA v3 (08/11/2026) – rascunho para revisão dos autores; números conferidos com os arquivos de resultados; não submeter sem revisão*

@p1
**Revisão das features de EEG consideradas para diagnóstico de TDAH: uma investigação de biomarcadores**

@p2
Murilo Miguel Narciso (murilo.narciso.102681@ga.ita.br)

@p3
Sarah Negreiros de Carvalho Leite (sarah.leite@gp.ita.br)

@p4
Instituto Tecnológico de Aeronáutica (ITA), Departamento de Engenharia Eletrônica e Computação (EEC-T)

@p5
**Resumo: **Quando se procura no eletroencefalograma (EEG) um marcador do Transtorno do Déficit de Atenção e Hiperatividade (TDAH), é comum corrigir todas as características do sinal pelas mesmas variáveis, em geral idade e sexo. Essa prática supõe que fatores alheios ao diagnóstico, como idade, sexo, artefatos musculares e oculares e equipamento de registro, afetam todas as características da mesma forma. Este trabalho testa essa suposição em três bancos públicos (adhdata, crianças em tarefa; Healthy Brain Network, crianças em repouso; TDBRAIN, adultos em repouso), com um pipeline harmonizado e o mesmo critério estatístico. Para o expoente e o offset aperiódicos, a potência do pico alfa, a razão teta/beta e uma representação espacial-riemanniana, mediram-se três coisas: a associação com idade, sexo, duração, artefatos musculares e oculares e potência do registro nos controles; o desequilíbrio desses fatores entre os grupos; e a mudança do efeito de grupo sob oito variantes de processamento e sob ajuste pela potência e pelos artefatos. As características espectrais foram dominadas pela idade (ρ de −0,43 a −0,51 nos controles) e acompanharam a potência do registro mesmo com a idade retirada; a distância espacial não acompanhou a idade, mas acompanhou a forma como o registro foi feito: a potência total do sinal, o lote de gravação e o artefato ocular. Um mesmo fator só confundiu onde também diferia entre os grupos: no HBN, equilibrado na aquisição, o ajuste não moveu o expoente; no TDBRAIN, o efeito do offset de olhos fechados desapareceu; no adhdata, um terço do efeito do expoente acompanhou artefatos mais frequentes no grupo TDAH. Limpar artefatos moveu pouco os efeitos; a referência e a faixa do ajuste moveram mais. O resultado é um tabuleiro de variáveis que indica, para cada característica e cada banco, de quais fatores desconfiar, e mostra que o componente aperiódico é o mais difícil de usar como marcador.

@p6
**Palavras-chave: **eletroencefalografia; TDAH; atividade aperiódica; confundidores; artefatos; reprodutibilidade.

@keep p7 [titulo] 1 INTRODUÇÃO

@p8
O diagnóstico do Transtorno do Déficit de Atenção e Hiperatividade (TDAH) é clínico, e há décadas se procura no eletroencefalograma (EEG) um marcador que o apoie [1]. A primeira candidata, a razão teta/beta (TBR), mostrou como essa busca pode falhar: a meta-análise de Arns et al. encontrou um efeito que diminuía ao longo dos anos, porque a TBR dos controles mudava entre os estudos, e Saad et al. argumentaram que o marcador é intrinsecamente sensível à maturação: ao longo do desenvolvimento, a potência teta cai e a beta sobe, de modo que a TBR diminui com a idade e um mesmo ponto de corte não vale para crianças de idades diferentes [2, 3]. A atenção passou depois ao componente aperiódico do espectro, separado do periódico pelo método de Donoghue et al. [4] e proposto como índice do balanço entre excitação e inibição [5]. A revisão sistemática do próprio Donoghue, porém, encontrou no TDAH relatos de expoente menor, maior e igual, e atribuiu a divergência à idade, à medicação e à condição de registro [6]. Em paralelo, representações espaciais e riemannianas da matriz de covariância, eficazes em interfaces cérebro-computador [7, 8], passaram a ser usadas para classificar o transtorno [9, 10].

@p9
Cada uma dessas famílias de características tem uma relação diferente com a idade, com os artefatos e com o equipamento. O expoente cai com a idade na infância e no envelhecimento [11–14]; a atividade muscular contamina o espectro acima de cerca de 20 Hz mesmo em sujeitos relaxados [15, 16]; os movimentos oculares agem nas frequências baixas [17]; e o local e o equipamento de coleta deslocam tanto o espectro quanto a covariância [18, 19]. Mesmo assim, os estudos de classificação, quando corrigem, costumam corrigir todas as características pelas mesmas variáveis, em geral idade e sexo, como se esses fatores agissem da mesma forma sobre todas elas.

@p10
Há ainda uma confusão de conceitos sobre o que é um confundidor. Confundir, aqui, é criar ou esconder uma diferença entre TDAH e controles que não vem do diagnóstico. Para isso, um fator precisa ter duas propriedades ao mesmo tempo, que são os critérios clássicos de confundimento, revistos por Greenland, Pearl e Robins com diagramas causais [20]: acompanhar a característica e diferir entre os grupos. A idade é o exemplo completo: o expoente cai com a idade e, se um grupo for mais novo que o outro, parte da diferença de expoente entre os grupos é só idade. Quando apenas uma das propriedades vale, não há confundimento. A idade acompanha o expoente, mas, se os dois grupos tiverem a mesma distribuição de idade, ela só acrescenta variação dentro de cada grupo, isto é, ruído, sem criar diferença entre eles; o sexo pode estar desequilibrado entre os grupos, mas, se não tiver relação com a frequência do pico alfa, é inofensivo para essa característica. Por fim, nem todo fator que tem as duas propriedades deve ser removido. Se as crianças com TDAH se mexem mais porque isso faz parte do quadro clínico, o artefato muscular está no caminho entre o diagnóstico e o sinal, e ajustar por ele retira parte do efeito que se quer medir, o viés de sobreajuste [21]; a ressonância magnética funcional enfrentou o mesmo dilema com o movimento da cabeça [22]. Por isso as duas propriedades precisam ser medidas separadamente, e o efeito, relatado com e sem o ajuste.

@p11
Este trabalho mede essas duas metades para cada característica, em três bancos públicos independentes, e acrescenta uma terceira dimensão: quanto o efeito de grupo muda com as escolhas de limpeza e de processamento. O produto é um tabuleiro de variáveis, que diz, para cada característica e cada banco, de quais fatores ela depende. A pergunta é anterior à classificação: antes de combinar características num modelo, é preciso saber o que cada uma carrega além do diagnóstico.

@keep p12 [titulo] 2 MÉTODOS

@p13
As Figuras 1 e 3 resumem o método. Os três bancos passaram pelo mesmo pré-processamento e, de cada registro, foram extraídas quatro famílias de características: aperiódica, periódica, TBR e espacial (Figura 1). Para cada característica e cada banco, foram feitas três medidas: a associação com fatores alheios ao diagnóstico nos controles (A), o desequilíbrio desses fatores entre TDAH e controles (B) e a mudança do efeito de grupo com o processamento e com o ajuste (C), cada uma com as variações que a Figura 3 enumera. As subseções seguintes detalham cada etapa nessa ordem.

@keep p14 [legenda] Figura 1 – Do registro às características: pipeline comum aos três bancos

@keep p15 [imagem] 

@keep p16 [fonte] Fonte: Elaborada pelos autores: Montagem, filtros, reamostragem e referência com MNE-Python [32], na ordem do 

@keep p17 [titulo] 2.1 Bancos de dados

@p18
Foram usados três bancos públicos (Tabela 1), escolhidos por cobrirem situações complementares: crianças em tarefa, crianças e adolescentes em repouso e adultos em repouso. Os três foram levados à mesma montagem de 19 canais do sistema 10–20, para que o mesmo pipeline e as mesmas definições de característica valessem para todos. O adhdata [23] já é registrado nesses 19 canais. Do TDBRAIN [25], registrado com 26 canais, foram mantidos os 19 do 10–20. O HBN [24] é registrado com a rede EGI HydroCel de 128 eletrodos mais o Cz de referência; dele foram mantidos os 19 eletrodos nas posições equivalentes às do 10–20, e o Cz, por ser a referência nativa, é um canal plano até a re-referência. O adhdata não informa idade nem sexo por participante, e o seu grupo TDAH havia usado metilfenidato. No HBN, os controles são participantes da mesma coorte de encaminhamento que não receberam diagnóstico no consenso clínico, o que aproxima o contraste do diagnóstico diferencial. O TDBRAIN compara pacientes de uma clínica com voluntários saudáveis, o desenho que Ransohoff e Feinstein mostraram inflar a acurácia aparente de um teste diagnóstico [26].

@keep p19 [legenda] Tabela 1 – Bancos de dados utilizados: canais nativos e usados após a harmonização, idade, condição e referênc

@keep p20 [fonte] Fonte: Nasrabadi et al. (2020), Alexander et al. (2017) e van Dijk et al. (2022); contagens após o pipeline; c

@keep p21 [titulo] 2.2 Processamento e características

@p22
Todos os registros passaram pelo mesmo pipeline: 19 canais, filtro passa-altas de 0,5 Hz, filtro *notch* na frequência da rede, reamostragem para 128 Hz e referência média comum (CAR). O espectro de cada canal foi estimado pelo método de Welch em épocas de 4 s e resumido pela mediana entre épocas. O FOOOF 1.1.1, no modo sem *knee* entre 1 e 40 Hz [4], forneceu o expoente χ (mediana dos canais, global e na região posterior), o offset e o pico alfa (frequência, IAF, e potência acima da reta aperiódica). A TBR clássica foi calculada na região frontocentral. Para a família espacial, a matriz de covariância de 1 a 40 Hz, estimada em 15 épocas fixas de 4 s, foi normalizada pelo traço, regularizada e projetada no espaço tangente log-euclidiano [7, 27]. As características seguem as definições dos trabalhos que as propuseram, para que os resultados possam ser comparados com a literatura. A Figura 2 ilustra o sinal em cada etapa.

@keep p23 [legenda] Figura 2 – O sinal em cada etapa do pré-processamento, do EEG bruto ao espectro

@keep p24 [imagem] 

@keep p25 [fonte] Fonte: Elaborada pelos autores: esquema desenhado pelos autores e ilustrado com o modelo generativo de imagem 

@p26
Formalmente, para cada canal *c*, o espectro mediano entre épocas, *P**_{c}*(*f*), foi decomposto em uma parte aperiódica e uma soma de picos [4]:

@keep p27 [equacao] log10 Pc(f) = bc − χc · log10 f + Σk Gk(f),	(1)

@p28
em que *b**_{c}* é o *offset* (a altura da reta no eixo log), *χ**_{c}* é o expoente (a inclinação da reta) e cada *G**_{k}* é uma gaussiana que descreve um pico periódico acima da reta. O ajuste é feito canal a canal, e os valores do participante são as medianas entre os 19 canais, *χ* = med*_{c}*(*χ**_{c}*) e *b* = med*_{c}*(*b**_{c}*); a mediana impede que um canal com ajuste ruim arraste o valor da pessoa. O expoente posterior é a mediana de *χ**_{c}* em P3, Pz, P4, O1 e O2. O pico alfa vem de um segundo ajuste de (1) ao espectro mediano desses cinco canais: a IAF é a frequência central do maior pico entre 7 e 14 Hz e a potência do alfa é a altura desse pico acima da reta aperiódica, em log_{10} de potência; sem pico nessa faixa, a potência do alfa é zero e a IAF é omissa. A TBR é

@keep p29 [equacao] TBR = log10 [ P(4–8 Hz) / P(13–30 Hz) ],	(2)

@p30
com *P* a potência média do espectro mediano em Fz, F3, F4, C3, C4 e Cz. Para a família espacial, a covariância *C* de cada participante (19 × 19, de 1 a 40 Hz) foi regularizada como *C*′ = 0,9 *C*/tr(*C*) + 0,1 *I*/19, o que remove a escala global do sinal e garante matrizes positivas definidas, e levada ao espaço tangente pelo logaritmo matricial, *L* = log(*C*′) [7, 27]. O vetor de cada participante é o triângulo superior de *L* − *L̄*_{controles}, com os elementos fora da diagonal multiplicados por √2, e a distância espacial usada no tabuleiro é a norma euclidiana desse vetor, isto é, quanto a covariância da pessoa se afasta da covariância média dos controles.

@keep p31 [titulo] 2.3 As duas metades de um confundidor

@p32
A Figura 3 mostra as três perguntas feitas a cada característica, em cada banco, e as variações testadas em cada uma; os parágrafos seguintes definem as medidas.

@keep p33 [legenda] Figura 3 – As três medidas e as variações testadas para cada característica e cada banco

@keep p34 [imagem] 

@keep p35 [fonte] Fonte: Elaborada pelos autores: (A) Correlação de Spearman nos controles com seis fatores. (B) Diferença padro

@p36
**Associação nos controles.** Para cada característica e cada banco, calculou-se a correlação de Spearman (ρ) com idade, sexo, duração do registro, índice muscular (potência de 30 a 45 Hz em T7 e T8), índice ocular (razão da potência de 1 a 4 Hz entre Fp1/Fp2 e O1/O2) e potência total (logaritmo do traço da covariância sem normalização), apenas nos controles, para que o diagnóstico não entre na associação. A família espacial foi resumida pela distância de cada pessoa ao centro dos controles no espaço tangente. Com 179 controles, |ρ| a partir de 0,19 já difere de zero a p < 0,01, e |ρ| = 0,4 corresponde a 16 % de variância compartilhada; por isso, as faixas de leitura foram |ρ| < 0,2, 0,2 a 0,4 e ≥ 0,4. Nos 60 controles do adhdata, o IC de cada ρ tem meia-largura de cerca de 0,26.

@p37
**Desequilíbrio entre os grupos.** Para os mesmos fatores, calculou-se a diferença padronizada entre TDAH e controles e, para os índices de aquisição, o g de Hedges ajustado por idade, idade², sexo, idade × sexo e logaritmo da duração (C0), com IC de 95 % por *bootstrap* estratificado (300 reamostragens) [28].

@p38
A diferença padronizada de um fator *x* entre os grupos é

@keep p39 [equacao] d = (x̄TDAH − x̄controle) / sp,	(3)

@p40
com *s**_{p}* o desvio-padrão combinado dos dois grupos. O g de Hedges ajustado generaliza essa medida para quando se mantêm fixas as covariáveis. Para a característica *v* do participante *i*, ajusta-se por mínimos quadrados a regressão linear

@keep p41 [equacao] vi = β0 + β1 · grupoi + Σj γj zij + εi,	(4)

@p42
em que grupo*_{i}* vale 1 para TDAH e 0 para controle, e *z**_{ij}* são as covariáveis do conjunto C0: idade padronizada, idade², sexo (1 para feminino), idade × sexo e logaritmo da duração do registro, padronizado. O coeficiente *β*_{1} é a diferença entre TDAH e controle entre participantes de mesma idade, sexo e duração; dividido pelo desvio-padrão residual da regressão, *σ**̂**_{ε}*, ele vira um efeito em unidades de desvio-padrão, e a correção de Hedges para amostras pequenas [28] completa a medida:

@keep p43 [equacao] g = J · β̂1 / σ̂ε,      J = 1 − 3 / (4n − 9),	(5)

@p44
com *n* o número de participantes. Sem covariáveis, (5) reduz-se à diferença das médias dividida pelo desvio-padrão combinado, como em (3); com centenas de participantes, *J* fica próximo de 1. O IC de 95 % é dado pelos percentis 2,5 e 97,5 de 300 reamostragens *bootstrap*, sorteadas com reposição dentro de cada grupo e com a regressão (4) refeita em cada uma. Nas análises de sensibilidade, a potência total e os índices muscular e ocular entram em (4) como novos termos *z**_{ij}*, e o g recalculado é comparado ao g com C0.

@p45
**Sensibilidade do efeito.** O efeito de grupo de cada característica é o g de Hedges ajustado por C0. Ele foi recalculado (i) em oito variantes de processamento: rejeição das épocas com mais músculo ou com mais olho (acima do percentil 90 do próprio registro), análise de componentes independentes com rotulagem pelo ICLabel e remoção dos componentes oculares [29, 30], só os 13 canais internos, só os primeiros 60 s, a referência nativa e o ajuste aperiódico entre 1 e 30 Hz; e (ii) com a potência total, e depois com ela e os índices muscular e ocular, acrescentados às covariáveis. Uma mudança de 0,20 no g é do tamanho do próprio efeito do expoente no HBN, e as faixas de leitura foram |Δg| < 0,10, 0,10 a 0,20 e ≥ 0,20. Para a família espacial, mediu-se a AUC do espaço tangente por participante, com validação cruzada estratificada (5 dobras × 5 repetições) e regressão logística regularizada.

@p46
**Critério de leitura.** Os valores do tabuleiro não foram corrigidos para comparações múltiplas [31] e são lidos como um mapa descritivo. Só é tratado como resultado o que se repete entre bancos, estados de registro ou variantes; o que aparece num banco só é tratado como linha de ação.

@keep p47 [titulo] 2.4 Implementação e reprodutibilidade

@p48
O processamento foi implementado em Python 3.11 com bibliotecas abertas e consolidadas, de modo que cada etapa corresponde a uma ferramenta já usada e citada na literatura (Tabela 2): MNE-Python para filtros, reamostragem, espectro, referências e ICA [32]; FOOOF 1.1.1, hoje distribuído como *specparam*, para a parametrização do espectro [4]; mne-icalabel para a rotulagem das componentes [29, 30]; e NumPy, SciPy e scikit-learn para o espaço tangente, a estatística e a classificação [33–35]. O espaço tangente foi escrito diretamente em NumPy, seguindo Barachant et al. [7], para que a centragem fosse ajustada só nos participantes de treino de cada dobra. A ordem das etapas segue a do pipeline PREP, com a rede elétrica removida antes da referência média [36], mas sem a sua detecção robusta de canais ruins; o relato segue, quando aplicável, as recomendações do comitê COBIDAS-MEEG [37]. A referência REST segue Yao [38]; o espectro, o método de Welch [39]; e a ICA, o algoritmo *infomax* estendido [40].

@p49
O código de análise soma 44 arquivos e cerca de 7.800 linhas de código (10.200 com comentários e linhas em branco): 2.500 de pré-processamento e extração de características, 1.400 de replicação das alegações publicadas, 1.700 de desenhos de teste e estatística, 800 de classificação e protocolo de avaliação, 1.000 de figuras e 360 de montagem do documento. As replicações reimplementam, a partir dos textos, as análises de Tröndle et al., Cellier et al., Donoghue et al. e TaghiBeyglou et al. [9, 13, 14, 41], com os mesmos parâmetros publicados sempre que descritos. O pré-registro e os módulos de extração foram congelados com SHA-256 antes das análises, e os hashes conferem com os arquivos usados; os módulos de extração têm testes automatizados com espectros sintéticos. O pré-processamento gerou uma única vez o conjunto de 19 canais a 128 Hz, que é a entrada fixa de todas as análises. As extrações rodaram em duas máquinas, de 8 e de 32 núcleos; a reextração dos três bancos nas nove versões do tabuleiro levou cerca de 2 h 20 min na estação local. Os dados por participante permanecem sob os termos de uso de cada banco. [A confirmar: repositório público do código.]

@keep p50 [legenda] Tabela 2 – Etapas do pipeline, ferramentas, parâmetros e referências.

@keep p51 [fonte] Fonte: elaborada pelos autores; parâmetros conferidos no código.

@keep p52 [titulo] 3 RESULTADOS E DISCUSSÃO

@keep p53 [titulo] 3.1 O tabuleiro: cada família tem a sua impressão digital

@p54
A Figura 4 resume, para cada característica e cada banco, as duas primeiras perguntas do trabalho: de que fatores ela depende nos controles e quanto o seu efeito de grupo muda com cada escolha de processamento. A última linha mostra o efeito de grupo na base. O número em cada casa dá o tamanho e a direção: o sinal importa, porque um fator que empurra a característica no mesmo sentido do efeito de grupo pode criá-lo, e um que empurra no sentido oposto pode escondê-lo.

@p55
**As características espectrais são dominadas pela idade.** Nos controles dos dois bancos de repouso, a idade acompanha o expoente (ρ = −0,43 no HBN e −0,45 no TDBRAIN), o offset (−0,47 e −0,50) e a TBR (−0,51 e −0,43). É um dos achados mais consistentes da literatura aperiódica, descrito no envelhecimento [11], na infância [12] e em amostras grandes dos 3 aos 24 anos [13, 14]; para a TBR, é a razão pela qual o efeito da meta-análise de Arns et al. diminuía ao longo dos anos [2]. Uma correção por idade, porém, não se transfere entre bancos: a inclinação etária do expoente nos controles do HBN foi 3,4 vezes maior que a publicada por Cellier et al. [13].

@keep p56 [legenda] Figura 4 – Tabuleiro de variáveis: quanto cada característica depende da amostra, dos artefatos e do processam

@keep p57 [imagem] 

@keep p58 [fonte] Fonte: Elaborada pelos autores. Cada coluna é uma característica (χ: expoente aperiódico, global ou da região 

@p59
A distância espacial não acompanha a idade, mas acompanha a forma como o registro foi feito. A distância de cada pessoa ao centro dos controles no espaço tangente tem ρ = −0,01 com a idade no HBN e −0,03 no TDBRAIN, e também não tem relação em forma de U (|ρ| ≤ 0,06 com o afastamento da idade mediana, nos quatro blocos). Ela acompanha, em compensação, a potência total do registro nos três bancos (ρ de +0,23 a +0,59), e o índice ocular no adhdata (−0,40). A geometria riemanniana capta tudo o que muda a distribuição do sinal entre os eletrodos [7], o que inclui o equipamento, o contato e o olho. No TDBRAIN, onde os grupos foram registrados em lotes diferentes, o espaço tangente separa adultos com e sem TDAH com AUC de 0,85 e, dentro dos controles saudáveis, prevê o lote de gravação com AUC de 0,86 a 0,90 (análise pré-registrada anterior); no HBN, com a aquisição equilibrada entre os grupos, fica no acaso (0,44 a 0,46). É o mesmo fenômeno em que o local de coleta previu o diagnóstico na competição ADHD-200 [42] e o hospital de origem previu a doença em radiografias de tórax [43].

@p60
A idade de um, o registro de outro. O que confunde o expoente (a idade) não confunde a família espacial, e o que confunde a família espacial (a aquisição: equipamento, lote de gravação, contato dos eletrodos e artefato ocular) atinge o expoente por outra via, como se vê a seguir. Corrigir todas as características pelo mesmo conjunto de covariáveis corrige o que não precisa e deixa passar o que importa.

@keep p61 [titulo] 3.2 As duas metades de um confundidor

@p62
A Figura 5 resume o raciocínio da Introdução com as três situações encontradas nos bancos: o mesmo fator pode ser inofensivo, confundidor ou parte do próprio quadro, conforme difira ou não entre os grupos e conforme pertença ou não ao transtorno. A Figura 6 mostra os dados.

@keep p63 [legenda] Figura 5 – Quando um fator confunde: três situações encontradas nos bancos

@keep p64 [imagem] 

@keep p65 [fonte] Fonte: Elaborada pelos autores: esquema desenhado pelos autores e ilustrado com o modelo generativo de imagem 

@p66
A Figura 6 junta as duas metades. O painel A mostra quanto cada fator difere entre TDAH e controles em cada banco; o painel B, quanto o efeito de três características muda quando a potência total e os índices de artefato entram nas covariáveis. Cada ponto do painel B é o g da equação (5), recalculado com três conjuntos de covariáveis em (4): C0; C0 mais a potência total; e C0 mais a potência total e os índices muscular e ocular. Um ponto que não se move indica que o fator não explica a diferença entre os grupos; um que caminha para o zero indica que parte da diferença acompanhava o fator; um que se afasta do zero indica que o fator a escondia.

@p67
Cada banco é uma situação diferente. No HBN, os grupos diferem pouco na idade (diferença padronizada de +0,20) e no sexo (−0,27) e quase nada nos índices de aquisição (−0,07 no músculo, +0,04 no olho e −0,04 na potência total, todos com IC incluindo zero depois do ajuste por C0). Ajustar pela potência e pelos artefatos não moveu nenhum efeito: o expoente ficou entre −0,22 e −0,25 de olhos fechados, e a potência do alfa, em torno de −0,25 de olhos fechados e −0,30 de olhos abertos. No TDBRAIN, os grupos diferem em quase tudo: duração (−0,83), sexo (−0,40), idade (−0,33), músculo (g = −0,30; IC −0,54 a −0,04) e potência total (−0,24; IC −0,41 a −0,10). Ali, o efeito do offset de olhos fechados passou de −0,19 para +0,01 com a potência total nas covariáveis. No adhdata, as crianças com TDAH têm mais atividade muscular (g = +0,56; IC 0,20 a 0,94) e ocular (+0,62; IC 0,28 a 0,97) durante a tarefa, e o ajuste reduziu o expoente de −0,91 para −0,63, ainda com IC excluindo zero; o expoente posterior e a TBR caíram para dentro do acaso.

@p68
**Confundidor**** ou parte do quadro?** No adhdata, o artefato pode não ser um confundidor no sentido estrito. Se as crianças com TDAH se movem mais durante a tarefa, a tensão muscular e o movimento ocular fazem parte da manifestação do transtorno, e removê-los do efeito é remover parte do que se quer medir, o viés de sobreajuste [21]. A ressonância magnética funcional enfrentou o mesmo dilema com o movimento da cabeça, que produz correlações espúrias, mas sistemáticas [22]. Nenhuma das versões é a verdadeira; por isso o efeito deve ser relatado com e sem o ajuste. O próprio ajuste pode criar efeitos: no TDBRAIN de olhos abertos, ele tornou o expoente negativo e significativo (de −0,21 para −0,31), sinal de que o músculo e o olho também podem esconder uma diferença.

@keep p69 [legenda] Figura 6 – As duas metades de um confundidor: desequilíbrio entre os grupos e efeito do ajuste

@keep p70 [imagem] 

@keep p71 [fonte] Fonte: Elaborada pelos autores. (A) Diferença padronizada (TDAH − controle, sem ajuste) de cada fator em cada 

@keep p72 [titulo] 3.3 A potência total, um confundidor ambíguo

@p73
A linha da potência total é a mais vermelha da família aperiódica. Ela mistura três origens: a amplitude do próprio sinal cerebral, que diminui ao longo do desenvolvimento [13]; a condução pelo crânio, cuja espessura se relaciona com a magnitude do alfa medida no escalpo [44]; e a aquisição, com impedância, ganho e artefatos. Um ganho puro, que multiplica o espectro inteiro pelo mesmo fator, só desloca para cima a reta ajustada pelo FOOOF: muda o offset, que é o deslocamento de banda larga do espectro [4], e deixa o expoente intacto (Figura 7). A correlação do offset com a potência total (+0,62 a +0,76, nos cinco blocos) é, assim, quase uma consequência da definição. A do expoente (+0,31 a +0,51) diz mais: as diferenças de potência entre as pessoas concentram-se nas frequências baixas, as mesmas que a idade, o olho e o movimento aumentam.

@keep p74 [legenda] Figura 7 – Por que a potência total afeta o offset e o expoente de formas diferentes

@keep p75 [imagem] 

@keep p76 [fonte] Fonte: Elaborada pelos autores: esquema desenhado pelos autores e ilustrado com o modelo generativo de imagem 

@p77
Uma objeção previsível é que expoente e potência andariam juntos só porque ambos dependem da idade. Retirada a idade, a correlação do expoente com a potência cai pouco no HBN (de 0,46 para 0,41 de olhos fechados) e no TDBRAIN de olhos abertos (de 0,47 para 0,46), e só se reduz no TDBRAIN de olhos fechados (de 0,31 para 0,19); a do offset não muda. A potência do registro é, portanto, um fator próprio. O componente periódico quase sempre escapa dela, porque a potência do pico é medida acima da reta aperiódica e a IAF é uma posição no eixo das frequências: a potência do alfa é a única característica cujo efeito no HBN manteve o IC excluindo zero em todas as variantes, nos dois estados de registro (g de −0,22 a −0,32).

@keep p78 [titulo] 3.4 Limpar muda pouco; processar muda muito

@p79
Rejeitar as épocas mais contaminadas por músculo ou por olho, remover os componentes oculares ou usar só 60 s moveu o g em menos de 0,10 na maior parte das casas. O expoente do adhdata ficou entre −0,85 e −0,97 com qualquer limpeza, e o do HBN de olhos fechados, entre −0,20 e −0,26. A ressalva é de desenho: como o espectro é a mediana entre épocas, descartar os 10 % piores tem efeito limitado por construção, e a contaminação muscular de fundo, presente em todas as épocas [15, 16], não sai com a rejeição. As limpezas que agem no sinal inteiro moveram mais: no TDBRAIN, a ICA ocular e a retirada das bordas reduziram a AUC espacial de 0,85 para 0,81 e 0,80, sinal de que parte do ganho espacial estava no olho e nas bordas [17, 45].

@keep p80 [legenda] Figura 8 – Limpar muda pouco o efeito de grupo; escolher a referência e a faixa do ajuste muda muito

@keep p81 [imagem] 

@keep p82 [fonte] Fonte: Elaborada pelos autores: esquema desenhado pelos autores e ilustrado com o modelo generativo de imagem 

@p83
Duas escolhas que quase nunca são relatadas moveram mais que qualquer limpeza (Figura 8). Restringir o ajuste aperiódico a 1–30 Hz reduziu o efeito do expoente no adhdata de −0,91 para −0,74 e no HBN de olhos abertos de −0,11 para −0,02, como previsto pela literatura metodológica [46]. A referência nativa criou no TDBRAIN um efeito positivo que não existe com CAR: o expoente de olhos fechados passou de +0,06 para +0,39 (IC 0,08 a 0,62). Essa referência, nas mastoides, também aumentou os índices de artefato do banco (o muscular em 0,38 de olhos fechados e o ocular em 0,32 de olhos abertos, em log10), o que sugere que parte do efeito é artefato somado a todos os canais [47]. Com três escolhas razoáveis de processamento, o TDBRAIN produz três conclusões diferentes sobre o expoente no TDAH adulto.

@p84
Os próprios índices de artefato não medem a mesma coisa em todo lugar. O índice muscular correlaciona-se negativamente com o expoente nos adultos (−0,35), como se espera do músculo, que eleva a ponta alta do espectro, mas positivamente no HBN (+0,11), onde acompanha sobretudo a potência de banda larga (+0,57 com o offset). O índice ocular muda de sinal até dentro do TDBRAIN: −0,14 com o expoente de olhos fechados e +0,54 de olhos abertos, quando mede sobretudo piscadas. Antes de entrar num modelo como covariáveis, os índices de artefato precisam ser validados em cada banco e estado.

@keep p85 [titulo] 3.5 O que se repete e o que é de um banco só

@p86
Separar os dois dá a cada afirmação o seu peso. Repetem-se em todos os bancos em que puderam ser testados: a dependência das características espectrais em relação à idade; a ligação do offset e do expoente com a potência do registro, que sobrevive à retirada da idade; a ligação da família espacial com a potência total; a ausência de relação entre a distância espacial e a idade; e o sentido negativo da potência do alfa no TDAH nos quatro blocos de repouso (HBN −0,25 e −0,30; TDBRAIN −0,17 e −0,20, com IC incluindo zero num banco em que o menor efeito detectável é de cerca de 0,36). Aparecem num banco só e são linhas de ação: o efeito grande do expoente no adhdata, a discriminação espacial do TDBRAIN e o efeito criado pela referência nativa nos adultos.

@p87
Os efeitos de grupo são pequenos onde a aquisição está equilibrada. No HBN, o expoente tem g = −0,22 (IC −0,33 a −0,10) de olhos fechados e não difere de zero de olhos abertos; com os grupos pareados por idade e sexo, a sua AUC fica em 0,55, no limite do nulo, o que concorda com o registered report de Panda et al. no mesmo banco [48]. Uma AUC de 0,55 equivale a d ≈ 0,18: o pareamento reduziu o poder, não necessariamente o efeito. A ausência de um efeito replicável, nesse quadro, é ausência de evidência, e não evidência de ausência [49]: com 179 controles, o HBN só detecta com segurança efeitos de g ≈ 0,22 ou maiores, o tamanho do único sinal estável do estudo.

@keep p88 [titulo] 3.6 Limitações

@p89
O adhdata não tem idade nem sexo por participante, o que impede separar o seu efeito das diferenças de condição e de amostra, e o seu grupo TDAH foi exposto ao metilfenidato; se o fármaco reduz o ruído 1/f, como relataram Pertermann et al. [50], ele atenuaria a diferença do expoente, e não a criaria. A medicação não está mensurada de forma útil em nenhum banco. O TDBRAIN tem o lote de gravação associado ao diagnóstico. A ICA muscular não foi avaliada, porque o ICLabel quase não rotula músculo a 128 Hz. Alguns números de aquisição do TDBRAIN vêm de uma análise pré-registrada anterior, e estão marcados. O cofre de confirmação do HBN permanece fechado, e os resultados devem ser confirmados nele.

@keep p90 [titulo] 4 CONCLUSÃO

@p91
Corrigir todas as características de EEG pelo mesmo conjunto de covariáveis parte de uma suposição que os três bancos refutam: a de que elas têm os mesmos confundidores. As características espectrais são dominadas pela idade e pela potência do registro; a família espacial não acompanha a idade, mas acompanha o equipamento, o lote de gravação e os artefatos do registro; e o componente periódico, medido acima da reta aperiódica, é o que menos depende dos dois.

@p92
Um fator só confunde onde também difere entre os grupos. O mesmo índice muscular é inofensivo no HBN, onde os grupos são equilibrados, e pesa no TDBRAIN e no adhdata, onde não são; e, no adhdata, pode ser parte do próprio quadro de TDAH. Por isso, o efeito deve ser relatado com e sem o ajuste, e a escolha entre eles é uma decisão causal, e não estatística.

@p93
O componente aperiódico é a característica mais difícil de usar como marcador: depende da idade, da potência do registro, dos artefatos, da referência e da faixa do ajuste, e o sinal do seu efeito muda com o banco e com o processamento. O offset, o único parâmetro com efeito nos três bancos, mede sobretudo a potência de banda larga. Limpar artefatos muda pouco o efeito de grupo; escolher a referência e a faixa do ajuste muda muito, e essas escolhas precisam ser fixadas e relatadas.

@p94
O tabuleiro de variáveis é a contribuição deste trabalho para os seguintes. Ele indica, para cada característica e cada banco, de quais fatores desconfiar. Um classificador generalizável precisará usar essa informação na sua estrutura, ajustando cada família pelas covariáveis que a confundem, em vez de corrigir todas da mesma forma.

@keep p95 [titulo] Declarações

@p96
Uso de ferramentas de inteligência artificial. Um modelo de linguagem (Claude, Anthropic) foi usado como apoio na triagem e na extração da revisão de escopo, com protocolo escrito e verificação por trecho literal, na escrita do código de análise, na conferência das referências e na redação, e um modelo generativo de imagem (Gemini 3 Pro Image, Google) ilustrou as Figuras 1, 2, 3, 5, 7 e 8 a partir de esquemas desenhados pelos autores. Os autores revisaram o conteúdo e assumem a responsabilidade por ele. Ética. Só foram usados bancos públicos e desidentificados, sob os termos de uso de cada um (HBN e TDBRAIN); cada banco tem a aprovação ética descrita pelos seus autores. Dados e código. Os bancos estão disponíveis nas fontes originais; os pré-registros, os hashes e o código de análise acompanham o projeto.

@keep p97 [titulo] REFERÊNCIAS

@p98
[1] FARAONE, S. V. et al. The World Federation of ADHD International Consensus Statement: 208 evidence-based conclusions about the disorder. **Neuroscience & Biobehavioral Reviews**, v. 128, p. 789–818, 2021. DOI: 10.1016/j.neubiorev.2021.01.022.

@p99
[2] ARNS, M.; CONNERS, C. K.; KRAEMER, H. C. A decade of EEG theta/beta ratio research in ADHD: a meta-analysis. **Journal of Attention Disorders**, v. 17, n. 5, p. 374–383, 2013. DOI: 10.1177/1087054712460087.

@p100
[3] SAAD, J. F. et al. Is the theta/beta EEG marker for ADHD inherently flawed? **Journal of Attention Disorders**, v. 22, n. 9, p. 815–826, 2018. DOI: 10.1177/1087054715578270.

@p101
[4] DONOGHUE, T. et al. Parameterizing neural power spectra into periodic and aperiodic components. **Nature Neuroscience**, v. 23, n. 12, p. 1655–1665, 2020. DOI: 10.1038/s41593-020-00744-x.

@p102
[5] GAO, R.; PETERSON, E. J.; VOYTEK, B. Inferring synaptic excitation/inhibition balance from field potentials. **NeuroImage**, v. 158, p. 70–78, 2017. DOI: 10.1016/j.neuroimage.2017.06.078.

@p103
[6] DONOGHUE, T. A systematic review of aperiodic neural activity in clinical investigations. **European Journal of Neuroscience**, v. 62, n. 7, art. e70255, 2025. DOI: 10.1111/ejn.70255.

@p104
[7] BARACHANT, A. et al. Multiclass brain–computer interface classification by Riemannian geometry. **IEEE Transactions on Biomedical Engineering**, v. 59, n. 4, p. 920–928, 2012. DOI: 10.1109/TBME.2011.2172210.

@p105
[8] CONGEDO, M.; BARACHANT, A.; BHATIA, R. Riemannian geometry for EEG-based brain-computer interfaces; a primer and a review. **Brain-Computer Interfaces**, v. 4, n. 3, p. 155–174, 2017. DOI: 10.1080/2326263X.2017.1297192.

@p106
[9] TAGHIBEYGLOU, B. et al. ADHD diagnosis in children using common spatial pattern and nonlinear analysis of filter banked EEG. In: IRANIAN CONFERENCE ON ELECTRICAL ENGINEERING (ICEE), 28., 2020. **Proceedings** [...]. [S. l.]: IEEE, 2020. p. 1–5. DOI: 10.1109/ICEE50131.2020.9260711.

@p107
[10] DUFFY, F. H. et al. A unique pattern of cortical connectivity characterizes patients with attention deficit disorders: a large electroencephalographic coherence study. **BMC Medicine**, v. 15, n. 1, art. 51, 2017. DOI: 10.1186/s12916-017-0805-9.

@p108
[11] VOYTEK, B. et al. Age-related changes in 1/f neural electrophysiological noise. **The Journal of Neuroscience**, v. 35, n. 38, p. 13257–13265, 2015. DOI: 10.1523/JNEUROSCI.2332-14.2015.

@p109
[12] HILL, A. T. et al. Periodic and aperiodic neural activity displays age-dependent changes across early-to-middle childhood. **Developmental Cognitive Neuroscience**, v. 54, art. 101076, 2022. DOI: 10.1016/j.dcn.2022.101076.

@p110
[13] CELLIER, D. et al. The development of theta and alpha neural oscillations from ages 3 to 24 years. **Developmental Cognitive Neuroscience**, v. 50, art. 100969, 2021. DOI: 10.1016/j.dcn.2021.100969.

@p111
[14] TRÖNDLE, M. et al. Decomposing the role of alpha oscillations during brain maturation. **eLife**, v. 11, art. e77571, 2022. DOI: 10.7554/eLife.77571.

@p112
[15] WHITHAM, E. M. et al. Scalp electrical recording during paralysis: quantitative evidence that EEG frequencies above 20 Hz are contaminated by EMG. **Clinical Neurophysiology**, v. 118, n. 8, p. 1877–1888, 2007. DOI: 10.1016/j.clinph.2007.04.027.

@p113
[16] MUTHUKUMARASWAMY, S. D. High-frequency brain activity and muscle artifacts in MEG/EEG: a review and recommendations. **Frontiers in Human Neuroscience**, v. 7, art. 138, 2013. DOI: 10.3389/fnhum.2013.00138.

@p114
[17] PLÖCHL, M.; OSSANDÓN, J. P.; KÖNIG, P. Combining EEG and eye tracking: identification, characterization, and correction of eye movement artifacts in electroencephalographic data. **Frontiers in Human Neuroscience**, v. 6, art. 278, 2012. DOI: 10.3389/fnhum.2012.00278.

@p115
[18] LI, M. et al. Harmonized-Multinational qEEG norms (HarMNqEEG). **NeuroImage**, v. 256, art. 119190, 2022. DOI: 10.1016/j.neuroimage.2022.119190.

@p116
[19] JARAMILLO-JIMENEZ, A. et al. ComBat models for harmonization of resting-state EEG features in multisite studies. **Clinical Neurophysiology**, v. 167, p. 241–253, 2024. DOI: 10.1016/j.clinph.2024.09.019.

@p117
[20] GREENLAND, S.; PEARL, J.; ROBINS, J. M. Causal diagrams for epidemiologic research. **Epidemiology**, v. 10, n. 1, p. 37–48, 1999. DOI: 10.1097/00001648-199901000-00008.

@p118
[21] SCHISTERMAN, E. F.; COLE, S. R.; PLATT, R. W. Overadjustment bias and unnecessary adjustment in epidemiologic studies. **Epidemiology**, v. 20, n. 4, p. 488–495, 2009. DOI: 10.1097/EDE.0b013e3181a819a1.

@p119
[22] POWER, J. D. et al. Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. **NeuroImage**, v. 59, n. 3, p. 2142–2154, 2012. DOI: 10.1016/j.neuroimage.2011.10.018.

@p120
[23] MOTI NASRABADI, A. et al. **EEG data for ADHD / control children**. [S. l.]: IEEE Dataport, 2020. [Conjunto de dados]. DOI: 10.21227/rzfh-zn36. Disponível em: https://ieee-dataport.org/open-access/eeg-data-adhd-control-children. Acesso em: 8 nov. 2026.

@p121
[24] ALEXANDER, L. M. et al. An open resource for transdiagnostic research in pediatric mental health and learning disorders. **Scientific Data**, v. 4, art. 170181, 2017. DOI: 10.1038/sdata.2017.181.

@p122
[25] VAN DIJK, H. et al. The two decades brainclinics research archive for insights in neurophysiology (TDBRAIN) database. **Scientific Data**, v. 9, art. 333, 2022. DOI: 10.1038/s41597-022-01409-z.

@p123
[26] RANSOHOFF, D. F.; FEINSTEIN, A. R. Problems of spectrum and bias in evaluating the efficacy of diagnostic tests. **New England Journal of Medicine**, v. 299, n. 17, p. 926–930, 1978. DOI: 10.1056/NEJM197810262991705.

@p124
[27] ARSIGNY, V. et al. Geometric means in a novel vector space structure on symmetric positive-definite matrices. **SIAM Journal on Matrix Analysis and Applications**, v. 29, n. 1, p. 328–347, 2007. DOI: 10.1137/050637996.

@p125
[28] HEDGES, L. V. Distribution theory for Glass's estimator of effect size and related estimators. **Journal of Educational Statistics**, v. 6, n. 2, p. 107–128, 1981. DOI: 10.3102/10769986006002107.

@p126
[29] PION-TONACHINI, L.; KREUTZ-DELGADO, K.; MAKEIG, S. ICLabel: an automated electroencephalographic independent component classifier, dataset, and website. **NeuroImage**, v. 198, p. 181–197, 2019. DOI: 10.1016/j.neuroimage.2019.05.026.

@p127
[30] LI, A. et al. MNE-ICALabel: automatically annotating ICA components with ICLabel in Python. **Journal of Open Source Software**, v. 7, n. 76, art. 4484, 2022. DOI: 10.21105/joss.04484.

@p128
[31] BENJAMINI, Y.; HOCHBERG, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. **Journal of the Royal Statistical Society: Series B**, v. 57, n. 1, p. 289–300, 1995. DOI: 10.1111/j.2517-6161.1995.tb02031.x.

@p129
[32] GRAMFORT, A. et al. MEG and EEG data analysis with MNE-Python. **Frontiers in Neuroscience**, v. 7, art. 267, 2013. DOI: 10.3389/fnins.2013.00267.

@p130
[33] HARRIS, C. R. et al. Array programming with NumPy. **Nature**, v. 585, n. 7825, p. 357–362, 2020. DOI: 10.1038/s41586-020-2649-2.

@p131
[34] VIRTANEN, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. **Nature** **Methods**, v. 17, n. 3, p. 261–272, 2020. DOI: 10.1038/s41592-019-0686-2.

@p132
[35] PEDREGOSA, F. et al. Scikit-learn: machine learning in Python. **Journal of Machine Learning Research**, v. 12, p. 2825–2830, 2011.

@p133
[36] BIGDELY-SHAMLO, N. et al. The PREP pipeline: standardized preprocessing for large-scale EEG analysis. **Frontiers in Neuroinformatics**, v. 9, art. 16, 2015. DOI: 10.3389/fninf.2015.00016.

@p134
[37] PERNET, C. et al. Issues and recommendations from the OHBM COBIDAS MEEG committee for reproducible EEG and MEG research. **Nature Neuroscience**, v. 23, n. 12, p. 1473–1483, 2020. DOI: 10.1038/s41593-020-00709-0.

@p135
[38] YAO, D. A method to standardize a reference of scalp EEG recordings to a point at infinity. **Physiological Measurement**, v. 22, n. 4, p. 693–711, 2001. DOI: 10.1088/0967-3334/22/4/305.

@p136
[39] WELCH, P. D. The use of fast Fourier transform for the estimation of power spectra: a method based on time averaging over short, modified periodograms. **IEEE Transactions on Audio and Electroacoustics**, v. 15, n. 2, p. 70–73, 1967. DOI: 10.1109/TAU.1967.1161901.

@p137
[40] LEE, T.-W.; GIROLAMI, M.; SEJNOWSKI, T. J. Independent component analysis using an extended infomax algorithm for mixed subgaussian and supergaussian sources. **Neural Computation**, v. 11, n. 2, p. 417–441, 1999. DOI: 10.1162/089976699300016719.

@p138
[41] DONOGHUE, T.; DOMINGUEZ, J.; VOYTEK, B. Electrophysiological frequency band ratio measures conflate periodic and aperiodic neural activity. **eNeuro**, v. 7, n. 6, ENEURO.0192-20.2020, 2020. DOI: 10.1523/ENEURO.0192-20.2020.

@p139
[42] BROWN, M. R. G. et al. ADHD-200 Global Competition: diagnosing ADHD using personal characteristic data can outperform resting state fMRI measurements. **Frontiers in Systems Neuroscience**, v. 6, art. 69, 2012. DOI: 10.3389/fnsys.2012.00069.

@p140
[43] ZECH, J. R. et al. Variable generalization performance of a deep learning model to detect pneumonia in chest radiographs: a cross-sectional study. **PLOS Medicine**, v. 15, n. 11, art. e1002683, 2018. DOI: 10.1371/journal.pmed.1002683.

@p141
[44] HAGEMANN, D. et al. Skull thickness and magnitude of EEG alpha activity. **Clinical Neurophysiology**, v. 119, n. 6, p. 1271–1280, 2008. DOI: 10.1016/j.clinph.2008.02.010.

@p142
[45] GONCHAROVA, I. I. et al. EMG contamination of EEG: spectral and topographical characteristics. **Clinical Neurophysiology**, v. 114, n. 9, p. 1580–1593, 2003. DOI: 10.1016/S1388-2457(03)00093-2.

@p143
[46] GERSTER, M. et al. Separating neural oscillations from aperiodic 1/f activity: challenges and recommendations. **Neuroinformatics**, v. 20, n. 4, p. 991–1012, 2022. DOI: 10.1007/s12021-022-09581-8.

@p144
[47] YAO, D. et al. Which reference should we use for EEG and ERP practice? **Brain Topography**, v. 32, n. 4, p. 530–549, 2019. DOI: 10.1007/s10548-019-00707-x.

@p145
[48] PANDA, E. J. et al. Relating ADHD to the aperiodic component of the EEG signal in both males and females across development. **Developmental Cognitive Neuroscience**, v. 79, art. 101732, 2026. DOI: 10.1016/j.dcn.2026.101732.

@p146
[49] ALTMAN, D. G.; BLAND, J. M. Statistics notes: absence of evidence is not evidence of absence. **BMJ**, v. 311, n. 7003, p. 485, 1995. DOI: 10.1136/bmj.311.7003.485.

@p147
[50] PERTERMANN, M. et al. The modulation of neural noise underlies the effectiveness of methylphenidate treatment in attention-deficit/hyperactivity disorder. **Biological Psychiatry: Cognitive Neuroscience and Neuroimaging**, v. 4, n. 8, p. 743–750, 2019. DOI: 10.1016/j.bpsc.2019.03.011.
