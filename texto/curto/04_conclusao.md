# 4 CONCLUSÃO

Corrigir todas as características de EEG pelo mesmo conjunto de covariáveis parte de uma suposição que os três bancos refutam: a de que elas têm os mesmos confundidores. As características espectrais são dominadas pela idade e pela potência do registro; a família espacial não acompanha a idade, mas carrega a aquisição; e o componente periódico, medido acima da reta aperiódica, é o que menos depende dos dois.

Um fator só confunde onde também difere entre os grupos. O mesmo índice muscular é inofensivo no HBN, onde os grupos são equilibrados, e pesa no TDBRAIN e no adhdata, onde não são; e, no adhdata, pode ser parte do próprio quadro de TDAH. Por isso, o efeito deve ser relatado com e sem o ajuste, e a escolha entre eles é uma decisão causal, e não estatística.

O componente aperiódico é a característica mais difícil de usar como marcador: depende da idade, da potência do registro, dos artefatos, da referência e da faixa do ajuste, e o sinal do seu efeito muda com o banco e com o processamento. O offset, o único parâmetro com efeito nos três bancos, mede sobretudo a potência de banda larga. Limpar artefatos muda pouco o efeito de grupo; escolher a referência e a faixa do ajuste muda muito, e essas escolhas precisam ser fixadas e relatadas.

O tabuleiro de variáveis é a contribuição deste trabalho para os seguintes. Ele indica, para cada característica e cada banco, de quais fatores desconfiar. Um classificador generalizável precisará usar essa informação na sua estrutura, ajustando cada família pelas covariáveis que a confundem, em vez de corrigir todas da mesma forma.
