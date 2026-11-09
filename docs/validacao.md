# Validação

O código foi conferido contra os resultados do artigo antes da publicação.

- **Pré-processamento.** Para cinco participantes de cada bloco do HBN e do TDBRAIN, `scripts/00_preparar.py --validar` refaz o pacote a partir do registro bruto: mesmo número de amostras e mesmos trechos, correlação de Pearson acima de 0,9999 em todos os canais com sinal e erro máximo da ordem do passo de quantização (menos de 0,06 % do desvio-padrão do canal). O arquivo bruto do adhdata não estava disponível nesta conferência.
- **Características e ensaios.** A partir do pacote do adhdata (121 participantes, todas as variantes de processamento, inclusive a ICA), os efeitos de grupo, as AUC e as tabelas do tabuleiro, das duas metades, do denominador fixo, do fluxo, da distância espacial sem o próprio participante e da impressão digital nos controles reproduzem os valores publicados (diferença absoluta menor que 0,0001). As tabelas que dependem só de agregados (tabuleiro numérico, escore de cuidado) e a de lote do TDBRAIN são idênticas às publicadas. As linhas do HBN e do TDBRAIN das demais tabelas vêm da execução original.
- **Seus dados.** O caminho do pacote às tabelas foi rodado com um bloco de nome livre (`--blocos meu_banco`) e reproduziu os mesmos números, e `tests/test_leitor_generico.py` grava arquivos EDF sintéticos, lê com o leitor genérico e confere canais, taxa, remoção da rede e o formato do pacote.
- **Figuras e documento.** As figuras de `figuras/` reproduzem as do artigo a partir das tabelas agregadas, e `scripts/04_documento.py` remonta o texto do artigo com a mesma numeração de citações.
- **Testes.** `python -m pytest tests` roda quatro testes com sinais sintéticos.
