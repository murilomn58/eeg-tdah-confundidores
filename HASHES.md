# Impressões digitais do código

O pré-registro e os módulos de extração foram congelados com SHA-256 antes das análises. O código deste repositório é uma reescrita compacta desses módulos, validada contra os resultados do artigo (ver [docs/validacao.md](docs/validacao.md)); por isso os hashes do congelamento não se aplicam a estes arquivos. Três arquivos congelados foram editados depois do congelamento, e o desvio foi declarado no pré-registro.

A tabela dá o SHA-256 de cada arquivo `.py` desta versão (1.0.0), como impressão digital. Para conferir:

```bash
sha256sum eegtdah/*.py scripts/*.py scripts/figuras_ilustradas/*.py tests/*.py texto/curto/v4/*.py
```

| Arquivo | SHA-256 |
|---|---|
| `eegtdah/__init__.py` | `46c663463329acd7039bdef0a9145293aad32e04d7f0dd06ca2719c99f0eb014` |
| `eegtdah/caracteristicas.py` | `7d3c1cb69ca68f95fef6205896fc15e1d80b5a87c8a69d1af9cd966c01b251ec` |
| `eegtdah/config.py` | `6eeeda40279b34ef747b2a0226b8c9b96f0a585f2776c06e1ebca4dfbc7522a6` |
| `eegtdah/documento.py` | `7844e2bcd95c24f549c22ccd0848bedc3808eede3ad680127e76a872883ec697` |
| `eegtdah/ensaios.py` | `f340a5af7f294bd422898c55581ad4d9fd99808b3f9089d0b71868dfd35d8e59` |
| `eegtdah/estatistica.py` | `8f45f654af3554876e3be43d9a8742b2bf769ba17aab3ab71e29bbe678d3f00e` |
| `eegtdah/figuras.py` | `67953f0bd84519aaadd4e5bf7ee844468f0022b25e6a5e2ab6c45afbab4b49f3` |
| `eegtdah/preparar.py` | `ce5da8ccd32e23abba2407980d658235874700fbdce68f26cfcbdcfb91949a55` |
| `eegtdah/replicacoes.py` | `b05060f6c5c26e3b9d30060da0151a70cc544c59636271e3681c56f373c2f1f6` |
| `eegtdah/sinais.py` | `b743717f588c1c6a98a61f71855a483ecf3d39a6a1e871dad267aacc038adc65` |
| `scripts/00_preparar.py` | `68dd54528f365a65133b63f794f2a431f682c4934ef6f6dfd6f52ce22d621120` |
| `scripts/01_caracteristicas.py` | `e14b4693bf9db539e35cf25828b29f1d6e8aaaf3b157c6203a70f4d5305aeb80` |
| `scripts/02_ensaios.py` | `bee02d66f612b8b33492eaecc4914d94172beb733cee64397a4596cbe3d90c24` |
| `scripts/03_figuras.py` | `b45c5db3abd618d45fde180e5e58cac19935c2c72ff58cbe5b24f9660b39419c` |
| `scripts/04_documento.py` | `3d5aaea89b234106efcdec8ba44cd8560e5406c54a0f12a2905afe5364cf2da6` |
| `scripts/contar.py` | `2036ae6ed02e82115eaaea51bb163f1e104c010570ed0364bef7c48df912c14d` |
| `scripts/figuras_ilustradas/corrigir_fig5.py` | `957407ef31d5e450d816586f9f57ff82701d23b3be910a7fc3527322f09bc8e8` |
| `scripts/verificar_publicavel.py` | `307585d510e055a056fc8d62e8122956d4e15a6baf5391ef16cbcfa0cca6689b` |
| `tests/test_basico.py` | `e7e1d4f560f5e254108ef3e0db2f3113e5b4ecdf5f62916877b12d7abaa7ae19` |
| `tests/test_leitor_generico.py` | `0aabdd76323757a46ec3e16b4235738cf2463966976bb7adf8a168372af9af25` |
| `texto/curto/v4/montar_v4.py` | `7a45518b1ff913de8bbe4c28a6932187533b6e21efdd61fbef007506483d4ff7` |
