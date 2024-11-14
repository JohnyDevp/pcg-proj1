# PCG projekt 1

-   autor: xholan11

## Měření výkonu (čas / 100 kroků simulace)

### Průběžné

N Time

|   N   | CPU [s]  | Step 0 [s] | Step 1 [s] | Step 2 [s] |
| :---: | -------- | ---------- | ---------- | ---------- |
| 4096  | 0.492139 | 0.179809s  | 0.120320s  | 0.097942s  |
| 8192  | 1.471328 | 0.362047s  | 0.242244s  | 0.194971s  |
| 12288 | 2.478942 | 0.542957s  | 0.363048s  | 0.291997s  |
| 16384 | 3.386801 | 0.724977s  | 0.484454s  | 0.389014s  |
| 20480 | 5.059240 | 0.904686s  | 0.605983s  | 0.486072s  |
| 24576 | 7.112179 | 1.086111s  | 0.728141s  | 0.583040s  |
| 28672 | 9.892856 | 1.268451s  | 0.848479s  | 0.680194s  |
| 32768 | 12.59829 | 1.448769s  | 0.970600s  | 0.777323s  |
| 36864 | 15.54297 | 1.631258s  | 1.092052s  | 0.874449s  |
| 40960 | 19.36099 | 1.812445s  | 1.213129s  | 0.971612s  |
| 45056 | 23.48723 | 1.993434s  | 1.333900s  | 1.068326s  |
| 49152 | 27.69359 | 2.174691s  | 1.455448s  | 1.165578s  |
| 53248 | 32.63063 | 2.355250s  | 1.577456s  | 1.262630s  |
| 57344 | 37.43660 | 4.004539s  | 2.698421s  | 2.320884s  |
| 61440 | 42.85863 | 4.305737s  | 2.905465s  | 2.487418s  |
| 65536 | 49.46104 | 4.600666s  | 3.106917s  | 2.653649s  |
| 69632 | 55.14939 | 4.886697s  | 3.300869s  | 2.819035s  |
| 73728 | 62.04446 | 5.176333s  | 3.494316s  | 2.984891s  |
| 77824 | 69.26138 | 5.463727s  | 3.688949s  | 3.150837s  |
| 81920 | 76.60071 | 5.750830s  | 3.882830s  | 3.316388s  |

### Závěrečné

|   N    | CPU [s]  |  GPU [s]  | Zrychlení | Propustnost [GiB/s] |   Výkon [GFLOPS]   |
| :----: | :------: | :-------: | :-------: | :-----------------: | :----------------: |
|  1024  |  1.0928  | 0.029182s |  37.448   |   0.129345703125    |      109.285       |
|  2048  |  0.5958  | 0.054568s |  10.918   |    0.11357421875    | 221.77457284275442 |
|  4096  |  0.6652  | 0.108550s |   6.128   |     0.105859375     | 446.8835633825511  |
|  8192  |  1.6599  | 0.208529s |   7.960   |   0.104716796875    | 897.0887730363756  |
| 16384  |  3.3655  | 0.407879s |   8.251   |    0.09912109375    | 1784.9229580266401 |
| 32768  | 12.7233  | 0.813445s |  15.641   |   0.098037109375    | 3570.8719969978556 |
| 65536  | 48.9732  | 2.720007s |  18.005   |   0.057099609375    | 4181.721532690339  |
| 131072 | 195.9965 | 8.046670s |  24.357   |   0.038525390625    | 4708.5613117180765 |

## Otázky

### Krok 0: Základní implementace

**Vyskytla se nějaká anomálie v naměřených časech? Pokud ano, vysvětlete:**
Ano. Můžeme pozorovat silně nelineární nárust doby výpočtu pro počty částic N=53248 a N=57344. Ten bude v tomto případě způsobem větším počtem bloků, než je počet dostupných SM procesorů. Tím dojde k zatížení některých SM procesorů vícekrát (zatímco ostatní budou stát) a tím se prodlouží doba výpočtu.

### Krok 1: Sloučení kernelů

**Došlo ke zrychlení?**
Ano

**Popište hlavní důvody:**
Došlo ke snížení celkového počtu iterací, kdy v kroku 1 byl proveden dvojnásobný počet těchto iterací.
Zároveň byl optimalizován přístup do paměti při aktualizaci polohy jednotlivých částic a upravení výpočtu pro minimalizaci nutných výpočetních operací.

### Krok 2: Sdílená paměť

**Došlo ke zrychlení?**
Ano

**Popište hlavní důvody:**
Načtení dat vlákny do sdílené paměti pro každý blok snižuje počet přístupů do globální paměti (v předchozích krocích byla data stále načítána z globální paměti).

### Krok 5: Měření výkonu

**Jakých jste dosáhli výsledků?**
Došlo k velkému zrychlení, což bylo očekávané při správném využití GPU. Největší zrychlení (37x) se projevilo v měřeních při počtu částic N=1024. Největší propustnosti paměti pak při N= a největší výkon byl dosažen s N= .

**Lze v datech pozorovat nějaké anomálie?**
Ano. Jednotlivá zrychlení GPU proti CPU s rostoucím počtem částic N nejdříve klesají a následně začnou stoupat. Toto chování je velmi pravděpodobně způsobené opět různých počtem bloků, kdy při drobném přesáhnutí maximálního počtu SM procesorů dojde k delší době výpočtu,
