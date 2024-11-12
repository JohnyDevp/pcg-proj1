# PCG projekt 1

-   autor: xlogin00

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

|   N    | CPU [s]  |  GPU [s]  | Zrychlení | Propustnost [GiB/s] | Výkon [GFLOPS] |
| :----: | :------: | :-------: | :-------: | :-----------------: | :------------: |
|  1024  |  1.0928  | 0.029182s |  37.448   |                     |                |
|  2048  |  0.5958  | 0.054568s |  10.918   |                     |                |
|  4096  |  0.6652  | 0.108550s |   6.128   |                     |                |
|  8192  |  1.6599  | 0.208529s |   7.960   |                     |                |
| 16384  |  3.3655  | 0.407879s |   8.251   |                     |                |
| 32768  | 12.7233  | 0.813445s |  15.641   |                     |                |
| 65536  | 48.9732  | 2.720007s |  18.005   |                     |                |
| 131072 | 195.9965 | 8.046670s |  24.357   |                     |                |

## Otázky

### Krok 0: Základní implementace

**Vyskytla se nějaká anomálie v naměřených časech? Pokud ano, vysvětlete:**
Za anomálii můžeme považovat téměř dvouvteřinový rozdíl ve výpočtech s počty částic 69632 a 73728. Ačkoli je nejvíc k více než jednovteřinovým nárůstům dochází v časech vícekrát, tento rozdíl je nejdramatičtější. Všechny rozdíly budou pravděpodobně způsobeny výpadky paměti, kdy pro zmíněné počty částic dojde bude velikost bloků taková, že částice načtené do CACHE v ní nevydrží do další iterace.

### Krok 1: Sloučení kernelů

**Došlo ke zrychlení?**
Ano
**Popište hlavní důvody:**
Došlo ke snížení celkového počtu iterací, kdy v kroku 1 byl proveden dvojnásobný počet těchto iterací.
Zároveň byl optimalizován přístup do paměti při aktualizaci polohy jednotlivých částic a upravení výpočtu, pro minimalizaci nutných výpočetních operací.

### Krok 2: Sdílená paměť

**Došlo ke zrychlení?**

**Popište hlavní důvody:**

### Krok 5: Měření výkonu

**Jakých jste dosáhli výsledků?**

**Lze v datech pozorovat nějaké anomálie?**
