# Semana 12: Análisis metataxonómico de comunidades bacterianas con datos de secuenciación Illumina

## Logro de la sesión:

Al finalizar la sesión, el estudiante utiliza herramientas bioinformáticas para procesar lecturas Illumina del gen 16S rRNA, inferir variantes de secuencia de amplicón (ASVs), asignarles taxonomía y comparar la composición y la diversidad de comunidades bacterianas acuáticas a lo largo de un gradiente altitudinal.

## Estructura de la práctica:

1. Acceso al servidor de cómputo
2. Preparación del entorno y de la base de datos
3. Importación y control de calidad de las lecturas
4. Recorte de los primers con cutadapt
5. Inferencia de ASVs con DADA2
6. Asignación taxonómica con SILVA 138.2
7. Filtrado de ASVs y composición taxonómica
8. Árbol filogenético y exportación de los datos
9. Análisis de diversidad en MicrobiomeAnalyst

> **Cómo está organizada esta práctica:** se analizan los datos públicos de Romero *et al.* (2021, *PLoS ONE* 16: e0250401; BioProject PRJNA646070): muestras de agua superficial del río Rímac tomadas en enero de 2020 en 13 localidades de seis subcuencas, desde Chicla (3 441 m s. n. m.) hasta el puente Gambetta, en el Callao (8 m s. n. m.). Se secuenció la región V3-V4 del gen 16S rRNA (primers 341F/805R) en un MiSeq, con lecturas pareadas de 2 × 251 nt. El set tiene **24 muestras**: una por localidad (01 a 13) y una réplica de 11 de ellas (14 a 24); Chacahuaro y Universitaria no tienen réplica.
>
> Todos los grupos trabajan con el mismo set de datos. El denoising con DADA2 (sección 5) tarda más de seis hora con los datos completos, por lo que **no se ejecuta en clase**: se copian sus salidas desde la carpeta compartida.

## Flujo de trabajo:

### Procesamiento de las lecturas en el servidor (secciones 2 a 8):

```mermaid
flowchart LR
    subgraph Q2["QIIME 2 (servidor)"]
        direction LR
        A["Lecturas pareadas FASTQ (manifest.txt)"] --> B["tools import (demuxed_seqs.qza)"]
        B --> C["cutadapt (recorte de primers 341F/805R)"]
        C --> D["DADA2 (ASVs, unión de pares, quimeras)"]
        D --> E["table.qza / representatives.qza"]
        E --> F["classify-sklearn (SILVA 138.2 V3-V4)"]
        F --> G["taxa filter (sin Unassigned, mitocondrias ni cloroplastos)"]
        G --> H["taxa barplot"]
        G --> I["MAFFT + FastTree (rooted_tree)"]
        G --> J["export + biom add-metadata (feature-table-tax.biom)"]
    end
```

### Análisis de diversidad en línea (sección 9):

```mermaid
flowchart LR
    subgraph MA["MicrobiomeAnalyst (en línea)"]
        direction LR
        A["feature-table-tax.biom + metadata_MA.txt + rooted_tree.nwk"] --> B["Data Integrity Check"]
        B --> C["Data Filtering"]
        C --> D["Normalization"]
        D --> E["Composición (barras apiladas)"]
        D --> F["Diversidad alfa (Chao1, Shannon)"]
        D --> G["Diversidad beta (PCoA, NMDS, Bray-Curtis)"]
    end
```

## Programas requeridos:

### Programas de acceso al servidor:

PuTTY v0.79 https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html
   - **Descripción:** cliente SSH para establecer conexiones de línea de comandos con el servidor.

WinSCP v6.1 https://winscp.net/eng/download.php
   - **Descripción:** cliente SFTP con interfaz gráfica para transferir archivos entre su computadora y el servidor. En esta práctica se usa para descargar las visualizaciones (`.qzv`) y los archivos que se cargan en MicrobiomeAnalyst.

### Programas bioinformáticos:

QIIME 2 amplicon v2024.10 https://qiime2.org/
   - **Descripción:** QIIME 2 (*Quantitative Insights Into Microbial Ecology*) es una plataforma para analizar datos de amplicones de principio a fin. Trabaja con **artefactos** (`.qza`), que guardan los datos junto con su tipo semántico, y **visualizaciones** (`.qzv`). Cada archivo registra su procedencia: qué comando, con qué parámetros y a partir de qué archivos se generó. En esta práctica se usan los plugins `cutadapt` (recorte de primers), `dada2` (inferencia de ASVs), `feature-classifier` (asignación taxonómica), `taxa` (filtrado y gráficos de barras) y `phylogeny` (MAFFT y FastTree).

biom-format v2 https://biom-format.org/
   - **Descripción:** formato estándar para tablas de abundancia de microbiomas (*Biological Observation Matrix*). El comando `biom`, instalado en el entorno de QIIME 2, se usa para añadir la taxonomía a la tabla de abundancias.

### Bases de datos:

Clasificador Naive Bayes SILVA 138.2 V3-V4 (341F/805R)
   - **Descripción:** SILVA es la base de datos de referencia de genes de ARNr. El clasificador se entrenó con scikit-learn sobre SILVA 138.2 (SSU no redundante al 99 %), curado con RESCRIPt y recortado *in silico* a la región que amplifican los primers 341F/805R. Como solo conoce la región secuenciada, clasifica más rápido, usa menos memoria y suele resolver mejor el género que un clasificador de secuencias completas.

### Herramientas bioinformáticas en línea:

QIIME 2 View https://view.qiime2.org/
   - **Descripción:** visor de los archivos `.qzv` en el navegador. Los datos no se suben a ningún servidor: se procesan en su computadora.

MicrobiomeAnalyst v2.0 https://www.microbiomeanalyst.ca/
   - **Descripción:** plataforma web para el análisis estadístico y visual de datos de microbioma. Su módulo *Marker Data Profiling* (MDP) recibe una tabla de abundancias con taxonomía, la metadata y, opcionalmente, un árbol filogenético, y calcula composición, diversidad alfa y beta, y abundancia diferencial.

## Metodología:

## 1. Acceso al servidor de cómputo:

### Abrir el programa PuTTY, colocar el hostname ( 10.142.250.66 ) y port ( 22 ), y dar clic en Open:

<img width="500" alt="image" src="https://github.com/user-attachments/assets/92f89dbb-1a21-411d-adb5-38fe486a5567" />

**Figura 1.** Ventana de configuración de PuTTY. En *Host Name* se coloca la dirección IP del servidor y en *Port* el puerto 22 (SSH); el tipo de conexión debe ser SSH.

> **Comentario (figura 1):** antes de abrir la sesión, vaya a **Connection** y en *Seconds between keepalives* escriba **30**. PuTTY enviará una señal cada 30 segundos y la conexión no se cortará por inactividad durante los pasos largos. Si la conexión se cae a mitad de un comando, ese comando se interrumpe y hay que volver a ejecutarlo.

### En la terminal abierta, escribir su usuario y contraseña correspondiente para tener acceso al servidor de cómputo:

<img width="700" alt="image" src="https://github.com/user-attachments/assets/4d246e93-c59c-4749-a2dd-03db25c53654" />

**Figura 2.** Terminal de PuTTY al iniciar la sesión. Al escribir la contraseña no se muestran caracteres en pantalla; es el comportamiento normal. Una vez dentro, la línea de comandos indica su usuario y el nombre del servidor.

### Abrir el programa WinSCP, colocar el hostname ( 10.142.250.66 ) y port ( 22 ), escribir su usuario y contraseña correspondiente para tener acceso al servidor de cómputo, y hacer clic en Login:

<img width="500" alt="image" src="https://github.com/user-attachments/assets/ef4dc253-ce4a-417d-b761-39692d2a011a" />

**Figura 3.** Ventana de inicio de sesión de WinSCP, con el protocolo SFTP, la dirección del servidor, el puerto 22, el usuario y la contraseña.

<img width="500" alt="image" src="https://github.com/user-attachments/assets/577debca-6085-47c5-9bbd-73688bfa8bb0" />

**Figura 4.** Sesión de WinSCP abierta. El panel izquierdo muestra las carpetas de su computadora y el derecho, las del servidor; los archivos se transfieren arrastrándolos de un panel al otro.

### Uso compartido del servidor

El servidor tiene 64 hilos y 180 GB de RAM, y lo usan **6 grupos a la vez**. Cada grupo usa como máximo **10 hilos** (6 × 10 = 60; quedan 4 para el sistema). No use valores mayores que los indicados en los comandos.

| Sección | Paso | Tiempo aproximado (10 hilos) |
|---|---|---|
| 3 | Importación y reporte de calidad | 5–10 min |
| 4 | cutadapt y reporte de calidad | 5–15 min |
| 5 | DADA2 | más de 1 h: **no se ejecuta**, se copian las salidas |
| 6 | Clasificación taxonómica | 5–20 min |
| 8 | Árbol filogenético (MAFFT + FastTree) | 5–15 min |

> **Comentario:** los tiempos dependen de la carga del servidor. Si un paso no termina a tiempo, copie su resultado desde la carpeta compartida `/data/2026_1/metagenomics/illumina/qiime2/`.

## 2. Preparación del entorno y de la base de datos

### Crear la estructura de carpetas de trabajo

```bash
cd

mkdir metataxonomic

cd metataxonomic

mkdir tmp
```

> **Comentario:**
> - `cd`: sin argumentos, lleva a su carpeta personal (`/home/<usuario>`).
> - `mkdir metataxonomic`: crea la carpeta de trabajo de la práctica.
> - `cd metataxonomic`: entra a esa carpeta.
> - `mkdir tmp`: crea una carpeta para los archivos temporales de QIIME 2.
> - Estos comandos no muestran nada en pantalla cuando funcionan correctamente.

```bash
export TMPDIR=$HOME/metataxonomic/tmp
```

> **Comentario:** QIIME 2 escribe archivos temporales grandes durante cada análisis. La variable `TMPDIR` los dirige a la carpeta `tmp` del usuario y no a `/tmp`, que comparten los 6 grupos y puede llenarse. `$HOME` se reemplaza por la carpeta personal de cada usuario. La variable se pierde al cerrar la terminal: si abre una nueva, ejecute de nuevo esta línea.

### Copiar el clasificador taxonómico

```bash
mkdir database

cd database

cp /data/2026_1/metagenomics/illumina/database/silva138.2_V3V4/silva-138.2-ssu-nr99-V3V4-341f-805r-classifier.qza .
```

> **Comentario:**
> - `mkdir database` y `cd database`: crea la carpeta de la base de datos y entra a ella.
> - `cp ... .`: copia el clasificador a la carpeta actual. El punto final (`.`) significa "aquí".
> - El clasificador se entrenó en el mismo entorno de QIIME 2 (2024.10) que se usa en la práctica. Esto es obligatorio: un clasificador entrenado con otra versión de scikit-learn falla al cargarse.

> **Nota sobre los nombres de los filos:** SILVA 138.2 usa la nomenclatura actual (*Pseudomonadota*, *Bacillota*, *Actinomycetota*, *Campylobacterota*). El artículo de Romero *et al.* usó SILVA 138, con nombres en parte distintos (*Proteobacteria*, *Firmicutes*). Téngalo en cuenta al comparar sus resultados con el artículo.

## 3. Importación y control de calidad de las lecturas

### Obtener los archivos manifest.txt y metadata.txt

```bash
cd ~/metataxonomic

mkdir illumina

cd illumina

conda activate qiime2-amplicon-2024.10
```

> **Comentario:**
> - `cd ~/metataxonomic`: vuelve a la carpeta de trabajo (`~` es su carpeta personal).
> - `mkdir illumina` y `cd illumina`: crea la carpeta del análisis y entra a ella. Todos los comandos de las secciones 3 a 8 se ejecutan aquí.
> - `conda activate qiime2-amplicon-2024.10`: activa el entorno donde está instalado QIIME 2. El nombre del entorno activo aparece entre paréntesis al inicio de la línea de comandos. Si abre una terminal nueva, actívelo otra vez.

```bash
qiime --version

q2cli version 2024.10.1
Run `qiime info` for more version details.
```

> **Comentario:** muestra la versión de QIIME 2. Anótela para la sección de metodología de la bitácora.

```bash
cp /data/2026_1/metagenomics/illumina/qiime2/manifest.txt .

cp /data/2026_1/metagenomics/illumina/qiime2/metadata.txt .
```

> **Comentario:**
> - `manifest.txt`: indica a QIIME 2 qué par de archivos FASTQ (forward y reverse) corresponde a cada muestra, con su ruta absoluta.
> - `metadata.txt`: contiene las variables de cada muestra: la subcuenca (`cuenca`) y la altitud en m s. n. m. (`altitud`).
> - Los archivos se copian ya preparados: al pegar texto en `nano` desde PuTTY, los tabuladores suelen convertirse en espacios y QIIME 2 no puede leer el archivo.

```bash
cat -A manifest.txt | head -n 3

sample-id^Iforward-absolute-filepath^Ireverse-absolute-filepath$
01.Chicla^I/data/2026_1/metagenomics/illumina/raw_data/SRR12217142_1.fastq.gz^I/data/2026_1/metagenomics/illumina/raw_data/SRR12217142_2.fastq.gz$
02.San_Mateo^I/data/2026_1/metagenomics/illumina/raw_data/SRR12217141_1.fastq.gz^I/data/2026_1/metagenomics/illumina/raw_data/SRR12217141_2.fastq.gz$
```

> **Comentario:**
> - **Comando:** `cat -A` muestra el archivo con sus caracteres invisibles: cada tabulador aparece como `^I` y cada final de línea como `$`. `| head -n 3` muestra solo las 3 primeras líneas.
> - **Resultado:** el encabezado tiene las tres columnas que exige el formato `PairedEndFastqManifestPhred33V2`: `sample-id`, `forward-absolute-filepath` y `reverse-absolute-filepath`, separadas por `^I`. Cada muestra tiene dos archivos: `_1.fastq.gz` (lectura forward) y `_2.fastq.gz` (lectura reverse). El nombre `SRR...` es el código de acceso de la corrida en el SRA del NCBI.

```bash
cat -A metadata.txt | head -n 3

Sample ID^Icuenca^Ialtitud$
01.Chicla^IUpper_Rimac^I3441$
02.San_Mateo^IParac-Upper_Rimac^I3422$
```

> **Comentario:** la metadata tiene tres columnas separadas por tabuladores. QIIME 2 acepta varios nombres para la primera columna (`Sample ID`, `sample-id`, `#SampleID`); lo importante es que los identificadores coincidan exactamente con los del manifest.

Correspondencia entre muestras, corridas del SRA y variables (contenido de `manifest.txt` y `metadata.txt`):

| Muestra | Réplica | Corrida SRA (muestra / réplica) | Subcuenca (`cuenca`) | Altitud (m) |
|---|---|---|---|---|
| 01.Chicla | 14.Chicla | SRR12217142 / SRR13224180 | Upper_Rimac | 3 441 |
| 02.San_Mateo | 15.San_Mateo | SRR12217141 / SRR13224179 | Parac-Upper_Rimac | 3 422 |
| 03.Tamboraque | 16.Tamboraque | SRR12217133 / SRR13224178 | Parac-Upper_Rimac | 3 023 |
| 04.Huanchor | 17.Huanchor | SRR12217132 / SRR13224177 | Santa_Eulalia-Parac | 2 868 |
| 05.Chacahuaro | – | SRR12217130 | Santa_Eulalia-Parac | 2 526 |
| 06.Santa_Eulalia | 18.Santa_Eulalia | SRR12217131 / SRR13224176 | Santa_Eulalia | 948 |
| 07.Chaclacayo | 19.Chaclacayo | SRR12217129 / SRR13224175 | Jicamarca-Santa_Eulalia | 657 |
| 08.Huachipa | 20.Huachipa | SRR12217128 / SRR13224174 | Jicamarca-Santa_Eulalia | 396 |
| 09.Libertadores | 21.Libertadores | SRR12217127 / SRR13224173 | Lower_Rimac | 213 |
| 10.Nuevo | 22.Nuevo | SRR12217126 / SRR13224172 | Lower_Rimac | 200 |
| 11.Universitaria | – | SRR13224171 | Lower_Rimac | 69 |
| 12.Faucett | 23.Faucett | SRR12217139 / SRR12217140 | Lower_Rimac | 33 |
| 13.Gambetta | 24.Gambetta | SRR12217138 / SRR13224170 | Lower_Rimac | 8 |

### Importar las lecturas a QIIME 2

```bash
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.txt --input-format PairedEndFastqManifestPhred33V2 --output-path demuxed_seqs.qza

Imported manifest.txt as PairedEndFastqManifestPhred33V2 to demuxed_seqs.qza
```

> **Comentario (comando):**
> - `--type 'SampleData[PairedEndSequencesWithQuality]'`: tipo semántico del artefacto: lecturas pareadas, con calidad, separadas por muestra.
> - `--input-path manifest.txt`: archivo con la ruta de los FASTQ de cada muestra.
> - `--input-format PairedEndFastqManifestPhred33V2`: formato del manifest. *PairedEnd*, lecturas pareadas; *Phred33*, la calidad de cada base está codificada como un carácter ASCII con desplazamiento 33 (el estándar de Illumina actual); *V2*, manifest separado por tabuladores.
> - `--output-path demuxed_seqs.qza`: artefacto de salida con las 24 muestras.
>
> **Comentario (resultado):** QIIME 2 confirma que importó el manifest y creó `demuxed_seqs.qza`. Si alguna ruta del manifest no existe, el comando se detiene e indica qué archivo de qué muestra no encontró.

```bash
qiime tools peek demuxed_seqs.qza

UUID:        0672b33e-4680-4064-8ac8-229fd820335f
Type:        SampleData[PairedEndSequencesWithQuality]
Data format: SingleLanePerSamplePairedEndFastqDirFmt
```

> **Comentario:**
> - **Comando:** `qiime tools peek` muestra la información básica de un artefacto sin abrirlo.
> - **Resultado:** `UUID` es el identificador único del artefacto (el suyo será distinto: cada ejecución genera uno nuevo); `Type`, su tipo semántico; y `Data format`, el formato interno de los datos: un par de archivos FASTQ por muestra, de un solo carril de secuenciación.

### Crear la carpeta de visualizaciones y generar el reporte de calidad

```bash
mkdir visualization
```

> **Comentario:** crea la carpeta donde se guardarán todas las visualizaciones (`.qzv`) de la práctica.

```bash
qiime demux summarize --i-data demuxed_seqs.qza --o-visualization visualization/seqs_quality.qzv

Saved Visualization to: visualization/seqs_quality.qzv
```

> **Comentario:**
> - `--i-data demuxed_seqs.qza`: lecturas importadas.
> - `--o-visualization visualization/seqs_quality.qzv`: reporte de salida.
> - El gráfico de calidad se calcula sobre una submuestra aleatoria de 10 000 lecturas (opción `--p-n`, que se deja en su valor por defecto).

### Descargar con WinSCP el archivo `seqs_quality.qzv` y abrirlo en QIIME 2 View (https://view.qiime2.org/)

> **[Captura pendiente: pestaña *Overview* de `seqs_quality.qzv`]**

**Figura 5.** Pestaña *Overview* de `seqs_quality.qzv`: resumen del número de lecturas por muestra y tabla con los conteos de lecturas forward y reverse de cada muestra.

> **Comentario (figura 5):**
> - Las 24 muestras suman **5 865 440 pares de lecturas**. El número de lecturas forward y reverse es idéntico en cada muestra, como corresponde a datos pareados.
> - La muestra con menos lecturas es **11.Universitaria (105 680)** y la que tiene más, **21.Libertadores (291 994)**; el resto tiene entre 175 381 y 287 437. Son los mismos valores que reporta el artículo (promedio de 244 393 lecturas por muestra, rango de 105 680 a 287 437).
> - Universitaria tiene menos de la mitad de lecturas que las demás. Se volverá a ver en todos los pasos siguientes.

> **[Captura pendiente: pestaña *Interactive Quality Plot* de `seqs_quality.qzv`]**

**Figura 6.** Pestaña *Interactive Quality Plot*: diagramas de cajas de la calidad Phred en cada posición de las lecturas forward (izquierda) y reverse (derecha).

> **Comentario (figura 6):**
> - Cada caja resume la calidad de las 10 000 lecturas submuestreadas en una posición: la línea central es la mediana y los extremos, los percentiles 25 y 75. Una calidad Phred de 30 equivale a una probabilidad de error de 1 en 1 000; una de 20, a 1 en 100.
> - Las lecturas miden **251 nt** (secuenciación de 2 × 251 en MiSeq).
> - **Forward:** la mediana se mantiene en 37–38 hasta cerca del final; en las últimas posiciones baja y en la posición 251 cae a 27.
> - **Reverse:** la mediana se mantiene en 37–38 hasta la posición 220, baja a 32–33 entre las posiciones 230 y 240 y cae a 13 en la última posición.
> - Es el patrón típico de Illumina: la calidad disminuye hacia el extremo 3′ de cada lectura, y más en la reverse.

> **Punto de control:** Anote el número de lecturas de cada muestra y la posición en la que empieza a caer la calidad de cada lectura. Los usará para explicar los resultados de DADA2.

## 4. Recorte de los primers con cutadapt

Las lecturas empiezan con la secuencia del primer, que no es información biológica de la muestra: proviene del oligonucleótido sintético usado en la PCR. Además, los primers tienen posiciones degeneradas que DADA2 interpretaría como variación real. Por eso se eliminan antes del denoising.

```bash
qiime cutadapt trim-paired --i-demultiplexed-sequences demuxed_seqs.qza --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACHVGGGTATCTAATCC --p-minimum-length 100 --p-discard-untrimmed --p-cores 10 --o-trimmed-sequences demuxed_seqs_trimmed.qza

Saved SampleData[PairedEndSequencesWithQuality] to: demuxed_seqs_trimmed.qza
```

> **Comentario (comando):**
> - `--i-demultiplexed-sequences demuxed_seqs.qza`: lecturas pareadas de entrada.
> - `--p-front-f CCTACGGGNGGCWGCAG`: primer forward 341F. *front* indica que se busca en el extremo 5′ de la lectura forward, y se elimina junto con todo lo que esté antes de él. Las letras N, W, H y V son códigos IUPAC de bases degeneradas (N = cualquier base; W = A o T; H = A, C o T; V = A, C o G).
> - `--p-front-r GACTACHVGGGTATCTAATCC`: primer reverse 805R, buscado en el extremo 5′ de la lectura reverse.
> - `--p-minimum-length 100`: descarta las lecturas que quedan con menos de 100 nt después del recorte.
> - `--p-discard-untrimmed`: descarta las lecturas en las que no se encontró el primer (productos inespecíficos o artefactos).
> - `--p-cores 10`: número de hilos.
> - `--o-trimmed-sequences demuxed_seqs_trimmed.qza`: lecturas recortadas.
> - Cutadapt acepta por defecto hasta un 10 % de diferencias entre el primer y la lectura (`--p-error-rate 0.1`).
> - En el protocolo Illumina 16S, la lectura empieza en el primer específico y no en la cola del adaptador; por eso basta con buscar el primer.

```bash
qiime demux summarize --i-data demuxed_seqs_trimmed.qza --o-visualization visualization/seqs_quality_trimmed.qzv

Saved Visualization to: visualization/seqs_quality_trimmed.qzv
```

> **Comentario:** genera el mismo reporte de calidad, ahora con las lecturas recortadas, para comparar el antes y el después.

### Descargar `seqs_quality_trimmed.qzv` y compararlo con `seqs_quality.qzv` en QIIME 2 View

> **[Captura pendiente: pestaña *Interactive Quality Plot* de `seqs_quality_trimmed.qzv`]**

**Figura 7.** Calidad por posición de las lecturas recortadas. Las lecturas son más cortas que las originales y la caja de la última posición tiene muy pocas lecturas.

> **Comentario (figura 7):**
> - Se conservan **5 793 891 de 5 865 440 pares (98,8 %)**. En ninguna muestra se pierde más del 2 %: los primers estaban en casi todas las lecturas, y `--p-discard-untrimmed` casi no descartó nada. Si se perdiera más del 20–30 %, habría que revisar la secuencia de los primers.
> - Tras el recorte, la lectura forward mide **~234 nt** (251 − 17 nt del primer 341F) y la reverse, **~230 nt** (251 − 21 nt del primer 805R).
> - **Solapamiento esperado:** el amplicón V3-V4 sin primers mide unos 400–430 nt. Las dos lecturas suman 464 nt, así que se superponen solo en unos 35–60 nt. DADA2 necesita al menos 12 nt de solapamiento para unirlas, y en esa zona están las bases de peor calidad: el final de cada lectura.

> **Punto de control:** ¿Cuántas lecturas se perdieron en cada muestra con el recorte? Con la longitud del amplicón y de las lecturas recortadas, calcule cuántos nucleótidos se superponen. ¿Qué pasaría si se truncaran las lecturas en las posiciones 200 y 200?

## 5. Inferencia de ASVs con DADA2

DADA2 aprende de los propios datos un modelo de los errores de la corrida y lo usa para corregir cada lectura. El resultado son **variantes de secuencia de amplicón (ASVs)**: secuencias exactas, que pueden diferir de otras en un solo nucleótido, en lugar de los OTUs agrupados al 97 % que se usaban antes.

### Comando usado (no se ejecuta en clase)

```bash
qiime dada2 denoise-paired --i-demultiplexed-seqs demuxed_seqs_trimmed.qza --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads 20 --o-table table.qza --o-representative-sequences representatives.qza --o-denoising-stats denoising_stats.qza
```

> **Comentario (comando):**
> - `--i-demultiplexed-seqs demuxed_seqs_trimmed.qza`: lecturas sin primers.
> - `--p-trunc-len-f 0` y `--p-trunc-len-r 0`: posición en la que se corta cada lectura antes del análisis. Con `0` no se truncan.
> - `--p-n-threads 20`: número de hilos.
> - `--o-table table.qza`: tabla de abundancias (cuántas lecturas tiene cada ASV en cada muestra).
> - `--o-representative-sequences representatives.qza`: la secuencia de cada ASV.
> - `--o-denoising-stats denoising_stats.qza`: cuántas lecturas sobreviven en cada etapa.
> - **Opciones que no se escriben** y quedan en su valor por defecto:
>   - `--p-max-ee-f 2` y `--p-max-ee-r 2`: descarta las lecturas con más de 2 errores esperados, calculados a partir de la calidad de cada base. Es el mismo filtro (`maxEE = 2`) que usaron Romero *et al.*
>   - `--p-trunc-q 2`: corta la lectura en la primera base con calidad ≤ 2.
>   - `--p-min-overlap 12`: solapamiento mínimo para unir las lecturas forward y reverse.
>   - `--p-pooling-method independent`: cada muestra se procesa por separado.
>   - `--p-chimera-method consensus`: método de detección de quimeras.
>   - `--p-n-reads-learn 1000000`: número de lecturas usadas para aprender el modelo de errores. Es la etapa más lenta.
> - Las etapas de DADA2 son: filtrado por calidad → aprendizaje del modelo de errores → denoising de cada muestra → unión de pares → eliminación de quimeras.

### Copiar las salidas de DADA2 a la carpeta de trabajo

```bash
cp /data/2026_1/metagenomics/illumina/qiime2/table.qza .

cp /data/2026_1/metagenomics/illumina/qiime2/representatives.qza .

cp /data/2026_1/metagenomics/illumina/qiime2/denoising_stats.qza .
```

> **Comentario:** copia a la carpeta actual (`~/metataxonomic/illumina`) las tres salidas de DADA2, generadas con el comando anterior.

### Generar los reportes del denoising

```bash
qiime metadata tabulate --m-input-file denoising_stats.qza --o-visualization visualization/denoising_stats.qzv

Saved Visualization to: visualization/denoising_stats.qzv
```

> **Comentario:**
> - `--m-input-file denoising_stats.qza`: estadísticas de DADA2. `metadata tabulate` convierte en una tabla interactiva cualquier artefacto que pueda leerse como metadata.
> - `--o-visualization visualization/denoising_stats.qzv`: tabla de salida.

```bash
qiime feature-table summarize --i-table table.qza --m-sample-metadata-file metadata.txt --o-visualization visualization/table_denoised.qzv

Saved Visualization to: visualization/table_denoised.qzv
```

> **Comentario:**
> - `--i-table table.qza`: tabla de abundancias de DADA2.
> - `--m-sample-metadata-file metadata.txt`: metadata de las muestras; en la pestaña *Interactive Sample Detail* permite ver cuántas muestras de cada subcuenca se conservarían con una profundidad de muestreo dada.
> - `--o-visualization visualization/table_denoised.qzv`: reporte de salida.

```bash
qiime feature-table tabulate-seqs --i-data representatives.qza --o-visualization visualization/representatives.qzv

Saved Visualization to: visualization/representatives.qzv
```

> **Comentario:**
> - `--i-data representatives.qza`: secuencias de los ASVs.
> - `--o-visualization visualization/representatives.qzv`: tabla con la secuencia y la longitud de cada ASV.

### Descargar los tres reportes y abrirlos en QIIME 2 View

> **[Captura pendiente: `denoising_stats.qzv`]**

**Figura 8.** Tabla de `denoising_stats.qzv`: número de lecturas de cada muestra al final de cada etapa de DADA2 y porcentaje respecto a la entrada.

> **Comentario (figura 8):**
> - Columnas: `input` (pares de entrada), `filtered` (pasan el filtro de calidad), `denoised` (quedan tras el denoising), `merged` (forward y reverse se unieron), `non-chimeric` (no son quiméricas) y los porcentajes de `filtered`, `merged` y `non-chimeric` respecto a la entrada. Las primeras filas:
>
> | Muestra | input | filtered | denoised | merged | non-chimeric | % non-chimeric |
> |---|---|---|---|---|---|---|
> | 01.Chicla | 269 928 | 191 317 | 180 865 | 121 477 | 97 653 | 36,18 |
> | 02.San_Mateo | 280 021 | 199 324 | 187 378 | 127 848 | 103 896 | 37,10 |
> | 03.Tamboraque | 268 032 | 190 905 | 179 133 | 118 605 | 88 390 | 32,98 |
> | 11.Universitaria | 103 981 | 69 887 | 65 344 | 45 056 | 35 456 | 34,10 |
>
> - **Filtrado:** pasa el **67–75 %** de los pares (72 % en promedio). Como no se truncaron las lecturas, las colas de baja calidad acumulan errores esperados y muchas lecturas superan el límite de `maxEE = 2`.
> - **Unión de pares:** es la etapa con **mayor pérdida**: solo se une el 67 % de los pares que llegan a ella, y queda el 38–50 % de la entrada. El solapamiento es corto (35–60 nt) y cae justo en el extremo de peor calidad de las dos lecturas: las diferencias en esa zona impiden unirlas.
> - **Quimeras:** se elimina alrededor del 23 % de las lecturas unidas.
> - **Resultado final:** se conserva el **30–41 %** de los pares de cada muestra (35 % en promedio), con un promedio de **83 966 lecturas por muestra**; 11.Universitaria es la de menor profundidad (35 456).
> - Los valores de `input` pueden diferir en unas pocas lecturas de los de `seqs_quality_trimmed.qzv`: las salidas que se copiaron se generaron con una versión anterior del recorte, que incluía la cola del adaptador.

> **[Captura pendiente: pestaña *Overview* de `table_denoised.qzv`]**

**Figura 9.** Pestaña *Overview* de `table_denoised.qzv`: número de muestras, número de ASVs (*features*), frecuencia total y distribución de lecturas por muestra y por ASV.

> **Comentario (figura 9):**
> - La tabla tiene **24 muestras y 13 298 ASVs**, con **2 015 186 lecturas** en total.
> - El histograma de frecuencia por ASV está muy sesgado: unos pocos ASVs concentran la mayoría de las lecturas y la mayor parte aparece pocas veces. Es lo habitual en comunidades microbianas.
> - Como referencia, Romero *et al.* obtuvieron 92 064 lecturas finales por muestra en promedio y 23 682 ASVs: procesaron los datos con DADA2 en R, juntaron las lecturas de cada muestra con su réplica y usaron otros parámetros de truncado.

> **[Captura pendiente: `representatives.qzv`]**

**Figura 10.** Tabla de `representatives.qzv`: estadísticas de longitud de los ASVs y lista de secuencias.

> **Comentario (figura 10):**
> - Los 13 298 ASVs miden entre **235 y 452 nt**, con una mediana de **422 nt**. La distribución tiene dos grupos, uno de ~402–405 nt y otro de ~422–428 nt: la longitud de la región V3-V4 varía entre linajes bacterianos.
> - Los ASVs mucho más cortos que el resto (menos de ~380 nt) son sospechosos: pueden ser amplificaciones inespecíficas o artefactos. Se revisan con la taxonomía (sección 6).
> - Al hacer clic en una secuencia se puede enviar directamente a BLAST del NCBI.

> **Punto de control:** ¿En qué etapa de DADA2 se pierden más lecturas? Relacione esa pérdida con la calidad de las lecturas (figura 6) y con el solapamiento que calculó en la sección 4. ¿Qué parámetro de DADA2 cambiaría para recuperar más lecturas, y qué riesgo tendría ese cambio?

## 6. Asignación taxonómica con SILVA 138.2

```bash
qiime feature-classifier classify-sklearn --i-classifier ~/metataxonomic/database/silva-138.2-ssu-nr99-V3V4-341f-805r-classifier.qza --i-reads representatives.qza --p-n-jobs 3 --p-reads-per-batch 2000 --o-classification taxa.qza

Saved FeatureData[Taxonomy] to: taxa.qza
```

> **Comentario (comando):**
> - `--i-classifier`: clasificador entrenado (sección 2).
> - `--i-reads representatives.qza`: secuencias de los ASVs.
> - `--p-n-jobs 3`: número de procesos en paralelo. Cada proceso carga su propia copia del clasificador (~6 GB de RAM); con 3 procesos por grupo y 6 grupos, el servidor usa unos 110 GB de sus 180 GB. Por eso no se usan más.
> - `--p-reads-per-batch 2000`: número de secuencias que procesa cada proceso por vez; valores menores reducen la memoria que se usa.
> - `--o-classification taxa.qza`: taxonomía asignada a cada ASV.
> - **Opción por defecto:** `--p-confidence 0.7`. El clasificador predice el linaje nivel por nivel, desde dominio hasta género, y solo asigna un nivel si la confianza es de al menos 0,7; si no la alcanza, la clasificación se detiene en el nivel anterior.

```bash
qiime metadata tabulate --m-input-file taxa.qza --o-visualization visualization/taxa.qzv

Saved Visualization to: visualization/taxa.qzv
```

> **Comentario:** convierte la taxonomía en una tabla interactiva, con una fila por ASV.

### Descargar `taxa.qzv` y abrirlo en QIIME 2 View

> **[Captura pendiente: `taxa.qzv`]**

**Figura 11.** Tabla de `taxa.qzv`: identificador de cada ASV (`Feature ID`), su linaje (`Taxon`) y la confianza de la asignación (`Confidence`).

> **Comentario (figura 11):**
> - `Feature ID` es un código MD5 calculado a partir de la secuencia del ASV: dos ASVs idénticos tienen el mismo código en cualquier estudio.
> - `Taxon` lista los niveles separados por punto y coma, con un prefijo que indica el rango: `d__` dominio, `p__` filo, `c__` clase, `o__` orden, `f__` familia, `g__` género y `s__` especie.
> - El clasificador se entrenó sin etiquetas de especie, por lo que `s__` siempre queda vacío: con un fragmento de ~420 nt del 16S rara vez se puede distinguir una especie de otra.
> - Use el buscador de la tabla para ubicar géneros de interés (por ejemplo, *Arcobacter* o *Flavobacterium*).

### Exportar la taxonomía y cuantificar la resolución de la asignación

```bash
qiime tools export --input-path taxa.qza --output-path exported_taxonomy

Exported taxa.qza as TSVTaxonomyDirectoryFormat to directory exported_taxonomy
```

> **Comentario:**
> - `--input-path taxa.qza`: artefacto a exportar.
> - `--output-path exported_taxonomy`: carpeta de salida, que contendrá el archivo `taxonomy.tsv`.
> - `qiime tools export` extrae los datos de un artefacto en un formato estándar que otros programas pueden leer.

```bash
head -n 3 exported_taxonomy/taxonomy.tsv

Feature ID      Taxon   Confidence
dcd4030bd55a4bbc6d5d935b8ba2797d        d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Arcobacteraceae;g__Arcobacter;s__    0.9999999995602877
6db4c5fe1f48ff886ffeb57d1dfabb2d        d__Bacteria;p__Fusobacteriota;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;g__Hypnocyclicus;s__       0.9999999999862439
```

> **Comentario:**
> - **Comando:** `head -n 3` muestra las primeras 3 líneas: el encabezado y dos ASVs.
> - **Resultado:** el primer ASV es un *Arcobacter* (Campylobacterota) y el segundo, un *Hypnocyclicus* (Fusobacteriota), ambos con una confianza cercana a 1. Son los dos géneros que el artículo destaca: *Arcobacter* en el Rímac bajo y *Hypnocyclicus* en la parte alta.

```bash
tail -n +2 exported_taxonomy/taxonomy.tsv | wc -l

13298
```

> **Comentario:** `tail -n +2` muestra el archivo desde la línea 2 (sin el encabezado) y `wc -l` cuenta las líneas: hay **13 298 ASVs** clasificados, uno por línea.

```bash
grep -c "p__" exported_taxonomy/taxonomy.tsv

13235
```

```bash
grep -c "f__" exported_taxonomy/taxonomy.tsv

12678
```

```bash
grep -c "g__" exported_taxonomy/taxonomy.tsv

11804
```

> **Comentario:**
> - **Comando:** `grep -c` cuenta las líneas que contienen el texto indicado. Cuando el clasificador no alcanza la confianza en un nivel, el linaje termina en el nivel anterior y el prefijo no aparece; por eso contar `p__`, `f__` o `g__` equivale a contar los ASVs clasificados hasta ese rango.
> - **Resultado:** el **99,5 %** de los ASVs se clasificó hasta filo (13 235), el **95,3 %** hasta familia (12 678) y el **88,8 %** hasta género (11 804). Romero *et al.* reportan 95–99 % de asignación hasta familia y 84 % hasta género con SILVA 138 completo.

```bash
grep -c "Unassigned" exported_taxonomy/taxonomy.tsv

1
```

```bash
grep -c "d__Archaea" exported_taxonomy/taxonomy.tsv

38
```

```bash
grep -c "Mitochondria" exported_taxonomy/taxonomy.tsv

34
```

```bash
grep -c "Chloroplast" exported_taxonomy/taxonomy.tsv

72
```

> **Comentario:**
> - **Resultado:** 1 ASV no se pudo asignar ni siquiera a dominio (`Unassigned`) y 38 son arqueas: los primers 341F/805R amplifican sobre todo bacterias. Hay **34 ASVs de mitocondrias y 72 de cloroplastos**: no son bacterias de vida libre, sino el gen 16S de los organelos de eucariotas (algas, restos vegetales) presentes en el agua. Se eliminan en la sección 7.
> - Además, 62 ASVs se clasificaron solo hasta dominio (13 298 − 13 235 − 1 = 62).

> **Punto de control:** ¿Qué porcentaje de los ASVs llega a género? ¿Por qué la resolución disminuye de filo a género? ¿Qué ventaja tiene un clasificador recortado a la región V3-V4 frente a uno de secuencias completas?

## 7. Filtrado de ASVs y composición taxonómica

### Filtrar los ASVs según su clasificación

```bash
qiime taxa filter-table --i-table table.qza --i-taxonomy taxa.qza --p-include p__ --p-exclude mitochondria,chloroplast --o-filtered-table table_f.qza

Saved FeatureTable[Frequency] to: table_f.qza
```

> **Comentario:**
> - `--i-table table.qza`: tabla de abundancias sin filtrar.
> - `--i-taxonomy taxa.qza`: taxonomía de cada ASV.
> - `--p-include p__`: conserva solo los ASVs cuya taxonomía contiene `p__`, es decir, los clasificados al menos hasta filo. Elimina el ASV `Unassigned` y los 62 que solo llegan a dominio.
> - `--p-exclude mitochondria,chloroplast`: elimina los ASVs cuya taxonomía contiene alguno de estos términos (separados por coma, sin distinguir mayúsculas de minúsculas).
> - `--o-filtered-table table_f.qza`: tabla filtrada. Quedan **13 129 ASVs** (13 298 − 1 − 62 − 34 − 72).

```bash
qiime taxa filter-seqs --i-sequences representatives.qza --i-taxonomy taxa.qza --p-include p__ --p-exclude mitochondria,chloroplast --o-filtered-sequences representatives_f.qza

Saved FeatureData[Sequence] to: representatives_f.qza
```

> **Comentario:**
> - `--i-sequences representatives.qza`: secuencias de los ASVs sin filtrar.
> - `--i-taxonomy`, `--p-include` y `--p-exclude`: los mismos criterios del comando anterior, para que la tabla y las secuencias tengan exactamente los mismos ASVs.
> - `--o-filtered-sequences representatives_f.qza`: secuencias filtradas, que se usan para construir el árbol filogenético.

### Gráfico de barras de la composición taxonómica

```bash
qiime taxa barplot --i-table table_f.qza --i-taxonomy taxa.qza --m-metadata-file metadata.txt --o-visualization visualization/taxa_barplot.qzv

Saved Visualization to: visualization/taxa_barplot.qzv
```

> **Comentario:**
> - `--i-table table_f.qza`: tabla filtrada.
> - `--i-taxonomy taxa.qza`: taxonomía de cada ASV.
> - `--m-metadata-file metadata.txt`: metadata; permite ordenar las barras por subcuenca o por altitud.
> - `--o-visualization visualization/taxa_barplot.qzv`: gráfico de barras interactivo.

### Descargar `taxa_barplot.qzv` y abrirlo en QIIME 2 View

> **[Captura pendiente: `taxa_barplot.qzv` en el nivel 2 (filo), con las muestras ordenadas por altitud]**

**Figura 12.** Gráfico de barras apiladas de `taxa_barplot.qzv` en el nivel taxonómico 2 (filo), con las muestras ordenadas por altitud de mayor a menor. Cada barra es una muestra y cada color, un filo; la altura de cada segmento es su abundancia relativa.

> **Comentario (figura 12):**
> - En *Taxonomic Level* se elige el rango (1 = dominio … 6 = género; el nivel 7 no agrega información, porque el clasificador no asigna especie). En *Sort Samples By* se ordenan las muestras según una columna de la metadata. *Download* permite descargar la figura y el CSV con las abundancias de cada nivel.
> - Abundancia relativa promedio de los filos principales en cada subcuenca (en %):
>
> | Filo | Upper Rímac | Parac-Upper Rímac | Santa Eulalia-Parac | Santa Eulalia | Jicamarca-Santa Eulalia | Lower Rímac |
> |---|---|---|---|---|---|---|
> | *Pseudomonadota* | 36,2 | 32,3 | 30,2 | 40,7 | 50,0 | 40,9 |
> | *Bacteroidota* | 26,6 | 34,2 | 32,1 | 33,8 | 30,4 | 15,6 |
> | *Campylobacterota* | 8,0 | 5,6 | 3,0 | 3,5 | 5,5 | 30,1 |
> | *Fusobacteriota* | 19,3 | 15,9 | 15,2 | 1,6 | 2,0 | 8,1 |
> | *Bacillota* | 4,5 | 4,6 | 1,9 | 10,1 | 4,7 | 2,8 |
>
> - Cinco filos reúnen alrededor del 93 % de las lecturas en todo el río. Los mismos cinco, con sus nombres anteriores (*Proteobacteria*, *Bacteroidota*, *Campylobacterota*, *Fusobacteriota* y *Firmicutes*), son los que reporta el artículo.
> - El cambio más marcado está en el **Rímac bajo**: *Campylobacterota* sube a 30 % (frente a 3–8 % en las demás subcuencas) y *Bacteroidota* baja a la mitad.
> - *Fusobacteriota* abunda en las tres subcuencas altas (15–19 %) y cae bruscamente desde Santa Eulalia (menos de 2 %). El artículo describe este mismo descenso de Fusobacteriia entre Chacahuaro y Huachipa.
> - Cambie al nivel 3 (clase): *Gammaproteobacteria*, *Bacteroidia*, *Campylobacteria*, *Fusobacteriia*, *Alphaproteobacteria* y *Clostridia* son las clases más frecuentes, como en la figura 4 del artículo.

> **Punto de control:** Cambie al nivel 6 (género). ¿Qué géneros dominan en las subcuencas altas y cuáles en el Rímac bajo? ¿Coinciden con los que reporta el artículo? Busque en la literatura qué tipo de ambiente o de contaminación indica cada uno. Este análisis forma parte de los resultados de la bitácora.

## 8. Árbol filogenético y exportación de los datos

### Construir el árbol filogenético de los ASVs

```bash
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences representatives_f.qza --p-n-threads 10 --o-alignment aligned_representatives.qza --o-masked-alignment masked_aligned_representatives.qza --o-tree unrooted_tree.qza --o-rooted-tree rooted_tree.qza

Saved FeatureData[AlignedSequence] to: aligned_representatives.qza
Saved FeatureData[AlignedSequence] to: masked_aligned_representatives.qza
Saved Phylogeny[Unrooted] to: unrooted_tree.qza
Saved Phylogeny[Rooted] to: rooted_tree.qza
```

> **Comentario (comando):**
> - `--i-sequences representatives_f.qza`: secuencias de los ASVs filtrados.
> - `--p-n-threads 10`: número de hilos.
> - `--o-alignment aligned_representatives.qza`: alineamiento múltiple de las secuencias, obtenido con MAFFT.
> - `--o-masked-alignment masked_aligned_representatives.qza`: el mismo alineamiento sin las columnas muy variables o con muchos huecos, que agregan ruido a la filogenia.
> - `--o-tree unrooted_tree.qza`: árbol sin raíz inferido con FastTree (máxima verosimilitud aproximada).
> - `--o-rooted-tree rooted_tree.qza`: el mismo árbol, enraizado en el punto medio de la rama más larga.
>
> **Comentario (resultado):** se generan cuatro artefactos, uno por cada salida. El árbol enraizado permite calcular métricas que consideran el parentesco entre los ASVs: UniFrac (diversidad beta) y la diversidad filogenética de Faith (diversidad alfa).

### Exportar la tabla de abundancias y el árbol

```bash
qiime tools export --input-path table_f.qza --output-path exported_table

Exported table_f.qza as BIOMV210DirFmt to directory exported_table
```

> **Comentario:** exporta la tabla filtrada en formato BIOM (versión 2.1). La carpeta `exported_table` contiene el archivo `feature-table.biom`.

```bash
qiime tools export --input-path rooted_tree.qza --output-path exported_tree

Exported rooted_tree.qza as NewickDirectoryFormat to directory exported_tree
```

> **Comentario:** exporta el árbol enraizado en formato Newick, un formato de texto que representa el árbol con paréntesis y longitudes de rama. La carpeta `exported_tree` contiene el archivo `tree.nwk`.

### Ajustar el encabezado de la taxonomía y generar el archivo BIOM con taxonomía

```bash
sed -i '1s/.*/#OTUID\ttaxonomy\tconfidence/' exported_taxonomy/taxonomy.tsv
```

> **Comentario:**
> - `sed`: editor de texto en la línea de comandos.
> - `-i`: modifica el archivo directamente.
> - `1s/.*/…/`: en la línea 1, reemplaza (`s`) la línea completa (`.*`) por el texto indicado.
> - `#OTUID\ttaxonomy\tconfidence`: el nuevo encabezado; `\t` es un tabulador.
> - El encabezado original (`Feature ID  Taxon  Confidence`) se cambia porque el comando `biom` espera que las columnas se llamen `#OTUID`, `taxonomy` y `confidence`.

```bash
head -n 2 exported_taxonomy/taxonomy.tsv

#OTUID  taxonomy        confidence
dcd4030bd55a4bbc6d5d935b8ba2797d        d__Bacteria;p__Campylobacterota;c__Campylobacteria;o__Campylobacterales;f__Arcobacteraceae;g__Arcobacter;s__    0.9999999995602877
```

> **Comentario:** el encabezado cambió y el resto del archivo quedó igual.

```bash
biom add-metadata -i exported_table/feature-table.biom -o feature-table-tax.biom --observation-metadata-fp exported_taxonomy/taxonomy.tsv --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy
```

> **Comentario:**
> - `-i exported_table/feature-table.biom`: tabla de abundancias exportada.
> - `-o feature-table-tax.biom`: tabla de salida, con la taxonomía incorporada.
> - `--observation-metadata-fp exported_taxonomy/taxonomy.tsv`: archivo con la información de cada ASV. En el formato BIOM, los ASVs son las *observations* y las muestras, las *samples*.
> - `--observation-header OTUID,taxonomy,confidence`: nombres de las columnas de ese archivo, en orden.
> - `--sc-separated taxonomy`: indica que la columna `taxonomy` contiene niveles separados por punto y coma (*semicolon-separated*), para guardarlos como una lista ordenada de rangos.
> - Solo se incorpora la taxonomía de los ASVs que están en la tabla filtrada; el comando no muestra nada en pantalla.

```bash
biom summarize-table -i feature-table-tax.biom | head -n 12
```

> **Comentario:** resume la tabla BIOM. Verifique que indique **24 muestras** (`Num samples`), **13 129 ASVs** (`Num observations`) y **2 009 113 lecturas** (`Total count`), y que en `Observation Metadata Categories` aparezcan `taxonomy` y `confidence`. Es la comprobación de que la tabla está completa y lleva la taxonomía antes de cargarla en MicrobiomeAnalyst. En el resumen por muestra, el mínimo (35 404 lecturas) corresponde a 11.Universitaria y el máximo (113 259), a 04.Huanchor.

### Preparar la metadata y el árbol para MicrobiomeAnalyst

```bash
sed '1s/^[^\t]*/#NAME/' metadata.txt > metadata_MA.txt
```

> **Comentario:**
> - `1s/^[^\t]*/#NAME/`: en la línea 1, reemplaza el texto desde el inicio de la línea (`^`) hasta el primer tabulador (`[^\t]*`) por `#NAME`. Es decir, cambia el nombre de la primera columna.
> - `> metadata_MA.txt`: guarda el resultado en un archivo nuevo. Sin la opción `-i`, `metadata.txt` no se modifica y se sigue usando en QIIME 2.

```bash
head -n 3 metadata_MA.txt

#NAME   cuenca  altitud
01.Chicla       Upper_Rimac     3441
02.San_Mateo    Parac-Upper_Rimac       3422
```

> **Comentario:** MicrobiomeAnalyst exige que la primera columna de la metadata se llame `#NAME`; las demás columnas no cambian.

```bash
cp exported_tree/tree.nwk rooted_tree.nwk
```

> **Comentario:** copia el árbol a la carpeta de trabajo con un nombre que indica que es el árbol enraizado.

### Descargar con WinSCP los archivos `feature-table-tax.biom`, `metadata_MA.txt` y `rooted_tree.nwk`

## 9. Análisis de diversidad en MicrobiomeAnalyst

> **Importante:** use el servidor oficial, https://www.microbiomeanalyst.ca. El servidor de desarrollo (`dev.microbiomeanalyst.ca`) puede dar errores sin mensaje en el paso de normalización.

### Ir a MicrobiomeAnalyst, entrar al módulo *Marker Data Profiling* y cargar los tres archivos

> **[Captura pendiente: página de carga de datos de *Marker Data Profiling*]**

**Figura 13.** Página de carga de datos del módulo *Marker Data Profiling*, con el formato BIOM seleccionado, la tabla `feature-table-tax.biom`, la metadata `metadata_MA.txt`, el árbol `rooted_tree.nwk` y las etiquetas taxonómicas de SILVA.

> **Comentario (figura 13):**
> - Seleccione el formato **BIOM** y cargue `feature-table-tax.biom`: la taxonomía viaja dentro del archivo, por lo que no hace falta una tabla de taxonomía aparte.
> - Cargue `metadata_MA.txt` como archivo de metadata y `rooted_tree.nwk` como árbol filogenético (opcional; habilita UniFrac y la diversidad filogenética).
> - En el tipo de etiquetas taxonómicas, indique **SILVA**.

### Revisar la integridad de los datos (*Data Integrity Check*)

> **[Captura pendiente: *Data Integrity Check*]**

**Figura 14.** Página *Data Integrity Check* de MicrobiomeAnalyst: opciones de filtrado por defecto y resumen de los datos cargados (*Microbiome data overview*).

> **Comentario (figura 14):**
> - **Default Filtering:** MicrobiomeAnalyst ofrece eliminar los ASVs constantes (*Constant features*) y los *singletons*, ya sea los que aparecen en una sola muestra (*One sample occurrence*) o los que tienen una sola lectura en total (*One total count*). Se deja sin marcar *Constant features* y en *Singleton* se elige **None**, para no perder ASVs raros antes del análisis de diversidad alfa.
> - **Resumen de los datos:**
>   - `Data type`: *OTU abundance table*. MicrobiomeAnalyst trata los ASVs como OTUs.
>   - `File format`: *biom*.
>   - `Sample names match (metadata vs. OTU table)`: **Yes**. Los nombres de las muestras coinciden entre la tabla y la metadata; si apareciera *No*, revise que la primera columna de la metadata se llame `#NAME` y que los nombres sean idénticos.
>   - `Normalized counts detected`: **No**. La tabla tiene conteos crudos, que es lo que se espera.
>   - `OTU annotation`: *SILVA*.
>   - `OTU number`: **13129/13129**. Los mismos 13 129 ASVs de la tabla filtrada en QIIME 2 (sección 7).
>   - `Number of experimental factors`: **2** (uno discreto, `cuenca`, y uno continuo, `altitud`).
>   - `Total read counts`: **2 009 113** y `Average counts per sample`: **83 713**. Coinciden con el resumen de `biom summarize-table`.
> - En el panel izquierdo, *Lib Size View* descarga el gráfico del tamaño de biblioteca (lecturas por muestra).

### Filtrar los datos (*Data Filtering*)

> **[Captura pendiente: *Data Filtering*]**

**Figura 15.** Página *Data Filtering*, con el filtro de conteos bajos (*Low count filter*) y el filtro de baja varianza (*Low variance filter*).

> **Comentario (figura 15):**
> - **Low count filter:** elimina los ASVs con muy pocas lecturas en pocas muestras, que suelen ser errores de secuenciación o contaminación de bajo nivel. *Minimum count* es el conteo mínimo y *Prevalence in samples (%)*, el porcentaje de muestras en que debe alcanzarse (también se puede filtrar por la media o la mediana).
> - **Low variance filter:** elimina el porcentaje indicado de ASVs con menor variación entre muestras (medida como rango intercuartílico, desviación estándar o coeficiente de variación), que no aportan a la comparación entre grupos.
> - Un filtro se desactiva llevando su barra al extremo izquierdo (valor 0). En la práctica se deja **Minimum count = 0** y **Percentage to remove = 0**, es decir, sin filtrar: la diversidad alfa depende de los ASVs raros, y el filtrado ya se hizo en QIIME 2 con criterios taxonómicos.

### Normalizar los datos (*Normalization*)

> **[Captura pendiente: *Data Normalization*]**

**Figura 16.** Página *Data Normalization*, con las opciones de rarefacción (*Data rarefying*), escalado (*Data scaling*) y transformación (*Data transformation*).

> **Comentario (figura 16):**
> - **Data rarefying:** submuestrea todas las muestras al mismo número de lecturas. La barra propone por defecto el tamaño de la biblioteca más pequeña (**35 404**, la de 11.Universitaria).
> - **Data scaling:** *TSS* (proporciones; cada muestra suma 1), *CSS* (escalado de suma acumulada, de metagenomeSeq) o *UQ* (normalización por el cuartil superior).
> - **Data transformation:** *RLE*, *TMM* o *CLR* (razón logarítmica centrada).
> - Se puede rarefaccionar y luego escalar **o** transformar, pero no escalar y transformar a la vez: los datos escalados o transformados ya no son conteos.
> - MicrobiomeAnalyst recomienda rarefaccionar cuando el tamaño de las bibliotecas varía más de 10 veces. Aquí la muestra más grande (113 259) tiene 3,2 veces más lecturas que la más pequeña (35 404), por lo que se deja **Do not rarefy my data**, **Do not scale my data** y **Do not transform my data**.

### Composición taxonómica

> **[Captura pendiente: gráfico de barras apiladas a nivel de género, agrupado por `cuenca`]**

**Figura 17.** Gráfico de barras apiladas (*Stacked Bar/Area Plot*) de la abundancia relativa de los géneros más abundantes en cada subcuenca.

> **Comentario (figura 17):** en *Taxonomic level* elija **Genus** y en *Experimental factor*, **cuenca**. Puede mostrar cada muestra o agrupar las muestras por subcuenca, y limitar el número de géneros que se dibujan (el resto se agrupa como *Others*).

### Diversidad alfa

> **[Captura pendiente: diversidad alfa con el índice de Chao1 a nivel de género]**

**Figura 18.** Diagramas de cajas del índice de Chao1 por subcuenca, a nivel de género, con el resultado de la prueba estadística.

> **Comentario (figura 18):**
> - La diversidad alfa describe la diversidad **dentro** de cada muestra. Se calculan dos índices:
>   - **Chao1:** estima la **riqueza** (número de taxones), incluidos los no observados, a partir de los taxones con una o dos lecturas.
>   - **Shannon:** combina riqueza y **equidad**; es mayor cuando hay muchos taxones con abundancias parecidas y menor cuando uno o pocos dominan.
> - Elija el índice (*Diversity measure*), el nivel taxonómico (*Taxonomic level*: **Family** y luego **Genus**) y el factor (*Experimental factor*: **cuenca**). Repita el análisis para las cuatro combinaciones.
> - Para la prueba estadística use la opción paramétrica (ANOVA) o la no paramétrica (Kruskal-Wallis). Con 2 a 4 muestras por subcuenca, la prueba no paramétrica es la más prudente.
> - Romero *et al.* calcularon Chao1 y Shannon a nivel de ASV y encontraron valores significativamente menores en el Rímac bajo (figura 2 del artículo).

### Diversidad beta

> **[Captura pendiente: PCoA con la distancia de Bray-Curtis a nivel de género]**

**Figura 19.** Análisis de coordenadas principales (PCoA) basado en la distancia de Bray-Curtis, con las muestras coloreadas por subcuenca y el resultado de la PERMANOVA.

> **[Captura pendiente: NMDS con la distancia de Bray-Curtis a nivel de género]**

**Figura 20.** Escalamiento multidimensional no métrico (NMDS) basado en la distancia de Bray-Curtis, con las muestras coloreadas por subcuenca.

> **Comentario (figuras 19 y 20):**
> - La diversidad beta compara la composición **entre** muestras. Cada punto es una muestra: cuanto más cerca están dos puntos, más se parecen sus comunidades.
> - **Distancia de Bray-Curtis:** mide la disimilitud entre dos muestras a partir de las abundancias de los taxones compartidos; va de 0 (composición idéntica) a 1 (ningún taxón compartido). Con el árbol cargado también puede usar **UniFrac**, que considera además el parentesco filogenético entre los ASVs.
> - **PCoA:** proyecta las distancias en ejes ordenados según la variación que explican; el porcentaje de cada eje indica cuánta variación de la matriz de distancias representa.
> - **NMDS:** ordena las muestras para conservar el **orden** de las distancias, no sus valores. Su calidad se mide con el **stress**: por debajo de 0,1 la representación es buena; por encima de 0,2, poco confiable.
> - **PERMANOVA:** prueba si la composición difiere entre los grupos del factor (aquí, las subcuencas). Reporte el valor de F, el R² (proporción de la variación explicada por el factor) y el valor de *p*.
> - Repita ambos métodos a nivel de **Family** y de **Genus**.
> - Romero *et al.* usaron un MDS con la distancia de Bray-Curtis sobre datos normalizados con VST (DESeq2) y encontraron que las muestras del Rímac bajo se agrupan entre sí y se separan de las demás subcuencas (figura 3 del artículo).

> **Punto de control:** ¿Qué subcuenca tiene la menor riqueza y la menor diversidad? ¿Cambia la conclusión entre familia y género, o entre Chao1 y Shannon? En la ordenación, ¿las réplicas quedan junto a su muestra? ¿Las muestras se separan por subcuenca o siguen un gradiente de altitud?

### Estructura de carpetas esperada al finalizar

```bash
~/metataxonomic/
├── tmp/                                    # Temporales de QIIME 2 (TMPDIR)
├── database/
│   └── silva-138.2-ssu-nr99-V3V4-341f-805r-classifier.qza   # Clasificador SILVA 138.2 V3-V4
└── illumina/
    ├── manifest.txt                        # Rutas de los FASTQ de cada muestra
    ├── metadata.txt                        # Metadata para QIIME 2
    ├── metadata_MA.txt                     # Metadata para MicrobiomeAnalyst (#NAME)
    ├── demuxed_seqs.qza                    # Lecturas importadas
    ├── demuxed_seqs_trimmed.qza            # Lecturas sin primers (cutadapt)
    ├── table.qza                           # Tabla de abundancias de DADA2 (copiada)
    ├── representatives.qza                 # Secuencias de los ASVs (copiadas)
    ├── denoising_stats.qza                 # Estadísticas de DADA2 (copiadas)
    ├── taxa.qza                            # Taxonomía de los ASVs
    ├── table_f.qza                         # Tabla filtrada
    ├── representatives_f.qza               # Secuencias filtradas
    ├── aligned_representatives.qza         # Alineamiento (MAFFT)
    ├── masked_aligned_representatives.qza  # Alineamiento enmascarado
    ├── unrooted_tree.qza                   # Árbol sin raíz (FastTree)
    ├── rooted_tree.qza                     # Árbol enraizado
    ├── rooted_tree.nwk                     # Árbol enraizado (Newick, para MicrobiomeAnalyst)
    ├── feature-table-tax.biom              # Tabla con taxonomía (para MicrobiomeAnalyst)
    ├── exported_table/                     # feature-table.biom
    ├── exported_taxonomy/                  # taxonomy.tsv
    ├── exported_tree/                      # tree.nwk
    └── visualization/
        ├── seqs_quality.qzv
        ├── seqs_quality_trimmed.qzv
        ├── denoising_stats.qzv
        ├── table_denoised.qzv
        ├── representatives.qzv
        ├── taxa.qzv
        └── taxa_barplot.qzv
```

### Bitácora bioinformática:

Debe incluir las siguientes secciones. Cada comando debe presentarse con su resultado, seguido de su explicación; cada tabla, figura o captura debe estar numerada, tener una leyenda y estar explicada en el texto:

1. **Carátula** (usar la carátula del modelo de bitácora, con los 6 integrantes del grupo)
2. **Título**
3. **Objetivo de la práctica**
4. **Metodología:** flujograma de los análisis realizados (lecturas FASTQ → importación → cutadapt → DADA2 → clasificación con SILVA 138.2 → filtrado → gráfico de barras → árbol filogenético → exportación a BIOM → MicrobiomeAnalyst: integridad, filtrado, normalización, composición, diversidad alfa y beta), con los archivos de entrada y salida de cada etapa
5. **Metodología:** estructura de las carpetas
6. **Metodología:** versión de cada programa y base de datos (QIIME 2, SILVA y MicrobiomeAnalyst, con la fecha del análisis) y parámetros principales (primers y opciones de cutadapt, parámetros de DADA2, umbral de confianza del clasificador, criterios de filtrado de ASVs, y opciones de filtrado, normalización, índices, nivel taxonómico, distancia y prueba estadística en MicrobiomeAnalyst)
7. **Resultados:** número de lecturas por muestra en los datos crudos y después del recorte de primers, porcentaje conservado y calidad por posición de las lecturas forward y reverse
8. **Resultados:** estadísticas del denoising (lecturas que se conservan en cada etapa de DADA2), número de ASVs, distribución de su longitud y profundidad final por muestra
9. **Resultados:** resolución de la asignación taxonómica (porcentaje de ASVs clasificados hasta filo, familia y género) y número de ASVs eliminados en el filtrado (sin asignar, solo dominio, mitocondrias y cloroplastos)
10. **Resultados:** composición taxonómica a nivel de filo y de género por subcuenca (gráfico de barras e interpretación de los géneros dominantes)
11. **Resultados:** diversidad alfa (Chao1 y Shannon) a nivel de familia y de género, con la prueba estadística
12. **Resultados:** diversidad beta (PCoA y NMDS con la distancia de Bray-Curtis) a nivel de familia y de género, con la PERMANOVA y el stress del NMDS
13. **Discusión:**
    - ¿Qué índice de diversidad alfa y qué nivel taxonómico representan mejor la diversidad de cada muestra? Justifique con sus resultados.
    - ¿Qué método de ordenación representa mejor la similitud o disimilitud entre las muestras? Considere la subcuenca, la altitud y la posición de las réplicas.
    - ¿Sus resultados coinciden con los de Romero *et al.* (2021)? Explique las diferencias según el método: procesamiento (QIIME 2 frente a DADA2 en R, réplicas analizadas por separado o juntas), versión de SILVA y nomenclatura de los filos, clasificador y resolución hasta especie, y normalización usada en la diversidad beta.
