---
title: "TP 14: cDNASeq"
output: github_document
---

# [cite_start]Introducción teórica [cite: 2]

[cite_start]El cáncer de mama triple negativo (TNBC, por sus siglas en inglés) es un subtipo de cáncer de mama que no expresa receptores de progesterona, estrógeno ni HER2, por lo tanto, no tiene blancos terapéuticos disponibles para tratar con las terapias existentes en la clínica[cite: 3]. [cite_start]Debido a esto, las únicas alternativas disponibles como primera línea de tratamiento son la cirugía, los rayos y la quimioterapia[cite: 4]. [cite_start]Desafortunadamente, aproximadamente el 35% de las pacientes con TNBC desarrollarán metástasis por falta de efectividad en dichos tratamientos[cite: 5]. [cite_start]De este modo, encontrar una terapia eficiente resulta una necesidad médica no satisfecha, resaltando la urgencia de encontrar nuevos marcadores de pronóstico y/o terapias específicas[cite: 6].

[cite_start]En el laboratorio estudiamos a la familia de factores de transcripción RUNX (RUNX1, -2 y -3), la cual regula numerosos procesos del desarrollo incluyendo el crecimiento celular y diferenciación[cite: 7]. [cite_start]Estas proteínas necesitan unirse a su cofactor, CBFβ, para regular la expresión génica ya que la unión aumenta 40 veces la afinidad de RUNX al ADN[cite: 8]. [cite_start]Los tres están involucrados tanto en procesos oncogénicos como de supresión tumoral, dependiendo del origen del tumor[cite: 9]. [cite_start]En TNBC, la expresión de RUNX1 correlaciona con un peor pronóstico y la evidencia sugiere que RUNX2 podría actuar como un oncogén, mientras que RUNX3 podría actuar como un gen supresor de tumores[cite: 10].

[cite_start]En el TP vamos a analizar datos provenientes de un cDNAseq realizado con la línea celular MDA-MB-468 que es una línea de TNBC derivada de un adenocarcinoma pulmonar (es decir, provienen de un tumor mamario que hizo metástasis en el pulmón)[cite: 11]. [cite_start]Las células fueron tratadas por 24 horas con un inhibidor comercial de la interacción entre RUNX y CBFβ (AI-10-104) o con el vehículo en el que está disuelta la droga, que en este caso fue dimetilsulfóxido (DMSO), como control[cite: 12]. [cite_start]El experimento se repitió por triplicado[cite: 13].

Pasado el tiempo de tratamiento, se levantaron las células, se extrajo el ARN y se procedió a armar la librería usando el kit PCR-cDNA Barcoding Kit (SQK-PCB111.24). [cite_start]Una vez armada la librería, esta se secuenció usando un PromethION 2 Solo (Oxford Nanopore Technologies)[cite: 14].

# [cite_start]Flujo de trabajo in vitro [cite: 15]

[cite_start]El preparado de la librería consiste, básicamente en cuatro pasos[cite: 16]:

### [cite_start]1- Retrotranscripción [cite: 17]
[cite_start]Al ARN extraído se le liga un adaptador (tiene una cola de poli-T) que va a asegurarnos que los ARN poli-adenilados se retrotranscriban por completo[cite: 18]. [cite_start]Posteriormente, nos quitamos de encima posibles contaminaciones con ADN usando bolitas (beads) magnéticas que tienen afinidad por el ARN[cite: 19]. [cite_start]Luego, usando un primer cuya secuencia es complementaria al extremo 3’ del adaptador recientemente ligado, realizamos un único ciclo de retrotranscripción[cite: 20]. [cite_start]Cuando la retrotranscriptasa termina de sintetizar la hebra complementaria, deja una colita de citocinas en el extremo 5’ que es aprovechada por otro primer para, junto con la enzima, seguir sintetizando un fragmento más[cite: 21]. [cite_start]De este modo, generamos un híbrido ARN-cDNA, donde el cDNA tiene en sus extremos, secuencias conocidas que serán de utilidad en el próximo paso[cite: 22].

### [cite_start]2- Barcoding [cite: 23]
[cite_start]A partir de estas secuencias conocidas en los extremos del cDNA sintetizado, vamos a agregarles un “código de barras” por muestra mediante PCR[cite: 24]. [cite_start]Es decir, a partir de ahora, los cDNAs tendrán una secuencia única por muestra[cite: 25]. [cite_start]Esto permitirá a la hora de secuenciar, saber de qué muestra provino el cDNA secuenciado[cite: 26]. [cite_start]Este paso se repite por 14 ciclos para amplificar los fragmentos, finalmente obteniendo muchas copias de cDNA doble cadena con una etiqueta por muestra[cite: 27]. [cite_start]Luego de esto, para quitarnos de encima el ARN presente, usamos otras beads magnéticas que ahora son afines al ADN[cite: 28]. [cite_start]Se evalúa la calidad del cDNA mediante una electroforesis en gel de agarosa y la cantidad mediante espectrofotometría[cite: 29].

### [cite_start]3- Mezcla de muestras (pooleo) [cite: 30]
[cite_start]Lo que hacemos, básicamente, es mezclar en un eppendorf cantidades equimolares de cada muestra a secuenciar[cite: 31]. [cite_start]Para ello hay que saber, aproximadamente, el tamaño de los transcritos que hay en las muestras[cite: 32]. [cite_start]Por ejemplo, en humanos el tamaño promedio de los ARNm es de 1,5 kilobases[cite: 33]. [cite_start]Una vez hecho esto, tenemos la biblioteca propiamente dicha[cite: 34].

### [cite_start]4- Agregado de adaptadores [cite: 35]
[cite_start]Previo a secuenciar, el último paso consiste en agregarle adaptadores a la biblioteca[cite: 36]. [cite_start]Estos adaptadores son proteínas motoras que reconocen una secuencia dentro del fragmento que agregamos durante el etiquetado[cite: 37]. [cite_start]Estas proteínas no sólo van a ser de guía para cada cDNA hasta el poro y así ser secuenciado, si no que también tienen actividad helicasa: Cuando la proteína se une al poro, el cDNA doble cadena comienza a “ser succionado” por éste y la actividad helicasa despliega el dímero para que la hebra pueda ser secuenciada[cite: 38]. [cite_start]La hebra complementaria queda desplazada con muy pocas chances de ser secuenciada también, dado a que se pliega, formando estructuras secundarias que le proporcionan inestabilidad[cite: 39].

# [cite_start]Diferencias con secuenciación por Illumina [cite: 40]

[cite_start]La diferencia clave está en la longitud de las lecturas[cite: 41]. [cite_start]Mediante secuenciación por nanoporos, las lecturas serán de cDNA o ARN de principio a fin (siempre y cuando no fallen las enzimas durante la RT y/o la PCR)[cite: 42]. [cite_start]Mientras que, para secuenciar por Illumina es necesario fragmentar el ARN en pedazos más pequeños (típicamente 150 pb) y luego usar algoritmos probabilísticos para determinar de qué transcripto provino)[cite: 43].

[cite_start]Obtener lecturas completas permite, además de evaluar la expresión diferencial de genes, evaluar el uso diferencial de isoformas para un mismo gen[cite: 44]. [cite_start]Esto, por lecturas cortas se puede hacer también, sin embargo no es ideal por la forma en la que se secuencia[cite: 45].

# [cite_start]Flujo de trabajo in silico [cite: 46]

## [cite_start]Previo al TP [cite: 47]

### [cite_start]1- Llamado de bases y de-multiplexado [cite: 48]
[cite_start]A partir de los archivos obtenidos del secuenciador, se hizo el llamado de bases (basecalling) y de-multiplexado de las muestras usando el software Dorado[cite: 49]. [cite_start]En criollo, tradujimos los cambios de corriente a través del poro en bases nitrogenadas y luego reconocimos las secuencias de código de barra para cada lectura, las removimos y guardamos las lecturas en archivos .fastq en una carpeta diferente por muestra[cite: 50].

### [cite_start]2- Control de calidad y filtrado [cite: 51]
[cite_start]Una vez obtenidos los archivos, evaluamos la calidad global de la secuenciación utilizando NanoComp[cite: 52]:

```bash
# NanoComp -o ./ --barcoded -f jpg --fastq *.fastq
[cite_start]``` [cite: 53]

[cite_start]Obtuvimos muchas lecturas cortas que no esperábamos y que al intentar alinear contra el transcriptoma en el paso siguiente, no fue exitoso[cite: 54]. [cite_start]Por lo tanto, realizamos un paso de filtrado de lecturas usando SeqKit[cite: 55]:

```bash
# seqkit seq input.fastq --min-len 250 --min-qual 8 > output.fastq
[cite_start]``` [cite: 56]

[cite_start]Acá le pedimos que nos genere un nuevo archivo .fastq sin las lecturas que sean de una longitud menor a 250 bases y una calidad de la lectura promedio menor a 8 (aceptable para lecturas largas)[cite: 57]. Ahora, al evaluar las lecturas por NanoComp:

[cite_start]![Histograma de longitudes](https://path-to-image.png) [cite: 58]

[cite_start]En el eje Y tenemos la cantidad de lecturas y en el eje X la longitud de las mismas, en escala logarítmica[cite: 59]. [cite_start]Como pueden ver, ahora hay un máximo de lecturas cercano a 1000 bases cuando antes había también un pico en tamaños de lecturas menor a 250 bases[cite: 60].

### [cite_start]3- Alineamiento contra transcriptoma [cite: 61]
[cite_start]Una vez tenemos las lecturas filtradas, usamos minimap2 para alinearlas al transcriptoma humano[cite: 62]:

```bash
# minimap2 -ax map-ont referencia.fa input.fastq > output.sam
[cite_start]``` [cite: 63]

[cite_start]Acá le pedimos que nos genere un archivo .sam (a) y que vamos a usar un protocolo de alineamiento pre-establecido (x) que es map-ont[cite: 64]. [cite_start]Este protocolo considera que la referencia y el archivo .fastq tienen secuencias continuas, es decir, que pueden alinear directamente y no hay que tener en cuenta intrones[cite: 65]. [cite_start]Además, considera que son lecturas largas y por lo tanto, modifica el algoritmo de puntuación para decidir si un alineamiento es bueno o malo respecto al algoritmo para lecturas cortas[cite: 66].

[cite_start]Evaluamos las métricas del mapeo usando samtools con el siguiente comando[cite: 67]:

```bash
# samtools flagstat archivo.sam
[cite_start]``` [cite: 68]

A continuación les muestro uno de los resultados tipo que obtuve:

[cite_start]![Métricas de mapeo](https://path-to-image.png) [cite: 69]

En el archivo tenemos 4.121.356 lecturas totales. [cite_start]De estas, sólo 1.093.269 alinean una única vez contra la referencia[cite: 70]. 3.006.408 lecturas se alinean en múltiples lugares. [cite_start]21.679 son lecturas divididas[cite: 71]. [cite_start]Del total de lecturas, se mapearon un 99.34% lo cual es un porcentaje muy alto y bueno[cite: 72].

### [cite_start]4- Cuantificación de lecturas mapeadas [cite: 73]
[cite_start]Por último, las lecturas mapeadas se cuantificaron usando Salmon[cite: 74]:

```bash
# salmon quant -t referencia.fa -l U -a input.sam --noErrorModel --ont -o dir/
[cite_start]``` [cite: 75]

[cite_start]Acá le dimos la referencia con el argumento -t; con el argumento -l U le decimos que la biblioteca que armamos no tiene direccionalidad (es decir que puede secuenciarse la hebra 3’->5’ como la complementaria; con el argumento -a le decimos que nuestro input es un archivo .sam ; --noErrorModel y --ont son dos argumentos para decirle al programa que las lecturas mapeadas provienen de una secuenciación de lecturas largas por nanoporos y no de Illumina, por lo tanto que ajuste el sistema de cuantificación para este modelo[cite: 76].

# [cite_start]Paso a paso del TP [cite: 78]

[cite_start]Primero **importantísimo** es que en **la terminal** escribas este código para sincronizar los datos entre tu computadora y la computadora maestra[cite: 79]:

```bash
rsync -rvP rsync://10.41.115.18/datos_genomica/14_TP_RNAseq_largas/ /media/libre/datos_genomica/14_TP_RNAseq_largas/
[cite_start]``` [cite: 80]

[cite_start]Después activamos el **entorno de conda** que contiene los programas y librerías para el TP y nos movemos a la carpeta **Documentos**[cite: 81].

```bash
conda activate nf
[cite_start]``` [cite: 82]

```bash
cd ~/Documentos
[cite_start]``` [cite: 83]

[cite_start]Creamos una **nueva carpeta** para guardar los archivos a generar y nos movemos a esta[cite: 84].

```bash
mkdir TP_cDNASeq
[cite_start]``` [cite: 85]

```bash
cd TP_cDNASeq
[cite_start]``` [cite: 86]

[cite_start]Abrimos el script de RStudio **desde la terminal**, con el comando[cite: 87]:

```bash
rstudio /media/libre/datos_genomica/14_TP_RNAseq_largas/cDNASeq/TP_GA.R
[cite_start]``` [cite: 88]

[cite_start]Vamos a seguir el TP desde ahí[cite: 89].

# [cite_start]CUESTIONARIO [cite: 90]

[cite_start]1- ¿Por qué hablamos de cDNA-seq y no de RNA-seq? [cite: 91]
[cite_start]2- Teniendo en cuenta cómo es la preparación de las bibliotecas para secuenciar por lecturas cortas y largas, ¿cómo esperarías que, luego de un alineamiento, te den las métricas usando samtools? [cite: 92]
[cite_start]3- ¿Qué tiene un archivo .gtf? [cite: 93]
[cite_start]4- ¿Qué es una matriz de conteos de expresión? [cite: 94]
[cite_start]5- ¿Qué significa TPM? [cite: 95]
6- ¿Qué es un gráfico de PCA? ¿Qué explican los ejes? [cite_start]¿Cómo le explicarías a tu director de tesis cómo dio el gráfico que hiciste en el TP? [cite: 96]
[cite_start]7- ¿Por qué se tiene en cuenta el p-valor ajustado en vez del p-valor? [cite: 97]
[cite_start]8- ¿Cómo se interpreta el log2FC? [cite: 98]
9- ¿Cuáles genes están más alterados en su expresión? [cite_start]Mencione 3 regulados positivamente y negativamente. [cite: 99]
[cite_start]10- ¿Hay más genes regulados negativa o positivamente? [cite: 100]
