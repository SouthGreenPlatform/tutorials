---
layout: page
title: "FROGs en Ligne de Commande"
permalink: /bioanalysis/frogsCL/
tags: [ frogs, swarm, OTU ]
description: comment lancer FROGS en ligne de commande sur le cluster IRD
---
<table class="table-contact">
<tr>
<td><img width="100%" src="{{ site.url }}/images/trainings-linux.png" alt="" />
</td>
<td>
<h1> FROGs en Ligne de Commande</h1><br />
This page describes how to launch FROGs in command line from IRD cluster. From fastq files to OTU and Phyloseq analysis.
</td>
</tr>
</table>


We need, in this tutorial:
* a compressed directory with fastq files
* a database used for the taxonomic assignation


### Author(s)

| Authors  | Julie ORJUELA  |
| :------------- | :------------- |
| Research Unit | UMR IPME-DIADE-BOREA   |
| Institut |  IRD |

### Keywords
flash, swarm, blast, vsearch, metabarcoding, 16S, 18S, ITS

### Files format
fastq, OTU tables

### Date
22/03/2019


## On crée un dossier ou vous voulez.

`mkdir /home/orjuela/TEST-FROGS/fromGitExemple`


### 1. Préparation de fastq

* Tous les fichiers fastq.gz (R1 et R2) seront mis dans un dossier qu'il faudra apres compresser en .tar.gz
Pour compresser le fichier il faut:


- se deplacer dans le dossier des fastq.gz
`cd test_dataset2/`

- Compresser 
`tar zcvf test_dataset.tar.gz * `

- vous obtenez un fichier test_dataset.tar.gz que vous pouvez deplacer avec
`mv test_dataset2.tar.gz .. `

* Vérifier que le dossier compressé a tous les fichiers et qu'il n'y a pas de sous-dossier.

Pour observer les fichiers sans le décompresser utilise :

`tar -tf test_dataset.tar.gz`

Vous devez avoir que la liste des fichiers fastq.gz sans sous-dossier.

{% highlight bash %}

splA_01_R1.fastq.gz
splA_01_R2.fastq.gz
splA_02_R1.fastq.gz
splA_02_R2.fastq.gz
splA_03_R1.fastq.gz
splA_03_R2.fastq.gz

{% endhighlight %}

Notes :

Documentation compression : https://openclassrooms.com/fr/courses/43538-reprenez-le-controle-a-laide-de-linux/41346-archiver-et-compresser

extraire : `tar zxvf`

create:   `tar zcvf`

examiner: `tar -tf`

### 2. Préparation d'un fichier tabulé "sample_metadata.tsv" qui R utilise

exemple 1 :

{% highlight bash %}

	Color	ID
splA_01 red	rep1
splA_02	blue	rep2
splA_03	green	rep1
{% endhighlight %}

exemple2:

{% highlight bash %}
Sample	Cell	Origin	Repetition	Color
17MET040	Cell1	SolNu	R1	red
17MET041	Cell1	SolNu	R2	red
17MET042	Cell1	SolNu	R3	red
17MET037	Cell1	Spergul	R1	green
17MET038	Cell1	Spergul	R2	green
17MET039	Cell1	Spergul	R3	green
17MET052	Cell2	Atriplex	R1	gray
17MET035	Cell2	Atriplex	R2	gray
17MET036	Cell2	Atriplex	R3	gray
17MET049	Cell2	SolNu	R1	pink
17MET050	Cell2	SolNu	R2	pink
17MET051	Cell2	SolNu	R3	pink
17MET046	Cell4	SolNu	R1	blue
17MET047	Cell4	SolNu	R2	blue
17MET048	Cell4	SolNu	R3	blue
17MET043	Cell4	Viperine	R1	orange
17MET044	Cell4	Viperine	R2	orange
17MET045	Cell4	Viperine	R3	orange
{% endhighlight %}


### 3. Connaitre le path de la base de données pour les etapes d'assignation tax

sur le cluster ird

`/usr/local/frogs_databases-2.01/silva_123_16S/silva_123_16S.fasta`


### 4. Visualiser/modifier le script avant de le lancer :

Ouvrir run_frogs_pipeline.sh dans un editeur.

Vous pouvez modifier les lignes 3 et 4 du script pour ajouter le chemin vers les fichiers sample_metadata et la base de données pour l'assignation taxonomique

`samplefile="/home/orjuela/TEST-FROGS/fromGitExemple/sample_metadata.tsv"`

`db="/usr/local/frogs_databases-2.01/silva_123_16S/silva_123_16S.fasta" `

+Le reste on ne touche pas sauf si vous savez ce que vous faites.+


### 5. Lancer le script /home/orjuela/scripts/run_frogs_pipeline.sh

Pour lancer le script place vous dans l'endroit ou vous voulez avoir les résultats :`/home/orjuela/TEST-FROGS/fromGitExemple`

Attention: les amorces doivent etre ecrit en 5'-3' 

`qsub -q bioinfo.q -N frogsCL -b yes -V -cwd -pe ompi 4 'bash /home/orjuela/scripts/run_frogs_pipeline.sh 380 460 GGCGVACGGGTGAGTAA GTGCCAGCNGCNGCGG 250 250 420 OUTPUT /home/orjuela/TEST-FROGS/fromGitExemple/test_dataset.tar.gz'`

`bash /home/orjuela/scripts/run_frogs_pipeline.sh`

{% highlight bash %}
1<minAmpliconSize>
2<maxAmpliconSize>
3<fivePrimPrimer>
4<threePrimPrimer>
5<R1size>
6<R2size>
7<expectedAmpliconSize>
8<out_dir>
9<datasetTarGz>
{% endhighlight %}

- Si tout se passe bien vous verrez ça:
{% highlight bash %}
380
460
GGCGVACGGGTGAGTAA
GTGCCAGCNGCNGCGG
250
250
420
OUTPUT
/home/orjuela/TEST-FROGS/fromGitExemple/test_dataset.tar.gz
Step preprocess ven. sept. 21 11:49:56 CEST 2018
Step clustering ven. sept. 21 11:52:29 CEST 2018
Step remove_chimera ven. sept. 21 11:52:44 CEST 2018
Step filters ven. sept. 21 11:54:30 CEST 2018
Step affiliation_OTU ven. sept. 21 11:54:33 CEST 2018 ...
{% endlight bash %}

Votre dossier OUTPUT doit rassembler à ça

{% highlight bash %}
orjuela@MPLCLTLP0157:~/Documents/tools/FROGS/test/OUT$ ll
total 91524
drwxr-xr-x 2 orjuela orjuela    53248 juin  15 15:23 ./
drwxr-xr-x 4 orjuela orjuela     4096 juil. 12 14:40 ../
-rw-r--r-- 1 orjuela orjuela 39183498 juin  15 15:12 01-prepro.fasta
-rw-r--r-- 1 orjuela orjuela    34922 juin  15 15:12 01-prepro.html
-rw-r--r-- 1 orjuela orjuela    19178 juin  15 15:12 01-prepro.log
-rw-r--r-- 1 orjuela orjuela  4668203 juin  15 15:12 01-prepro.tsv
-rw-r--r-- 1 orjuela orjuela  4941307 juin  15 15:13 02-clustering.biom
-rw-r--r-- 1 orjuela orjuela  4325539 juin  15 15:13 02-clustering_compo.tsv
-rw-r--r-- 1 orjuela orjuela 17010668 juin  15 15:13 02-clustering.fasta
-rw-r--r-- 1 orjuela orjuela     3278 juin  15 15:13 02-clustering.log
-rw-r--r-- 1 orjuela orjuela  1610390 juin  15 15:14 03-chimera.biom
-rw-r--r-- 1 orjuela orjuela  5455138 juin  15 15:14 03-chimera.fasta
-rw-r--r-- 1 orjuela orjuela    13943 juin  15 15:14 03-chimera.html
-rw-r--r-- 1 orjuela orjuela    81978 juin  15 15:14 03-chimera.log
-rw-r--r-- 1 orjuela orjuela   989852 juin  15 15:14 04-affiliation.biom
-rw-r--r-- 1 orjuela orjuela    15831 juin  15 15:14 04-affiliation.html
-rw-r--r-- 1 orjuela orjuela     1824 juin  15 15:14 04-affiliation.log
-rw-r--r-- 1 orjuela orjuela   224506 juin  15 15:14 04-filters.biom
-rw-r--r-- 1 orjuela orjuela   307008 juin  15 15:14 04-filters.excluded
-rw-r--r-- 1 orjuela orjuela   661810 juin  15 15:14 04-filters.fasta
-rw-r--r-- 1 orjuela orjuela   130418 juin  15 15:14 04-filters.html
-rw-r--r-- 1 orjuela orjuela     1478 juin  15 15:14 04-filters.log
-rw-r--r-- 1 orjuela orjuela   187007 juin  15 15:14 05-clustersStat.html
-rw-r--r-- 1 orjuela orjuela      933 juin  15 15:14 05-clustersStat.log
-rw-r--r-- 1 orjuela orjuela   248668 juin  15 15:14 06-affiliationsStat.html
-rw-r--r-- 1 orjuela orjuela     1170 juin  15 15:14 06-affiliationsStat.log
-rw-r--r-- 1 orjuela orjuela     1120 juin  15 15:14 07-biom2tsv.log
-rw-r--r-- 1 orjuela orjuela    23252 juin  15 15:14 07-biom2tsv.multi
-rw-r--r-- 1 orjuela orjuela  1138093 juin  15 15:14 07-biom2tsv.tsv
-rw-r--r-- 1 orjuela orjuela   154897 juin  15 15:14 08-affiliation_multihit.tsv
-rw-r--r-- 1 orjuela orjuela   844425 juin  15 15:14 08-affiliation_std.biom
-rw-r--r-- 1 orjuela orjuela      338 juin  15 15:14 08-biom2stdbiom.log
-rw-r--r-- 1 orjuela orjuela  1124530 juin  15 15:14 09-tsv2biom.biom
-rw-r--r-- 1 orjuela orjuela   641114 juin  15 15:14 09-tsv2biom.fasta
-rw-r--r-- 1 orjuela orjuela      882 juin  15 15:14 09-tsv2biom.log
-rw-r--r-- 1 orjuela orjuela   182059 juin  15 15:17 10a-tree.html
-rw-r--r-- 1 orjuela orjuela     1346 juin  15 15:17 10a-tree.log
-rw-r--r-- 1 orjuela orjuela    64971 juin  15 15:17 10a-tree.nwk
-rw-r--r-- 1 orjuela orjuela   180238 juin  15 15:18 10b-tree.html
-rw-r--r-- 1 orjuela orjuela     1082 juin  15 15:18 10b-tree.log
-rw-r--r-- 1 orjuela orjuela    67318 juin  15 15:18 10b-tree.nwk
-rw-r--r-- 1 orjuela orjuela  1267549 juin  15 15:18 11-phylo_import.html
-rw-r--r-- 1 orjuela orjuela     1295 juin  15 15:18 11-phylo_import.log
-rw-r--r-- 1 orjuela orjuela    69734 juin  15 15:18 11-phylo_import.Rdata
-rw-r--r-- 1 orjuela orjuela  4033991 juin  15 15:19 12-phylo_composition.html
-rw-r--r-- 1 orjuela orjuela     1027 juin  15 15:19 12-phylo_composition.log
-rw-r--r-- 1 orjuela orjuela  1214435 juin  15 15:22 13-phylo_alpha_div.html
-rw-r--r-- 1 orjuela orjuela     1077 juin  15 15:22 13-phylo_alpha_div.log
-rw-r--r-- 1 orjuela orjuela      223 juin  15 15:19 13-phylo_alpha_div.tsv
-rw-r--r-- 1 orjuela orjuela   789726 juin  15 15:22 14-phylo_beta_div.html
-rw-r--r-- 1 orjuela orjuela     1010 juin  15 15:22 14-phylo_beta_div.log
-rw-r--r-- 1 orjuela orjuela      971 juin  15 15:23 16-phylo_clustering.log
-rw-r--r-- 1 orjuela orjuela   865732 juin  15 15:23 16-phylo_clutering.html
-rw-r--r-- 1 orjuela orjuela   748887 juin  15 15:23 17-phylo_manova.html
-rw-r--r-- 1 orjuela orjuela      951 juin  15 15:23 17-phylo_manova.log
-rw-r--r-- 1 orjuela orjuela       67 juin  15 15:22 Jaccard_binary.tsv
-rw-r--r-- 1 orjuela orjuela       67 juin  15 15:22 Unifrac.tsv
{% endlight bash %}

Rapatrier les dossier OUTPUT dans votre machine local et visualiser les html.

Utilisez Fillezilla par exemple




