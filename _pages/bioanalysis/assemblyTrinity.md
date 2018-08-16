---
layout: page
title: "AssemblyTrinity"
permalink: /bioanalysis/trinity/
tags: [ RNASeq, trinity, transcriptomics ]
description: De novo transcriptome assembly and functional annotation of transcrits
---
<table class="table-contact">
<tr>
<td><img width="70%" src="{{ site.url }}/images/trainings-trinity.png" alt="" />
</td>
<td>
<h3> De novo transcriptome assembly and functional annotation of transcrits </h3><br />

</td>
</tr>
</table>

| Name | Transcriptome Assembly and Funtional Annotation |
| :------------- | :------------- | :------------- | :------------- |
| Description | This page describes a serie of tools and linux commands used to manipulate raw data (fastq file) for transcriptome assembly  and funtional annotation of transcrits using Trinity and Trinonate. |
| Authors | Julie Orjuela (julie.orjuela_at_ird.fr)  |
| Research Unit | UMR BOREA IPME DIADE |
| Institut |  IRD |
| Creation Date | 10/08/2018 |
| Last Modified Date | 10/08/2018 |

We need, in this tutorial:
* A directory with fastq files
* Samples information (biological replicates?).

### Keywords
Trinity, assembly, de novo, normalisation, RNAseq, transcriptomics

### Files format
fastq, sam, bam

### Date
10/08/2018

***

## Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->

- [1. Checking quality control and cleaning reads](#preAssembly)
  - [1.1. Quality control using `fastqc`](#fastqc)
  - [1.2. Quality trimming and adapter removal using `Trimmomatic`](#trimmomatic)
  - [1.3. Normalization using `Trinity`](#ReadsNormalisation)
  - [1.4. Removing Ribosomal RNA using `sortmerna`](#SORTMERNA)
  
  
- [2. Generating a Trinity de novo RNA-Seq assembly](#trinity)
  - [2.1. Evaluating the quality of the assembly](#assemblyQuality)
  - [2.2. Quantifying transcript expression levels](#QuantifyingTranscriptExpression)
  - [2.3. Identifying differentially expressed (DE) transcripts](#DEtranscrits)
  
  
- [3. Functional annotation of transcripts using `Trinotate` and predicting coding regions using `TransDecoder`](#trinonate)
  - [3.1. Examining functional enrichments for DE transcripts using GOseq](#GO)
  - [3.2. Interactively Exploring annotations and expression data via TrinotateWeb](#trinonateWeb)


<!-- /TOC -->

***
In this section, $shortName it is the sample name.

<a name="preAssembly"></a>
## 1. Checking quality control and cleaning reads

<a name="fastqc"></a>
### 1.1. Quality control of reads using `fastqc`
Follow fastqc protocol from bioanalysis/polymorphism/#fastq-info (creer lien)

<a name="trimmomatic"></a>
### 1.2. Quality trimming and adapter removal using `Trimmomatic`

{% highlight bash %}
java -Xmx4G -jar $path_to_trimmomatic/trimmomatic-0.33.jar PE -phred33 -threads 16 \
-trimlog logfile_$shortName $shortName_R1.fastq.gz $shortName_R2.fastq.gz \
$path_to_trimmomatic_results/$shortName_R1.PairedTrimmed.fastq.gz \
$path_to_trimmomatic_results/$shortName_R1.PairedUntrimmed.fastq.gz \
$path_to_trimmomatic_results/$shortName_R2.PairedTrimmed.fastq.gz \
$path_to_trimmomatic_results/$shortName_R2.PairedUntrimmed.fastq.gz \
ILLUMINACLIP:"$path_to_trimmomatic_adapters":2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50"
{% endhighlight %}


<a name="ReadsNormalisation"></a>
### 1.3. Normalization using `Trinity`

{% highlight bash %}
perl $path_to_trinity/util/insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
--left $path_to_trimmomatic_results/$shortName_R1.PairedTrimmed.fastq.gz \
--right $path_to_trimmomatic_results/$shortName_R2.PairedTrimmed.fastq.gz \
--pairs_together --PARALLEL_STATS --CPU 8 --output $path_to_normalized_data/
{% endhighlight %}

<a name="SORTMERNA"></a>
### 1.4. Removing Ribosomal RNA using `sortmerna`

##### i. indexing the rRNA databases

{% highlight bash %}
echo "indexing sortmerna" 
path_to_sortmerna="/usr/local/sortmerna-2.1/";
$path_to_sortmerna/indexdb_rna --ref \
$path_to_sortmerna/rRNA_databases/silva-bac-16s-id90.fasta,$path_to_sortmerna/index/silva-bac-16s-db:\
$path_to_sortmerna/rRNA_databases/silva-bac-23s-id98.fasta,$path_to_sortmerna/index/silva-bac-23s-db:\
$path_to_sortmerna/rRNA_databases/silva-arc-16s-id95.fasta,$path_to_sortmerna/index/silva-arc-16s-db:\
$path_to_sortmerna/rRNA_databases/silva-arc-23s-id98.fasta,$path_to_sortmerna/index/silva-arc-23s-db:\
$path_to_sortmerna/rRNA_databases/silva-euk-18s-id95.fasta,$path_to_sortmerna/index/silva-euk-18s-db:\
$path_to_sortmerna/rRNA_databases/silva-euk-28s-id98.fasta,$path_to_sortmerna/index/silva-euk-28s:\
$path_to_sortmerna/rRNA_databases/rfam-5s-database-id98.fasta,$path_to_sortmerna/index/rfam-5s-db:\
$path_to_sortmerna/rRNA_databases/rfam-5.8s-database-id98.fasta,$path_to_sortmerna/index/rfam-5.8s-db 
echo "done"
{% endhighlight %}

##### ii.  merging reads
{% highlight bash %}
cd $path_sortmerna_results/
echo "=> Starting Interleaving of reads .."
zcat $path_to_normalized_data/$shortName_R1.fastq.gz | perl -pe 's/\n/\t/ if $. %4' - > $path_sortmerna_results/TMP_$shortName_R1.fastq
zcat $path_to_normalized_data/$shortName_R2.fastq.gz | perl -pe 's/\n/\t/ if $. %4' - > $path_sortmerna_results/TMP_$shortName_R2.fastq
echo "   Interleaving R1 and R2 .."
paste -d '\n' $path_sortmerna_results/TMP_$shortName_R1.fastq $path_sortmerna_results/TMP_$shortName_R2.fastq |\
tr "\t" "\n" > $path_sortmerna_results/$shortName.interleaved.fastq
echo "   Removing temporal files  .."
rm $path_sortmerna_results/TMP_$shortName_R1.fastq $path_sortmerna_results/TMP_$shortName_R2.fastq
echo "   Interleaving was done."
{% endhighlight %}

##### iii. Filtering out rRNA from reads
{% highlight bash %}
$path_to_sortmerna/sortmerna --fastx -a 8 --log --paired_out -e 0.1 --id 0.97 --coverage 0.97 --otu_map\
--ref $path_to_sortmerna/rRNA_databases/silva-bac-16s-id90.fasta,$path_to_sortmerna/index/silva-bac-16s-db:\
$path_to_sortmerna/rRNA_databases/silva-bac-23s-id98.fasta,$path_to_sortmerna/index/silva-bac-23s-db:\
$path_to_sortmerna/rRNA_databases/silva-arc-16s-id95.fasta,$path_to_sortmerna//index/silva-arc-16s-db: \
$path_to_sortmerna/rRNA_databases/silva-arc-23s-id98.fasta,$path_to_sortmerna//index/silva-arc-23s-db:\
$path_to_sortmerna/rRNA_databases/silva-euk-18s-id95.fasta,$path_to_sortmerna//index/silva-euk-18s-db:\
$path_to_sortmerna/rRNA_databases/silva-euk-28s-id98.fasta,$path_to_sortmerna//index/silva-euk-28s:\
$path_to_sortmerna/rRNA_databases/rfam-5s-database-id98.fasta,$path_to_sortmerna//index/rfam-5s-db:\
$path_to_sortmerna/rRNA_databases/rfam-5.8s-database-id98.fasta,$path_to_sortmerna//index/rfam-5.8s-db \
--reads $path_sortmerna_results/$shortName.interleaved.fastq \
--other $path_sortmerna_results/$shortName.sortmerna.mRNA \
--aligned $path_sortmerna_results/$shortName.sortmerna.aligned -v
{% endhighlight %}


##### iv.  unmerging reads
{% highlight bash %}
echo "=> Starting un-interleave .."
echo "   Processing R1 .. "
perl -pe 's/\n/\t/ if $. %4' $path_sortmerna_results/$shortName.sortmerna.mRNA.fastq | awk 'NR%2 {print}' | tr "\t" "\n" >| $path_sortmerna_results/$shortName_R1.sortmerna.mRNA.fastq
echo "   Processing R2 .."
perl -pe 's/\n/\t/ if $. %4' $path_sortmerna_results/$shortName.sortmerna.mRNA.fastq | awk '(NR+1)%2 {print}'| tr "\t" "\n" >| $path_sortmerna_results/$shortName_R2.sortmerna.mRNA.fastq
echo "   Un-interleaving was done."

{% endhighlight %}
