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
| Description | This page describes a serie of tools and linux commands used to manipulate fastq files for transcriptome assembly and funtional annotation of transcrits using Trinity and Trinotate. |
| Authors | Julie Orjuela (julie.orjuela_at_ird.fr)  |
| Research Unit | UMR BOREA IPME DIADE |
| Institut |  IRD |
| Creation Date | 10/08/2018 |
| Last Modified Date | 13/02/2019 |

We need, in this tutorial:
* A directory with fastq files
* Samples information (biological replicates?).

### Keywords
Trinity, assembly, de novo, normalisation, RNAseq, transcriptomics

### Files format
fastq, sam, bam

***

## Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->

- [1. Checking quality control and cleaning reads](#preAssembly)
  - [1.1. Quality control using `fastqc`](#fastqc)
  - [1.2. Quality trimming and adapter removal using `Trimmomatic`](#trimmomatic)
  - [1.3. Removing Ribosomal RNA using `sortmerna`](#SORTMERNA)
  - [1.4. Normalization using `Trinity`](#ReadsNormalisation)

  
  
- [2. Generating a Trinity de novo RNA-Seq assembly](#trinity)
  - [2.1. Evaluating the quality of the assembly](#assemblyQuality)
  - [2.2. Identifying differentially expressed (DE) transcripts](#DEtranscrits)
  
  
- [3. Functional annotation of transcripts using `Trinotate` and predicting coding regions using `TransDecoder`](#trinotate)


<!-- /TOC -->

***
In this section, $shortName it is the sample name.

<a name="preAssembly"></a>
## 1. Checking quality control and cleaning reads

<a name="fastqc"></a>
### 1.1. Quality control of reads using `fastqc`
Follow fastqc protocol [here](https://southgreenplatform.github.io/tutorials//bioanalysis/rnaSeq/#fastqc)

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

If you are sure of quality reads and parameters, you can directly run trimmomatic and assembly of reads usign Trinity.
Similar for normalisation.

<a name="SORTMERNA"></a>
### 1.3. Removing Ribosomal RNA using `sortmerna`

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
zcat $path_to_trimmomatic_results/$shortName_R1.PairedTrimmed.fastq.gz | perl -pe 's/\n/\t/ if $. %4' - > $path_sortmerna_results/TMP_$shortName_R1.fastq
zcat $path_to_trimmomatic_results/$shortName_R2.PairedTrimmed.fastq.gz | perl -pe 's/\n/\t/ if $. %4' - > $path_sortmerna_results/TMP_$shortName_R2.fastq
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

<a name="ReadsNormalisation"></a>
### 1.4. Normalisation using `Trinity` 

If you don't have biological replicates, you can directly done a alone normalisation of reads by sample.

{% highlight bash %}
perl $path_to_trinity/util/insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
--left $shortName_R1.sortmerna.mRNA.fastq \
--right $shortName_R2.sortmerna.mRNA.fastq \
--pairs_together --PARALLEL_STATS --CPU 8 --output $path_to_normalized_data/
{% endhighlight %}

If biological replicates, you can run trinity assembly with option `--normalize_by_read_set` in section 2 or give R1 and R2 reads for each condition to `insilico_read_normalization.pl` :

{% highlight bash %}
perl $path_to_trinity/util/insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
--left $shortNameReplique1_R1.fastq.gz,\
$shortNameReplique2_R1.sortmerna.mRNA.fastq.gz, \
$shortNameReplique3_R1.sortmerna.mRNA.fastq.gz  \
--right $shortNameReplique1_R2.sortmerna.mRNA.fastq.gz,\
$shortNameReplique2_R2.sortmerna.mRNA.fastq.gz, \
$shortNameReplique3_R2.sortmerna.mRNA.fastq.gz \
--pairs_together --PARALLEL_STATS --CPU 8 --output $path_to_normalized_data/
{% endhighlight %}

<a name="trinity"></a>
## 2. Generating a Trinity de novo RNA-Seq assembly
You can assembly reads from one sample :
{% highlight bash %}
$R1=$path_sortmerna_results/$shortName_R1.sortmerna.mRNA.fastq
$R2=$path_sortmerna_results/$shortName_R2.sortmerna.mRNA.fastq
Trinity --seqType fq --left $R1 --right $R2 --max_memory 50G --CPU 8 --output trinity_OUT
{% endhighlight %}

If you want assembly reads using the whole of samples of a specie (several tissues of a specie without biological replicates) OR
if you have biological replicates in your experiment and you want to obtain a transcriptome by condition :
{% highlight bash %}
Trinity --seqType fq --max_memory 80G --CPU 8 --normalize_by_read_set --samples_file samples.txt --output trinity_OUT 
{% endhighlight %}

Remember that is possible run trimmomatic, normalisation and assembly in one command line :
{% highlight bash %}
Trinity --seqType fq --max_memory 50G --CPU 4 --samples_file sample.txt --trimmomatic --quality_trimming_params "ILLUMINACLIP:illumina.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 --normalize_by_read_set
{% endhighlight %}

Samples.txt file exemple (tabulated file)
{% highlight bash %}
condA\tcondA_rep1\tcondA_rep1_R1.fastq.gz\tcondA_rep1_R2.fastq.gz\n
condA\tcondA_rep2\tcondA_rep2_R1.fastq.gz\tcondA_rep2_R2.fastq.gz\n
condB\tcondB_rep1\tcondB_rep1_R1.fastq.gz\tcondB_rep1_R2.fastq.gz\n
condB\tcondB_rep2\tcondB_rep2_R1.fastq.gz\tcondB_rep2_R2.fastq.gz\n
{% endhighlight %}



<a name="assemblyQuality"></a>
### 2.1. Evaluating the quality of the assembly


#### Assembly metrics

{% highlight bash %}
$path_to_trinity/util/TrinityStats.pl Trinity.fasta
{% endhighlight %}


#### Reads mapping back rate :

A typical ‘good’ assembly has ~80 % reads mapping to the assembly and \~80% are properly paired

- Alignment methods : bowtie2 -RSEM, kallisto, salmon `--est_method` 

{% highlight bash %}
perl $path_to_trinity/util/align_and_estimate_abundance.pl \
--transcripts Trinity.fasta \
--seqType fq \
--left $R1 --right $R2\
--est_method RSEM --aln_method bowtie2 \
--trinity_mode --prep_reference \
--output_dir outdir

OR

perl $path_to_trinity/util/align_and_estimate_abundance.pl \
--transcripts Trinity.fasta \
--seqType fq \
--samples_file samples.txt \
--est_method salmon
--trinity_mode --prep_reference
--output_dir outdir

{% endhighlight %}

We suggest visualise mapping back using IGV. Recovery BAM and Trinity.fasta files and import it in IGV browser. You must to index BAMs files before. Use `samtools index BAM` to do it.

If you don't have replicates and you want only mapping reads agains transcriptome obtained by trinity use :

{% highlight bash %}
$path_to_trinity/util/bowtie_PE_separate_then_join.pl --seqType fq --left left.fq --right right.fq --target Trinity.fasta --aligner bowtie -- -p 4 --all --best --strata -m 300
{% endhighlight %}

To get alignment statistics, run the following:
{% highlight bash %}
$path_to_trinity/util/SAM_nameSorted_to_uniq_count_stats.pl bowtie_out/bowtie_out.nameSorted.bam
{% endhighlight %}

- Expression matrix construction

{% highlight bash %}
$path_to_trinity/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix Trinity_trans\
--name_sample_by_basedir\
cond_A_rep1/abundance.tsv\
cond_A_rep2/abundance.tsv\
cond_B_rep1/abundance.tsv\
cond_B_rep2/abundance.tsv
{% endhighlight %}

You have to obtain two matrices: The firts one containing the estimated counts, and the second one containing the TPM expression values that are cross-sample normalized using the TMM method `Trinity_trans.TMM.EXPR.matrix`. TMM normalization assumes that most transcripts are not differentially expressed, and linearly scales the expression values of samples to better enforce this property.

- Compute N50 based on the top-most highly expressed transcripts (Ex50)

{% highlight bash %}
$path_to_trinity/util/misc/contig_ExN50_statistic.pl Trinity_trans.TMM.EXPR.matrix Trinity.fasta > ExN50.stats
{% endhighlight %}

{% highlight bash %}
$path_to_trinity/util/misc/contig_ExN50_statistic.pl Trinity_trans.TMM.EXPR.matrix Trinity.fasta | tee ExN50.stats
{% endhighlight %}

Plotting ExN50

{% highlight bash %}
% /usr/local/trinityrnaseq-2.5.1/util/misc/plot_ExN50_statistic.Rscript ExN50.stats
{% endhighlight %}

If you want to know, how many transcripts correspond to the Ex 90 peak, you could:
{% highlight bash %}
cat transcripts.TMM.EXPR.matrix.E-inputs |  egrep -v ^\# | awk '$1 <= 90' | wc -l
{% endhighlight %}

https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats


#### Tools to evaluate transcriptomes

To avoid redundant transcripts, we kept the longest isoform for each “gene” identified by TRINITY (unigene) using the `get_longest_isoform_seq_per_trinity_gene.pl` utility in TRINITY:

{% highlight bash %}
$path_to_trinity/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta > Trinity.longest.fasta
{% endhighlight %}

- Validation using Transrate

{% highlight bash %}
$path_to_transrate/transrate --assembly  Trinity.fasta --left $R1 --right $R2  --output transrate_outdir
{% endhighlight %}

- Validation using BUSCO

{% highlight bash %}
BUSCOPathDB="/home/orjuela/BUSCO_DB/actinopterygii_odb9"
python $path_to_busco/scripts/run_BUSCO.py -i Trinity.fasta -o outputBusco -l $BUSCOPathDB -m transcriptome -c 8
{% endhighlight %}

- Validation using BLASTX

First, we downloaded and indexed the database:

{% highlight bash %}
wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
$path_to_ncbi-blast+/makeblastdb -in uniprot_sprot.fasta -dbtype prot -out SwissProt_no_seqids
{% endhighlight %}

Then, we ran BLASTX to get the top match hit:

{% highlight bash %}
$path_to_ncbi-blast+/blastx -db SwissProt_no_seqids -query Trinity.longest.fasta \
-num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-20 > SwissProt_1E20_TrinityLongest_blastx.outfmt6
{% endhighlight %}

Finally, we examined the percent of alignment coverage:

{% highlight bash %}
$path_to_trinity/util/misc/blast_outfmt6_group_segments.pl SwissProt_1E20_TrinityLongest_blastx.outfmt6 \
> SwissProt_1E20_TrinityLongest_blastx.outfmt6.grouped
{% endhighlight %}

{% highlight bash %}
$path_to_trinity/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl SwissProt_1E20_TrinityLongest_blastx.outfmt6.grouped \
> SwissProt_1E20_TrinityLongest_blastx.outfmt6.grouped.output
{% endhighlight %}

If you generate assemblies at a range of different read depths up to and including your assembly leveraging all available reads, you can perform this full-length transcript analysis separately for each of your assemblies, and then plot the number of full-length transcripts vs. number of input RNA-Seq fragments.

<a name="DEtranscrits"></a>
### 2.2 Identifying differentially expressed (DE) transcripts

{% highlight bash %}
 $path_to_trinity/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix Trinity.isoform.counts.matrix \
--samples_file samples.txt \
--method DESeq2 \
--output DESeq2_trans
{% endhighlight %}

- Extracting differentially expressed transcripts and generating heatmaps

Extract those differentially expressed (DE) transcripts that are at least 4-fold  (C is set to 2^(2) ) differentially expressed at a significance of <= 0.001 (-P 1e-3) in any of the pairwise sample comparisons

{% highlight bash %}
cd DESeq2_trans/
$path_to_trinity/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix Trinity.isoform.TMM.EXPR.matrix \
--samples samples.txt -P 1e-3 -C 2 
{% endhighlight %}

- Extract transcript clusters by expression profile by cutting the dendrogram

Extract clusters of transcripts with similar expression profiles by cutting the transcript cluster dendrogram at a given percent of its height (ex. 60%), like so:

{% highlight bash %}
$path_to_trinity/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl \
--Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData
{% endhighlight %}

- Run the DE analysis at the gene level
{% highlight bash %}
$path_to_trinity/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix Trinity.gene.counts.matrix \
--samples_file samples.txt \
--method DESeq2 \
--output DESeq2_gene
{% endhighlight %}

https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-FAQ#ques_why_so_many_transcripts

- Most downstream analyses should be applied to the entire set of assembled transcripts, including functional annotation and differential expression analysis.

If you decide that you want to filter transcripts to exclude those that are lowly expressed, you can use the following script:

{% highlight bash %}
$path_to_trinity/util/filter_low_expr_transcripts.pl
{% endhighlight %}


<a name="trinotate"></a>

## 3. Functional annotation of transcripts using `Trinotate` and predicting coding regions using `TransDecoder`

Transcrits assembled using Trinity can be easily annotate using trinotate https://github.com/Trinotate/Trinotate.github.io/wiki.

Trinotate use different methods for functional annotation including homology search to known sequence data (BLAST+/SwissProt), protein domain identification (HMMER/PFAM), protein signal peptide and transmembrane domain prediction (signalP/tmHMM), and take advantage from annotation databases (eggNOG/GO/Kegg). These data are integrated into a SQLite database which allows to create an annotation report for a transcriptome.

Two bash scripts were created to obtain the whole of files obligatories to build a Sqlite database and create reports. 

The fist one, `trinotate-JAv1.0.sh` https://github.com/julieaorjuela/scripts/blob/master/trinotate-JAv1.0.sh, needs as input a repertory containing the fasta files you want to annotate. It generates three repertories : Trinonate, sh, and trash and a submitQsub.sge file that launch every fasta analysis in job array mode. The bash repertory contains scripts created automatically for every fasta file, the Trinotate repertory contains annotation results and the trash contains the log files for every step in the process.

{% highlight bash %}
bash ~/scripts/trinotate-JAv1.0.sh -f /repertory/containing/fastaFiles/
{% endhighlight %}

{% highlight bash %}
qsub /repertory/containing/fastaFiles/jobArray-Trinotate/submitQsub.sge
{% endhighlight %}

To understand steps run by trinotate-JAv1.0.sh we can view a script generated from HNglobal fasta file as exemple :

{% highlight bash %}
more jobArray-Trinotate/sh/2_Trinotate.sh 
{% endhighlight %}

{% highlight bash %}

# Charging modules
module load bioinfo/Trinotate/3.0.1
module load bioinfo/TransDecoder/3.0.0
module load bioinfo/hmmer/3.1b2
module load bioinfo/diamond/0.7.11
 
# Defining scratch and destination repertories\n
pathToScratch="/scratch/orjuela/Trinotate_$JOB_ID.$SGE_TASK_ID/"
pathToDest="/repertory/containing/fastaFiles/jobArray-Trinotate/Trinotate"
mkdir -p $pathToScratch
mkdir -p $pathToScratch/DB
 
# Copie du fichier Trinity.fasta vers la partition /scratch du noeud
scp /repertory/containing/fastaFiles/HNglobal*.fasta $pathToScratch/
 
# Copie des bases uniprot_sprot*, Pfam-A.hmm*, et uniref90.fasta.dmnd vers la partition /scratch du noeud
scp /usr/local/Trinotate-3.0.1/uniprot_sprot.dmnd /usr/local/Trinotate-3.0.1/uniprot_sprot.pep /usr/local/Trinotate-3.0.1/uniprot_sprot.pep.phr /usr/local/Trinotate-3.0.1/uniprot_sprot.pep.pin /usr/local/Trinota
te-3.0.1/uniprot_sprot.pep.psq $pathToScratch/DB/
scp /data/projects/banks//uniref90.fasta.dmnd $pathToScratch/DB/
scp /usr/local/Trinotate-3.0.1/Pfam-A.hmm /usr/local/Trinotate-3.0.1/Pfam-A.hmm.h3f /usr/local/Trinotate-3.0.1/Pfam-A.hmm.h3i /usr/local/Trinotate-3.0.1/Pfam-A.hmm.h3m /usr/local/Trinotate-3.0.1/Pfam-A.hmm.h3p $
pathToScratch/DB/
 
cd $pathToScratch/ 
mkdir $pathToScratch/results_HNglobal
cd $pathToScratch/results_HNglobal/
 
# Running tool
 
# Calculing trinity_component_distribution
perl /usr/local/trinityrnaseq-2.5.1/util/misc/trinity_component_distribution.pl $pathToScratch/HNglobal.fasta 
cmd=" perl /usr/local/trinityrnaseq-2.5.1/util/misc/trinity_component_distribution.pl $pathToScratch/HNglobal.fasta "
echo "commande executee: $cmd"
 
# 1 getting gene to trans map
perl /usr/local/trinityrnaseq-2.5.1/util/support_scripts/get_Trinity_gene_to_trans_map.pl $pathToScratch/HNglobal.fasta > $pathToScratch/results_HNglobal/HNglobal.fasta_gene_trans_map 
cmd=" perl /usr/local/trinityrnaseq-2.5.1/util/support_scripts/get_Trinity_gene_to_trans_map.pl $pathToScratch/HNglobal.fasta \> $pathToScratch/results_HNglobal/HNglobal.fasta_gene_trans_map "
echo "commande executee: $cmd"
 
# 2 generation of peptide file
 
# 2.1 generation of longestOrf
TransDecoder.LongOrfs -t $pathToScratch/HNglobal.fasta --gene_trans_map $pathToScratch/results_HNglobal/HNglobal.fasta_gene_trans_map -m 50 
cmd=" TransDecoder.LongOrfs -t $pathToScratch/HNglobal.fasta --gene_trans_map $pathToScratch/results_HNglobal/HNglobal.fasta_gene_trans_map -m 50  "
echo "commande executee: $cmd"
 
# 2.2a recherche d’identité parmis les longorfs hmmscan
hmmscan --cpu 10 --domtblout pfam_longorfs.domtblout $pathToScratch/DB//Pfam-A.hmm $pathToScratch/results_HNglobal/HNglobal.fasta.transdecoder_dir/longest_orfs.pep 
cmd=" hmmscan --cpu 10 --domtblout pfam_longorfs.domtblout $pathToScratch/DB//Pfam-A.hmm $pathToScratch/results_HNglobal/HNglobal.fasta.transdecoder_dir/longest_orfs.pep  "
echo "commande executee: $cmd"
 
# 2.2b recherche d’identité parmis les longorfs diamond
 /usr/local/diamond-0.8.29/diamond blastp --query $pathToScratch/results_HNglobal/HNglobal.fasta.transdecoder_dir/longest_orfs.pep --db $pathToScratch/DB//uniprot_sprot --out diamP_uniprot_longorfs.outfmt6 --out
fmt 6 --max-target-seqs 1 
cmd=" /usr/local/diamond-0.8.29/diamond blastp --query $pathToScratch/results_HNglobal/HNglobal.fasta.transdecoder_dir/longest_orfs.pep --threads 10 --db $pathToScratch/DB//uniprot_sprot.pep --out diamP_uniprot_
longorfs.outfmt6 --outfmt 6 --max-target-seqs 1  "
echo "commande executee: $cmd"
 
# #2.3 Prediction peptides
 TransDecoder.Predict --cpu 10 -t $pathToScratch/HNglobal.fasta --retain_pfam_hits $pathToScratch/results_HNglobal/pfam_longorfs.domtblout --retain_blastp_hits $pathToScratch/results_HNglobal/diamP_uniprot_longo
rfs.outfmt6 
cmd=" TransDecoder.Predict --cpu 10 -t $pathToScratch/HNglobal.fasta --retain_pfam_hits pfam_longorfs.domtblout --retain_blastp_hits diamP_uniprot_longorfs.outfmt6  "
echo "commande executee: $cmd"
 
# 3 Recherche de similarité en utilisant Diamond
 
# blastp diamP_uniprott 
 /usr/local/diamond-0.8.29/diamond blastp --query $pathToScratch/results_HNglobal/HNglobal.fasta.transdecoder.pep --threads 10 --db $pathToScratch/DB//uniprot_sprot --out $pathToScratch/results_HNglobal/diamP_un
iprot.outfmt6 --outfmt 6 --max-target-seqs 1 --more-sensitive 
cmd=" /usr/local/diamond-0.8.29/diamond blastp --query $pathToScratch/results_HNglobal/HNglobal.fasta.transdecoder.pep --threads 10 --db $pathToScratch/DB//uniprot_sprot --out $pathToScratch/results_HNglobal/dia
mP_uniprot.outfmt6 --outfmt 6 --max-target-seqs 1 --more-sensitive  "
echo "commande executee: $cmd"
{% endhighlight %}


### Building a Sqlite database and report

The second bash script, `build_Sqlite_trinotate_database_and_report-JAv1.2.0.sh`https://github.com/julieaorjuela/scripts/blob/master/build_Sqlite_trinotate_database_and_report-JAv1.2.0.sh, needs as input the assembled transcrits and the repertory containing the whole of results obtained by trinotate-JAv1.0.sh in the last step.

{% highlight bash %}
qsub -q bioinfo.q -N reportTrinonate -V -b yes -cwd 'bash ~/scripts/build_Sqlite_trinotate_database_and_report-JAv1.2.0.sh -f /repertory/containing/fastaFiles/longestAGglobal-Trinity.fasta -r /repertory/containing/fastaFiles/jobArray-Trinotate/Trinotate/results_longestAGglobal-Trinity/'
{% endhighlight %}
