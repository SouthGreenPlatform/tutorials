---
layout: page
title: "RNASeq"
permalink: /bioanalysis/rnaSeq/
tags: [ rnaSeq, ]
description: RNA Seq differential expression
---
<table class="table-contact">
<tr>
<td><img width="100%" src="{{ site.url }}/images/trainings-linux.png" alt="" />
</td>
<td>
<h1> RNA Seq differential expression</h1><br />
This page describes a serie of tools and linux commands used to manipulate raw data (fastq file) for differential expression.
</td>
</tr>
</table>


We need, in this tutorial:
* a directory with fastq files
* a reference file used for the mapping step.


### Author(s)

| Authors  | Gaetan Droc  |
| :------------- | :------------- |
| Research Unit | UMR DIADE   |
| Institut |  <img src="http://s3f-haiti.cirad.fr/var/projets/storage/images/media/media_s3f_haiti/logo_agap/54594-1-fre-FR/logo_agap_medium.jpg" width="20%"> |

### Keywords
fastqc, cutadapt, hitsat2, stringtie, ballgown, bedtools

### Files format
fastq, sam, bam, bed

### Date
16/03/2017



***

## Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Performing a quality control check with `fastqc`](#fastqc)
- [Using `cutadapt` to remove adapters and to trim reads based on quality](#cutadapt)
- [Mapping reads with `hisat2`](#hisat2)
- [Convert and sort SAM to BAM with `samtools`](#samtools)
- [Transcript assembly and quantification with `StringTie`](#stringtie)


<!-- /TOC -->

***
<a name="fastqc"></a>
## Performing a quality control check with `fastqc`

`FastQC` perform some simple quality control checks to ensure that the raw data looks good and there are no problems or biases in data which may affect how user can usefully use it.

http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

### `fastqc` command

For one file

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ mkdir FastQC
[droc@cc2-login RNASeq_fastq_MGX]$ qsub -b y -q normal.q -N fastqc -V fastqc --format fastq --outdir FastQC Leaf_CGATGT_L006_R2_017.fastq
{% endhighlight %}

Loop on every file for a directory and print command on terminal

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ for i in *fastq; do echo qsub -b y -q normal.q -N fastqc -V fastqc --format fastq --outdir ~/work/FastQC $i;done
qsub -b y -q normal.q -N fastqc -V fastqc --format fastq --outdir /homedir/droc/work/FastQC Leaf_CGATGT_L006_R1_001.fastq
qsub -b y -q normal.q -N fastqc -V fastqc --format fastq --outdir /homedir/droc/work/FastQC Leaf_CGATGT_L006_R1_002.fastq
qsub -b y -q normal.q -N fastqc -V fastqc --format fastq --outdir /homedir/droc/work/FastQC Leaf_CGATGT_L006_R1_003.fastq
qsub -b y -q normal.q -N fastqc -V fastqc --format fastq --outdir /homedir/droc/work/FastQC Leaf_CGATGT_L006_R1_004.fastq
qsub -b y -q normal.q -N fastqc -V fastqc --format fastq --outdir /homedir/droc/work/FastQC Leaf_CGATGT_L006_R1_005.fastq
qsub -b y -q normal.q -N fastqc -V fastqc --format fastq --outdir /homedir/droc/work/FastQC Leaf_CGATGT_L006_R1_006.fastq
{% endhighlight %}

Loop on every file for a directory and write command on a file (cmd_fastqc.sh)

`[droc@cc2-login RNASeq_fastq_MGX]$ for i in *fastq; do echo qsub -b y -q normal.q -N fastqc -V fastqc --format fastq --outdir ~/work/FastQC $i >> cmd_fastqc.sh;done`

`[droc@cc2-login RNASeq_fastq_MGX]$ ./cmd_fastqc.sh `

HTML Report

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ cd FastQC
Leaf_CGATGT_L006_R2_017_fastqc.html  Leaf_CGATGT_L006_R2_017_fastqc.zip
{% endhighlight %}

<a name="cutadapt"></a>
## Using `cutadapt` to remove adapters and to trim reads based on quality

`Cutadapt` is a tool specifically designed to remove adapters from NGS data. 
https://code.google.com/p/cutadapt/

### `cutadadapt` command

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ module load bioinfo/cutadapt/1.8.1
[droc@cc2-login RNASeq_fastq_MGX]$ for i in *fastq; do echo qsub -q normal.q -b yes -V -N CUTADAPT cutadapt -a AGATCGGAAGAGCG -O 10 -q 30,30 -f fastq -m 30 -o /work/droc/sugarcane/cutadapt/$i.cutadapt.fastq /work/NGSwaiting4newNAS/sugarcane/BackupCarine/RNASeq/RNASeq_fastq_MGX/$i >> cmd_cutadapt.sh ;done
./cmd_cutadapt.sh
{% endhighlight %}

-q 30, 30 : by default, only the 3’ end of each read is quality-trimmed. If you want to trim the 5’ end as well, use the -q option with two comma-separated cutoffs

<a name="hisat2"></a>
## Mapping reads with `hisat2`

https://ccb.jhu.edu/software/hisat2/index.shtml

There are several steps involved in mapping sequence reads.

### Creating an index of the reference genome if necessary

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ module load bioinfo/hisat2/2.0.5
[droc@cc2-login RNASeq_fastq_MGX]$ qsub -V -b y -q normal.q -N index hisat2-build reference_genome reference_genome
{% endhighlight %}

### Performing mapping

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ qsub -V -b y -q bigmem.q -N hisat2 hisat2 -x reference_genome -1 Leaf_CGATGT_L006_R1_001.fastq -2  Leaf_CGATGT_L006_R2_001.fastq -S Leaf_CGATGT_L006.sam
{% endhighlight %}

Here is a description for the contents of the SAM file: https://samtools.github.io/hts-specs/SAMv1.pdf

Do the same thing for all the library. For this you can use and adapt this Perl script



{% highlight perl %}
[droc@cc2-login RNASeq_fastq_MGX]$ nedit fastq2tab.pl&

#!/usr/bin/perl
use File::Basename;
# Get the directory of fastq
my $pwd = shift;

# List all of the file
open(IN,"ls $pwd/*fastq|");
my %fastq;
while(<IN>){
    chomp;
    my $file = fileparse($_);
    # In this example, $file = Leaf_CGATGT_L006_R1_006.fastq
    # Split the name by _
    my ($tissu,$library,$strand,$number) = (split(/\_/,$file))[0,2,3,4];
    # Create an uniq identifier $tag
    my $tag = join("_",$tissu,$library,$number);
    # Create a hash with 2 keys, the uniq name ($tag) and strand (R1 or R2)
    $fastq{$tag}{$strand} = $_;
}
# Next foreach
foreach my $tag (keys %fastq) {
    if ($fastq{$tag}{R1} && $fastq{$tag}{R2}) {
        print "qsub -V -b y -q bigmem.q -N hisat2 hisat2 -x reference_genome -1 $fastq{$tag}{R1} -2  $fastq{$tag}{R2} -S $tag.sam\n";
    }
}

[droc@cc2-login RNASeq_fastq_MGX]$ perl fastq2tab.pl /work/NGSwaiting4newNAS/sugarcane/BackupCarine/RNASeq/RNASeq_fastq_MGX/ > cmd_hisat2.sh
[droc@cc2-login RNASeq_fastq_MGX]$ ./cmd_hisat2.sh
{% endhighlight %}

<a name="samtools"></a>
## Convert and sort SAM to BAM with `samtools`
http://samtools.sourceforge.net/


{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ module load bioinfo/samtools/1.3
[droc@cc2-login RNASeq_fastq_MGX]$ qsub -V -q normal.q -N samtools -b y "samtools view -bS Leaf_CGATGT_L006.sam -o Leaf_CGATGT_L006.bam"
[droc@cc2-login RNASeq_fastq_MGX]$ qsub -V -q normal.q -N samtools -b y "samtools sort Leaf_CGATGT_L006.bam -o Leaf_CGATGT_L006_sorted.bam"
{% endhighlight %}

To run on batch, you can use the command for

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ for i in *sam; do echo qsub -V -q normal.q -N samtools -b y samtools view -bS $i -o $i.bam >> cmd_samtools.sh;done;
[droc@cc2-login RNASeq_fastq_MGX]$ rm -f *sam
[droc@cc2-login RNASeq_fastq_MGX]$ for i in *bam; do echo qsub -V -q normal.q -N samtools -b y samtools sort $i -o $i.sorted.bam >> cmd_samtools_sorted.sh;done;
{% endhighlight %}

<a name="stringtie"></a>
## Transcript assembly and quantification with `StringTie`
https://ccb.jhu.edu/software/stringtie/

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ module load bioinfo/stringtie/1.2.1
[droc@cc2-login RNASeq_fastq_MGX]$ qsub -V -b y -q normal.q -N stringtie stringtie  -p 16 -o Leaf_CGATGT_L006.gtf Leaf_CGATGT_L006_sorted.bam
{% endhighlight %}

For all sample

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ for i in *sorted.bam; do echo qsub -V -b y -q normal.q -N stringtie stringtie  -p 16 -o $i.gtf $i >> cmd_stringtie.sh;done
{% endhighlight %}

### Merge transcripts from all sample

{% highlight ruby %}
[droc@cc2-login RNASeq_fastq_MGX]$ qsub -V -b y -q normal.q -N merge stringtie --merge -p 8 -o stringtie_merged.gtf mergelist.txt
{% endhighlight %}

Here mergelist.txt is a text file that has the names of the gene transfer format (GTF) files created in the previous step, with each file name on a single line.

