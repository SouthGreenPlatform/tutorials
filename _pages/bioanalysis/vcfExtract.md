---
layout: page
title: "vcfExtract"
permalink: /bioanalysis/vcfExtract/
tags: [ polymorphism, SNP, indel, vcf ]
description: Commands to manipulate VCF files
---
<table class="table-contact">
<tr>
<td><img width="100%" src="{{ site.url }}/images/trainings-linux.png" alt="" />
</td>
<td>
<h1> Commands to manipulate VCF files</h1><br />
This page describes a serie of tools and linux commands used to manipulate VCF files.
</td>
</tr>
</table>


We need, in this tutorial:
* 1 vcf file
* GATK tools
* bcftools

### Author(s)

| Authors  | Christine Tranchant-Dubreuil  |
| :------------- | :------------- |
| Research Unit | UMR DIADE   |
| Institut |  <img src="https://www.ird.fr/extension/ird/design/ird/images/picto/logo_ird.png" width="20%"> |


#### _Keywords_ : `gatk`,`bcftools`

#### _Files format_ : vcf

#### _Date_ : 10/03/2017


***

## Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Extracting list of samples from a vcf file with `grep`and `cut`commands](#sample-list)
- [Extracting a subset of samples from a multigenome vcf file with `GATK selectVariants`](#sample-from-vcf-gatk)
- [Extracting a subset of samples from a multigenome vcf file with `bcftools`](#sample-from-vcf-bcftools)
- [Calculating the nucleotide diversity from a vcf file with `vcftools`](#calculating-pi)

<!-- /TOC -->

***

<a name="sample-list"></a>
## Extracting list of samples from a vcf file with `grep`and `cut`commands

one line with all samples

{% highlight ruby %}
grep "#CHROM" output | cut -f 10-
{% endhighlight %}

one line by sample

{% highlight ruby %}
grep "#CHROM" output | cut -f 10- | xargs -n 1
{% endhighlight %}

<a name="sample-from-vcf-gatk"></a>
## Extracting a subset of samples from a multigenome vcf file with `GATK selectVariants`

#### Select two samples out of a vcf with many samples

{% highlight ruby %}
java -Xmx12g -jar /usr/local/gatk-3.6/GenomeAnalysisTK.jar -T SelectVariants -R reference.fa -V inputFileName.vcf -o outputFilename.vcf -sn sample1 -sn sample2
{% endhighlight %}

###### Rk : if you get the following error message "_Fasta dict file ... for reference ... does not exist_", please see https://www.broadinstitute.org/gatk/guide/article?id=1601 for help creating in.


#### Select genotypes from a file containing a list of samples to include

{% highlight ruby %}
java -Xmx12g -jar /usr/local/gatk-3.6/GenomeAnalysisTK.jar -T SelectVariants -R reference.fa -V inputFileName.vcf -o outputFileName.vcf --sample_file barthii.only.RG.list  --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES
{% endhighlight %}

#### Select genotypes from a file containing a list of samples to exclude

{% highlight ruby %}
java -Xmx12g -jar /usr/local/gatk-3.6/GenomeAnalysisTK.jar -T SelectVariants -R reference.fa -V inputFileName.vcf -o outputFileName.vcf --exclude_sample_file barthii.only.RG.list  --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES
{% endhighlight %}

###### Rk : if you get the following error message : "_Bad input: Samples entered on command line (through -sf or -sn)) that are not present in the VCF_", run with --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES


<a name="sample-from-vcf-bcftools"></a>
## Extracting a subset of samples from a multigenome vcf file with `bcftools`
#### Select genotypes from a file containing a list of samples to include with `bcftools`

{% highlight ruby %}
bcftools view -S barthii.only.RG.list inputFileName.vcf --force-samples -o outputFilename.vcf`
{% endhighlight %}


<a name="calculating-pi"></a>
## Calculating the nucleotide diversity from a vcf file with `vcftools`

{% highlight ruby %}
vcftools --vcf inputFilename.vcf  --out outputFilename.PI  --window-pi 100000 --remove-filtered-all
{% endhighlight %}

{% highlight ruby %}
grep "PI" OgOb-all-MSU7-CHR2.GATKSV.VCFTOOLS.stats-100000.windowed.pi -v | awk '{ sum+=$5; print $5,"; ",sum , "* ", NR ; } END { print "PI average :", sum / NR; }'
{% endhighlight %}
