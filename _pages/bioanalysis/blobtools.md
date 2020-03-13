---
layout: page
title: "Qiime2 en ligne de Commande"
permalink: /bioanalysis/qiime2/
tags: [ dada2, qiime2 ]
description: Metabarcoding using Qiime2 and DADA2
---
<table class="table-contact">
<tr>
<td><img width="100%" src="{{ site.url }}/images/qiime2.png" alt="" />
</td>
<td>
<h1> Metabarcoding using Qiime2 and DADA2</h1><br />
This page describes how to analyse metabarcoding data using Qiime2 and DADA2. From fastq files to SVs table obtention and Phyloseq analysis.
</td>
</tr>
</table>



### Author(s)

| Authors  | Julie Orjuela |
| :------------- | :------------- |
| Research Unit | UMR IPME-DIADE   |
| Institut |  IRD |

### Keywords
contamination, blobtools, taxonomy

### Files format
fasta, diamond

### Date
13/03/2020

### Checking contaminations with bloobtools in your assemblies

#### diamond in assemblies against a uniprot database

{% highlight bash %}
[orjuela@node25 bloobtools]$ more diamond.slurm 
#!/bin/bash
#SBATCH --export=ALL
#SBATCH -J diamondblastX
#SBATCH -n 8                         # coeurs
#SBATCH --mem-per-cpu=5GB           # mémoire
#SBATCH -t 2-00:00                  # durée job (D-HH:MM)
#SBATCH -o diamondblastX.%N.%j.out        # STDOUT
#SBATCH -e diamondblastX.%N.%j.err        # STDERR
#SBATCH --partition supermem	     # partition

#modules 
module load bioinfo/diamond/0.9.30

# diamond blastX contre la bd reformate uniprot voir https://github.com/blobtoolkit/blobtools2#additional-dependencies partie 4
diamond blastx \
 --query /scratch/orjuela/bloobtools/mil_best_assembly.fasta \
 --db /scratch/orjuela/bloobtools/reference_proteomes.dmnd \
 --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
 --sensitive \
 --max-target-seqs 1 \
 --evalue 1e-25 \
 --threads 8 
 --out diamond 
 #pd: prochaine fois ajouter --out sinon il ecrit dans le log.out

{% endhighlight %}

