---
layout: page
title: "blobtools en ligne de commande"
permalink: /bioanalysis/blobtools
tags: [blobtools]
description: Checking contaminations with blobtools in your assemblies
---
<table class="table-contact">
<tr>
<td><img width="100%" src="{{ site.url }}/images/blob.png" alt="" />
</td>
<td>
<h1> Checking contaminations with blobtools in your assemblies</h1><br />
bloobtools is nice!
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
 # --out sinon il ecrit dans le log.out
{% endhighlight %}

#### minimap2 assembly vs your fastq ONT.

{% highlight bash %}
#!/bin/bash
#SBATCH --export=ALL
#SBATCH -J minimap2
#SBATCH -n 10                         # coeurs
#SBATCH --mem-per-cpu=40GB           # mémoire
#SBATCH -t 10-00:00                  # durée job (D-HH:MM)
#SBATCH -o minimap2.%N.%j.out        # STDOUT
#SBATCH -e minimap2.%N.%j.err        # STDERR
#SBATCH --partition supermem	     # partition

#modules
module load bioinfo/minimap2/2.16 
module load bioinfo/samtools/1.9

# mapping
echo "minimap2 -x map-ont -m 250 --MD -t 20 --no-long-join -r 50 -a /scratch/orjuela/bloobtools/mil_best_assembly.fasta ALL-ONT-fastqpass-HAC.fastq | samtools sort -@ 4 -O BAM -o minimap2_assembly_vs_fastqONT.ba
m"

minimap2 -x map-ont -m 250 --MD -t 20 --no-long-join -r 50 -a /scratch/orjuela/bloobtools/mil_best_assembly.fasta /scratch/orjuela/bloobtools/ALL-ONT-fastqpass-HAC.fastq | samtools sort -@ 4 -O BAM -o minimap2_a
ssembly_vs_fastqONT.bam

# echo "END"
{% endhighlight %}


#### Lauch bloobtools

{% highlight bash %}
# chargement des modules
module load bioinfo /samtools/1.9
module load system/singularity/3.3.0

# sort and index your bam
samtools sort -o minimap2_assembly_vs_fastqONT_sorted.bam minimap2_assembly_vs_fastqONT.bam
samtools index minimap2_assembly_vs_fastqONT_sorted.bam      

# run singularity blobtools v1
mkdir blob1; cd blob1
singularity run /data3/projects/containers/CULEBRONT/bloobtools-v1.1.1.simg blobtools create -i ../mil_best_assembly.fasta -b ../minimap2_assembly_vs_fastqONT_sorted.bam -t ../diamondblastX.csv -o milbloob1 --names /usr/local/blobtools-blobtools_v1.1.1/data/names.dmp --nodes /usr/local/blobtools-blobtools_v1.1.1/data/nodes.dmp
singularity run /data3/projects/containers/CULEBRONT/bloobtools-v1.1.1.simg blobtools view -i milbloob1.blobDB.json --cov -o mil
singularity run /data3/projects/containers/CULEBRONT/bloobtools-v1.1.1.simg plot -i milbloob1.blobDB.json
{% endhighlight %}

