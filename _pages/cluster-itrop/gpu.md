---
layout: page
title: "Utilisation du noeud GPU"
permalink: /cluster-itrop/gpu/
tags: [linux, HPC, cluster, GPU ]
description:  Use of GPU node for i-Trop cluster
---

| Description | Know how to use GPU node in I-Trop cluster |
| :------------- | :------------- | :------------- | :------------- |
| Author | Julie ORJUELA (julie.orjuela_at_ird.fr) and Aurore COMTE (aurore.comte_at_ird.fr) |
| Creation date |27/01/2020 |
| modification date | 27/01/2020 |


-----------------------


### Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->
* [Objective](#part-1)
* [Launch jobs in GPU node with Slurm](#part-2)
* [Resources supervision with nvidia](#part-3)
* [Liens](#liens)
* [License](#license)


-----------------------
<a name="part-1"></a>
## Objectives

Know how to launch a Slurm job in GPU node in I-Trop Cluster and monitoring this jobs in GPU

-------------------------------------------------------------------------------------

<a name="part-2"></a>

## Launch jobs in GPU node with Slurm 

Node GPU in I-trop cluster has 8 graphic cards RTX2080, each with 124G de RAM. In total this node has 24 threads.
We recommend to basecaller a data set using a graphic card to obtain results in only one folder. If you split data you can enjoy of the whole of graphic cards but your data results will be in several folders. In each results folder, reads can be share names. So, you can lost information if you decide to merge it.

Create a sbatch script to allocate ressources. Here, slurm script `lauchGuppyGPU.sbash` takes 4 threads for lauch guppy-gpu, partition `-p gpu`. If you are using i-Trop GPU you are into `gpu_group` so, give this parametter to slurm whit `-A ` option.

Guppy is a data processing toolkit that contains the Oxford Nanopore Technologies’ basecalling algorithms, and several bioinformatic post-processing features.

Basecalling with guppy can be launch using gyppy-gpu tool. In guppy commande you have to specify data containig raw read files (fast5) (-i), the output repertory to write fastq files (-o), How many worker threads you are using	–cpu_threads_per_caller	(-c) and the number of parallel basecallers to create	(-num_callers), we recommend to compress the fastq output (-compress_fastq)


{% highlight bash %} 
#!/bin/bash
#SBATCH -J Basecalling
#SBATCH -p gpu
#SBATCH -A gpu_group
#SBATCH -c 4
INPUT=$1
OUTPUT=$2
CUDA=$3

#loading modules
module load bioinfo/guppy-gpu/3.2.4

#running basecalling
guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i ${INPUT} -r -s ${OUTPUT} --num_callers 4 --gpu_runners_per_device 8 --qscore_filtering --min_qscore 7 -x cuda:${CUDA}
The following command allocate computing resources ( nodes, memory, cores) and immediately launch the command on each allocate resource.

Launch lauchGuppyGPU.sbash script:

{% highlight bash %}$ sbatch lauchGuppyGPU.sbash {% endhighlight %} 
  
Note:
Beside the path of our fast5 files (-i), the basecaller requires an output path (-s) and a config file or the flowcell/kit combination. In order to get a list of possible flowcell/kit combinations and config files, we use:

{% highlight bash %}$ guppy_basecaller --print_workflows {% endhighlight %}
  
  
  
  <a name="part-3"></a>
## Resources supervision with nvidia
  
  
### Liens
<a name="liens"></a>

* Cours liés : [Slurm Trainings](https://southgreenplatform.github.io/tutorials//cluster-itrop/Slurm/)


-----------------------

### License
<a name="license"></a>

<div>
The resource material is licensed under the Creative Commons Attribution 4.0 International License (<a href="http://creativecommons.org/licenses/by-nc-sa/4.0/">here</a>).
<center><img width="25%" class="img-responsive" src="http://creativecommons.org.nz/wp-content/uploads/2012/05/by-nc-sa1.png"/>
</center>
</div>
