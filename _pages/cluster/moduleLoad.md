---
layout: page
title: "SGE"
permalink: /cluster/moduleLoad/
tags: [ HPC, SGE, cluster, module load, module unload ]
description: Build own module load
---
<table class="table-contact">
<tr>
<td><img width="100%" src="{{ site.url }}/images/trainings-linux.png" alt="" />
</td>
<td>
<h1> Build own module load</h1><br />
This page describes how you can build own module load.
</td>
</tr>
</table>



| Authors  | Sébastien RAVEL  |
| :------------- | :------------- |
| Research Unit | <img src="http://printemps-baillarguet.e-monsite.com/medias/images/bgpi.jpg" width="10%">    |
| Institut | <img src="{{ site.url }}/images/logo-cirad.png" width="10%">  |

#### _Keywords_ : `module load, private, own, cc2-login`

#### _Date_ : 12/06/2017


## Create and use the personal module directory

{% highlight ruby %}
mkdir $HOME/privatemodules
{% endhighlight %}

add permanently into the .bashrc

{% highlight ruby %}
#use own module
module use --append $HOME/privatemodules
{% endhighlight %}

## Module command

* module avail
* module help
* module whatis

## Example of file module load (name toggleDev)

Into __$HOME/privatemodules__ add new file __toggleDev__

{% highlight ruby %}
vim $HOME/privatemodules/toggleDev
{% endhighlight %}

{% highlight ruby %}
#%Module1.0
##

## Required internal variables
set     prefix       $env(HOME)/TOGGLE-DEV/
set     version      "TOGGLE-DEV"

if {![file exists $prefix]} {
	puts stderr "\t[module-info name] Load Error: $prefix does not exist"
	break
	exit 1
}

## List conflicting modules here
conflict toggleMerge

## List prerequisite modules here
set		fullname	TOGGLE-DEV
set		externalurl	"https://github.com/SouthGreenPlatform/TOGGLE-DEV\n"
set		description	"A framework to quickly build pipelines and to perform large-scale NGS analysis"

## Required for "module help ..."

proc ModulesHelp { } {
  global description externalurl
  puts stderr "Description - $description"
  puts stderr "More Docs   - $externalurl"
}

## Required for "module display ..." and SWDB
module-whatis   "loads the [module-info name] environment"

## Software-specific settings exported to user environment
module load system/java/jre-1.8.111
module load system/perl/5.24.0
module load bioinfo/FastQC/0.11.5
module load bioinfo/bwa/0.7.12
module load bioinfo/picard-tools/2.5.0
module load bioinfo/samtools/1.3.1
module load bioinfo/gatk/3.6
module load bioinfo/cutadapt/1.10
module load system/libgtextutils/0.7
module load bioinfo/fastx_toolkit/0.0.14
module load bioinfo/tophat/2.1.1
module load bioinfo/bowtie/1.1.2
module load bioinfo/bowtie2/2.2.9
module load bioinfo/cufflinks/2.2.1
module load bioinfo/tgicl_linux/1.0
module load bioinfo/trinityrnaseq/2.2.0
module load bioinfo/stacks/1.43
module load bioinfo/snpEff/4.3
module load bioinfo/ngsutils/0.5.9

# Path
prepend-path PATH $prefix
prepend-path PERL5LIB $prefix/modules:$prefix/test/pipelines
prepend-path TOGGLE_PATH $prefix

{% endhighlight %}
