---
layout: page
title: "SGE"
permalink: /cluster/SGE/
tags: [ HPC, SGE, cluster ]
description: How use SGE on HPC
---
<table class="table-contact">
<tr>
<td><img width="100%" src="{{ site.url }}/images/trainings-linux.png" alt="" />
</td>
<td>
<h1> SGE Job submission for cc2 login cluster</h1><br />
This page describes how you can submit job with SGE system.
</td>
</tr>
</table>



| Authors  | Sébastien RAVEL  |
| :------------- | :------------- |
| Research Unit | <img src="http://printemps-baillarguet.e-monsite.com/medias/images/bgpi.jpg" width="10%">    |
| Institut | <img src="{{ site.url }}/images/logo-cirad.png" width="10%">  |


#### _Keywords_ : `qsub, qrsh, job, cc2-login`

#### _Date_ : 02/06/2017

<a name="summary"></a>
## Summary

- [Summary](#summary)
- [How to run correctly Jobs](#how-to-run-correctly-jobs)
- [1 - Use Module load](#1-use-module-load)
	- [A - Load modules by default on connection](#a-load-modules-by-default-on-connection)
	- [B - Loading before job submission (before qsub)](#b-loading-before-job-submission-before-qsub)
	- [C - Loading in Job Script](#c-loading-in-job-script)
- [2 - Run Job](#2-run-job)
	- [A - qsub mode](#a-qsub-mode)
	- [B - qrsh mode (for test script)](#b-qrsh-mode-for-test-script)
- [3 - Knowing / asking for resources](#3-knowing-asking-for-resources)
- [4 - Jobs Resources](#4-jobs-resources)
	- [A - Request more RAM](#a-request-more-ram)
	- [B - Request more Threads](#b-request-more-threads)
- [5 - Get information about running jobs](#5-get-information-about-running-jobs)
- [6 - Delete job](#6-delete-job)
- [7 - More infos](#7-more-infos)


<a name="how-to-run-correctly-jobs"></a>
## How to run correctly Jobs

**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Not starting a program on the master node of the cluster: (cf: cc2-admin)**

<a name="1-use-module-load"></a>
## 1 - Use Module load

On the cluster the tools do not load by default.
Each user must "load" programs to use them.
This system allows you to have different versions of the same tools without creating conflicts.
The disadvantage is that if you forget to load the module the program does not work ...
Here are 3 methods to manage modules:

<a name="a-load-modules-by-default-on-connection"></a>
### A - Load modules by default on connection

If you often use the same program you can add to the .bashrc file
For example the program GIT or python2.7 that I often use:

{% highlight ruby %}
gedit /homedir/{YourName}/.bashrc
# Copy/Paste the two lines:
# add module load permanently
module load system/git/2.8.3 system/python/2.7.9
# Save file, quit gedit and actualise connection with: (or kill and open new connection):
source ~/.bashrc
{% endhighlight %}

<a name="b-loading-before-job-submission-before-qsub"></a>
### B - Loading before job submission (before qsub)

This method loads modules before qsub.
It will be **MANDATORY** to put the **-V** parameter in the qsub to pass the modules to the compute node.

{% highlight ruby %}
##### load modules
module load mpi/openmpi/1.6.5 compiler/gcc/4.9.2 bioinfo/RAxML/8.2.4
##### run job
qsub -V -N raxmlALL -cwd -q long.q -pe parallel_smp 24 ./run_raxml.sh
{% endhighlight %}

<a name="c-loading-in-job-script"></a>
### C - Loading in Job Script

The best way for the job to work. (We always forget to load before qsub ....)
It's simple enough to have the modules loaded by the compute node.
It is enough to add in the script.sh the module load before the control of the program.
Example in file **__run_raxml.sh__**:

{% highlight ruby %}
module load mpi/openmpi/1.6.5 compiler/gcc/4.9.2 bioinfo/RAxML/8.1.17
raxmlHPC-PTHREADS -T 24 -n 2241Ortho-82souches -f d -m GTRGAMMA -p 12345 -s /work/sravel/phylogenomique3/raxmlConcat/2241Ortho-82souches.fasta
{% endhighlight %}

{% highlight ruby %}
##### run job
qsub -V -N raxmlALL -cwd -q long.q -pe parallel_smp 24 ./run_raxml.sh
{% endhighlight %}

PS: leave the -V argument in your job, because sometimes it is necessary despite the loading of the modules in the script.

<a name="2-run-job"></a>
## 2 - Run Job

<a name="a-qsub-mode"></a>
### A - qsub mode

To submit a job you must use the command **qsub**

{% highlight ruby %}
# For normal job
qsub -V -N NomJob -cwd -q long.q script.sh

# example if the program command is in a script.sh (here run_Lorma.sh)
qsub -V -N Lorma -cwd -q long.q run_Lorma.sh

# OR you can run the command directly (to avoid because no trace of the parameters)
qsub -V -N Lorma -cwd -q long.q /work/sravel/lorma.sh -s /work/sravel/MinION/Minion_Sanger/eBSMYV_SeqSangerBAC29H14.fasta
{% endhighlight %}

With the following arguments:

| Arguments submission Jobs ||
| :------------- | :------------- |
| **-i file** | Use file as standard input for this job |
| **-o file** | Set the job's standard output to a file (which should be displayed in the terminal) |
| **-e file** | Set the standard error output of the job to the file (when returns error) |
| **-N name** | Name the job for output files (replaces STDIN by default) |
| **-cwd** | Uses the current working directory for input and output, rather than /homedir/username/ |
| **-q queue** | Specifies a queue |
| **-M my address@work** | Receive by e-mail job info:
||**-m beas**: Allows you to select events to receive:    |
||-b: warned at first   |
||-e: end of the job    |
||-a: job interrupted   |
||-s: job suspended   |
| **-l mem_free=nG** | run job with "n" Go de RAM (see below) |
| **-pe parallel_smp n** | run job with "n" threads (see below) |

<a name="b-qrsh-mode-for-test-script"></a>
### B - qrsh mode (for test script)

There is a method for not running tests on the master node:
Request an **interactive job**
This type of job brings advantages but also disadvantages:
* **Benefits**: allows you to test scripts to debug programs
* **Disadvantages**: the job and kill if the terminal and closed

It is therefore necessary to use it ONLY for the debug.

The command is similar to the * qrsh * and takes at least the -q argument:

{% highlight ruby %}
# For interactive job
qrsh -q normal.q
{% endhighlight %}

<a name="3-knowing-asking-for-resources"></a>
## 3 - Knowing / asking for resources

In computer science there are 2 types of resources to make a program:
* CPU / Thread
* RAM

The CPU is the number of cores of a machine, and each core divides into threads.
A program uses the CPU when it needs to do a lot of calculations.

RAM, express in Go or To corresponds to the active memory of the computer. It preloads information to perform the calculation more quickly.
(The access time is much faster than reading from disk)
A program that needs to load files in memory uses more RAM

Once this theory is understood, one can see in the programs parameters such as
* CPU
* Thread
* Java -Xmx18G
* ...

These are the famous parameters to increase the computing power.

On the Cluster the available resources are quite important:

| queue | NB Threads | NB RAM |
| :------------- | :------------- | :------------- |
|  long.q and normal.q | 48 | 200Go |
|  bigmem.q | 96 | 2.6To |

<a name="4-jobs-resources"></a>
## 4 - Jobs Resources

By default a job uses 1 thread and 10GB of RAM

{% highlight ruby %}
# for njob with default parameters
qsub -V -N NomJob -cwd -q long.q script.sh
{% endhighlight %}

To boost the job you have to ask for more resources.

<a name="a-request-more-ram"></a>
### A - Request more RAM

{% highlight ruby %}
# for a job with 20Go of RAM:
qsub -V -N NomJob -cwd -q long.q -l mem_free=20G script.sh
{% endhighlight %}

<a name="b-request-more-threads"></a>
### B - Request more Threads

{% highlight ruby %}
# for a parallele job with 24 threads:
qsub -V -N NomJob -cwd -q long.q -pe parallel_smp 24 script.sh
{% endhighlight %}

<a name="5-get-information-about-running-jobs"></a>
## 5 - Get information about running jobs

| Jobs status||
| :------------- |:-------------|
| **qstat** | Displays the status of all jobs |
| **qstat -f** | Displays the status of all queues (long list) |
| **qstat -u "*"** | Displays the status of all jobs belonging to all users |
| **qstat -g c** | Resources available |
| **qstat -j jobid** | Displays the status of a particular job (jobid = 1st qstat column) |

<a name="6-delete-job"></a>
## 6 - Delete job


{% highlight ruby %}
# delete one job
qdel jobID

# delete  multiple jobs:
qdel echo `seq -f "%.0f" 876775 876778`

# where 876775 is the first job to delete and 876778 the last
{% endhighlight %}

<a name="7-more-infos"></a>
## 7 - More infos

https://doc.cc.in2p3.fr/en:ge_submit_a_job_qsub
