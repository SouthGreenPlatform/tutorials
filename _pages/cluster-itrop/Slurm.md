---
layout: page
title: "Utilisation de slurm"
permalink: /cluster-itrop/Slurm/
tags: [ linux, HPC, cluster, OS ]
description:  Use of Slurm for i-Trop cluster
---

| Description | Know how to use Slurm|
| :------------- | :------------- | :------------- | :------------- |
| Author | Ndomassi TANDO (ndomassi.tando@ird.fr)  |
| Creation date |08/11/2019 |
| modification date | 08/11/2019 |


-----------------------


### Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->
* [Objective](#part-1)
* [Launch jobs with Slurm](#part-2)
* [Resources supervision with Slurm](#part-3)
* [Liens](#liens)
* [License](#license)


-----------------------
<a name="part-1"></a>
## Objectives

Know how to launch different type of jobs with Slurm.

Jobs monitoring with Slurm

-------------------------------------------------------------------------------------

<a name="part-2"></a>
## Launch jobs with Slurm:
       
### Launch commands from the master node

The following command allocate computing resources ( nodes, memory, cores) and immediately launch the command on each allocate resource.


    {% highlight bash %}$ srun + command{% endhighlight %} 
    
  Example:
   
    {% highlight bash %}$ srun hostname{% endhighlight %} 
  
  Allow to obtain the name of the computing resource used   
  
### Reserve computing resources to launch Slurm commands

We use:

{% highlight bash %}$ salloc{% endhighlight %}

This sommand allows to reserve one or several computing resources and work on the master node at the same time.

Commands on the computing resources can be launched later with the command `srun + arguments`.

When you use this command, it is important to precise a reservation time with the option --time


Example:We reserve 2 nodes( option -N) at the same time for 5 minutes and later we run the hostname command thanks to srun

{% highlight bash %}$ salloc --time=05:00 -N 2
$ srun hostname{% endhighlight %}

We obtain:

{% highlight bash %}$[tando@master0 ~]$ srun hostname
node21.alineos.net
node14.alineos.net
 {% endhighlight %}

### Connect to a node in interactive mode and launch commands:

To connect to a node in interactive mode for X minutes , use the following command:


{% highlight bash %}$ srun --time=X:00 --pty bash -i{% endhighlight %}

Then you can launch on this node without using the `srun` prefix 



### Main options for Slurm:

`salloc`, `srun` or `sbatch` can be used with the following options:

| actions | Option 
| :------------- | :------------- | 
|Choose a  partition |	-p [queue]| 
| Number of nodes to use | -N [min[-max]]|
| Number of cpus to use| -n [count]| 
| Time limitation|-t [min] ou -t [days-hh:mm:ss] |
| Precise a output file| -o [file_name] |
| Precise a error file| -e [file_name] |
| Combine  STDOUT et STDERR files| utiliser -o sans -e|
| Copy the environnement|	--export=[ALL , NONE , variables]|
|Type of notifications to send|	--mail-type=[events]|
|Send a mail|--mail-user=[address]|
|Job Name|--job-name=[name]|
| Relaunch job in case of problem|--requeue|
| Set the working dir|--workdir=[dir_name] |
| Memory size |--mem=[mem][M,G,T] ou-mem-per-cpu=[mem][M,G,T]|
| Charge to a account|	--account=[account]|
|Task per node|--tasks-per-node=[count]|
| cpus per task| --cpus-per-task=[count]|
|Job dependency|	--depend=[state:job_id]|
| Job host preference| --nodelist=[nodes] ET/OU --exclude=[nodes]|
| Job arrays|	--array=[array_spec]|
| Begin Time|--begin=YYYY-MM-DD[THH:MM[:SS]]|

### Launch jobs via a script

The batch mode allows to launch an analysis by following the steps described into a script. 

Slurm allows to use different types of  scripts such as bash, perl or python.

Slurm allocates the desired computing resources and launch analyses on these resources in background

To be interpreted by Slurm, the script should contain a specific header with all the keyword `#BATCH` to precise the Slurm options.

Slurm script example:

{% highlight bash %}#!/bin/bash
## Define the job name
#SBATCH --job-name=test
## Define the output file
#SBATCH --output=res.txt
## Define the number of tasks
#SBATCH --ntasks=1
## Define the execution time limit
#SBATCH --time=10:00
## Define 100Mo of memory per cpu
#SBATCH --mem-per-cpu=100
sleep 180 #launch a 3 minutes sleep {% endhighlight %}

To launch an analysis use the following command:

{% highlight bash %}$ sbatch script.sh{% endhighlight %}

With `script.sh` the name of the script to use

### Environment variables:

       SLURM_JOB_ID		The ID of the job allocation.
       SLURM_JOB_NAME		The name of the job.
       SLURM_JOB_NODELIST	List of nodes allocated to the job.
       SLURM_JOB_NUM_NODES	Number of nodes allocated to the job.
       SLURM_NTASKS		Number of CPU tasks in this job.
       SLURM_SUBMIT_DIR	The directory from which sbatch was invoked.

### Delete a job

{% highlight bash %}$ scancel <job_id>{% endhighlight %}

With `<job_id>`: the job number

----------------------------------------------------------------------------------------------

<a name="part-3"></a>
## Monitor resources:


### Get jobs infos :

use the command

{% highlight bash %}$ squeue {% endhighlight %}

To refresh the infos every 5 seconds

{% highlight bash %}$ squeue -i 5 {% endhighlight %}

Infos on a particular job:

{% highlight bash %}$ scontrol show job <job_id>{% endhighlight %}
  
With `<job_id>`: the job number  

Infos on the jobs of a particular user

{% highlight bash %}$ squeue -u <user> {% endhighlight %}

With `<user>`: the user login

More infos on jobs:

{% highlight bash %}$ sacct --format=JobID,elapsed,ncpus,ntasks,state,node{% endhighlight %}

Infos on resources used by a finished job 

{% highlight bash %}$ seff <job_id>  {% endhighlight %}

With `<job_id>`: the job number

You can add the following command at the end of your script to get infos of the jobs in your output file.

{% highlight bash %}$ seff $SLURM_JOB_ID{% endhighlight %}

### Get infos on partition:

Type the following command:

{% highlight bash %}$ sinfo {% endhighlight %}

It gives infos on partitions and nodes


To get more informations:

{% highlight bash %}$ scontrol show partitions {% endhighlight %}

 `scontrol show` can be used with nodes, user, account etc...

Konow the time limit for each partition:

{% highlight bash %}$sinfo -o "%10P %.11L %.11l"{% endhighlight %}

### Get infos on nodes

Type the command:

{% highlight bash %}$ sinfo -N -l {% endhighlight %}

Several states are possible:

- alloc : The node is fully uses

- mix : The node is partially used

- idle : No job is runing on the node

- drain : node is finishing the received jobs but it doesn't accept new ones ( when the node will be stopped for maintenance)

To obtain more informations :

{% highlight bash %}$ scontrol show nodes {% endhighlight %}

 
-----------------------

### Liens
<a name="liens"></a>

* Cours liés : [HPC Trainings](https://southgreenplatform.github.io/trainings/HPC/)


-----------------------

### License
<a name="license"></a>

<div>
The resource material is licensed under the Creative Commons Attribution 4.0 International License (<a href="http://creativecommons.org/licenses/by-nc-sa/4.0/">here</a>).
<center><img width="25%" class="img-responsive" src="http://creativecommons.org.nz/wp-content/uploads/2012/05/by-nc-sa1.png"/>
</center>
</div>
                  
 
