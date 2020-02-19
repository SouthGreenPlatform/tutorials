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
| modification date | 19/02/2020 |


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

{% highlight bash %}$ salloc --time=05:00 -N 2 -p short
$ srun hostname{% endhighlight %}

We obtain:

{% highlight bash %}$[tando@master0 ~]$ srun hostname
node21.alineos.net
node14.alineos.net
 {% endhighlight %}

### Connect to a node in interactive mode and launch commands:

To connect to a node in interactive mode for X minutes , use the following command:


{% highlight bash %}$ srun -p short --time=X:00 --pty bash -i{% endhighlight %}

Then you can launch on this node without using the `srun` prefix 


### Connect to a node in interactive mode with x11 support:

The x11 support allows you to launch graphical software within  a node.

You first have to connect to the bioinfo-master.ird.fr with the -X option:

{% highlight bash %}$ ssh -X login@bioinfo-master.ird.fr {% endhighlight %}

Then you can launch this command with the `--x11` option

{% highlight bash %}$ srun -p short --x11 --pty bash -i{% endhighlight %}

### Partitions available:

Depending on the type of jobs you want to launch you have the choice between several partitions.
 
 The partitions can be considered job queues, each of which has an assortment of constraints such as job size limit, job time limit, users permitted to use it, etc. 
 
 Priority-ordered jobs are allocated nodes within a partition until the resources (nodes, processors, memory, etc.) within that partition are exhausted.
 
 Here are the available partitions:

 | partition | role  | nodes list | Number of Cores | Ram on nodes
| :------------- | :------------- | :------------- |:------------- |:------------- |
|short|	Short Jobs < 1 day (higher priority,interactive jobs)| node0,node1,node2,node13,node14| 12 cores | 48 to 64 GB
| normal | job of maximum 7 days| node0,node1,node2,node5,node13,node14,node15,node16,node17,node18,node19,node20,node22,node23,node24 | 12 to 24 cores| 64 to 96GB|
| long| <7 days< long jobs< 45 days| node3,node8,node9,node10,node11,node12|12 to 24 cores| 48 GB|
| highmem| jobs with more memory needs |node4, node7,node17,node21| 12 to 24 cores| 144 GB|
| supermem| jobs with much more memory needs|  node25| 40 cores | 1 TB|
| gpu |Need of analyses on GPU cores| node26| 24 cpus and 8 GPUS cores | 192 GB|


Note that the gpu node access is restricted, a request access should be done here: [request access to gpu](https://itrop-glpi.ird.fr/plugins/formcreator/front/formdisplay.php?id=15)

### Main options for Slurm:

`salloc`, `srun` or `sbatch` can be used with the following options:

| actions | Slurm Options | SGE options
| :------------- | :------------- | :------------- |
|Choose a  partition |	-p [queue]|  -q [queue]|
| Number of nodes to use | -N [min[-max]]|N/A |
| Number of tasks to use| -n [count]| -pe [PE] [count]|
| Time limitation|-t [min] ou -t [days-hh:mm:ss] |-l h_rt=[seconds]|
| Precise a output file| -o [file_name] |  -o [file_name] |
| Precise a error file| -e [file_name] | -e [file_name] |
| Combine  STDOUT et STDERR files| utiliser -o sans -e| -j yes |
| Copy the environnement|	--export=[ALL , NONE , variables]|
|Type of notifications to send|	--mail-type=[events]| -V|
|Send a mail|--mail-user=[address]| -M [address] |
|Job Name|--job-name=[name]| -N [name] |
| Relaunch job in case of problem|--requeue| -r [yes,no]
| Set the working dir|--workdir=[dir_name] | -wd [directory] |
| Memory size |--mem=[mem][M,G,T] ou-mem-per-cpu=[mem][M,G,T]| -l mem_free=[memory][K,M,G]|
| Charge to a account|	--account=[account]| -A [account] |
|Tasks per node|--tasks-per-node=[count]| (Fixed allocation_rule in PE)|
| cpus per task| --cpus-per-task=[count]| N/A|
|Job dependency|	--depend=[state:job_id]| -hold_jid [job_id , job_name]|
| Job host preference| --nodelist=[nodes] ET/OU --exclude=[nodes]| -q [queue]@[node] OR -q 
[queue]@@[hostgroup] |
| Job arrays|	--array=[array_spec]| -t [array_spec]|
| Begin Time|--begin=YYYY-MM-DD[THH:MM[:SS]]| -a [YYMMDDhhmm] |

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

### Submit  a array job

{% highlight bash %}#!/bin/bash
#SBATCH --partition=short      ### Partition
#SBATCH --job-name=ArrayJob    ### Job Name
#SBATCH --time=00:10:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --ntasks=1             ### Number of tasks per array job
#SBATCH --array=0-19%4           ### Array index from 0  to 19 with 4 runnings jobs

 
echo "I am Slurm job ${SLURM_JOB_ID}, array job ${SLURM_ARRAY_JOB_ID}, and array task ${SLURM_ARRAY_TASK_ID}."{% endhighlight %}

You  have to  use the `$SBATCH --array`option to  define the range

The variable `${SLURM_JOB_ID}` precise the job id

`${SLURM_ARRAY_JOB_ID}`precise the id of the job array

`${SLURM_ARRAY_TASK_ID}` precise the number of the job array task.

The script should give a answer like:

{% highlight bash %}$ sbatch array.srun
Submitted batch job 20303
$ cat slurm-20303_1.out
I am Slurm job 20305, array job 20303, and array task 1.
$ cat slurm-20303_19.out
I am Slurm job 20323, array job 20303, and array task 19.{% endhighlight %}

### Submit a R job
You can use the same syntax than before for Slurm.
You just have to launch your R script with the R script command


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
Rscript script.R #launch the R script script.R {% endhighlight %}

### Submit a job with several command in parallel at the same time

You have to use the options `--ntasks` and `--cpus-per-task` 

Example:

{% highlight bash %}#!/bin/bash

#SBATCH --ntasks=2
#SBATCH --cpu-per-task=2

srun --ntasks=1 sleep 10 & 
srun --ntasks=1 sleep 12 &
wait{% endhighlight %}

In this example, we use 2 tasks with 2 cpus allocated per task that is to say 4 cpus allocated for this job.

For each task a sleep is launched at the same time.

Notice the use of srun to launch a parallelised command and the  `&` to launch the command in background

The `wait` is needed here to ask the job to wait for the end of each command before stopping


### Submit an OpenMP job:

A OpenMP job is a job using several cpus on the same single node. Therefore the number of nodes will always be one.

This will work with a program compiled with openMP

{% highlight bash %}#!/bin/bash
#SBATCH --partition=short   ### Partition
#SBATCH --job-name=HelloOMP ### Job Name
#SBATCH --time=00:10:00     ### WallTime
#SBATCH --nodes=1           ### Number of Nodes
#SBATCH --ntasks-per-node=1 ### Number of tasks (MPI processes)
#SBATCH --cpus-per-task=28  ### Number of threads per task (OMP threads)
#SBATCH --account=hpcrcf    ### Account used for job submission
 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
 
./hello_omp{% endhighlight %}


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
                  
 
