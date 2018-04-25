---
layout: page
title: "Linux Tutorial"
permalink: /toolbox/linux/
tags: [ linux, survival guide ]
description: Linux page
---

| Description |  Linux tutorial |
| :------------- | :------------- | :------------- | :------------- |
| Authors | christine Tranchant-Dubreuil (christine.tranchant@ird.fr)  |
| Creation Date | 26/02/2018 |
| Last Modified Date | 25/03/2018 |


-----------------------

### Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->
- [Moving around the Filesystem  and manipulating files/folders](#filesystem)
  - [Printing the full path of the current directory `pwd`](#pwd)
  - [Listing files in a directory  `ls`](#ls)
  - [Moving in the file tree  `cd`](#cd)
  - [Making directories `mkdir directory_name`](#mkdir)
  - [Copying files `cp`](#cp)
  - [Moving files `mv`](#mv)
  - [Removing files and directories using `rm`and `rmdir`](#rm)
- [Displaying the contents of a file on the screen](#display)
  - [clear (clear screen) `clear`](#clear)
  - [Display the content of file using `cat`](#cat)
  - [Display the content of file using `less`](#less)
  - [Display the begin of a file using `head`](#head)
  - [Display the end of a file using `tail`](#tail)
- [Searching the contents of a file](#searching)
  - [Searching word in a file using `grep`](#grep)
  - [Count for word, line count in a file using `wc`](#wc)
- [TIPS](#tips):
  - [Creating and extracting a tar gz archive using `tar`](#tar)
  - [Compressing and extracting files using `gzip`](#gzip)
  - [Knowing how much space a file or directory is using on a disk with `du`](#du)
  - [Creating a file shortcut/ a symbolic link with `ln -s`](#ln)
  - [Downloading a file over HTTP with `wget`](#wget)
  
- [License](#license) 

-----------------------
<a name="filesystem"></a>
### Moving around the Filesystem  and manipulating files/folders

<a name="pwd"></a>
##### Printing the name/full path of the current directory `pwd`

{% highlight bash %}
[tranchant@master0 ~]$ pwd
/home/tranchant
{% endhighlight %}

<a name="ls"></a>
##### Listing files in a directory `ls`
`ls`: list all files in the current directory

{% highlight bash %}
[tranchant@master0 ~]$ ls 
AIRAIN
All-EST-coffea.fasta
all-gene.gff3.100000.gene-density
all_no_TE.gff3
all_no_TE.gff3.100000.gene-density
all_no_TE.gff3.10000.gene-density
All-SNP.vcf.density100000.snpden
blast_batch.pl
blast.sh
circos
DEDUP2-1
DEDUP2-all
Desktop
FASTA-TRINITY
fastq-stat.txt
GBS
{% endhighlight %}

###### Listing files in a directory gived as argument `ls directory_name`

{% highlight bash %}
[tranchant@master0 ~]$ ls /home 
abate
abdelrah
adam
adeoti
adereeper
admin
agrondi
aichatou
{% endhighlight %}


##### `ls -l` : Display the long format listing of all files in the current directory 

{% highlight bash %}
[tranchant@master0 ~]$ ls -l
total 148272
drwxr-xr-x 10 tranchant ggr      4096 13 mars   2017 AIRAIN
-rwxr-xr-x  1 tranchant ggr  51128305 11 sept. 14:16 All-EST-coffea.fasta
-rw-r--r--  1 tranchant ggr     95117 24 févr.  2017 all-gene.gff3.100000.gene-density
-rwxr-xr-x  1 tranchant ggr  64221458 24 févr.  2017 all_no_TE.gff3
-rw-r--r--  1 tranchant ggr     93796 24 févr.  2017 all_no_TE.gff3.100000.gene-density
-rw-r--r--  1 tranchant ggr    889389 24 févr.  2017 all_no_TE.gff3.10000.gene-density
-rwxr-xr-x  1 tranchant root    90498  5 déc.  07:21 All-SNP.vcf.density100000.snpden
-rwxr-xr-x  1 tranchant ggr     10366 29 mars   2017 blast_batch.pl
-rwxr-xr-x  1 tranchant ggr       878 29 mars   2017 blast.sh
drwxr-xr-x  3 tranchant ggr      4096 19 janv. 15:02 circos
drwxr-sr-x  2 tranchant ggr      4096  8 juin   2017 DEDUP2-1
drwxr-sr-x  2 tranchant ggr      4096  8 juin   2017 DEDUP2-all
drwxr-xr-x  4 root      root     4096  9 juil.  2015 Desktop
drwxr-sr-x  2 tranchant ggr      8192  8 juin   2017 FASTA-TRINITY
-rw-r--r--  1 tranchant ggr      1056 30 nov.  15:31 fastq-stat.txt
drwxr-xr-x  2 tranchant ggr        48 10 août   2016 GBS
{% endhighlight %}

##### `ls -a` : Display all the files and directories (even hidden files)  

{% highlight bash %}
[tranchant@master0 ~]$ ls -a
.                                                localConfig.pm
..                                               .localConfig.pm.swp
AIRAIN                                           localConfig-USR.pm
All-EST-coffea.fasta                             LOG
all-gene.gff3.100000.gene-density                .ls_couleur
all_no_TE.gff3                                   .Mathematica
all_no_TE.gff3.100000.gene-density               .matplotlib
all_no_TE.gff3.10000.gene-density                .Megan.def
All-SNP.vcf.density100000.snpden                 Microsatellite_markers
.apollo                                          moduleLoad
.bash_history                                    .mozilla
.bash_logout                                     .oracle_jre_usage
.bash_profile                                    out.sites.pi
.bashrc                                          out.windowed.pi
blast_batch.pl                                   parse-go.pl
blast.sh                                         password.txt
.cache                                           path-new.txt
circos                                           path.sh
.config                                          PATH.sh
.cpan                                            path.txt
Data-TP                                          PERL
.dbus                                            perl5
{% endhighlight %}

-----------------------

<a name="cd"></a>
##### Moving in the file tree  `cd`

`cd DIRECTORY_NAME` :  change the current working directory to 'directory'
_change directory_

{% highlight bash %}
[tranchant@master0 ~]$ cd FASTQ
{% endhighlight %}

-----------------------

<a name="mkdir"></a>
##### Making directories `mkdir directory_name`

`mkdir directory_name`: make a directory in your current working directory type
_make directory_

{% highlight bash %}
[tranchant@master0 ~]$ mkdir results
{% endhighlight %}


-----------------------

<a name="cp"></a>
##### Copying files `cp`
 _copy_
  
* `cp file1 file2`: 	makes a copy of file1 and calls it file2

* `cp file_name directory_name`: copy the file  _file_name_  to the directory _directory_name', keeping the same name.

-----------------------

<a name="mv"></a>
##### Moving files `mv`
 _move_

* `mv file1 file2` :  moves (or renames) file1 to file2. It is used to rename a file, by moving the file to the same directory, but giving it a different name.
* `mv file_na
me directory` : To move the file _file_name_ from one directory to another (here _directory_name_). This has the effect of moving rather than copying the file, so you end up with only one file rather than two.

>Be careful : use preferentially cp command rather than mv command to move big files
-----------------------

<a name="rm"></a>
#####  Removing files and directories using `rm`and `rmdir`
_remove_

`rm file_name`: delete (remove) a file (here _file_name_)
To delete (remove) a file, use the rm command. As an example, we are going to create a copy of the science.txt file then delete it.

`rmdir directory_name` : remove a directory (make sure it is empty first because linux will not let remove a non-empty directory).


-----------------------
<a name="display"></a>
### Displaying the contents of a file on the screen

<a name="clear"></a>
##### clear (clear screen) `clear`

`clear` : clears the terminal. Before displaying files, it's possible to clear the terminal window with this command.

-----------------------

<a name="cat"></a>
##### Display the content of file using `cat`
`cat file1` : displays the contents of a file on the screen
_concatenate_

`cat file1 file2`

`cat *.fasta`

-----------------------

<a name="less"></a>
##### Display the content of file using `less`
`less file_name` : writes the contents of a file onto the screen a page at a time. Less is used in preference to cat for long files.

* Press the [space-bar] if you want to see another page
* Type [q] if you want to quit reading. 
* still in less, type a forward slash [/] followed by the word to search

 -----------------------

<a name="head"></a>
##### Display the begin of a file using `head`
`head file_name` : writes the first ten lines of a file to the screen.

{% highlight bash %}
[tranchant@master0 ~]$ head 
head -n 5
{% endhighlight %}


 -----------------------

<a name="tail"></a>
##### Display the end of a file using `head`
`tail file_name` : writes the last ten lines of a file to the screen.

{% highlight bash %}
[tranchant@master0 ~]$ tail
tail -n 5
{% endhighlight %}

-----------------------
<a name="searching"></a>
### Searching the contents of a file

<a name="grep"></a>
##### Searching word in a file using `grep`
`grep motif file_name` : searches files for specified words or patterns. First clear the screen, then type

> To search for a phrase or pattern, you must enclose it in single quotes (the apostrophe symbol). For example to search for spinning top, type

{% highlight bash %} 
# printed out each line containg the word science.

# The grep command is case sensitive; it distinguishes between Science and science.

#To ignore upper/lower case distinctions, use the -i option, i.e. type
[tranchant@master0 ~]$grep -i science science.txt

#
[tranchant@master0 ~]$grep -i 'spinning top' science.txt

{% endhighlight %}

> Some of the other options of grep are:

> -v display those lines that do NOT match 
> -n precede each matching line with the line number 
> -c print only the total count of matched lines 

> Try some of them and see the different results. Don't forget, you can use more than one option at a time. For example, the number of lines without the words science or Science is grep -ivc science science.txt


-----------------------

<a name="wc"></a>
##### Count for word, line count in a file using `wc`
`wc`:  short for word count
_(word count)_

-----------------------

<a name="tar"></a>
##### Creating and extracting a tar gz archive using `tar`

* To create a tar.gz archive from a given folder 

{% highlight bash %} 
# compress the contents of source-folder-name to a tar.gz archive named tar-archive-name.tar.gz
[tranchant@master0 ~]$tar -zcvf tar-archive-name.tar.gz source-folder-name 
{% endhighlight %}


* To extract a tar.gz compressed archive 
{% highlight bash %} 
# extract the archive to the folder tar-archive-name. 
tar -zxvf tar-archive-name.tar.gz 
{% endhighlight %}

-----------------------
<a name="tips"></a>

<a name="gzip"></a>
##### Compressing and extracting files using `gzip`
* `gzip {filename}`
* `gzip -d {.gz file}`

-----------------------

<a name="du"></a>
##### Knowing how much space a file or directory is using on a disk with `du`

`du [options] [file or dir]`

If you use it with no arguments you will the usage of all files and directories (recursively) of the working directory.

> -h : Shows the in human readable format
> -s : Sumarize, so it displays the total of each subdirectory and not for its contents

-----------------------

<a name="ln"></a>
##### Creating a file shortcut/ a symbolic link with `ln -s`

A symbolic link, also termed a soft link, is a special kind of file that points to another file, much like a shortcut in Windows or a Macintosh alias.

`ln -s source_file myfile`

* _source_file_ corresponds to the name of the existing file for which the symbolic link  is created (this file can be any existing file or directory across the file systems). 
*  _myfile_ is the name of the symbolic link. 

Note :
* If the source file is deleted or moved it to a different location, the symbolic file will not function properly and should be either deleted or moved it.

-----------------------

<a name="wget"></a>
##### Downloading a file over HTTP with `wget`

`wget http://website.com/files/file.zip` : downloads the file `http://website.com/files/file.zip` into the working directory.

### License
<a name="license"></a>

<div>
The resource material is licensed under the Creative Commons Attribution 4.0 International License (<a href="http://creativecommons.org/licenses/by-nc-sa/4.0/">here</a>).
<center><img width="25%" class="img-responsive" src="http://creativecommons.org.nz/wp-content/uploads/2012/05/by-nc-sa1.png"/>
</center> 
</div>
