---
layout: page
title: "R Survival Manuel"
permalink: /toolbox/RsurvivalManual/
tags: [ R, survival manual ]
description: R tutoroal
---

| Description |  A step by step guide to start on R |
| :------------- | :------------- | :------------- | :------------- |
| Authors | christine Tranchant-Dubreuil (christine.tranchant@ird.fr)  |
| Creation Date | 21/09/2018 |
| Last Modified Date | 21/09/2018 |



### Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->
- [Starting slowly on R](#start)
  - [Installing new packages `install.packages(packages)`](#install)
  - [Loading library`library`](#library)
  - [Setting a working directory](#setwd)
  - [Getting the current working directory `getwd`](#getwd)
  
- [Importing Data From `csv` File Into a `Dataframe`](#importCsv)
  - [Storing the contebnt of a file *allGenomeVersion.csv* into a dataframe *myGenome* - `read.csv2()` (french format)](#
  - [Displaying the whole dataframe](#printDF)
  - [Getting the type - `class(dataframe)`](#getType)
  - [Getting the structure of the dataframe (data type, number of levels) - `str(dataframe)`](#structureDF)
  - [Getting basic stats about te content of a dataframe - `summary(dataframe)`](#statDF)

- [License](#license) 

-----------------------
<a name="start"></a>
### Starting slowly on R

<a name="install"></a>
##### Installing New Packages/libraries - `install.packages(package) `

{% highlight bash %}
>install.packages('ggplot2')
{% endhighlight %}

<a name="library"></a>
##### Loading libraries  - `library(package)`

{% highlight bash %}
>library(ggplot2) # plot
>library(scales)  # plot scales
>library(stringr) # str_detect
>library(gridExtra) # grid`
{% endhighlight %}

<a name="setwd"></a>
#####  Setting a Working Directory `setwd(path)`

{% highlight bash %}
>setwd(/users/tranchan/riceAnalysis)
{% endhighlight %}

<a name="getwd"></a>
#####  Getting the current Working Directory `getwd()`

{% highlight bash %}
>getwd()
[1] "/Users/tranchan/Documents/Bioanalyse/panGenome"
{% endhighlight %}

<a name="importCsv"></a>
### Importing Data From `csv` File Into a `Dataframe`

<a name="readCsv"></a>
##### Storing a file *allGenomeVersion.csv* into the dataframe *myGenome* - `read.csv2()` (french format)

{% highlight bash %}
#File with 6 columns and 2010 lines
Name	Type	Length	%GC	Organism	Type
Chr01	N	32613412	43.53	dna:chromosome	v1
Chr01	N	39656875	43.05	dna:chromosome	v2
Chr02	N	29142869	43.19	dna:chromosome	v1
Chr02	N	34451706	42.81	dna:chromosome	v2
Chr03	N	33053699	43.20	dna:chromosome	v1
Chr03	N	35526804	43.23	dna:chromosome	v2
Chr04	N	26299011	43.54	dna:chromosome	v1
Chr04	N	31572834	43.19	dna:chromosome	v2
Chr05	N	23192814	42.91	dna:chromosome	v1
Chr05	N	26274045	42.82	dna:chromosome	v2
Chr06	N	24174485	42.90	dna:chromosome	v1
Chr06	N	29275208	42.90	dna:chromosome	v2
Chr07	N	21799424	43.09	dna:chromosome	v1
Chr07	N	26599614	42.82	dna:chromosome	v2
Chr08	N	20292731	42.97	dna:chromosome	v1
Chr08	N	25472747	42.64	dna:chromosome	v2
Chr09	N	17607432	43.51	dna:chromosome	v1
Chr09	N	21796211	42.91	dna:chromosome	v2
Chr10	N	16910673	43.10	dna:chromosome	v1
Chr10	N	22418184	42.96	dna:chromosome	v2
Chr11	N	20796451	42.59	dna:chromosome	v1
Chr11	N	26393634	42.28	dna:chromosome	v2
Chr12	N	19154523	42.47	dna:chromosome	v1
Chr12	N	25020143	42.38	dna:chromosome	v2
ADWL01002872.1	N	9871	37.93	dna:scaffold	v1
ADWL01002880.1	N	2202	52.91	dna:scaffold	v1
ADWL01002881.1	N	2284	63.09	dna:scaffold	v1
....

{% endhighlight %}

{% highlight bash %}
>myGenome <- read.csv2("allGenomeVersion.csv",header = TRUE)
{% endhighlight %}

<a name="printDF"></a>
### Displaying the whole dataframe

{% highlight bash %}
>myGenome

Name  Type  Length  X.GC  Organism  Type.1
Chr10	N	16910673	43.10	dna:chromosome	v1
Chr09	N	17607432	43.51	dna:chromosome	v1
Chr12	N	19154523	42.47	dna:chromosome	v1
Chr08	N	20292731	42.97	dna:chromosome	v1
Chr11	N	20796451	42.59	dna:chromosome	v1
Chr07	N	21799424	43.09	dna:chromosome	v1
....
{% endhighlight %}

<a name="getType"></a>
### Get the type - `class(dataframe)` 

{% highlight bash %}
>class(myGenome)
[1] "data.frame"
{% endhighlight %}


<a name="structureDF"></a>
### Get the structure of the dataframe (data type, number of levels) - `str(dataframe)`

{% highlight bash %}
>str(myGenome) 
'data.frame':	2010 obs. of  6 variables:
 $ Name    : Factor w/ 1998 levels "ADWL01002872.1",..: 1210 1209 1212 1208 1211 1207 1205 1206 1204 1202 ...
 $ Type    : Factor w/ 1 level "N": 1 1 1 1 1 1 1 1 1 1 ...
 $ Length  : int  16910673 17607432 19154523 20292731 20796451 21799424 23192814 24174485 26299011 29142869 ...
 $ X.GC    : Factor w/ 1191 levels "26.08","30.70",..: 512 547 454 499 466 511 493 492 550 520 ...
 $ Organism: Factor w/ 2 levels "dna:chromosome",..: 1 1 1 1 1 1 1 1 1 1 ...
 $ Type.1  : Factor w/ 2 levels "v1","v2": 1 1 1 1 1 1 1 1 1 1 ...
{% endhighlight %}

<a name="statDF"></a>
### Get basics stats about dataframe content - `summary(dataframe)`

{% highlight bash %}
>summary(myGenome)
      Name      Type         Length              X.GC                Organism    Type.1   
 Chr01  :   2   N:2010   Min.   :    1026   44.90  :   7   dna:chromosome:  24   v1:1951  
 Chr02  :   2            1st Qu.:    3229   41.27  :   6   dna:scaffold  :1986   v2:  59  
 Chr03  :   2            Median :    5742   41.73  :   6                                  
 Chr04  :   2            Mean   :  330219   42.82  :   6                                  
 Chr05  :   2            3rd Qu.:   13154   42.91  :   6                                  
 Chr06  :   2            Max.   :39656875   43.00  :   6                                  
 (Other):1998                               (Other):1973  
{% endhighlight %}

***


### License
<a name="license"></a>

<div>
The resource material is licensed under the Creative Commons Attribution 4.0 International License (<a href="http://creativecommons.org/licenses/by-nc-sa/4.0/">here</a>).
<center><img width="25%" class="img-responsive" src="http://creativecommons.org.nz/wp-content/uploads/2012/05/by-nc-sa1.png"/>
</center> 
</div>




