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
| Last Modified Date | 22/09/2018 |



### Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->
- [Starting slowly on R](#start)
  - [Installing new packages `install.packages(packages)`](#install)
  - [Loading library`library`](#library)
  - [Setting a working directory](#setwd)
  - [Getting the current working directory `getwd`](#getwd)
  
- [Importing Data From `csv` File Into a `Dataframe`](#importCsv)
  - [Storing the content of a file *allGenomeVersion.csv* into a dataframe *myGenome* - `read.csv2()` (french format)](#readCsv)
  - [Displaying the whole dataframe](#printDF)
  - [Getting the type - `class(dataframe)`](#getType)
  - [Getting the structure of the dataframe (data type, number of levels) - `str(dataframe)`](#structureDF)
  - [Getting basic stats about te content of a dataframe - `summary(dataframe)`](#statDF)

- [Displaying basic informations about the dataframe structure](#info)
  - [Column names - `names(dataframe)` or  `colnames(dataframe)`](#colnames)
  - [Lines and columns number - `dim(dataframe)`](#dim)
  - [Lines number - `nrow(dataframe)`](#nrow)
  - [Columns number - `ncol(dataframe)`](#ncol)

- [Displaying the dataframe content](#print)
  - [Printing the first lines - `head(dataframe)`](#head)
  - [Printing the last lines - `last(dataframe)`](#last)
  - [Printing a column - `dataframe[colNum|colName]` or  `dataframe$colName`](#df0)
  - [Printing one whole line - `dataframe[lineNum,]`](#df)
  - [Printing one line, one column - `dataframe[lineNum,colNum]`](#df2)
  - [Printing some lines - `dataframe[LineStart:LineEnd,]`](#df3)
  - [Printing one column of some lines - `dataframe[c(line1,line2,line3),colNum]`](#df4)
  
- [Manipulating the content of a dataframe](#manipulating)
  - [Adding a new column](#add)
  - [Extracting unique values of a column - `unique(dataframe$colName)`](#unique)
  - [Extracting one part of a dataframe - `subset(dataframe)`](#extract)
  - [Calculating a sum - `sum(dataframe)` with filtering on an other column](#cal)
  - [Ordering dataframe on one column](#order)

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

-----------------------




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
##### Displaying the whole dataframe

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
##### Get the type - `class(dataframe)` 

{% highlight bash %}
>class(myGenome)
[1] "data.frame"
{% endhighlight %}


<a name="structureDF"></a>
##### Get the structure of the dataframe (data type, number of levels) - `str(dataframe)`

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
##### Get basics stats about dataframe content - `summary(dataframe)`

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

-----------------------




<a name="info"></a>
### Displaying basic informations about the dataframe structure

<a name="colnames"></a>
##### Column names - `names(dataframe)` or  `colnames(dataframe)`
{% highlight bash %}
> names(myGenome)
[1] "Name"     "Type"     "Length"   "X.GC"     "Organism" "Type.1"  

> colnames(myGenome)
[1] "Name"     "Type"     "Length"   "X.GC"     "Organism" "Type.1"  
{% endhighlight %}

<a name="dim"></a>
##### Lines and columns number - `dim(dataframe)`
{% highlight bash %}
> dim(myGenome) 
[1] 2010    6
{% endhighlight %}

<a name="nrow"></a>
##### Lines number - `nrow(dataframe)`
{% highlight bash %}
> nrow(myGenome)
[1] 2010
{% endhighlight %}

<a name="ncol"></a>
##### Columns number - `ncol(dataframe)`]  
{% highlight bash %}
> ncol(myGenome)
[1] 6
{% endhighlight %}  



-----------------------

<a name="print"></a>
### Displaying and manipulating the dataframe content

<a name="head"></a>
##### Printing the first lines - `head(dataframe)`
{% highlight bash %}
> head(myGenome)
   Name Type   Length  X.GC       Organism Type.1
1 Chr10    N 16910673 43.10 dna:chromosome     v1
2 Chr09    N 17607432 43.51 dna:chromosome     v1
3 Chr12    N 19154523 42.47 dna:chromosome     v1
4 Chr08    N 20292731 42.97 dna:chromosome     v1
5 Chr11    N 20796451 42.59 dna:chromosome     v1
6 Chr07    N 21799424 43.09 dna:chromosome     v1

> head(myGenome,n=2)
   Name Type   Length  X.GC       Organism Type.1
1 Chr10    N 16910673 43.10 dna:chromosome     v1
2 Chr09    N 17607432 43.51 dna:chromosome     v1
{% endhighlight %}

<a name="last"></a>
##### Printing the last lines - `last(dataframe)`
{% highlight bash %}
> tail(myGenome)
            Name Type Length  X.GC     Organism Type.1
2005 ChrUN-Ctg78    N  21053 43.05 dna:scaffold     v2
2006 ChrUN-Ctg79    N  19363 40.63 dna:scaffold     v2
2007 ChrUN-Ctg80    N  15054 41.18 dna:scaffold     v2
2008 ChrUN-Ctg81    N  13998 42.98 dna:scaffold     v2
2009 ChrUN-Ctg82    N  13471 42.29 dna:scaffold     v2
2010 ChrUN-Ctg83    N  10534 42.37 dna:scaffold     v2

> tail(myGenome,n=1)
            Name Type Length  X.GC     Organism Type.1
2010 ChrUN-Ctg83    N  10534 42.37 dna:scaffold     v2
{% endhighlight %}


<a name="df0"></a>
##### Printing a column - `dataframe[colNum|colName]` or  `dataframe$colName`

* output result (line)

{% highlight bash %}
> head(myGenome[,"Length"])
[1] 16910673 17607432 19154523 20292731 20796451 21799424

{% endhighlight %}

{% highlight bash %}
> head(myGenome$Length) 
[1] 16910673 17607432 19154523 20292731 20796451 21799424

{% endhighlight %}

* output result (column)

{% highlight bash %}
> head(myGenome[1])
   Name
1 Chr10
2 Chr09
3 Chr12
4 Chr08
5 Chr11
6 Chr07
{% endhighlight %}

{% highlight bash %}
> head(myGenome['Length'])
    Length
1 16910673
2 17607432
3 19154523
4 20292731
5 20796451
6 21799424
{% endhighlight %}


<a name="df"></a>
##### Printing one whole line - `dataframe[lineNum,] `
Line 2
{% highlight bash %}
> myGenome[2,]
   Name Type   Length  X.GC       Organism Type.1
2 Chr09    N 17607432 43.51 dna:chromosome     v1
{% endhighlight %}

<a name="df2"></a>
##### Printing one line, one column - `dataframe[lineNum,colNum]`
Line 2, column 6 then line 2, column 3
{% highlight bash %}
> myGenome[2,6]
[1] v1
Levels: v1 v2

> myGenome[2,3]
[1] 17607432
{% endhighlight %}

<a name="df3"></a>
##### Printing some lines - `dataframe[LineStart:LineEnd,]`
{% highlight bash %}
> myGenome[1:3,]
   Name Type   Length  X.GC       Organism Type.1
1 Chr10    N 16910673 43.10 dna:chromosome     v1
2 Chr09    N 17607432 43.51 dna:chromosome     v1
3 Chr12    N 19154523 42.47 dna:chromosome     v1

> myGenome[10:15,]
    Name Type   Length  X.GC       Organism Type.1
10 Chr02    N 29142869 43.19 dna:chromosome     v1
11 Chr01    N 32613412 43.53 dna:chromosome     v1
12 Chr03    N 33053699 43.20 dna:chromosome     v1
13 Chr01    N 39656875 43.05 dna:chromosome     v2
14 Chr02    N 34451706 42.81 dna:chromosome     v2
15 Chr03    N 35526804 43.23 dna:chromosome     v2
{% endhighlight %}

<a name="df4"></a>
##### Printing one column of some lines - `dataframe[c(line1,line2,line3),colNum]`
Print the complete lines 1, 3, 7, 6 then just the column 3 of the same lines

{% highlight bash %}
> myGenome[c(1,3,7,6),]
   Name Type   Length  X.GC       Organism Type.1
1 Chr10    N 16910673 43.10 dna:chromosome     v1
3 Chr12    N 19154523 42.47 dna:chromosome     v1
7 Chr05    N 23192814 42.91 dna:chromosome     v1
6 Chr07    N 21799424 43.09 dna:chromosome     v1

> myGenome[c(1,3,7,6),3]
[1] 16910673 19154523 23192814 21799424

{% endhighlight %}

-----------------------



<a name="manipulating"></a>
### Manipulating the content of a dataframe


<a name="add"></a>
##### Adding a new column 

* A new collumn firstly filled with "NA"
{% highlight bash %}
> myGenome["mb"] <- NA 

> head(myGenome)
   Name Type   Length  X.GC       Organism Type.1 mb
1 Chr10    N 16910673 43.10 dna:chromosome     v1 NA
2 Chr09    N 17607432 43.51 dna:chromosome     v1 NA
3 Chr12    N 19154523 42.47 dna:chromosome     v1 NA
4 Chr08    N 20292731 42.97 dna:chromosome     v1 NA
5 Chr11    N 20796451 42.59 dna:chromosome     v1 NA
6 Chr07    N 21799424 43.09 dna:chromosome     v1 NA
{% endhighlight %}

* The new column receives the result of an operation
{% highlight bash %}
> myGenome$mb <- as.integer(myGenome$Length/1000000)

> head(myGenome)
   Name Type   Length  X.GC       Organism Type.1 mb
1 Chr10    N 16910673 43.10 dna:chromosome     v1 16
2 Chr09    N 17607432 43.51 dna:chromosome     v1 17
3 Chr12    N 19154523 42.47 dna:chromosome     v1 19
4 Chr08    N 20292731 42.97 dna:chromosome     v1 20
5 Chr11    N 20796451 42.59 dna:chromosome     v1 20
6 Chr07    N 21799424 43.09 dna:chromosome     v1 21
{% endhighlight %}

<a name="unique"></a>
##### Extracting unique values of a column - `unique(dataframe$colName)`
{% highlight bash %}
> myGenome$Organism
   [1] dna:chromosome dna:chromosome dna:chromosome dna:chromosome dna:chromosome dna:chromosome dna:chromosome
   [8] dna:chromosome dna:chromosome dna:chromosome dna:chromosome dna:chromosome dna:chromosome dna:chromosome
  [15] dna:chromosome dna:chromosome dna:chromosome dna:chromosome dna:chromosome dna:chromosome dna:chromosome
  [22] dna:chromosome dna:chromosome dna:chromosome dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
  [29] dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
  [36] dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
  [43] dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
  [50] dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
  [57] dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
  [64] dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
 ... 
 [981] dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
 [988] dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
 [995] dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold   dna:scaffold  
 [ reached getOption("max.print") -- omitted 1010 entries ]
Levels: dna:chromosome dna:scaffold

> unique(myRef$Organism)
[1] dna:chromosome dna:scaffold  
Levels: dna:chromosome dna:scaffold

> unique(myGenomef$Type.1)
[1] v1 v2
Levels: v1 v2
{% endhighlight %}

<a name="extract"></a>
##### Extracting one part of a dataframe into a new dataframe - `subset(dataframe)`
{% highlight bash %}
> myGenomeSubset <- subset(myGenome, Organism=="dna:chromosome")
> head(myrefGenome)

> myrefSubset <- subset(myGenome, Organism=="dna:chromosome")
> head(myrefGenome)
   Name Type   Length  X.GC       Organism Type.1 mb
1 Chr10    N 16910673 43.10 dna:chromosome     v1 16
2 Chr09    N 17607432 43.51 dna:chromosome     v1 17
3 Chr12    N 19154523 42.47 dna:chromosome     v1 19
4 Chr08    N 20292731 42.97 dna:chromosome     v1 20
5 Chr11    N 20796451 42.59 dna:chromosome     v1 20
6 Chr07    N 21799424 43.09 dna:chromosome     v1 21

> dim(myRGenome)
[1] 2010    7

> dim(myGenomeSubset)
[1] 24  7

{% endhighlight %}

<a name="cal"></a>
##### Calculating a sum - `sum(dataframe)` with filtering on an other column
{% highlight bash %}
> sum(myGenomeSubset$Length)
[1] 629495529

> sum(myGenomeSubset$Length[myGenomeSubset$Type.1=="v1"])
[1] 285037524

> sum(myGenomeSubset$Length[myGenomeSubset$Type.1=="v2"])
[1] 344458005

{% endhighlight %}

<a name="order"></a>
##### Ordering dataframe on one column
{% highlight bash %}
myGenomeOrdered <- myGenomeSubset[order(myGenomeSubset$Name),]

> head(myGenomeOrdered)
    Name Type   Length  X.GC       Organism Type.1 mb
11 Chr01    N 32613412 43.53 dna:chromosome     v1 32
13 Chr01    N 39656875 43.05 dna:chromosome     v2 39
10 Chr02    N 29142869 43.19 dna:chromosome     v1 29
14 Chr02    N 34451706 42.81 dna:chromosome     v2 34
12 Chr03    N 33053699 43.20 dna:chromosome     v1 33
15 Chr03    N 35526804 43.23 dna:chromosome     v2 35
{% endhighlight %}

-----------------------
<a name="barplot"></a>
### Creating a barplot graphic

size (Mb) per chromosome (for each genome version)

http://www.sthda.com/french/wiki/ggplot2-barplots-guide-de-demarrage-rapide-logiciel-r-et-visualisation-de-donnees#barplots-basiques

<a name="head"></a>
##### with `gggplot()` - `head(dataframe)`

* basic barplot with `stat = "identity"` and `by type - `fill` 

{% highlight bash %}
#basic barplot
p <- ggplot(data = myrefSubset, aes(x = Name, y=mb, fill=Type.1)) + 
      geom_bar(stat = "identity")
p 
{% endhighlight %} 

<img class="img-responsive" width="50%" src="{{ site.url }}/images/R.barplot-type.png" alt="barplot" />

* distinct plot - `position=position_dodge()`

{% highlight bash %}
#basic barplot
q <- ggplot(data = myrefSubset, aes(x = Name, y=mb, fill=Type.1)) + 
      geom_bar(stat = "identity", position=position_dodge())
#horizontal
q + coord_flip()
{% endhighlight %} 

* save into a file

{% highlight bash %}
jpeg(file="pseudomolOMAP.jpg");
#basic barplot
q <- ggplot(data = myrefSubset, aes(x = Name, y=mb, fill=Type.1)) + 
      geom_bar(stat = "identity", position=position_dodge())
#horizontal
q + coord_flip()
dev.off;
{% endhighlight %} 

-----------------------


<a name="print"></a>
### Displaying and manipulating the dataframe content

<a name="head"></a>
##### Printing the first lines - `head(dataframe)`
{% highlight bash %}

{% endhighlight %}

-----------------------

### License
<a name="license"></a>

<div>
The resource material is licensed under the Creative Commons Attribution 4.0 International License (<a href="http://creativecommons.org/licenses/by-nc-sa/4.0/">here</a>).
<center><img width="25%" class="img-responsive" src="http://creativecommons.org.nz/wp-content/uploads/2012/05/by-nc-sa1.png"/>
</center> 
</div>




