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
  
- [Import Data From `csv` File Into a `Dataframe`](#importCsv)
- [License](#license) 

-----------------------
<a name="start"></a>
### Starting slowly on R

<a name="install"></a>
##### Installing New Packages/libraries - `install.packages(package) `

{% highlight bash %}
install.packages('ggplot2')
{% endhighlight %}

<a name="library"></a>
##### Loading libraries  - `library(package)`

{% highlight bash %}
library(ggplot2) # plot
library(scales)  # plot scales
library(stringr) # str_detect
library(gridExtra) # grid`
{% endhighlight %}

<a name="setwd"></a>
#####  Setting a Working Directory `setwd(path)`

{% highlight bash %}
setwd(/users/tranchan/riceAnalysis)
{% endhighlight %}

<a name="getwd"></a>
#####  Getting the current Working Directory `getwd()`

{% highlight bash %}
getwd()
[1] "/Users/tranchan/Documents/Bioanalyse/panGenome"
{% endhighlight %}

<a name="importCsv"></a>
### Import Data From `csv` File Into a `Dataframe`

<a name="importCsv"></a>

### License
<a name="license"></a>

<div>
The resource material is licensed under the Creative Commons Attribution 4.0 International License (<a href="http://creativecommons.org/licenses/by-nc-sa/4.0/">here</a>).
<center><img width="25%" class="img-responsive" src="http://creativecommons.org.nz/wp-content/uploads/2012/05/by-nc-sa1.png"/>
</center> 
</div>




