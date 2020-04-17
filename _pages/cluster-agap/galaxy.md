--- 
layout: page
title: "Galaxy"
permalink: /cluster-agap/galaxy/
tags: [ ]
description: Galaxy Howto page
---

| Description | HowTos for Galaxy |
| :------------- | :------------- | :------------- | :------------- |
| Authors | Marilyne Summo (marilyne.summo@cirad.fr)  |
| Creation Date | 16/04/20 |
| Last Modified Date | 16/04/20  |


-----------------------

### Summary

<!-- TOC depthFrom:2 depthTo:2 withLinks:1 updateOnSave:1 orderedList:0 -->
* [How to: Register to Galaxy South Green](#register)
* [How to: Upload a file in my history](#upload)
* [How to: Upload big datasets](#bigdata)



* [Links](#links)
* [License](#license)


-----------------------

<a name="register"></a>
### Register

You can create your account directly by clicking on the register button on our Galaxy instance. 
South Green members or partner benefits from 25Go. Please contact admin-bioinfo@cirad.fr after creating your account.
<img width="50%" class="img-responsive" src="{{ site.url }}/images/galaxy_register.png"/>

-----------------------


<a name="upload"></a>
### How to : Upload a file in my history

##### There are several solutions to import a file

* Import a file stored locally on your computer by clicking on "choose a file"
* Import a file from a URL by copying the address in the "URL / Text" frame
* Copy the text of the file directly in the "URL / Text" frame
* Import a shared file in the “shared data”

<img width="50%" class="img-responsive" src="{{ site.url }}/images/galaxy_upload.png"/>

When loading a file, Galaxy can detect type automatically but you can also choose the type of your file (txt, fasta,…).
<img width="50%" class="img-responsive" src="{{ site.url }}/images/galaxy_filetype.png"/>

To import files from shared data go to "Shared Data" => " Data Libraries "
Select the files you want to import and click the " To history " button.

Monitoring of imports in your history:

* Blue: job has been submitted
* Yellow: the job is being processed
* Green: the job ended successfully
* Red: the job is in error

You can have as many historie as you want and switch between histories. However, we recommend that you organize your data as follows :
1 history = 1 analysis
and to name the history in a recognizable way.

-----------------------


<a name="bigdata"></a>
### How to : Upload big datasets

<img width="50%" class="img-responsive" src="{{ site.url }}/images/galaxy_loadbigdata.png"/>

#### Transfer your file to the HPC, using Filezilla or any other FTP client, into your personnal "User directory":

/work/GALAXY/galaxy/users_libraries/your.name@mail.com

If you do not have any directory, ask galaxy-dev-southgreen@cirad.fr

#### Add the file into your personnal "data library"

Shared Data => Data Libraries.

Then select the library corresponding to your name.

If you do not have any library, please contact: galaxy-dev-southgreen@cirad.fr

#### Import the file into one of your histories for analysis

In order to add data into your library, click on the icon as below. 

Then you can import data by browsing your personal "User directory" (corresponding to the directory /work/GALAXY/galaxy/users_libraries/your.name@mail.com)

<img width="50%" class="img-responsive" src="{{ site.url }}/images/galaxy_data_libraries.png"/>


-----------------------
<a name="howto-11"></a>
### How to : Cite the Itrop platform in your publications

Please just copy the following sentence:

`The authors acknowledge the IRD itrop HPC (South Green Platform) at IRD montpellier for providing HPC  resources that have contributed to the research results reported within this paper.
    URL: https://bioinfo.ird.fr/- http://www.southgreen.fr` 

-----------------------

### Links
<a name="links"></a>

* Related courses : [Linux for Dummies](https://southgreenplatform.github.io/trainings/linux/)
* Related courses : [HPC](https://southgreenplatform.github.io/trainings/HPC/)
* Tutorials : [Linux Command-Line Cheat Sheet](https://southgreenplatform.github.io/trainings/linux/linuxTuto/)

-----------------------

### License
<a name="license"></a>

<div>
The resource material is licensed under the Creative Commons Attribution 4.0 International License (<a href="http://creativecommons.org/licenses/by-nc-sa/4.0/">here</a>).
<center><img width="25%" class="img-responsive" src="http://creativecommons.org.nz/wp-content/uploads/2012/05/by-nc-sa1.png"/>
</center>
</div>
                  
 
