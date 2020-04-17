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
<img width="100%" class="img-responsive" src="{{ site.url }}/images/galaxy_filetype.png"/>


##### Open FileZilla and save the IRD cluster into the site manager

<img width="50%" class="img-responsive" src="{{ site.url }}/images/tpLinux/tp-filezilla1.png"/>

In the FileZilla menu, go to _File > Site Manager_. Then go through these 5 steps:

1. Click _New Site_.
2. Add a custom name for this site.
3. You have 3 possible choices:

  - bioinfo-nas2.ird.fr (nas2) to transfer  to /data/project
  - bioinfo-nas.ird.fr (nas) to transfer to /home/user, /data2/projects or /teams
  - bioinfo-nas3.ird.fr  (nas3)to transfer to /data3/project
  
<img width="50%" class="img-responsive" src="{{ site.url }}/images/transfert_cluster.png"/>

      
4. Set the Logon Type to "Normal" and insert your username and password used to connect on the IRD cluster
5. Choose port 22 and press the "Connect" button.


##### Transferring files

<img width="50%" class="img-responsive" src="{{ site.url }}/images/tpLinux/tp-filezilla2.png"/>

1. From your computer to the cluster : click and drag an text file item from the left local colum to the right remote column 
2. From the cluster to your computer : click and drag an text file item from he right remote column to the left local column




-----------------------


<a name="howto-2"></a>
### How to : Connect to the cluster via `ssh`

#### First connection:

Your password has to be changed at your first connection.

"Mot de passe UNIX (actuel)":   you are asked to type the password provided in the account creation email.

Then type your new password twice.

The session will be automatically closed.

You will need to open a new session with your new password.

#### From a windows computer:
 
In mobaXterm:
1. Click the session button, then click SSH.
* In the remote host text box, type: bioinfo-master.ird.fr
* Check the specify username box and enter your user name
2. In the console, enter the password when prompted.
Once you are successfully logged in, you will be use this console for the rest of the lecture. 

#### From a mac or a linux computer:

Open the terminal application and type the following command:

`ssh login@bioinfo-master.ird.fr`

with login: your cluster account


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
                  
 
