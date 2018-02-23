---
layout: page
title: "Docker for dumpies"
permalink: /toolbox/docker/
tags: [ git ]
description: Docker for dumpies
---
<table class="table-contact">
<tr>
<td><img width="100%" src="{{ site.url }}/images/trainings-linux.png" alt="" />
</td>
<td>
<h1> Docker For Dummies</h1><br />
This page describes describes the main commands you need in order to use Docker.
</td>
</tr>
</table>


### Author(s)

| Authors  | Christine Tranchant-Dubreuil  |
| :------------- | :------------- |
| Research Unit | UMR DIADE   |
| Institut |  <img src="https://www.ird.fr/extension/ird/design/ird/images/picto/logo_ird.png" width="20%"> |


<a name="module1-intro"></a>
##  Module 1 : Introduction Docker

### build
* Workflows Developer : production et compilation image (DEV)
* outil : Docker toolbox
> pour Mac et windows (pas demon natif comme sous linux), permet de faire tourner le demon docker qui va tourner dans une machine virtuelle linux (virtual box) <-> engine docker <=> docker for Mac and Windows now : add url)
* Se compose de :
> Docker Engine
> Docker Compose
> Docker Swarm (partie cluster)
> Docker Client
> Kitematic : outil graphique  (pas intéressant)


### ship
* registry services: transport et stockage des images produites pendant la phase de dev
* outil :
> ** Docker Hub **
 <=> "github" docker = plateforme mutualisée en ligne de reference https://hub.docker.com/
 stocke image sous cloud ou sur son propre réseau/registry

> ** Docker Trusted Registry **
avec interface (commerciale)
--> interface web
--> gestion utilisateur avec accès

> ** Docker Registry open source **
version gratuite sans interface


### run
en envt dev ou prod : gestion des dockers
* ** docker cloud ** (plateforme en ligne payante ) qui permet de deployer les conteners crées
* ** docker Universal Control Plane ** / Docker data center
> deploiement sur cloud, vm ou hybrides

Rq : toutes ces couches s appuient sur le Docker Engine


===> mettre image

***

<a name="module1-intro"></a>
##  Module 2 :

### Docker Engine

* Program that enables containers to be distribued and run
* Several natives fonctionalities
* keys :
> NAMESPACES : isolation ps et syst fichiers
> CGROUPS : mesurer et limiter utilisation ressources + donner accès à des periphériques
> Regles IPTABLES : communication entre containers sur meme hote + comm entre containers et exterieur

* Archi client / server
* Docker Containers and Images
* Registry and repository
image

* Docker orchestration

3 outils
- Docker Machine : provisionne hotes dockers et installer Docker engine dessus
- Docker Compose : pour creer et gerer applications multi containers
- Docker Engine swarm mode : pour creer des clusters dans la version Engine facilement en 3 commandes

<a name="module3-install"></a>
##  Module 3 : Installation Docker

### Docker Engine
* on linux : url / cf dia
* on mac (> yosemite) / window (>)
* community edition / entreprise edition

<a name="module4-image"></a>
##  Module 4 : Introduction image

<a name="module5-run"></a>
### Module 5 : Running and Managing containers

### Module 6 : Building images
image

3 ways
* commit changes from a container as a new image
* build from a DockerFile
* import a tarball into Docker as a standalone base layer

#### Committing changes in a container
pas automatisable!

#### build from a DockerFile
* programs to install
* base images to use **FROM**
* commands to run ** RUN cmde **

### Module 8 : volume

### Advanced docker topics

