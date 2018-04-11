---
layout: page
title: "Git for dummies"
permalink: /toolbox/git/
tags: [ git ]
description: Git for dummies
---
<table class="table-contact">
<tr>
<td><img width="100%" src="{{ site.url }}/images/trainings-linux.png" alt="" />
</td>
<td>
<h1> Git For Dummies</h1><br />
This page describes describes the main commands you need in order to use Git.
</td>
</tr>
</table>


### Author(s)

| Authors  | Christine Tranchant-Dubreuil  |
| :------------- | :------------- |
| Research Unit | UMR DIADE   |
| Institut |  <img src="https://www.ird.fr/extension/ird/design/ird/images/picto/logo_ird.png" width="20%"> |


### Keywords
git

### Date
10/03/2017

***

## Summary

- [Download the repository using the `git clone`command](#clone)
- [Update the downloaded repository using the `git pull` command](#pull)
- [Add a file, commit and pull with `git add`, `git commit` and `git pull`](#file-add)
- [Remove a file using `git rm`](#file-add)
- [Branching](#branch)

<a name="clone"></a>
## Download the repository using the `git clone`command

{% highlight ruby %}
 git clone https://github.com/SouthGreenPlatform/TOGGLE-DEV.git <directory name>
{% endhighlight %}

{% highlight ruby %}
#Example
[tranchant@master0 TOGGLE-ON-THE-FLY]$ git clone https://github.com/SouthGreenPlatform/TOGGLE-DEV.git .
Cloning into '.'...
remote: Counting objects: 7945, done.
remote: Compressing objects: 100% (124/124), done.
remote: Total 7945 (delta 78), reused 0 (delta 0), pack-reused 7820
Receiving objects: 100% (7945/7945), 170.06 MiB | 23.03 MiB/s, done.
Resolving deltas: 100% (5503/5503), done.
Checking out files: 100% (364/364), done.
{% endhighlight %}

<a name="pull"></a>
## Update the downloaded repository using the `git pull` command

update your copy of repository with the version on remote server

{% highlight ruby %}
 git pull https://github.com/SouthGreenPlatform/TOGGLE-DEV.git <branch name>
{% endhighlight %}`

{% highlight ruby %}
#Example
[tranchant@master0 TOGGLE-ON-THE-FLY]$ git pull https://github.com/SouthGreenPlatform/TOGGLE-DEV.git  master
From https://github.com/SouthGreenPlatform/TOGGLE-DEV
 * branch            master     -> FETCH_HEAD
Already up-to-date.
{% endhighlight %}


<a name="file-add"></a>
## Add a file, commit and pull with `git add`, `git commit` and `git pull`

Don't forget to pull to download the latest changes before pushing

#### To add a file (a change) to your local index with `git add`

{% highlight ruby %}
git add <filename>
{% endhighlight %}

#### To actually commit these changes with `git commit`

{% highlight ruby %}
git commit -m "message" <file name>
{% endhighlight %}

#### To send those changes to your remote repository with `git pull`
{% highlight ruby %}
git push https://github.com/SouthGreenPlatform/TOGGLE-DEV.git <branch_name>
{% endhighlight %}

Example

{% highlight ruby %}
[tranchant@master0 TOGGLE-ON-THE-FLY]$ git add update.txt

[tranchant@master0 TOGGLE-ON-THE-FLY]$ git commit -m "Adding update.txt file" update.txt
[master ebb0a1c] Adding update.txt file
 1 file changed, 0 insertions(+), 0 deletions(-)
 create mode 100644 update.txt

[tranchant@master0 TOGGLE-ON-THE-FLY]$ git push https://github.com/SouthGreenPlatform/TOGGLE-DEV.git master
Counting objects: 3, done.
Delta compression using up to 16 threads.
Compressing objects: 100% (2/2), done.
Writing objects: 100% (2/2), 271 bytes | 0 bytes/s, done.
Total 2 (delta 1), reused 0 (delta 0)
remote: Resolving deltas: 100% (1/1), completed with 1 local objects.
To https://github.com/SouthGreenPlatform/TOGGLE-DEV.git
   fec3a1f..ebb0a1c  master -> master
{% endhighlight %}

<a name="file-rm"></a>
## Remove a file using `git rm`

{% highlight ruby %}
git rm <file name>
{% endhighlight %}

{% highlight ruby %}
[tranchant@master0 TOGGLE-ON-THE-FLY]$ git rm update.txt
rm 'update.txt'

[tranchant@master0 TOGGLE-ON-THE-FLY]$ git commit -m "Remove update.txt file" update.txt
[master 9fa50b4] Remove update.txt file
 1 file changed, 0 insertions(+), 0 deletions(-)
 delete mode 100644 update.txt

[tranchant@master0 TOGGLE-ON-THE-FLY]$ git push https://github.com/SouthGreenPlatform/TOGGLE-DEV.git master
Counting objects: 3, done.
Delta compression using up to 16 threads.
Compressing objects: 100% (2/2), done.
Writing objects: 100% (2/2), 236 bytes | 0 bytes/s, done.
Total 2 (delta 1), reused 0 (delta 0)
remote: Resolving deltas: 100% (1/1), completed with 1 local objects.
To https://github.com/SouthGreenPlatform/TOGGLE-DEV.git
   ebb0a1c..9fa50b4  master -> master
{% endhighlight %}

<a name="branch"></a>
## Branching

Branches are used to develop new features or modify codes isolated from each other. The _master_ branch is the "default" branch when a repository is created. Use other branches for development and merge them back to the master branch.


#### View all branches that were ever checked out on your local copy using ` git branch `

{% highlight ruby %}
git branch
{% endhighlight %}

{% highlight ruby %}
[tranchant@master0 TOGGLE-ON-THE-FLY]$ git branch
* master
{% endhighlight %}
* indicates the branch used actually

#### View all distant branches using ` git branch `

{% highlight ruby %}
git branch -r
{% endhighlight %}

{% highlight ruby %}
[tranchant@master0 TOGGLE-ON-THE-FLY]$ git branch -r
  origin/HEAD -> origin/master
  origin/master
  origin/picardtools-samtofastq
  origin/samtoolsBlocks
  origin/structuralVariant
  origin/tgicl
  origin/transabyss
  origin/trinity
{% endhighlight %}

### Create your own branch on your local copy then transfer it on remote server

Create the branch
{% highlight ruby %}
git branch <branch name>
{% endhighlight %}

Move into this branch
{% highlight ruby %}
git checkout <branch name>
{% endhighlight %}

Commit the changes

{% highlight ruby %}
git commit -m "mon commentaire"
{% endhighlight %}

Push this local branch on the remote server

{% highlight ruby %}
 git push https://github.com/SouthGreenPlatform/TOGGLE-DEV.git nom_branche
{% endhighlight %}


### Get a distant branch on the local repository if the branch don't exist locally

#### Method 1

{% highlight ruby %}
git checkout <remote branch name>
git pull https://github.com/SouthGreenPlatform/TOGGLE-DEV.git nom_branche_distante `
{% endhighlight %}

#### Method 2

{% highlight ruby %}
git branch <remote branch name>
git pull https://github.com/SouthGreenPlatform/TOGGLE-DEV.git <remote branch name>
git checkout <remote branch name>
{% endhighlight %}


###  To merge another branch  (ex: samtoolsBlock)  into your active branch (e.g. master)


#### Move into the "active" branch (e.g. master)

{% highlight ruby %}
git checkout master
{% endhighlight %}

#### Update your local repository to the newest commit,

{% highlight ruby %}
git pull https://github.com/SouthGreenPlatform/TOGGLE-DEV.git master
{% endhighlight %}

#### Merging

{% highlight ruby %}
git merge samtoolsBlock
{% endhighlight %}`

#### Check and resolve the conflicts generated

You are responsible to merge those conflicts manually by editing the files shown by `git status`.

{% highlight ruby %}
 git status
{% endhighlight %}

#### Commit and push the changes and the merge on the distant server

{% highlight ruby %}
 git commit -m "Branch merging samtoolsBlock-master" -a
 git push https://github.com/SouthGreenPlatform/TOGGLE-DEV.git master `
{% endhighlight %}

### Remove a branche

#### on the remote server

{% highlight ruby %}
 git push https://github.com/SouthGreenPlatform/TOGGLE-DEV.git :nom_branche_a_suppr `
{% endhighlight %}`

#### on our local copy

{% highlight ruby %}
git branch nom-branche_a_suppr -d
{% endhighlight %}`

## Back to the change just before the last commit without losing the work done

{% highlight ruby %}
# create one branch
git branch readDir

# move on this branch
git checkout readDir

# push
git push https://github.com/SouthGreenPlatform/TOGGLE-DEV.git readDir `

# back to the former branch
git checkout dev

# Revert the commit (number given on the terminal)
git revert d10a97d

# Push
git push

# Back to the branch
git checkout readDir
{% endhighlight %}


##  +DIVERSES COMMANDES+

To get status

{% highlight ruby %}
git status
{% endhighlight %}

To get log

{% highlight ruby %}
git log
git log --graph --pretty=format:'%C(red)%h%Creset -%C(yellow)%d%Creset %s %C(green)(%cr) %C(bold blue)<%an>%Creset' --abbrev-commit
{% endhighlight %}

Pour mettre la commande du dessus en alias dans git (exemple avec git lg)
` git config --global alias.lg "log --color --graph --pretty=format:'%C(red)%h%Creset -%C(yellow)%d%Creset %s %C(green)(%cr) %C(bold blue)<%an>%Creset' --abbrev-commit" `

* pour virer un fichier du git quand il est trop gros (et qu'on ne peut plus pusher)

` git filter-branch --index-filter 'git rm --cached --ignore-unmatch DATA/expectedData/snpEffdata/MSU6.1/sequences.fa' --prune-empty --tag-name-filter cat -- --all `

## Docs:
> https://ccwiki.in2p3.fr/developpements:formation:git
> http://www.moussu.fr/git/
