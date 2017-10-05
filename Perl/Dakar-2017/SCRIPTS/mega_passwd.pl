#!/usr/bin/perl
#
use strict;
use warnings;

my $passwd_file="password.txt";
my %user;

#ouverture du fichieri de mot de passe
open(my $fin , "<", $passwd_file) or die "Je n'arrive pas a ouvrir le fichier $passwd_file\n";

#Lecture du fichier de mot de passe ligne par ligne
while (my $line = <$fin>)
{
   chomp($line); # Je retire le saut de ligne
   #print "----> $line \n\n";
   
   # Je stocke ma ligne dans un tableau en utilisant la fonction split
   # Le separateur de ligne dans mon fichier passwrd.txt est :
   my @tab = split /:/ , $line;
   #print "Ma ligne dans le tableau : @tab \n\n";
   #print "J affiche ma colonne 1 & 8 : $tab[0] ----- $tab[7]--- \n\n";

   # Insertion de l utilasateur dans le hash %user
   $user{$tab[0]}=$tab[7]; #cle = login / valeur = mdp
  
}

# affichage %user
foreach my $nom (keys %user)
{
   print "J affiche la cle de mon hash user : $nom ------ $user{$nom} \n\n";
}


#fermture du fichier password
close($fin);
