Version 6.2 : 14/01/2019
-----------
- Modification des fichiers gen/Makefile* :
==> Création du répertoire obj par le Makefile 
    de sorte à compenser que Github ne peut pas fournir de répertoire vide dans l'archive fournie au téléchargement.
- Modification du fichier inc/SOS.h :
==> Le nombre de couches du profil par défaut passe de 26 à 100 (constante SOS_OS_NT dans SOS.h)
    pour ne pas avoir des couches de trop grande épaisseur optique.

Version 6.1 : 05/09/2017
-----------
- Le code SOS passe sous licence « GNU GPL V3.0 or later »:
==> Ajout du fichier COPYING.txt au niveau du répertoire des codes sources (src), des éxécutables (exe) et de la documentation (doc).
==> Introduction du copyright dans le manuel utilisateur 
    (version française et anglaise, répertoire doc).
==> Introduction du copyright dans les codes sources (répertoire src).
- Modification d'une variable pour la lecture du fichier des angles (routines SOS_AEROSOLS.F, SOS_SURFACE.F et SOS.F)

Version 6.0 : 20/01/2016
-----------
- Adaptation du code pour une compilation sous linux / gfortran
==> Modification du fichier "Aerosols" de sortie de SOS_AEROSOLS
    et du fichier utilisateur de fonctions de phase externes. 
- Modification de la définition des granulométries LND 
en faveur de la définition usuelle des LND.
==> La valeur de la variance change en paramètre d'entrée.
Sig_new = Sig_old * ln(10) pour ln = log népérien.
Voir le manuel pour plus de précisions.
- Quelques corrections d'incohérence de sortie en erreur pour cause de
  paramètres non définis (mais inutiles pour le cas de simulation souhaité).
- Changement de la précision du fichier de profil : l'épaisseur optique 
qui était codée en F9.5 devient codée en F10.6 (pour régler des problèmes
de couche à épaisseur optique interprétée nulle, en milieu très ténu).
- Correction de quelques anomalies dans les fichiers de spectre d'indice de réfraction 
des modèles d'aérosols de Shettle&Fenn et de la WMO, ainsi que dans le fichier de 
granulométrie des modèles de Shettle&Fenn (particules "Large Rural" à 70% d'humidité relative).
- Correction du calcul de la réflexion solaire directe sur mer plate.


Version 5.2 : 23/10/2014
-----------
Ajout d'un makefile pour gfortran


Version 5.1 : 24/06/2010
-----------

- Les étiquettes de FORMAT en ,X, passent en ,1X, pour éviter des erreurs de compilation.
- Correction pour effacer le fichier FICOS_TMP à la fin des calculs de transmission.



Version 5.0 : 01/06/2010
-----------

- Complete reprise du parametrage : utilisation de "-keyword value".

- Introduction d'une fonction SOS_ANGLES.F qui gere la definition des angles de Gauss
  et l'utilisation d'angles utilisateurs. Elle gere également les ordres limites des developpements.
	
- Definition de tableaux (luminances, fonctions de phase) sur des domaines plus importants 
  que le domaine utilisé : permet un usage du code avec une possibilité de modification du 
  nombre d'angles utilisés sans recompilation (et l'ajout d'angles utilisateurs variables).
  
- Adaptation de la routine SOS_AEROSOLS.F pour permettre l'utilisation d'un fichier 
  de fonctions de phase externes (utile pour des particules non-spheriques).
  
- Adaptation de la routine SOS.F pour permettre le calcul des transmissions diffuses.
	    
- Tracabilité des constantes par copie du fichier SOS.h sur l'espace de compilation 
  et sur l'espace des resultats (modification du main_SOS.ksh).

- Mise a jour du Manuel Utilisateur + version anglaise.

       
######################################################       
       
Version 4.1 : 08/09/2008
-----------

- Intègre une prise en compte des modeles bimodaux (Dubovik) avec correctifs / version 4.0.

- Evolutions depuis la version 4.0 :
     - Adaptation de la routine SOS_AEROSOLS.F (VERSION:2.1):
        *  Correction de depassements de zone d'indentation pour compilation f77.
	
        * Suppression de la sortie des calculs pour fin de fichier de MIE atteinte.
          Car si c'est le cas, il s'agit d'une erreur : le fichier de MIE est 
          incomplet (probable manque d'espace disque à sa génération).
          --> Gestion du cas d'erreur.
	  
        * Correction pour le calcul des proportions de composants des modèles bimodaux
          Dubovik :
	  Cas du calcul à partir des rapports d'épaisseur optique des deux modes
          (fin et grossier) et de l'épaisseur optique totale pour une longueur
           d'onde de référence 
	  --> utilisation de l'indice de réfraction des particules pour la longueur 
	  d'onde de référence (uniquement).
	    
     - Evolution du script main_SOS.ksh pour l'appel de la fonction SOS_AEROSOLS.exe
       (gestion des arguments)


