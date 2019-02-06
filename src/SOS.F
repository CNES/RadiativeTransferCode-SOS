C******************************************************************************
C* Copyright 2017, Centre National d'Etudes Spatiales (CNES)
C* 
C* This file is part of SOS.
C* 
C* SOS is free software: you can redistribute it and/or modify
C* it under the terms of the GNU General Public License as published by
C* the Free Software Foundation, either version 3 of the License, or
C* (at your option) any later version.
C* 
C* SOS is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C* GNU General Public License for more details.
C* 
C* You should have received a copy of the GNU General Public License
C* along with SOS. If not, see <http://www.gnu.org/licenses/>.
C******************************************************************************


C******************************************************************************
C* FICHIER: SOS.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Module principal du code des OS
C*           - gere les entrees / sorties et options de traitements,
C*           - lit les parametres de granulometries aerosols,
C*           - lit le profil atmospherique et l'adapte a une eventuelle 
C*             troncature de la fonction de phase des aerosols,
C*           - resout l'equation de transfert radiatif par ordre 
C*             successif d'interaction (diffusion/reflexion),
C*           - recombine les parametres de Stokes a partir de leur 
C*             decomposition en series de Fourier.
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et developpement : CS SI, B.Lafrance, C. Casteras, C. Ruffel.
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0: Mise a niveau du code des OS du LOA
C*             2.0: Lecture des parametres aerosols et profils, 
C*                  Mise a jour des commentaires.
C*             2.1: 07/07/2008 
C*                  Modification du format de lecture des parametres WIND et TETAS
C*                  pour correction d'une erreur de décodage si valeur entière.
C*			WIND passe de F4.1 en F4
C*                      TETAS passe de F9.5 en F9
C*
C* MOD:VERSION:3.0: 24/11/2009 
C*       - Modification du passage des arguments : 
C*           - passage par couples (mot clef, valeur). 
C*           - adaptation de la gestion des cas d'erreurs. 
C*       - Conversion des messages d'erreurs en anglais.
C*
C*
C* MOD:VERSION:3.1: 22/12/2009 
C*                  Modification du format de lecture des parametres Alpha et Beta du modele Nadal
C*                  pour correction d'une erreur de décodage si valeur entière.
C*			ALPHA passe de F6.4 en F6
C*			BETA  passe de F4.1 en F4
C*
C*                  
C* MOD:VERSION:4.0: 10/01/2010 
C*
C*   - Suppression de la constante SOS_OS_FICANGLE (de SOS.h) 
C*     au profit d'un passage par argument du nom du fichier des angles (produit par SOS_ANGLES).
C*
C*   - Modification de la taille des tableaux de luminance :
C*     ==> Passage de SOS_OS_NBMU a SOS_OS_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'a LUM_NBMU egal au nombre effectif d'angles utilises.
C*
C*   - Suppression de la constante  SOS_OS_NB (de SOS.h)  
C*     au profit : d'un passage par argument de la valeur OS_NB (lue dans le fichier des angles)
C*
C*   - Modification de la taille des tableaux de decomposition en polynomes de Legendre :
C*     ==> Passage de SOS_OS_NB a SOS_OS_NB_MAX.
C*     ==> Boucle sur les ordre jusqu'a OS_NB. 
C* 
C*   - Gestion des angles : 
C*         * Suppression de la lecture d'un fichier d'angles de Gauss
C*         * Ajout d'un appel a la routine SOS_ANGLES pour calculer les angles de Gauss et leurs poids
C*           et pour ajouter des angles utilisateurs (de poids nuls) : angles utiles = LUM_NBMU.
C*         * Suppression de la lecture de SOS.Imus et SOS.thetas : la recuperation de ces parametres
C*           est faite via la lecture du contenu du fichier des angles.
C*         * Ajout des parametres d'entree optionnels SOS.ResFileUp.UserAng et SOS.ResFileDown.UserAng
C*           pour l'enregistrement des luminances de sortie uniquement pour les angles "utilisateur"
C*           --> Adaptation de la routine SOS_TRPHI_OPTION.                                  
C*
C*   - Modification de la taille des tableaux de luminance :
C*     ==> Passage de SOS_OS_NBMU a SOS_OS_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'a LUM_NBMU egal au nombre effectif d'angles utilises.
C*
C*   - Routine SOS_TRPHI_OPTION : 
C*         * Modification de la ligne de lecture des tableaux de
C*           la  matrice de reflexion a l'ordre IS de la decomposition en series de Fourier :
C*           ==> limitation aux bornes utiles -LUM_NBMU:LUM_NBMU
C*         * Passage du tableau IND_ANGOUT definissant les angles pour lesquels une sortie
C*           est requise dans les fichiers optionnels SOS.ResFileUp.UserAng 
C*           et SOS.ResFileDown.UserAng
C*
C*   - Adaptation des routines SOS_OS & SOS_TRPHI_OPTION au passage du parametre
C*     LUM_NBMU (nb d'angles utiles) et de OS_NB (ordre max des developpements en SF).
C*
C*   - Adaptation au passage du parametre FICTRANS (-SOS.Trans) pour sortir les valeurs
C*     de transmission.
C*
C* MOD:VERSION:4.1: 24/06/2010 
C*   -  Les etiquettes de FORMAT en ,X, passent en ,1X, pour eviter des erreurs de compilation.
C*   -  Correction pour effacer le fichier FICOS_TMP a la fin des calculs de transmission.
C*
C* MOD:VERSION:4.2: 08/12/2015
C*    - Suppression des formats à la lecture des paramètres d'entrée en faveur 
C*      de READ(PROG_ARG,*,ERR=NNN)
C*      ==> Permet une compilation sous gfortran et évite les mauvaises conversions 
C*          si la valeur d'entrée ne respecte pas strictement le format attendu 
C*          (par exemple, saisie d'un entier N au lieu d'un réel N. 
C*           aboutissant à une valeur interprétée fausse). 
C*
C*    - Changement de la manière de lire les fichiers formatés (Angles, ...)
C*      pour être compatible avec le compilateur gfortran.
C*
C*    - Ajustement des lignes utilisant SOS_PI pour éviter une troncature 
C*      au décodage de la constante.
C*
C*    - Modification de la récupération du nom du fichier optionnel des 
C*      angles utilisateur (FICANGLES_USER) pour compatibilité avec gfortran/linux
C*
C*    - Ajout d'un chiffre significatif de plus aux valeurs lues du profil pour 
C*      être cohérent avec le nouveau format d'enregistrement des fichiers de
C*      profil atmosphèrique.
C*      Etiquette 70 : Passage du format FORMAT(2X,I4,3(F9.5)) 
C*                     au format FORMAT(2X,I4,3(F10.6))
C*
C*
C* MOD:VERSION:4.3: 06/09/2017
C*     Changes to be compliant with a f77 compilation
C*         - Change a variable type to read the file FICANGLES  
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    SOS_LENDIR : Longueur des noms de repertoires.
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_LENKEYWORD : Taille maximale des Keywords pour le passage des arguments.
C    SOS_OS_NBMU_MAX : Nombre maximal d'angles positifs pour la resolution de l'ETR
C    SOS_OS_NB_MAX :  Valeur maximale pour l'ordre limite du developpement en 
C                     polynomes de Legendre.
C    SOS_OS_NT  : Nombre de couches du profil atmospherique.
C
C Constantes spécifiques :
C    LENLINE  : Longueur maximale d'une ligne de données lue dans un fichier
C----------------------------------------------------------------------------
#include "SOS.h"
#define LENLINE 1000

C==============================================================================
C PROGRAMME: SOS
C ==========      
C      Programme principal du code des OS : gere les E/S et options de traitements,
C      enchaine les differentes routines.
C
C Donnees en entree du programme 
C ------------------------------
C    -->Les parametres/arguments du programme
C
C        Les arguments sont passes sous la forme d'une liste de couples (-Keyword Value)
C
C        Definition des mots clef, statut obligatoire / option, format de la valeur associee :
C         -SOS.AngFile  ==> Obligatoire (chemin complet au fichier ANGLES_LUM) 
C             FICANGLES (CHARACTER*SOS_LENFIC2) : nom complet du fichier des angles utilises 
C                                                 pour les calculs de luminance et BRDF
C                 (repertoire + nom fichier + extension)
C
C         -SOS.AerFile ==> Obligatoire (chemin complet au fichier AEROSOLS)
C	      FICGRANU (CHARACTER*SOS_LENFIC2) : nom complet du fichier resultat de GRANU
C                                 (repertoire + nom fichier + extension)
C
C         -SOS.ProfileFile ==> Obligatoire (chemin complet au fichier PROFIL)
C             FICPROFIL(CHARACTER*SOS_LENFIC2) : nom complet du fichier resultat de PROFIL
C                                 (repertoire + nom fichier + extension)
C
C         -SOS.ResBin ==> Obligatoire (chemin complet au fichier resultat binaire des OS)
C             FICOS (CHARACTER*SOS_LENFIC2) : Nom complet du fichier resultat des OS.
C                               (repertoire + nom fichier + extension)
C       
C         -SOS.ResFileUp ==> Obligatoire (chemin complet au fichier resultat texte des OS)
C             FICUP (CHARACTER*SOS_LENFIC2) : Nom complet du fichier resultat pour le champ 
C                                             montant.
C 					     (repertoire + nom fichier + extension)	
C       
C         -SOS.ResFileDown ==> Obligatoire (chemin complet au fichier resultat texte des OS)
C             FICDOWN (CHARACTER*SOS_LENFIC2) : Nom complet du fichier resultat pour le champ 
C                                               descendant.
C 					        (repertoire + nom fichier + extension)	      
C
C         -SOS.Log ==> Obligatoire (0 ou chemin complet au fichier Trace)
C             FICLOG (CHARACTER*SOS_LENFIC2): nom complet du fichier trace des OS
C                                            (repertoire + nom fichier + extension)
C                                            ou =0, si pas de  trace	 
C
C         -SOS.Trans  ==> Optionnel (chemin complet)
C             FICTRANS (CHARACTER*SOS_LENFIC2): nom complet du fichier des valeurs de transmissions
C                                            (repertoire + nom fichier + extension)
C    
C         -SOS.ResFileUp.UserAng ==> Optionnel (chemin complet au fichier resultat texte des OS)
C                                    Necessite que -SOS.AngFile soit defini.
C             FICUP_ANGLES_USER (CHARACTER*SOS_LENFIC2) : Nom complet du fichier resultat pour 
C                                                         le champ  montant.
C 					                  (repertoire + nom fichier + extension)	
C       
C         -SOS.ResFileDown.UserAng ==> Optionnel (chemin complet au fichier resultat texte des OS)
C                                      Necessite que -SOS.AngFile soit defini.
C             FICDOWN_ANGLES_USER (CHARACTER*SOS_LENFIC2) : Nom complet du fichier resultat pour 
C                                                           le champ descendant.
C 					                    (repertoire + nom fichier + extension)     
C
C         -SOS.MDF ==> Obligatoire
C            PROG_ARG (CHARACTER*10)  : Facteur de depolarisation (moleculaire) : RON
C                                       (format de décodage : libre --> Float)
C
C         -SURF.Alb ==> Obligatoire
C            PROG_ARG (CHARACTER*10)  :	Albedo de surface RHO.
C                                       (format de décodage : libre --> Float)   	 	 	
C 	 	 	
C         -SURF.Type ==> Obligatoire
C            PROG_ARG (CHARACTER*10)  : Type de surface ISURF 
C                                       (format de décodage : libre --> Integer) 
C      		        0 :surface lambertienne d'albedo RHO.
C			1 :surface lambertienne + glitter.
C			2 :surface lambertienne + reflexion de Fresnel.
C		        3 :surface lambertienne + BRDF de Roujean.
C			4 :surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
C			5 :surface lambertienne + BRDF de Roujean + BPDF de Breon.
C			6 :surface lambertienne + BRDF de Roujean + BPDF de Nadal.		 	
C
C         -SURF.Ind ==> Associe a : -SURF.Type 1 2 4 5 ou 6
C            PROG_ARG (CHARACTER*10)  : Indice de refraction de la surface / air.
C                                       (format de décodage : libre --> Float)
C
C         -SURF.File ==> Associe a : -SURF.Type 1 3 4 5 ou 6 
C             FICSURF (CHARACTER*SOS_LENFIC2) : nom complet du fichier contenant la 
C                   matrice de reflexion de surface  (fichier Glitter ou autres BRDF/BPDF).
C                   (repertoire + nom fichier + extension)
C                    => DEFAULT si on utilise un fichier de la base des BRDF / BPDF.	 
C	 	 
C         -SURF.Dir ==> Associe a : -SURF.Type 1 3 4 5 ou 6 et -SURF.File "DEFAULT"
C             REPSURF (CHARACTER*SOS_LENDIR) : Repertorie de stockage des fichiers de 
C                   matrice de reflexion de surface (fichier Glitter ou autres BRDF/BPDF).
C		
C         -SURF.Glitter.Wind ==> Associe a : -SURF.Type 1
C            PROG_ARG (CHARACTER*10)  : Vitesse du vent m/s.	
C                                       (format de décodage : libre --> Float)
C		
C         -SURF.Roujean.K0 ==> Associe a : -SURF.Type >= 3
C            PROG_ARG (CHARACTER*10)  : Parametre k0 du modele de BRDF de Roujean. 
C                                       (format de décodage : libre --> Float)
C		
C         -SURF.Roujean.K1 ==> Associe a : -SURF.Type >= 3
C            PROG_ARG (CHARACTER*10)  : Parametre k1 du modele de BRDF de Roujean.  
C                                       (format de décodage : libre --> Float)
C		
C         -SURF.Roujean.K2 ==> Associe a : -SURF.Type >= 3
C C            PROG_ARG (CHARACTER*10)  : Parametre k2 du modele de BRDF de Roujean.  
C                                       (format de décodage : libre --> Float)
C
C         -SURF.Nadal.Alpha ==> Associe a : -SURF.Type 6
C            PROG_ARG (CHARACTER*10)  : Parametre alpha du modele de BPDF de Nadal.  
C                                       (format de décodage : libre --> Float)
C
C         -SURF.Nadal.Beta ==> Associe a : -SURF.Type 6
C            PROG_ARG (CHARACTER*10)  : Parametre beta du modele de BPDF de Nadal.   
C                                       (format de décodage : libre --> Float)
C 
C         -SOS.OutputLevel ==> optionnel 
C                              Si non defini : le code prend par defaut -1
C            PROG_ARG (CHARACTER*10)  :	 Niveau de sortie (IOUT):
C                                       (format de décodage : libre --> Integer) 
C			-1 =>  Sortie standard:
C				 Luminance emergeante au TOA
C				 + descendante au niveau du sol.
C			 n =>  Luminance montante et descendante
C			       pour le niveau n (0 <= n <= SOS_OS_NT).	     	 
C		 
C         -SOS.Ipolar ==> optionnel 
C                         Si non defini : le code prend par defaut 1 (polarisation incluse)
C            PROG_ARG (CHARACTER*10)  : Option de coupure de la polarisation.
C                                       (format de décodage : libre --> Integer)
C                       0 : simulation sans polarisation.
C			1 : simulation normale avec polarisation.
C
C         -SOS.IGmax ==> Obligatoire
C            PROG_ARG (CHARACTER*10)  : Ordre maximal d'interaction.
C                                       (format de décodage : libre --> Integer) 
C
C         -SOS.View ==> Obligatoire
C            PROG_ARG (CHARACTER*10)  : Indice de type de sortie (ITRPHI) :
C                                       (format de décodage : libre --> Integer) 
C		        1 =>  Plan de visee pour une valeur d'azimut.		           
C		        2 =>  Diagramme polaire en fonction de la direction
C			      d'azimut et de visee.
C
C         -SOS.View.Phi ==> Associe a : -SOS.View 1
C            PROG_ARG (CHARACTER*10)  : Azimut relatif en degres (pour ITRPHI=1). 
C                                       (format de décodage : libre --> Float)
C
C         -SOS.View.Dphi ==> Associe a : -SOS.View 2
C            PROG_ARG (CHARACTER*10)  : Pas sur l'azimut en degres (pour ITRPHI=2)
C                                       (format de décodage : libre --> Integer) 
C  
C
C    -->Fichiers
C       Le fichier des angles et parametres internes pour les limites des developpements.
C       Le fichier Matrice de reflexion de surface
C       Le fichier Aerosols
C       Le fichier Profil
C       Le fichier des transmissions
C
C
C Resultats fournis par le programme
C ----------------------------------
C   Fichier resultat des OS (decomposition en series de Fourier)
C   Fichier resultat pour le champ montant      | de base pour tous les angles utiles 
C   Fichier resultat pour le champ descendant   | + optionnellement juste pour les angles utilisateurs
C   Fichier resultat des transmissions directes et diffuses (option)
C   Fichier de trace (optionnel)
C              (No d'unite logique du fichier trace =99)
C
C   Si un ou des fichiers de meme noms existent deja, ils sont ecrases
C
C
C 
C Description des fichiers utilises
C ---------------------------------
C
C   -->Contenu du fichier BRDF/BPDF
C       Fichier binaire non formatte, cree par SOS_SURFACE
C       Il contient la matrice de reflexion de surface 
C       (fichier Glitter ou autres BRDF/BPDF):
C       c'est a dire les elements Pij de la matrice de reflexion :
C              P11(LUM_NBMU,LUM_NBMU) (reel),
C              P12(LUM_NBMU,LUM_NBMU) (reel),
C              P13(LUM_NBMU,LUM_NBMU) (reel),
C              P21(LUM_NBMU,LUM_NBMU) (reel),
C              P22(LUM_NBMU,LUM_NBMU) (reel),
C              P23(LUM_NBMU,LUM_NBMU) (reel),
C              P31(LUM_NBMU,LUM_NBMU) (reel),
C              P32(LUM_NBMU,LUM_NBMU) (reel),
C              P33(LUM_NBMU,LUM_NBMU) (reel)
C
C        --> Enregistrements sur le domaine limite aux angles utiles
C                                   --------------------------------
C
C   -->Contenu du fichier resultat des OS
C       Fichier binaire non formatte, cree par SOS_OS
C       Il contient les parametres de Stokes Q , U et I pour un ordre S
C       du developpement en series de Fourier, en fonction de l'angle mu
C       (valeurs en double precision).
C
C   -->Contenu du fichier temporaire FICOS_TMP
C       Identique au fichier de resultat des OS.
C       Fichier cree uniquement si l'option de sortie des transmissions est
C       activee. Fichier temporaire, non utilisee en tant que tel. Detruit en 
C       fin de traitement.
C
C   -->Contenu du fichier resultat pour le champ montant
C       Fichier Ascii, format d'un enregistrement 4(E13.6,4X) ou 5(E13.6,4X)
C       2 cas : ITRPHI = 1, 
C               on obtient le plan de visees pour PHIos .
C               Nombre d'enregistrements : 
C                   un pour chaque valeur de -LUM_NBMU a LUM_NBMU
C               Chaque enregistrement contient 4 valeurs :
C                   angle de visee(en degres) et parametres de Stokes I,Q,U 
C
C               ITRPHI = 2, 
C               on obtient un diagramme polaire en fonction de l'angle 
C               de visee (teta) et de l'angle azimutal (phi) 
C               Nombre d'enregistrements :
C                   un pour chaque valeur de 1 a LUM_NBMU
C               Chaque enregistrement contient 5 valeurs :
C                   coordonnees polaires X et Y et parametres de Stokes I,Q,U 
C
C   -->Contenu du fichier resultat pour le champ descendant
C      Fichier Ascii, format d'un enregistrement 4(E13.6,4X) ou 5(E13.6,4X)
C      Meme description que pour le fichier resultat du champ montant
C
C    -->Contenu du fichier des angles:
C       Fichier ascii formatte, cree par SOS_ANGLES.
C       Il contient l'information sur les angles a utiliser 
C       et sur les ordres des developpements adaptes : 
C		NB_TOTAL_ANGLES :    nb d'angles total a utiliser
C               NB_GAUSS_ANGLES :    nb d'angles de Gauss
C               ANGLES_USERFILE :  Nom du fichier des angles utilisateur
C	        SOLAR ZENITH ANGLE : Angle zenithal solaire (degres)
C		INTERNAL_IMUS :      Numero d'indice de table correspondant
C               INTERNAL_OS_NB :   Ordre des developpements des fonctions de phase 
C				   en Polynomes de Legendre a utiliser.
C               INTERNAL_OS_NS :     Ordre maximal a utiliser pour le developpement en
C                                    polynomes de Legendre des elements de la matrice de Fresnel 
C                                    et en series de Fourier des calculs de luminance.
C                                    (non utilise dans ce programme)
C               INTERNAL_OS_NM :     Ordre maximal du developpement en series de Fourier  
C                                    de la fonction G des calculs de matrice de reflexion.
C                                    (non utilise dans ce programme)
C               INDEX COS_ANGLE WEIGHT OUTPUT: lignes donnant le numero d'angle, son cosinus 
C                                              et son poids, s'il est donne en sortie "utilisateur"
C                                              au format I4,1X,2D21.14,1X,I4
C         
C   -->Contenu du fichier resultat des valeurs de transmissions
C      Fichier Ascii, format d'un enregistrement 
C
C          Solar Zenithal Angle : Valeur (pas de format specifie)
C          Direct transmission  TOA -> Surface : Valeur (pas de format specifie)
C          Diffuse transmittance : TOA -> Surface
C               Valeur : FORMAT(12h    thetas =,1x,f6.3,3x,12htd(thetas) =,1x,f7.4)  
C          Diffuse transmittance : surface -> TOA
C               1 ligne par angle de visee
C               Valeur : FORMAT(12h    thetav =,1x,f6.3,3x,12htd(thetav) =,1x,f7.4)  
C
C
C Common utilise:
C --------------
C  Aucun
C
C
C Cas d'erreur :
C ------------        
C     -Decodage des parametres incorrect
C     -Erreur dans un sous-programme
C     -Erreur a l'ouverture du fichier Trace
C     -Parametres non valables
C
C     Affichage d'un message a l'ecran, arret du programme et 
C     retour du status 1 au shell
C   
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel aux routines:
C     - SOS_OS
C     - SOS_TRPHI_OPTION
C
C============================================================================
      PROGRAM SOS
      
      IMPLICIT NONE

C* Definition des variables                         
C*-----------------------------------------------------------------
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
                                             !Cosinus des angles de Gauss (et utilisateur).
      DOUBLE PRECISION GA(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
                                             !Poids des angles de Gauss (et utilisateur).				
 
 
      DOUBLE PRECISION TETAS	! Angle zenithal solaire en degres (0 < tetas < 90).
      				! (angle complementaire : mus= -cos(tetas))

      DOUBLE PRECISION RHO	! Albedo de surface.
      
      DOUBLE PRECISION IND_SURF  ! Indice de refraction de la surface / air.
      
      DOUBLE PRECISION WIND      ! Vitesse du vent en m/s.
      
      DOUBLE PRECISION K0_ROUJEAN ! Parametre k0 du modele de BRDF de Roujean.
      DOUBLE PRECISION K1_ROUJEAN ! Parametre k1 du modele de BRDF de Roujean.
      DOUBLE PRECISION K2_ROUJEAN ! Parametre k2 du modele de BRDF de Roujean.
      
      DOUBLE PRECISION ALPHA_NADAL ! Parametre alpha du modele de BPDF de Nadal.
      DOUBLE PRECISION BETA_NADAL  ! Parametre beta du modele de BPDF de Nadal.
      
      DOUBLE PRECISION RON	! Facteur de depolarisation (moleculaire).
      
      DOUBLE PRECISION A	! Coefficient de troncature.
      DOUBLE PRECISION PIZ	! Albedo de simple diffusion.
      DOUBLE PRECISION PIZTR	! Albedo de simple diffusion apres troncature.

      DOUBLE PRECISION ALPHA(0:SOS_OS_NB_MAX)	! Parametres alpha(k).
      DOUBLE PRECISION BETA(0:SOS_OS_NB_MAX)	! Parametres beta(k).
      DOUBLE PRECISION GAMMA(0:SOS_OS_NB_MAX)	! Parametres gamma(k).
      DOUBLE PRECISION ZETA(0:SOS_OS_NB_MAX)	! Parametres zeta(k).

      DOUBLE PRECISION VA,VATR
      DOUBLE PRECISION VR

      DOUBLE PRECISION H(0:SOS_OS_NT)	   ! Epaisseur optique pour chaque niveau
      DOUBLE PRECISION HTR(0:SOS_OS_NT)	   ! Epaisseur optique tronquee pour chaque niveau
      DOUBLE PRECISION XDEL(0:SOS_OS_NT)   ! Pourcentage d'aerosols par niveau
      DOUBLE PRECISION YDEL(0:SOS_OS_NT)   ! Pourcentage de molecules par niveau

      DOUBLE PRECISION TTOT_TRONC ! Epaisseur optique totale au sol (tronquee).
      DOUBLE PRECISION TTOT_VRAI! Epaisseur optique totale au sol (non tronquee).
      DOUBLE PRECISION TAUOUT	! Epaisseur optique au niveau de
      				! sortie des resultats pour le champ montant (tronquee).
				
      DOUBLE PRECISION PHIOS	! Azimut relatif en degres (pour ITRPHI=1).  

      DOUBLE PRECISION TDIFMUS  !Transmission diffuse pour la direction d'incidence solaire
      				!de l'atmosphere equivalente
      DOUBLE PRECISION TDIR_TRONC  !Transmission directe pour la direction d'incidence solaire
      				   !de l'atmosphere equivalente
      DOUBLE PRECISION TDIR_VRAI   !Transmission directe pour la direction d'incidence solaire
      				   !de l'atmosphere vraie
				               
     
      INTEGER*4 N0	! Numero de mu de Gauss pour l'angle
      			! zenithal solaire.
			! (angle complementaire : mus = -rmu(N0))
			! (N0 > 0 => pour en tenir compte).

      INTEGER*4 ISURF	! Type de surface :
			  ! 0 :surface lambertienne d'albedo RHO.
			  ! 1 :surface lambertienne + glitter.
			  ! 2 :surface lambertienne + reflexion de Fresnel.
			  ! 3 :surface lambertienne + BRDF de Roujean.
			  ! 4 :surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
			  ! 5 :surface lambertienne + BRDF de Roujean + BPDF de Breon.
			  ! 6 :surface lambertienne + BRDF de Roujean + BPDF de Nadal.
				
      INTEGER*4 IOUT	! Type de sortie :
			!      -1 =>  Sortie standard:
			!             Luminance emergeante au TOA
			!           + descendante au niveau du sol.
			!       n =>  Luminance montante et descendante
			!             pour le niveau n (0 <= n <= SOS_OS_NT).

      INTEGER*4 IGMAX	! Ordre maximal de diffusion.
      
      INTEGER*4 IPOLAR	! Option de simulation sans polarisation.
      			! 0 : calculs sans polarisation.
			! 1 : calculs normaux avec polarisation.

      INTEGER*4 ITRPHI	! Indice de type de sortie.
			!   1 =>  Plan de visee pour une valeur d'azimut.
			!   2 =>  Diagramme polaire en fonction de la direction
			!         d'azimut et de visee.
			
      INTEGER*4 PAS_PHI	! Pas sur l'azimut en degres (pour ITRPHI=2)

      INTEGER*4 IGLI	! Indice precisant si les simulations
      			! introduisent le Glitter :
			!    1 => pour calcul du Glitter (il faut alors utiliser N0 > 0).
			!    Sinon => pas de simulation du Glitter.            
	
      INTEGER*4 IFRESNEL ! Indice precisant si les simulations
      			 ! introduisent la reflexion de Fresnel sur mer plate :
			 !    1 => pour un calcul de la reflexion de Fresnel.
			 !    Sinon => pas de prise en compte.            
	
      INTEGER*4 IROUJEAN ! Indice precisant si les simulations
      			 ! utilisent le modele de reflexion de Roujean :
			 !    1 => pour calcul de BRDF Roujean (il faut utiliser N0 > 0).
			 !    Sinon => pas de simulation de BRDF Roujean. 
			 
      INTEGER*4 IRONDEAUX ! Indice precisant si les simulations
      			  ! utilisent le modele de reflexion de Rondeaux - Herman :
			  !    1 => pour calcul de BPDF RH (il faut utiliser N0 > 0).
			  !    Sinon => pas de simulation de BPDF RH. 			 
			
      INTEGER*4 IBREON   ! Indice precisant si les simulations
      			 ! utilisent le modele de reflexion de Breon :
			 !    1 => pour calcul de BPDF Breon (il faut utiliser N0 > 0).
			 !    Sinon => pas de simulation de BPDF Breon. 

      INTEGER*4 INADAL   ! Indice precisant si les simulations
      		         ! utilisent le modele de reflexion de Nadal :
		         !    1 => pour calcul de BPDF Nadal (il faut utiliser N0 > 0).
			 !    Sinon => pas de simulation de BPDF Nadal. 				 
	
      INTEGER*4 IMAT_SURF ! Indice precisant si les simulations
      			  ! utilisent un modele de reflexion BRDF ou BPDF :
			  !    1 => pour calcul de BRDF ou BPDF (il faut utiliser N0 > 0).
			  !    Sinon => pas de simulation de BRDF ou BPDF. 	
			  			  			 											
      INTEGER*4 IBORM	! Ordre maximal du developpement en series de Fourier
      			! en cours de traitement :
			!   _ 2 pour une atmosphere purement moleculaire,
			!   _ NB en presence d'aerosols.

      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon

      INTEGER*4 K       ! Indice de tableau
      INTEGER*4 I,J   	! Indice de tableau
      
      INTEGER*4 LDIRSURF ! longueur de la chaine de caractere REPSURF.
      INTEGER*4 LFICSURF ! longueur de la chaine de caractere FICSURF.
      
      INTEGER*4 NBMU_GAUSS      ! Nombre d'angles de Gauss utilises       
      INTEGER*4 LUM_NBMU        ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB 	        ! Ordre des developpements des fonctions de phase 
			        ! en Polynomes de Legendre.      

      INTEGER*4 IND_ANGOUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
                                ! Indices des tableaux d'angles pour lesquels il doit y avoir 
      				! une sortie dans le fichier résultat des angles utilisateur 
				! (1 : pour une sortie requise)				

			   
      INTEGER*4 IARGC           !Type de la fonction IARGC
      INTEGER*4 NPARAM 		!Numero du parametre en entrée (type INTEGER*4 obligatoire
                                !pour le f77)
      
      INTEGER*4 IFIN            ! Nb de caractères significatifs d'une chaine. 

      INTEGER*4 NVAL 		! Integer pour la lecture de fichier ANGLES
     
      INTEGER*2 NBARGS		! Nombre d'arguments en entrée
      INTEGER*2 DEF_FICGRANU	! Indicateur d'affectation d'une valeur a FICGRANU
      INTEGER*2 DEF_FICPROFIL	! Indicateur d'affectation d'une valeur a FICPROFIL
      INTEGER*2 DEF_RON		! Indicateur d'affectation d'une valeur a RON
      INTEGER*2 DEF_RHO		! Indicateur d'affectation d'une valeur a RHO
      INTEGER*2 DEF_ISURF	! Indicateur d'affectation d'une valeur a ISURF
      INTEGER*2 DEF_IND_SURF	! Indicateur d'affectation d'une valeur a IND_SURF
      INTEGER*2 DEF_FICSURF	! Indicateur d'affectation d'une valeur a FICSURF
      INTEGER*2 DEF_REPSURF	! Indicateur d'affectation d'une valeur a REPSURF
      INTEGER*2 DEF_WIND	! Indicateur d'affectation d'une valeur a WIND
      INTEGER*2 DEF_K0_ROUJEAN	! Indicateur d'affectation d'une valeur a K0_ROUJEAN 
      INTEGER*2 DEF_K1_ROUJEAN	! Indicateur d'affectation d'une valeur a K1_ROUJEAN 
      INTEGER*2 DEF_K2_ROUJEAN	! Indicateur d'affectation d'une valeur a K2_ROUJEAN 
      INTEGER*2 DEF_ALPHA_NADAL	! Indicateur d'affectation d'une valeur a ALPHA_NADAL
      INTEGER*2 DEF_BETA_NADAL	! Indicateur d'affectation d'une valeur a BETA_NADAL
      INTEGER*2 DEF_IOUT	! Indicateur d'affectation d'une valeur a IOUT
      INTEGER*2 DEF_IGMAX	! Indicateur d'affectation d'une valeur a IGMAX
      INTEGER*2 DEF_IPOLAR	! Indicateur d'affectation d'une valeur a IPOLAR
      INTEGER*2 DEF_ITRPHI	! Indicateur d'affectation d'une valeur a ITRPHI
      INTEGER*2 DEF_PHIOS	! Indicateur d'affectation d'une valeur a PHIOS
      INTEGER*2 DEF_PAS_PHI	! Indicateur d'affectation d'une valeur a PAS_PHI
      INTEGER*2 DEF_FICOS	! Indicateur d'affectation d'une valeur a FICOS  
      INTEGER*2 DEF_FICUP	! Indicateur d'affectation d'une valeur a FICUP  
      INTEGER*2 DEF_FICDOWN	! Indicateur d'affectation d'une valeur a FICDOWN 
      INTEGER*2 DEF_FICUP_ANGLES_USER	! Indicateurs d'affectation d'une valeur 
      INTEGER*2 DEF_FICDOWN_ANGLES_USER	! a FICUP_ANGLES_USER  et FICDOWN_ANGLES_USER
      INTEGER*2 DEF_FICTRANS    ! Indicateur d'affectation d'une valeur a FICTRANS
      INTEGER*2 DEF_FICLOG	! Indicateur d'affectation d'une valeur a FICLOG  
      INTEGER*2 DEF_FICANGLES   ! Indicateur d'affectation d'une valeur a FICANGLES
      
      
      
      LOGICAL TRACE     ! = vrai,si ecriture dans le fichier trace
      LOGICAL LTA       ! = vrai,si epaisseur optique des aerosols nulle
      LOGICAL EX	! Verifie l'existence d'un fichier

      CHARACTER*SOS_LENFIC2 FICOS       ! Nom du fichier resultat des OS.
      CHARACTER*SOS_LENFIC2 FICOS_TMP   ! Nom du fichier temporaire resultat des OS 
      					! pour les calculs de transmission.
      CHARACTER*SOS_LENFIC2 FICUP	! Nom du fichier resultat pour le champ
      				        ! montant.					
      CHARACTER*SOS_LENFIC2 FICDOWN	! Nom du fichier resultat pour le champ
      				        ! descendant.
      CHARACTER*SOS_LENFIC2 FICUP_ANGLES_USER	! Nom du fichier resultat pour le champ
      				                ! montant (pour les angles utilisateur).					
      CHARACTER*SOS_LENFIC2 FICDOWN_ANGLES_USER	! Nom du fichier resultat pour le champ
      				                ! descendant (pour les angles utilisateur).			
      CHARACTER*SOS_LENFIC2 FICTRANS    ! Nom du fichier optionnel de sortie des transmissions.
      CHARACTER*SOS_LENDIR REPSURF	! Repertoire de stockage des matrice de reflexion de
      				        ! surface (fichier Glitter ou autres BRDF/BPDF).
      CHARACTER*SOS_LENFIC2 FICSURF	! Fichier contenant la matrice de reflexion de
      				        ! surface (fichier Glitter ou autres BRDF/BPDF).
      CHARACTER*SOS_LENFIC2 FICLOG      ! nom du fichier trace des OS
      CHARACTER*SOS_LENFIC2 FICANGLES           ! chemin et nom du fichier des angles effectifs    
      CHARACTER*SOS_LENFIC2 FICANGLES_USER      ! chemin et nom du fichier des angles "utilisateurs"
      
      CHARACTER*SOS_LENFIC2 FICGRANU	! nom du fichier d'entree AEROSOLS
      CHARACTER*SOS_LENFIC2 FICPROFIL	! nom du fichier d'entree PROFIL

      CHARACTER*SOS_LENFIC1 FICS     ! nom recompose du fichier SURFACE a partir 
      				     ! des arguments de simulation de la BRDF / BPDF.

      CHARACTER*SOS_LENKEYWORD KEYWORD	       	! Mot clef du passage d'arguments.

      CHARACTER*LENLINE LIGNE_TEXTE    ! Ligne de caractères pour décodage de paramètres 
      CHARACTER*1000 CVAL	       ! Chaine de caracteres pour la lecture de valeurs inutilisees
             
      CHARACTER*10 PROG_ARG		    ! Arguments du programme


 
C* Lecture des parametres/arguments d'entree et decodage 
C------------------------------------------------------    

C --- Initialisation des indicateurs de definition des parametres (0 pour indefinis)
      DEF_FICGRANU=0
      DEF_FICPROFIL=0
      DEF_RON=0
      DEF_RHO=0
      DEF_ISURF=0
      DEF_IND_SURF=0
      DEF_FICSURF=0
      DEF_REPSURF=0
      DEF_WIND=0
      DEF_K0_ROUJEAN=0
      DEF_K1_ROUJEAN=0
      DEF_K2_ROUJEAN=0
      DEF_ALPHA_NADAL=0
      DEF_BETA_NADAL=0
      DEF_IOUT=0
      DEF_IGMAX=0
      DEF_IPOLAR=0
      DEF_ITRPHI=0
      DEF_PHIOS=0
      DEF_PAS_PHI=0    
      DEF_FICOS=0  
      DEF_FICUP=0
      DEF_FICDOWN=0
      DEF_FICLOG=0
      DEF_FICANGLES=0
      DEF_FICUP_ANGLES_USER=0
      DEF_FICDOWN_ANGLES_USER=0
      DEF_FICTRANS=0
      

C --- Nombre d'arguments en entree
      NBARGS = IARGC()

C --- Boucle de lecture des couples d'arguments       
      DO NPARAM=1,NBARGS, 2
         CALL GETARG(NPARAM,KEYWORD)
	 
	 IF (KEYWORD.EQ."-SOS.AerFile") THEN
	     CALL GETARG (NPARAM+1,FICGRANU)	!Fichier AEROSOLS
	     DEF_FICGRANU=1
	 ENDIF	
	 
	 IF (KEYWORD.EQ."-SOS.ProfileFile") THEN
	     CALL GETARG (NPARAM+1,FICPROFIL)	!Fichier PROFIL
	     DEF_FICPROFIL=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-SOS.ResBin") THEN
	     CALL GETARG (NPARAM+1,FICOS)	!Nom du fichier OS binaire.
	     DEF_FICOS=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileUp") THEN
	     CALL GETARG (NPARAM+1,FICUP)	!Nom du fichier de sortie pour le champ UP.
	     DEF_FICUP=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileDown") THEN
	     CALL GETARG (NPARAM+1,FICDOWN)	!Nom du fichier de sortie pour le champ DOWN.
	     DEF_FICDOWN=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.Log") THEN
	     CALL GETARG (NPARAM+1,FICLOG)		!Nom du fichier Trace.
	     DEF_FICLOG=1
	 ENDIF	 
     
	 IF (KEYWORD.EQ."-SOS.AngFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES)		!Nom du fichier des angles (chemin complet)
	     DEF_FICANGLES=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-SOS.ResFileUp.UserAng") THEN
	     CALL GETARG (NPARAM+1,FICUP_ANGLES_USER)	!Nom du fichier de sortie pour le champ UP 
	     DEF_FICUP_ANGLES_USER=1			!aux angles "utilisateur".
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileDown.UserAng") THEN
	     CALL GETARG (NPARAM+1,FICDOWN_ANGLES_USER)	!Nom du fichier de sortie pour le champ DOWN.
	     DEF_FICDOWN_ANGLES_USER=1			!aux angles "utilisateur".
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.Trans") THEN
	     CALL GETARG (NPARAM+1,FICTRANS)	!Nom du fichier de sortie des valeurs de transmission.
	     DEF_FICTRANS=1			
	 ENDIF
         	 	 	
	 IF (KEYWORD.EQ."-SOS.MDF") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Facteur de depolarisarion moleculaire
	     READ(PROG_ARG,*,ERR=800) RON
	     DEF_RON=1
	 ENDIF 	
		
	 IF (KEYWORD.EQ."-SURF.Alb") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Albedo de surface
	     READ(PROG_ARG,*,ERR=801) RHO
	     DEF_RHO=1
	 ENDIF 	 	 
		
	 IF (KEYWORD.EQ."-SURF.Type") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Type de surface
	     READ(PROG_ARG,*,ERR=803) ISURF
	     DEF_ISURF=1
	 ENDIF 	 	 
		
	 IF (KEYWORD.EQ."-SURF.Ind") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Indice de refraction de la surface / air
	     READ(PROG_ARG,*,ERR=802) IND_SURF
	     DEF_IND_SURF=1
	 ENDIF 	 
	 	 
	 IF (KEYWORD.EQ."-SURF.File") THEN
	     CALL GETARG (NPARAM+1,FICSURF)	!Nom du fichier de Matrice de reflexion
	     DEF_FICSURF=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Dir") THEN
	     CALL GETARG (NPARAM+1,REPSURF)	!Repertoire de stockage des Matrices de BRDF/BPDF
	     DEF_REPSURF=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-SURF.Glitter.Wind") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Vitesse de vent
	     READ(PROG_ARG,*,ERR=804) WIND 
	     DEF_WIND=1
	 ENDIF

      	 IF (KEYWORD.EQ."-SURF.Roujean.K0") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Parametre k0 du modele de BRDF de Roujean
             READ(PROG_ARG,*,ERR=805) K0_ROUJEAN    
	     DEF_K0_ROUJEAN=1
	 ENDIF	
	 
	 IF (KEYWORD.EQ."-SURF.Roujean.K1") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Parametre k1 du modele de BRDF de Roujean
             READ(PROG_ARG,*,ERR=806) K1_ROUJEAN
	     DEF_K1_ROUJEAN=1
	 ENDIF	
	 
	 IF (KEYWORD.EQ."-SURF.Roujean.K2") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Parametre k2 du modele de BRDF de Roujean
             READ(PROG_ARG,*,ERR=807) K2_ROUJEAN 
	     DEF_K2_ROUJEAN=1
	 ENDIF	

	  	 
	 IF (KEYWORD.EQ."-SURF.Nadal.Alpha") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Parametre alpha du modele de BPDF de Nadal
             READ(PROG_ARG,*,ERR=808) ALPHA_NADAL  
	     DEF_ALPHA_NADAL=1
	 ENDIF	 

         IF (KEYWORD.EQ."-SURF.Nadal.Beta") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Parametre beta du modele de BPDF de Nadal
             READ(PROG_ARG,*,ERR=809) BETA_NADAL
             DEF_BETA_NADAL=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.OutputLevel") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Niveau de sortie
             READ(PROG_ARG,*,ERR=812) IOUT  
             DEF_IOUT=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.IGmax") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Ordre maximal d'interaction
             READ(PROG_ARG,*,ERR=813) IGMAX
             DEF_IGMAX=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.Ipolar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Option de coupure de la polarisation
             READ(PROG_ARG,*,ERR=814) IPOLAR
             DEF_IPOLAR=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.View") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Option de visee
             READ(PROG_ARG,*,ERR=815) ITRPHI
             DEF_ITRPHI=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.View.Phi") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Angle azimutal relatif 
             READ(PROG_ARG,*,ERR=816) PHIOS   
             DEF_PHIOS=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.View.Dphi") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)    		!Pas en angles azimutaux
             READ(PROG_ARG,*,ERR=817) PAS_PHI     
             DEF_PAS_PHI=1
         ENDIF	 
        	 	 	 	 	 	 	 	 	 		  
      ENDDO	!Fin de boucle de lecture des arguments



C --- Gestion de la coherence des arguments 

C     Argument FICANGLES obligatoire       
      IF (DEF_FICANGLES.EQ.0) GOTO 500   
	  
C     Fichier OS binaire : obligatoire           
      IF (DEF_FICOS.EQ.0) GOTO 600

C     Fichier de sortie pour le champ UP : obligatoire           
      IF (DEF_FICUP.EQ.0) GOTO 601
      
C     Fichier de sortie pour le champ DOWN : obligatoire           
      IF (DEF_FICDOWN.EQ.0) GOTO 602
      
C     Argument FICLOG obligatoire          
      IF (DEF_FICLOG.EQ.0) THEN
         GOTO 603
      ELSE !Determination de la requete de suivi des traitements
         IF ((FICLOG(1:1).EQ.'0').AND.(FICLOG(2:2).EQ.' ')) THEN
	    TRACE=.FALSE.
         ELSE
            TRACE=.TRUE.
         ENDIF
      ENDIF	 
                        
C     Fichier de granulometrie aerosol : obligatoire     
      IF (DEF_FICGRANU.EQ.0) GOTO 700
      
C     Fichier de profil d'epaisseurs optiques : obligatoire           
      IF (DEF_FICPROFIL.EQ.0) GOTO 701
           
C     Argument RON (facteur de depolarisation moleculaire) obligatoire, controle de valeur        
      IF (DEF_RON.EQ.0) THEN
          GOTO 702
      ELSE
          !Controle elementaire de la valeur du facteur de depolarisation moleculaire
          IF (RON.LT.0.D-00) GOTO 1006
      ENDIF
      
C     Argument RHO (albedo de surface) obligatoire, controle de sa valeur         
      IF (DEF_RHO.EQ.0) THEN
          GOTO 703
      ELSE
          !Controle elementaire de la valeur de l'albedo de surface
          IF (RHO.LT.0.D-00) GOTO 1002
      ENDIF

C     Argument ISURF (type de surface) obligatoire , controle de sa valeur        
      IF (DEF_ISURF.EQ.0) THEN
          GOTO 704
      ELSE
      
          !Controle du domaine couvert par ISURF
	  IF ((ISURF.LT.0).OR.(ISURF.GT.6)) GOTO 1003
	  
          !Argument IND_SURF (indice air / surface) a definir si ISURF = 1, 2, 4, 5 ou 6
	  IF (DEF_IND_SURF.EQ.0) THEN
	  
	      IF ((ISURF.EQ.0).OR.(ISURF.EQ.3)) THEN
	         !Pour le cas d'une reflexion lambertienne ou de Roujean : 
	         !on n'a pas besoin de definir IND_SURF : non utilise par le code SOS, 
		 !mais l'argument est attendu en entree de SOS_OS et SOS_TRPHI_OPTION
	         ! ==> on definit une valeur quelleconque
	         IND_SURF=-9999.0
	      ELSE
	         !Cas ISURF = 1, 2, 4, 5 ou 6 : l'indice IND_SURF doit etre defini
		 IF ( (ISURF.EQ.1).OR.(ISURF.EQ.2).OR.(ISURF.EQ.4)
     &                .OR.(ISURF.EQ.5).OR.(ISURF.EQ.6) ) GOTO 705
	     ENDIF
	     
	  ELSE !Cas IND_SURF est defini
	       !Controle elementaire de la valeur de l'indice de surface 
	       !si sa definition est requise (Cas ISURF = 1, 2, 4, 5 ou 6)
	       IF ( (ISURF.EQ.1).OR.(ISURF.EQ.2).OR.(ISURF.EQ.4)
     &                .OR.(ISURF.EQ.5).OR.(ISURF.EQ.6) ) THEN
  	            IF (IND_SURF.LT.1.D-00) GOTO 1005
               ENDIF
	  ENDIF !DEF_IND_SURF
	      
	  
	  !Argument FICSURF (fichier BRDF/BPDF surface) a definir si ISURF = 1, 3, 4, 5 ou 6	  
          IF (DEF_FICSURF.EQ.0) THEN 
	     
             IF ((ISURF.EQ.0).OR.(ISURF.EQ.2)) THEN
	         !Pour le cas d'une reflexion lambertienne ou de Fresnel : 
	         !on n'a pas besoin de definir FICSURF car il n'y a pas de fichier 
		 !BRDF/BPDF utilise, mais l'argument est attendu en entree de SOS_OS
	         ! ==> on definit une chaine de caracteres quelleconque
	         FICSURF="NOT_USED"
	     ELSE
	         !Cas ISURF = 1, 3, 4, 5 ou 6 : le fichier FICSURF doit etre defini
		 IF ( (ISURF.EQ.1).OR.(ISURF.EQ.3).OR.(ISURF.EQ.4)
     &                .OR.(ISURF.EQ.5).OR.(ISURF.EQ.6) ) GOTO 706
	     ENDIF
	  
	  ELSE
	     !Si FICSURF est defini et vaut 'DEFAULT', il faut aussi que REPSURF 
	     !(repertoire de stockage des fichiers BRDF/BPDF) soit egalement defini
	     !si ISURF = 1, 3, 4, 5 ou 6
	     IF ((FICSURF.EQ.'DEFAULT').AND.((ISURF.EQ.1).OR.(ISURF.EQ.3).
     &            OR.(ISURF.EQ.4).OR.(ISURF.EQ.5).OR.(ISURF.EQ.6))) THEN
		 IF (DEF_REPSURF.EQ.0) GOTO 707
	     ENDIF 
          ENDIF !DEF_FICSURF

	  !Controle WIND
	  IF (ISURF.EQ.1) THEN
	      IF (DEF_WIND.EQ.0) THEN
	         GOTO 708
	      ELSE !WIND est defini : controle elementaire de sa valeur
	         IF (WIND.LT.0.D-00) GOTO 1004   
	      ENDIF	 
	  ELSE
	      !Hormis pour le cas du glitter on n'a pas besoin de definir WIND : 
	      !mais l'argument est attendu en entree de SOS_TRPHI_OPTION
	      ! ==> on definit une valeur  quelleconque
	      WIND = -9999.0
	  ENDIF
	  
	  !Controle passage des arguments du modele de Roujean
	  IF (ISURF.GE.3) THEN
	     IF (DEF_K0_ROUJEAN.EQ.0) GOTO 709
	     IF (DEF_K1_ROUJEAN.EQ.0) GOTO 709
	     IF (DEF_K2_ROUJEAN.EQ.0) GOTO 709
	     
	     ! Controle de coherence par rapport au format F7.3
             ! Valeurs de K0, K1,K2 forcees au format F7.3
             ! -----------------------------------------	
             IF((K0_ROUJEAN.LT.0.).OR.(K0_ROUJEAN.GE.1000.)) GOTO 1013
	     IF((K1_ROUJEAN.LT.0.).OR.(K1_ROUJEAN.GE.1000.)) GOTO 1013
	     IF((K2_ROUJEAN.LT.0.).OR.(K2_ROUJEAN.GE.1000.)) GOTO 1013
	
	     K0_ROUJEAN = DNINT(K0_ROUJEAN*1000.D+00) / 1000.D+00
	     K1_ROUJEAN = DNINT(K1_ROUJEAN*1000.D+00) / 1000.D+00
	     K2_ROUJEAN = DNINT(K2_ROUJEAN*1000.D+00) / 1000.D+00	 
	  ELSE
	     !Hormis pour les cas ISURF >= 3 on n'a pas besoin de definir  
	     !K*_ROUJEAN mais l'argument est attendu en entree de SOS_TRPHI_OPTION
	     ! ==> on definit une valeur  quelleconque
	     K0_ROUJEAN = -9999.0
	     K1_ROUJEAN = -9999.0
	     K2_ROUJEAN = -9999.0
	  ENDIF
	  
	  !Controle passage des arguments du modele de Nadal
          IF (ISURF.EQ.6) THEN	
             IF (DEF_ALPHA_NADAL.EQ.0) GOTO 710
	     IF (DEF_BETA_NADAL.EQ.0)  GOTO 710
	     
	     ! Controle de coherence par rapport aux formats F6.4 et F4.1
	     ! Valeurs de alpha et beta forcees a ces formats
	     ! -----------------------------------------	
             IF((ALPHA_NADAL.LT.0.).OR.(ALPHA_NADAL.GE.10.)) GOTO 1014
	     IF((BETA_NADAL.LT.0.).OR.(BETA_NADAL.GE.100.)) GOTO 1014

	     ALPHA_NADAL = DNINT(ALPHA_NADAL*10000.D+00) / 10000.D+00
	     BETA_NADAL = DNINT(BETA_NADAL*10.D+00) / 10.D+00		
          ELSE
	     !Hormis pour les cas ISURF = 6 on n'a pas besoin de definir  
	     !ALPHA_NADAL et BETA_NADAL mais les arguments sont attendus
	     !en entree de SOS_TRPHI_OPTION
	     ! ==> on definit des valeurs  quelleconque
	     ALPHA_NADAL = -9999.0
	     BETA_NADAL = -9999.0	       
	  ENDIF 	  
	 	 
      ENDIF !DEF_ISURF


C     Argument TETAS obligatoire    

C     Niveau de sortie des simulations de luminance : optionnel  
C     Si non defini : la sortie est standard : 
C        Champ emergeant au TOA, champ descendant au SOL : IOUT = -1     
      IF (DEF_IOUT.EQ.0) IOUT=-1
      IF ((IOUT.LT.-1).OR.(IOUT.GT.SOS_OS_NT)) GOTO 1007
             
C     Option de coupure de la polarisation : optionnelle
C     Si non defini : la simulation prend en compte la polarisation : IPOLAR = 1        
      IF (DEF_IPOLAR.EQ.0) IPOLAR=1  
      IF ((IPOLAR.LT.0).OR.(IPOLAR.GT.1)) GOTO 1009            
      
C     Ordre maximal d'interaction : obligatoire        
      IF (DEF_IGMAX.EQ.0) THEN 
          GOTO 730
      ELSE !Controle de validite de la valeur du parametre
         IF (IGMAX.LT.1) GOTO 1008
      ENDIF
      
C     Option de visee : obligatoire  
      IF (DEF_ITRPHI.EQ.0) THEN
          GOTO 740
      ELSE
          IF ((ITRPHI.LT.1).OR.(ITRPHI.GT.2)) GOTO 1010
	  
          !Si choix ITRPHI=1 : il faut definir l'azimut relatif
          IF (ITRPHI.EQ.1) THEN
	      IF (DEF_PHIOS.EQ.0) THEN
	         GOTO 741
	      ELSE
	        !Pour le cas d'une definition du plan azimutal d'observation 
	        !on n'a pas besoin de definir PAS_PHI car il n'est pas utilise par le code SOS, 
	        !mais l'argument est attendu en entree SOS_TRPHI_OPTION
	        ! ==> on definit une valeur quelleconque (entiere)
	        PAS_PHI=-9999
	     ENDIF
	  ENDIF
	  
	  !Si choix ITRPHI=2 : il faut definir le pas en azimut
	  IF (ITRPHI.EQ.2) THEN
	     IF (DEF_PAS_PHI.EQ.0) THEN
	        GOTO 742
	     ELSE
	        !Pour le cas d'une definition du champ d'observation sur pas en azimut 
	        !on n'a pas besoin de definir PHIOS car il n'est pas utilise par le code SOS, 
	        !mais l'argument est attendu en entree SOS_TRPHI_OPTION
	        ! ==> on definit une valeur quelleconque
	        PHIOS=-9999.0
		
		!Controle elementaire du pas en phi 
		IF (PAS_PHI.LE.0) GOTO 1011
	     ENDIF
	  ENDIF
      ENDIF !DEF_ITRPHI
            

C     Argument optionnels
C     --> Sorties pour un fichier d'angles "utilisateur"            
      IF (DEF_FICUP_ANGLES_USER.EQ.0)   FICUP_ANGLES_USER="NO_OUTPUT"	           
      IF (DEF_FICDOWN_ANGLES_USER.EQ.0) FICDOWN_ANGLES_USER="NO_OUTPUT"

C     --> Sortie d'un fichier des transmissions     
      IF (DEF_FICTRANS.EQ.0) FICTRANS="NO_OUTPUT"
      
	 
C* Ouverture du fichier Trace 
C------------------------------------------------------
      IF (TRACE) THEN
         OPEN (99,FILE=FICLOG,ERR=911)
      ENDIF

C* Initialisation
C-----------------
      IER=0
      IGLI=0
      IFRESNEL=0
      IROUJEAN=0
      IRONDEAUX=0
      IBREON=0
      INADAL=0
      IMAT_SURF=0
      IF (ISURF.EQ.1) IGLI=1
      IF (ISURF.EQ.2) IFRESNEL=1
      IF (ISURF.GE.3) IROUJEAN=1
      IF (ISURF.EQ.4) IRONDEAUX=1
      IF (ISURF.EQ.5) IBREON=1
      IF (ISURF.EQ.6) INADAL=1
      IF ((IGLI.EQ.1).OR.(IROUJEAN.EQ.1).OR.(IRONDEAUX.EQ.1).
     &     OR.(IBREON.EQ.1).OR.(INADAL.EQ.1)) IMAT_SURF=1
     
      LTA=.TRUE.			!epaisseur optique des aerosols nulle
      VA=0.
      VATR=0.
      VR=0.




C* Lecture des tables d'angles et de poids utiles aux simulations
C* -----------------------------------------------------------------     
      OPEN(UNIT=10,FILE=FICANGLES,STATUS='OLD',ERR=950)
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) LUM_NBMU
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) NBMU_GAUSS
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      !On ne conserve que la partie après le séparateur ':'
      IFIN=INDEX(LIGNE_TEXTE,':')	
      LIGNE_TEXTE=LIGNE_TEXTE(IFIN+1:LENLINE)	
      !On recherche la fin de ligne utile ' ' pour réduire la taille de la chaine
      IFIN=INDEX(LIGNE_TEXTE,' ')
      FICANGLES_USER=LIGNE_TEXTE(1:IFIN)

      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) TETAS

      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) N0

      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NB

      READ(10,*,ERR=951) CVAL
      READ(10,*,ERR=951) CVAL
      READ(10,*,ERR=951) CVAL

      DO J=1,LUM_NBMU
         READ(10,*,ERR=951) NVAL, RMU(J), GA(J), IND_ANGOUT(J)
      ENDDO 


      DO J=1,LUM_NBMU
         RMU(-J) = -RMU(J)
	 GA(-J) =  GA(J)
	 IND_ANGOUT(-J) = IND_ANGOUT(J)
      ENDDO 
      
      !Controle d'applicabilite
      IF (ABS(N0).GT.LUM_NBMU) GOTO 1000
      IF ((TETAS.LT.0.D-00).OR.(TETAS.GE.90.D-00)) GOTO 1001
      IF ((ISURF.NE.0).AND.(N0.LT.0)) GOTO 1012
      
      
      IF (FICANGLES_USER.EQ.'NO_USER_ANGLES') THEN
	  FICUP_ANGLES_USER="NO_OUTPUT"
	  FICDOWN_ANGLES_USER="NO_OUTPUT"
      ENDIF
      
       
      IF (TRACE) THEN
         WRITE (99,*,ERR=960) "ANGLES UTILISES"
	 WRITE (99,*,ERR=960) "****************"
	 WRITE (99,*,ERR=960) " "
	 WRITE (99,*,ERR=960) "Nombre d'angles au total : ", LUM_NBMU
	 WRITE (99,*,ERR=960) "Nombre d'angles de Gauss : ", NBMU_GAUSS
	 WRITE (99,*,ERR=960) "incidence solaire : "
	 WRITE (99,*,ERR=960) "  --> Angle zenithal solaire : ", TETAS
	 WRITE (99,*,ERR=960) "  --> Indice dans la table d'angles : ", N0
	 
	 WRITE (99,*,ERR=960) " "
         IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	     WRITE(99,*,ERR=960)  
     &	       "Usage du fichier d'angles complementaires :",
     &         FICANGLES_USER
     	     WRITE(99,*,ERR=960) "Nb d'angles complementaires :",
     &         (LUM_NBMU-NBMU_GAUSS)
	 ELSE
	     WRITE(99,*,ERR=960)  
     &	       "Pas de fichier d'angles complementaires"
             IF (DEF_FICUP_ANGLES_USER.EQ.1) THEN
	         WRITE(99,*,ERR=960) 
     &		 'Pas de creation du fic. Up pour angles utilisateur'
	     ENDIF           
             IF (DEF_FICDOWN_ANGLES_USER.EQ.1) THEN
	         WRITE(99,*,ERR=960)  
     &		 'Pas de creation du fic. Down pour angles utilisateur'
	     ENDIF
	 ENDIF !Fin test Pas de fichier d'angles utilisateur
         WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,ERR=960) "Ordre max NB ", OS_NB
	 WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,err=960) "Liste des angles utilises :"
	 WRITE (99,*,err=960)"   Indice, mu, angle (deg.), poids, user"
	 DO 1200 J=-LUM_NBMU,LUM_NBMU
	    IF (J.EQ.0) GOTO 1200
	    IF (J.LT.0) THEN 
	        WRITE(99,223,ERR=960) J, RMU(J),
     &	              -DACOS(-RMU(J))*180./SOS_PI, GA(J),
     &                 IND_ANGOUT(J)
            ELSE
	        WRITE(99,223,ERR=960) J, RMU(J),
     &	               DACOS(RMU(J))*180./SOS_PI, GA(J), 
     &      	       IND_ANGOUT(J) 
	    ENDIF
1200     CONTINUE	
	 WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,ERR=960) " "  	
      ENDIF !Fichier Trace	 



C* Determination du fichier de BRDF ou BPDF de surface
C  ---------------------------------------------------

      IF (IMAT_SURF.EQ.1) THEN

C* Restitution du nom du fichier de Reflexion de surface
C  ----------------------------------------------------------	
         IF (FICSURF.EQ.'DEFAULT') THEN

	     IF (ISURF.EQ.1) THEN
	         CALL SOS_NOM_FICGLITTER(NBMU_GAUSS,FICANGLES_USER,
     &		                         IND_SURF,WIND,FICS)	        
	     ENDIF
	
	     IF (ISURF.EQ.3) THEN
	         CALL SOS_NOM_FICROUJ(NBMU_GAUSS,FICANGLES_USER,
     &                                K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                                FICS)    	        
	     ENDIF
	     
	     IF (ISURF.EQ.4) THEN
	         CALL SOS_NOM_FICRH(NBMU_GAUSS,FICANGLES_USER,IND_SURF,
     &	                            K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                              FICS)
	     ENDIF
	     
	     IF (ISURF.EQ.5) THEN
	         CALL SOS_NOM_FICBREON(
     &				NBMU_GAUSS,FICANGLES_USER,IND_SURF,
     &                          K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                          FICS)
	     ENDIF
	     
	     IF (ISURF.EQ.6) THEN
	         CALL SOS_NOM_FICNADAL(
     &                          NBMU_GAUSS,FICANGLES_USER,IND_SURF,
     &                          ALPHA_NADAL,BETA_NADAL,
     &	                        K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                          FICS)
	     ENDIF
	 
	     LDIRSURF = INDEX(REPSURF,' ')
             LDIRSURF = LDIRSURF - 1
             IF(LDIRSURF.LE.0) LDIRSURF = SOS_LENDIR
             FICSURF = REPSURF(1:LDIRSURF)//'/'//FICS
	 
	 ENDIF        
		
C* Controle de l'existence du fichier de Matrice de reflexion
C  ----------------------------------------------------------	
         LFICSURF=INDEX(FICSURF,' ')-1 
         INQUIRE(FILE=FICSURF(1:LFICSURF),EXIST=EX)	 	
         IF (EX.EQV..FALSE.) GOTO 1020   !Cas d'erreur  
	       
      ENDIF ! Fin du test pour l'utilisation du fichier de Matrice de reflexion



C* Lecture du fichier AEROSOLS
C  ---------------------------
      OPEN(1,FILE=FICGRANU,STATUS='OLD',ERR=913)
 
     
      ! Section efficace d'extinction
      READ(1,'(a)',ERR=923) LIGNE_TEXTE

      ! Section efficace de diffusion
      READ(1,'(a)',ERR=923) LIGNE_TEXTE

      ! Facteur d'asymétrie (non tronqué)
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
          
      ! Coefficient de troncature  
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) A
	 
      ! Albédo de simple diffusion (éventuellement tronqué) 
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*)  PIZTR

      !Lecture des informations de contenu des colonnes
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      
      DO K=0,OS_NB
         READ(1,60,ERR=923) ALPHA(K),BETA(K),GAMMA(K),ZETA(K)
      ENDDO
	 
	

      
C* Lecture du fichier PROFIL 
C  et stockage de l'epaisseur optique non tronquee
C  -----------------------------------------------

      OPEN(2,FILE=FICPROFIL,STATUS='OLD',ERR=912)
      
      DO I=0,SOS_OS_NT
         READ(2,70,err=922)J,H(I),XDEL(I),YDEL(I)
      ENDDO

      TTOT_VRAI=H(SOS_OS_NT)
      
C* Calcul de l'albedo de simple diffusion non tronquee 
C* a partir de la valeur tronquee
C  ---------------------------------------------------
      PIZ=PIZTR/(1+0.5*A*(PIZTR-1))
            
C* Adaptation du profil a la troncature si A non nul
C  -------------------------------------------------
      HTR(0)=H(0)
      
      IF(A.NE.0.) THEN
         DO I=1,SOS_OS_NT	
	    VA = XDEL(I)*(H(I)-H(I-1))
	    VATR = VA*(1-PIZ*0.5*A)
	    VR = YDEL(I)*(H(I)-H(I-1))
	   
	    HTR(I) = (VATR+VR) + HTR(I-1)
	    XDEL(I) = VATR/(VATR+VR)
	    YDEL(I) = VR/(VATR+VR)
         ENDDO
      ENDIF

      DO I=0,SOS_OS_NT
      	 IF(A.NE.0.) H(I)=HTR(I)	! Cas avec troncature
         XDEL(I)=XDEL(I)*PIZTR
      	 IF(XDEL(I).NE.0.) LTA=.FALSE.	!TA non nulle
      ENDDO
      
C* Ordre maximal de la decomposition en series de Fourier
C* (2 pour le rayleigh pur)
C---------------------------------------------------------
      IBORM=OS_NB    
      IF (LTA) IBORM=2
	
C* Appel des OS
C------------------
      CALL SOS_OS(LUM_NBMU,RMU,GA,OS_NB,
     &            FICSURF,FICOS,
     &            N0,TETAS,RHO,IMAT_SURF,IFRESNEL,IND_SURF,
     &            H,XDEL,YDEL,RON,
     &            ALPHA,BETA,GAMMA,ZETA,
     &            IOUT,IGMAX,IBORM,IPOLAR,TRACE,TDIFMUS,IER)
      IF (IER.NE.0) GOTO 995   ! cas d'erreur

C* Definition de l'epaisseur optique au niveau du profil 
C* pour la sortie du champ montant 
C-------------------------------------------------------
      IF (IOUT.EQ.-1) THEN
	       TAUOUT=H(0)		!Toa.
      ELSE
	       TAUOUT=H(IOUT)		!Niveau utilisateur.
      ENDIF

C* Epaisseur optique au niveau du sol (eventuellement tronquee)
      TTOT_TRONC=H(SOS_OS_NT)

C* Recombinaison sur l'azimut de visee
C------------------------------------- 

      IF (TRACE) THEN
          WRITE(99,*) 
	  WRITE(99,*) "************************************************"
	  WRITE(99,*) "Recombinaison du champ de rayonnement par azimut"
	  WRITE(99,*) "************************************************"
	  WRITE(99,*) 
	  IF (ITRPHI.EQ.1) THEN
	     WRITE(99,*) "Recombinaison pour l'azimut (degres) : ",
     &                   PHIOS
	  ELSE
	     WRITE(99,*) "Recombinaison sur les azimuts de 0 a 360 " //
     &                   "degres par pas (degres): ",PAS_PHI
	  ENDIF   
	  WRITE(99,*) 
	  WRITE(99,*) "Fichier Champ Up : ",FICUP
	  WRITE(99,*) "Fichier Champ Down : ",FICDOWN
	  WRITE(99,*) 
	  IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
	     WRITE(99,*) 
	     WRITE(99,*) "Fichier Champ Up pour angles utilisateur: ",
     &	                  FICUP_ANGLES_USER
	  ENDIF
	  IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
	     WRITE(99,*) 
	     WRITE(99,*)"Fichier Champ Down pour angles utilisateur: ",
     &	                 FICDOWN_ANGLES_USER
	  ENDIF
	  WRITE(99,*)
      ENDIF !TRACE
      
      CALL SOS_TRPHI_OPTION(LUM_NBMU,RMU,GA,
     &                      FICOS,TTOT_TRONC,TAUOUT,
     &                      IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &                      IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                      IRONDEAUX,IBREON,
     &                      INADAL,ALPHA_NADAL,BETA_NADAL,
     &                      ITRPHI,PHIOS,PAS_PHI,IPOLAR,
     &			    FICUP,FICDOWN,IND_ANGOUT,
     &                      FICUP_ANGLES_USER,FICDOWN_ANGLES_USER,
     &                      IER)
       IF (IER.NE.0) GOTO 996   ! cas d'erreur  



C* Calcul des transmissions directes et diffuses
C-----------------------------------------------

C     Ecriture des resultats dans un fichier FICOS_TMP detruit en fin de traitement
C     Respect des conditions pour simuler une transmission diffuse :
C        RHO = 0 : sol noir (pas de reflexion lambertienne)
C        IMAT_SURF = 0 : sol noir (pas de reflexion speculaire)
C        IFRESNEL = 0 : sol noir (pas de reflexion speculaire)
C        IOUT = -1 : sortie pour champ down = sol
C        IBORM = 0 : developpement limite a l'ordre IS=0
C        TRACE = .FALSE. : on ne complete pas le fichier TRACE 
C
C     Les transmissions diffuses en sortie de SOS_OS sont fournies pour 
C     l'atmosphere equivalente (suite a la troncature de la fonction de phase
C     des aerosols). Elles sont ramenees a la valeur correspondant a 
C     l'atmosphere vraie.

      IF (FICTRANS.NE.'NO_OUTPUT') THEN
          OPEN(UNIT=20,FILE=FICTRANS,FORM='FORMATTED',ERR=997)

C     --> Transmission directe pour epaisseur optique vraie (non tronquee) 	  	  
  
	  TDIR_TRONC = 
     &     DEXP(-TTOT_TRONC/DCOS(SOS_PI*TETAS/180.D+00))
	  TDIR_VRAI  = 
     &     DEXP(-TTOT_VRAI/DCOS(SOS_PI*TETAS/180.D+00))
	  
	  WRITE(20,*,err=998)'Solar Zenithal Angle  : ',TETAS 
	  WRITE(20,*,err=998)'Direct transmission  TOA -> surface :',
     &	                      TDIR_VRAI
	  WRITE(20,*,err=998) ' '

C     --> Transmission diffuse TOA --> Sol	 
	  FICOS_TMP="FICOS_TRANSM_TMP"
	  	 	  
	  CALL SOS_OS(LUM_NBMU,RMU,GA,OS_NB,FICSURF,FICOS_TMP,
     &                N0,TETAS,0.D+00,0,0,IND_SURF,
     &                H,XDEL,YDEL,RON,ALPHA,BETA,GAMMA,ZETA,
     &                -1,IGMAX,0,IPOLAR,.FALSE.,TDIFMUS,IER)
          IF (IER.NE.0) GOTO 995   ! cas d'erreur	  
	
	  WRITE(20,*,err=998)'Diffuse transmittance : TOA -> surface'
	  WRITE(20,101,err=998) 
     &          TETAS,(TDIFMUS + TDIR_TRONC - TDIR_VRAI)
	  WRITE(20,*,err=998) ' '
	  
C     --> Transmission diffuse Sol --> TOA : 
C         Boucle sur les valeurs IMUS : N0 = J	pour calcul TOA --> Sol (reciprocite)  
          WRITE(20,*,err=998)'Diffuse transmittance : surface -> TOA'
	  DO J=1,LUM_NBMU
	  
	     TDIR_TRONC = DEXP(-TTOT_TRONC/RMU(J))
	     TDIR_VRAI  = DEXP(-TTOT_VRAI/RMU(J))
	  
             CALL SOS_OS(LUM_NBMU,RMU,GA,OS_NB,FICSURF,FICOS_TMP,
     &                   J,TETAS,0.D+00,0,0,IND_SURF,
     &                   H,XDEL,YDEL,RON,ALPHA,BETA,GAMMA,ZETA,
     &                   -1,IGMAX,0,IPOLAR,.FALSE.,TDIFMUS,IER)
             IF (IER.NE.0) GOTO 995   ! cas d'erreur
	      
	     WRITE(20,201,err=998) 
     &             DACOS(RMU(J))*180.D+00/SOS_PI, 
     &	           (TDIFMUS + TDIR_TRONC - TDIR_VRAI) 
          ENDDO	  

C     --> Destruction du fichier temporaire
	  OPEN(UNIT=3,FILE=FICOS_TMP,FORM='UNFORMATTED') 
          CLOSE(3,STATUS='DELETE')    
	  
      ENDIF !Fin test si fichier transmission a creer
      
      
C* Fermeture fichier
C-------------------
      CLOSE(99)
      CLOSE(1)
      CLOSE(2)
      CLOSE(20)
      
C* Format 
C----------
  100 FORMAT(17hNB_TOTAL_ANGLES :,I4)
  110 FORMAT(17hNB_GAUSS_ANGLES :,I4)
  120 FORMAT(17hANGLES_USERFILE :,A)  
  200 FORMAT(16hINTERNAL_OS_NB :,I4)
  410 FORMAT(20hSOLAR ZENITH ANGLE :,F7.3)
  420 FORMAT(15hINTERNAL_IMUS :,I4)
  610 FORMAT(I4,1X,2D21.14,1X,I4)

101   FORMAT(12h    thetas =,1x,f6.3,3x,12htd(thetas) =,1x,f7.4)  
201   FORMAT(12h    thetav =,1x,f6.3,3x,12htd(thetav) =,1x,f7.4)     

C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  500 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.AngFile Value'
      CALL EXIT(1)
      GOTO 9999  
      
  600 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.ResBin Value'
      CALL EXIT(1)
      GOTO 9999  

  601 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.ResFileUp Value'
      CALL EXIT(1)
      GOTO 9999
      
  602 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.ResFileDown Value'
      CALL EXIT(1)
      GOTO 9999      
            
  603 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.Log Value'
      CALL EXIT(1)
      GOTO 9999      

  700 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.AerFile Value'
      CALL EXIT(1)
      GOTO 9999  
      
  701 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.ProfileFile Value'
      CALL EXIT(1)
      GOTO 9999  
      
  702 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.MDF Value (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999  

  703 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Alb Value (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999  

  704 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Type Value (format : I4)'
      CALL EXIT(1)
      GOTO 9999  

  705 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Ind Value (format : F5.3)'
      CALL EXIT(1)
      GOTO 9999  

  706 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.File Value'
      CALL EXIT(1)
      GOTO 9999                                

  707 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Dir Value'
      CALL EXIT(1)
      GOTO 9999      

  708 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Glitter.Wind Value (format : F4)'
      CALL EXIT(1)
      GOTO 9999   
           
  709 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Roujean.K0 Value (format : F7.3)',
     &           '     -SURF.Roujean.K1 Value (format : F7.3)',
     &           '     -SURF.Roujean.K2 Value (format : F7.3)'
      CALL EXIT(1)
      GOTO 9999        
           
  710 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Nadal.Alpha Value (format : F6.4)',
     &           '     -SURF.Nadal.Beta Value (format : F4.1)'
      CALL EXIT(1)
      GOTO 9999                                 
 
  730 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.IGmax Value (format I4)'
      CALL EXIT(1)
      GOTO 9999   
 
  740 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.View Value (format I4)'
      CALL EXIT(1)
      GOTO 9999    
 
  741 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Option "-SOS.View 1" '
      WRITE(6,*) '      requires associated (-Keyword Value): ',
     &           '      -SOS.View.Phi Value (format F9)'
      CALL EXIT(1)
      GOTO 9999    
 
  742 WRITE(6,*) '  ERROR on parameters for SOS'
      WRITE(6,*) '  --> Option "-SOS.View 2" '
      WRITE(6,*) '      requires associated (-Keyword Value): ',
     &           '      -SOS.View.Dphi Value (format I4)'
      CALL EXIT(1)
      GOTO 9999    

  799 WRITE(6,*) '  ERROR on -SOS.NbAngGauss argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for number of Gauss angles (I4).'
      CALL EXIT(1)
      GOTO 9999
      	                                             
  800 WRITE(6,*) '  ERROR on -SOS.MDF argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for molecular depolarisation factor (F9.5).'
      CALL EXIT(1)
      GOTO 9999
      
  801 WRITE(6,*) '  ERROR on -SURF.Alb argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for surface albedo (F9.5).'
      CALL EXIT(1)
      GOTO 9999   
      
  802 WRITE(6,*) '  ERROR on -SURF.Ind argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for atm. / surface refractive index (F5.3).'
      CALL EXIT(1)
      GOTO 9999            
      
  803 WRITE(6,*) '  ERROR on -SURF.Type argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for the option of surface type (I4).'
      CALL EXIT(1)
      GOTO 9999   

  804 WRITE(6,*) '  ERROR on -SURF.Glitter.Wind loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for the wind velocity (I4).'
      CALL EXIT(1)
      GOTO 9999   

  805 WRITE(6,*) '  ERROR on -SURF.Roujean.K0 argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) "      for the Roujean's model parameter K0 (F7.3)."
      CALL EXIT(1)
      GOTO 9999       

  806 WRITE(6,*) '  ERROR on -SURF.Roujean.K1 argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) "      for the Roujean's model parameter K1 (F7.3)."
      CALL EXIT(1)
      GOTO 9999

  807 WRITE(6,*) '  ERROR on -SURF.Roujean.K2 argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) "      for the Roujean's model parameter K2 (F7.3)."
      CALL EXIT(1)
      GOTO 9999

  808 WRITE(6,*) '  ERROR on -SURF.Nadal.Alpha argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) "      for the Nadal's model parameter Alpha (F6.4)."
      CALL EXIT(1)
      GOTO 9999

  809 WRITE(6,*) '  ERROR on -SURF.Nadal.Beta argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) "      for the Nadal's model parameter Beta (F4.1)."
      CALL EXIT(1)
      GOTO 9999                                            

  812 WRITE(6,*) '  ERROR on -SOS.OutputLevel argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for the output level number (I4).'
      CALL EXIT(1)
      GOTO 9999              

  813 WRITE(6,*) '  ERROR on -SOS.Igmax argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for the maximal order of interaction (I4).'
      CALL EXIT(1)
      GOTO 9999                

  814 WRITE(6,*) '  ERROR on -SOS.Ipolar argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for the polarization option (I4).'
      CALL EXIT(1)
      GOTO 9999      

  815 WRITE(6,*) '  ERROR on -SOS.View argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for the viewing representation option (I4).'
      CALL EXIT(1)
      GOTO 9999   

  816 WRITE(6,*) '  ERROR on -SOS.View.Phi argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for the  viewing azimutal angle (F9, °).'
      CALL EXIT(1)
      GOTO 9999   

  817 WRITE(6,*) '  ERROR on -SOS.View.Dphi argument loading for SOS'
      WRITE(6,*) '  --> Check required formats '
      WRITE(6,*) '      for the step on viewing azimutal angle (I4, °)'
      CALL EXIT(1)
      GOTO 9999         

  911 WRITE(6,*) '  ERROR on logfile opening for SOS'
      CALL EXIT(1)  
      GOTO 9999  

  950 WRITE(6,*) '  ERROR on angles file opening for SOS'
      CALL EXIT(1)
      GOTO 9999

  951 WRITE(6,*) '  ERROR on angles file reading for SOS'
      CALL EXIT(1)
      GOTO 9999
      
  960 WRITE(6,*) '  ERROR on logfile writing for SOS'
      IER=-1
      goto 9999
            
  912 WRITE(6,*) '  ERROR on PROFILE file opening for SOS'
      CALL EXIT(1)  
      GOTO 9999 
      
  913 WRITE(6,*) '  ERROR on AEROSOLS file opening for SOS'
      CALL EXIT(1)  
      GOTO 9999    
      
  922 WRITE(6,*) '  ERROR on PROFILE file reading for SOS'
      CALL EXIT(1) 
      GOTO 9999
      
  923 WRITE(6,*) '  ERROR on AEROSOLS file reading for SOS'
      CALL EXIT(1) 
      GOTO 9999
      
  995 WRITE(6,*) '  ERROR on subroutine SOS_OS'
      CALL EXIT(1)
      GOTO 9999
      
  996 WRITE(6,*) '  ERROR on subroutine SOS_TRPHI_OPTION'
      CALL EXIT(1)
      GOTO 9999

  997 WRITE(6,*) '  ERROR on transmission file opening for SOS'
      CALL EXIT(1)  
      GOTO 9999  

  998 WRITE(6,*) '  ERROR on transmission file writing for SOS'
      CALL EXIT(1)  
      GOTO 9999  
              
 1000 WRITE(6,*) '  ERROR on angles file value for IMUS' 
      WRITE(6,*) '  Out of available range.'
      WRITE(6,*) '  --> Defined value : ',N0
      WRITE(6,*) '      Have to be define between 1 and ',LUM_NBMU
      CALL EXIT(1)
      GOTO 9999
            
 1001 WRITE(6,*) '  ERROR on -SOS.Thetas value :' 
      WRITE(6,*) '  Out of available range.'
      CALL EXIT(1)
      GOTO 9999      

 1012 WRITE(6,*) '  ERROR on angles file value for IMUS' 
      WRITE(6,*) '  Out of available range.'
      WRITE(6,*) '  --> -SURF.Type > 0 requires IMUS > 0 '
      WRITE(6,*) '      IMUS corresponding to the number '
      WRITE(6,*) '      of the solar angle in the angles file.'
      CALL EXIT(1)
      GOTO 9999     
            
 1002 WRITE(6,*) '  ERROR on -SURF.Alb value : have to be positive.'
      CALL EXIT(1)
      GOTO 9999
      
 1003 WRITE(6,*) '  ERROR on -SOS.Type : out of available range.'
      CALL EXIT(1)
      GOTO 9999      
      
 1004 WRITE(6,*) '  ERROR on -SURF.Glitter.Wind value.'
      CALL EXIT(1)
      GOTO 9999  
           
 1005 WRITE(6,*) '  ERROR on -SURF.Ind value : have to be > 1'
      CALL EXIT(1)
      GOTO 9999   
         
 1006 WRITE(6,*) '  ERROR on -SOS.MDF value : have to be positive.'
      CALL EXIT(1)
      GOTO 9999
      
 1007 WRITE(6,*) '  ERROR on -SOS.OutputLevel value : '
      WRITE(6,*) '  --> Defined value : ',IOUT
      WRITE(6,*) '      Have to be defined between 0 (surface)'
      WRITE(6,*) '      and maximal profile level : ',SOS_OS_NT
      CALL EXIT(1)
      GOTO 9999    
                 
 1008 WRITE(6,*) '  ERROR on -SOS.Igmax value : '
      WRITE(6,*) '  --> maximal scattering order have to be > 0'
      CALL EXIT(1)
      GOTO 9999 

 1009 WRITE(6,*) '  ERROR on -SOS.Ipolar value : '
      WRITE(6,*) '  --> Simulation including polarization : 1'
      WRITE(6,*) '  --> Simulation neglecting polarization : 0'
      CALL EXIT(1)
      GOTO 9999 

 1010 WRITE(6,*) '  ERROR on -SOS.View value : '
      WRITE(6,*) '  --> Azimutal plane : 1'
      WRITE(6,*) '  --> Polar diagram  : 2'
      CALL EXIT(1)
      GOTO 9999             

 1011 WRITE(6,*) '  ERROR on -SOS.View.Dphi value : '
      WRITE(6,*) '  --> Step for azimut on polar diagram > 0'
      CALL EXIT(1)
      GOTO 9999                                
                       

 1013 WRITE(6,*) '  ERROR on a -SOS.Roujean.K* value : '
      WRITE(6,*) '  --> Inconsistent parameter format.'
      CALL EXIT(1)
      GOTO 9999                        
           
 1014 WRITE(6,*) '  ERROR on a -SOS.Nadal.* value : '
      WRITE(6,*) '  --> Inconsistent parameter format.'
      CALL EXIT(1)
      GOTO 9999      
    
 1020 WRITE(6,*) '  ERROR on -SURF.File value : file does not exist.'
      CALL EXIT(1)
      GOTO 9999         
        
                       
C* Format
C---------   
   60 FORMAT(E15.8,3(1X,E15.8))  
   70 FORMAT(2X,I4,3(F10.6))
  223 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X,I2)
    
 9999 STOP
      END	!FIN DU PROGRAMME SOS




