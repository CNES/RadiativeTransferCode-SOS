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
C* FICHIER: SOS_SURFACE.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module genere les matrices de reflexion pour differents types de
C*	 modeles.  
C*       
C*
C* AUTEUR: 
C*   Fondement algorithmique: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et developpement : CS SI, B.Lafrance, C.Ruffel.
C*
C* DATE: creation 1.0 : 30/04/2002
C*
C* MOD:VERSION:1.1: Mise a niveau des codes du LOA (codes relatifs au Glitter).
C*                  Developpements d'adaptation des codes aux autres modeles.
C* MOD:VERSION:2.0: Correction d'une erreur sur la polarisation de surface
C*                  pour la routine SOS_MAT_REFLEXION.
C* MOD:VERSION:2.1: 07/07/2008 
C*                  Modification du format de lecture du parametre WIND,
C*                  pour correction d'une erreur de décodage si valeur entière.
C*			WIND passe de F4.1 en F4
C* MOD:VERSION:3.0: 17/11/2009 
C*                  * Modification du passage des arguments : 
C*                      - passage par couples (mot clef, valeur). 
C*                      - adaptation de la gestion des cas d'erreurs. 
C*                  * Conversion des messages d'erreurs en anglais. 
C* MOD:VERSION:3.1: 22/12/2009 
C*                  Modification du format de lecture des parametres Alpha et Beta du modele Nadal
C*                  pour correction d'une erreur de décodage si valeur entière.
C*			ALPHA passe de F6.4 en F6
C*			BETA  passe de F4.1 en F4
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
C*   - Suppression des constantes SOS_OS_NB, SOS_OS_NS et SOS_OS_NM (de SOS.h)  
C*     au profit : d'un passage par argument des valeurs OS_NB, OS_NS, OS_NM (lus dans le fichier des angles)
C*
C*   - Modification de la taille des tableaux de decomposition en polynomes de Legendre :
C*     ==> Passage de SOS_OS_NB a SOS_OS_NB_MAX.
C*     ==> Passage de SOS_OS_NS a SOS_OS_NS_MAX.
C*     ==> Passage de SOS_OS_NM a SOS_OS_NM_MAX.
C*     ==> Boucle sur les ordre jusqu'a OS_NB, OS_NS ou OS_NM 
C*
C*   - Routines SOS_MISE_FORMAT & SOS_BPDF_AJOUT_BRDF : Modification de la  
C*     ligne d'enregistrement ou de lecture des tableaux de la  matrice de 
C*     reflexion a l'ordre IS de la decomposition en series de Fourier :
C*     ==> limitation aux bornes utiles -LUM_NBMU:LUM_NBMU
C*       Rem : pas de modification des bornes d'enregistrement des fichiers 
C*             temporaire RES_FRESNEL, SOS_MAT_REFLEX car il ne s'agit pas 
C*             de produits finaux.
C*
C*   - Adaptation de l'appel des routines de denomination des fichiers
C*     (SOS_NOM_FICGLITTER, SOS_NOM_FICROUJ, SOS_NOM_FICRH, SOS_NOM_FICBREON, SOS_NOM_FICNADAL)
C*     aux nouveaux parametres NBMU_GAUSS et FICANGLES_USER
C*
C*   - Adaptation de la routine SOS_BPDF_AJOUT_BRDF au passage du parametre
C*     LUM_NBMU (nb d'angles utiles) et de OS_NB
C*
C*   - Adaptation de l'appel des routines SOS_GLITTER, SOS_ROUJEAN, SOS_RONDEAUX_BREON, 
C*     SOS_NADAL et SOS_BPDF_AJOUT_BRDF aux passages des parametres
C*     LUM_NBMU (nb d'angles utiles), RMU (table des cosinus des angles) 
C*     et CHR (table des poids des angles), ainsi que OS_NB, OS_NS et OS_NM
C*
C* MOD:VERSION:4.1: 24/06/2010 
C*     Les etiquettes de FORMAT en ,X, passent en ,1X, pour eviter des erreurs de compilation.
C*   
C*   
C* MOD:VERSION:4.2: 08/12/2015
C*    - Suppression des formats à la lecture des paramètres d'entrée en faveur 
C*      de READ(PROG_ARG,*,ERR=NNN)
C*      ==> Permet une compilation sous gfortran et évite les mauvaises conversions 
C*          si la valeur d'entrée ne respecte pas strictement le format attendu 
C*          (par exemple, saisie d'un entier N au lieu d'un réel N. 
C*           aboutissant à une valeur interprétée fausse). 
C*
C*    - Ajustement des lignes utilisant SOS_PI pour éviter une troncature 
C*      au décodage de la constante.
C*
C*    - Changement de la manière de lire les fichiers formatés (Angles, ...)
C*      pour être compatible avec le compilateur gfortran.
C*
C*    - Modification de la récupération du nom du fichier optionnel des 
C*      angles utilisateur (FICANGLES_USER) pour compatibilité avec gfortran/linux
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
C    SOS_LENKEYWORD : Taille maximale des Keywords pour le passage des arguments
C    SOS_PI : Valeur de PI
C    SOS_OS_NBMU_MAX : Nombre maximal d'angles positifs.
C    SOS_OS_NB_MAX :  Valeur maximale pour l'ordre limite du developpement en 
C                     polynomes de Legendre.
C    SOS_OS_NS_MAX :  Valeur maximale pour l'ordre limite du developpement en polynomes 
C                     de Legendre des elements de la matrice de Fresnel 
C                     et en series de Fourier des luminances.
C    SOS_OS_NM_MAX :  Valeur maximale pour l'ordre limite  du developpement 
C                     en series de Fourier de la fonction G.
C
C Constantes spécifiques :
C    LENLINE  : Longueur maximale d'une ligne de données lue dans un fichier
C----------------------------------------------------------------------------
#include "SOS.h"
#define LENLINE 1000

C==============================================================================
C PROGRAMME: SOS_SURFACE
C ==========  
C      Ce programme calcule les matrice de reflexion (BRDF et/ou BPDF) 
C      de differents modeles:
C	- Reflexion du glitter,
C	- Modele de BRDF de Roujean pour simuler la reflexion en intensite 
C       de couverts vegetaux,
C	- Modele simplifie de BPDF de Rondeaux et Herman (RH) pour simuler la 
C	reflexion en polarisation des couverts vegetaux,
C	- Modele simplifie de BPDF de Breon, pour simuler la reflexion
C        en polarisation des sols nus,
C	- Modele de BPDF de Nadal dont l'interet est de permettre de simuler 
C	la polarisation de differents types de surface, d'une maniere 
C	simplifiee mais representative.
C
C       Aux 3 derniers modeles de BPDF (Rondeaux, Breon, Nadal), 
C       on ajoute la BRDF de Roujean.
C
C Donnees en entree du programme
C ------------------------------
C     -->Les parametres/arguments du programme
C
C        Les arguments sont passes sous la forme d'une liste de couples (-Keyword Value)
C
C        Definition des mots clef, statut obligatoire / option, format de la valeur associee :
C         -SURF.AngFile  ==> Obligatoire (chemin complet au fichier ANGLES_LUM) 
C             FICANGLES (CHARACTER*SOS_LENFIC2) : nom complet du fichier des angles utilises 
C                                                 pour les calculs de luminance et BRDF
C                 (repertoire + nom fichier + extension)
C
C         -SURF.Type ==> Obligatoire
C             PROG_ARG   (CHARACTER*10)  : Type de surface
C                             		   (format de décodage : libre --> Integer)
C      		        0 :surface lambertienne d'albedo RHO.
C			1 :surface lambertienne + glitter.
C			2 :surface lambertienne + reflexion de Fresnel.
C		        3 :surface lambertienne + BRDF de Roujean.
C			4 :surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
C			5 :surface lambertienne + BRDF de Roujean + BPDF de Breon.
C			6 :surface lambertienne + BRDF de Roujean + BPDF de Nadal.	
C
C         -SURF.Log ==> Obligatoire (0 ou chemin complet au fichier Trace)
C             FICLOG (CHARACTER*SOS_LENFIC2) : fichier trace SURFACE
C
C         -SURF.Dir ==> Obligatoire
C             DIRS   (CHARACTER*SOS_LENDIR)  : repertoire de stockage des fichiers SURFACE
C
C         -SURF.Ind ==> Associe a : -SURF.Type 1 2 4 5 ou 6
C      	      PROG_ARG   (CHARACTER*10)  : Indice de refraction du dioptre.
C                             		   (format de décodage : libre --> Integer)
C
C         -SURF.Glitter.Wind ==> Associe a : -SURF.Type 1
C      	      PROG_ARG   (CHARACTER*10)  : Vitesse du vent m/s.
C                             		   (format de décodage : libre --> Float)
C
C         -SURF.DirBRDFtoBPDFcalculation ==> Associe a : -SURF.Type 3 4 5 ou 6
C             DIRROUJEAN   (CHARACTER*SOS_LENDIR)  : repertoire de stockage des fichiers 
C                                                    de BRDF de Roujean, appliquees aux 
C                                                    modeles de BPDF.
C
C         -SURF.Roujean.K0 ==> Associe a : -SURF.Type 3 4 5 ou 6
C      	      PROG_ARG   (CHARACTER*10)  : Parametre K0 (modele de Roujean)
C                             		   (format de décodage : libre --> Float)
C
C         -SURF.Roujean.K1 ==> Associe a : -SURF.Type 3 4 5 ou 6
C      	      PROG_ARG   (CHARACTER*10)  : Parametre K1 (modele deRoujean)
C                             		   (format de décodage : libre --> Float)
C
C         -SURF.Roujean.K2 ==> Associe a : -SURF.Type 3 4 5 ou 6
C      	      PROG_ARG   (CHARACTER*10)  : Parametre K2 (modele deRoujean)
C                             		   (format de décodage : libre --> Float)
C
C         -SURF.Nadal.Alpha ==> Associe a : -SURF.Type 6
C      	      PROG_ARG   (CHARACTER*10)  : Parametre alpha (modele de Nadal)
C                             		   (format de décodage : libre --> Float)
C         -SURF.Nadal.Beta
C      	      PROG_ARG   (CHARACTER*10)  : Parametre beta (modele de Nadal)
C                             		   (format de décodage : libre --> Float)
C                                 
C
C    -->Fichiers
C       Le fichier Trace.
C       Le fichier resultat contenant l'expression, pour le code des OS, de la 
C       decomposition en series de Fourier de la matrice de reflexion de surface.
C       Des fichiers intermediaires (non conserves) : 
C           commentaires et description dans les routines appelees.
C       Le fichier des angles et parametres internes pour les limites des developpements.
C         
C       
C
C Resultats fournis par le programme
C ----------------------------------
C fichier BRDF / BPDF de la matrice de reflexion de surface
C pour une décomposition en séries de Fourier
C
C Variable d'environnement
C ------------------------
C   Aucune
C
C
C Description des fichiers utilises
C ---------------------------------
C
C   -->Contenu du fichier BRDF / BPDF resultat
C       Fichier binaire non formatte, cree par SOS_MISE_FORMAT
C       Il contient la matrice de reflexion de surface 
C       (fichier GLITTER ou autres BRDF / BPDF):
C       c'est a dire les elements Pij de la matrice de reflexion 
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
C        --> Enregistrement sur le domaine limite aux angles utiles
C                                  --------------------------------
C
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
C               INTERNAL_OS_NM :     Ordre maximal du developpement en series de Fourier  
C                                    de la fonction G des calculs de matrice de reflexion.
C               INDEX COS_ANGLE WEIGHT OUTPUT: lignes donnant le numero d'angle, son cosinus 
C                                              et son poids, s'il est donne en sortie "utilisateur"
C                                              au format I4,1X,2D21.14,1X,I4
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
C     - SOS_NOM_FICGLITTER   |
C     - SOS_NOM_FICROUJ      |
C     - SOS_NOM_FICRH        | -> Fichier SOS_NOM_FIC_SURFACE.F
C     - SOS_NOM_FICBREON     |
C     - SOS_NOM_FICNADAL     |
C     - SOS_GLITTER
C     - SOS_ROUJEAN
C     - SOS_RONDEAUX_BREON
C     - SOS_NADAL
C     - SOS_BPDF_AJOUT_BRDF
C
C============================================================================
      PROGRAM SOS_SURFACE
      
      IMPLICIT NONE
      
C* Definition des variables   
C*-------------------------------------------------------------------------- 

      CHARACTER*10 PROG_ARG		! Arguments du programme
         
      CHARACTER*SOS_LENDIR DIRS	! chemin d'acces au fichier SURFACE resultat .
      CHARACTER*SOS_LENDIR DIRROUJEAN   ! repertoire de stockage des fichiers de BRDF 
                                        ! de Roujean, appliquees aux modeles de BPDF.                                               
      CHARACTER*SOS_LENFIC1 FICS   	! Nom du fichier SURFACE (sans arborescence).
      CHARACTER*SOS_LENFIC2 FICLOG      ! nom du fichier trace SURFACE.
      CHARACTER*SOS_LENFIC2 FICANGLES   ! chemin et nom du fichier des angles effectifs
      CHARACTER*SOS_LENFIC2 FICANGLES_USER  ! chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*SOS_LENFIC2 FICGLITTER ! Nom complet du fichier Glitter.
      CHARACTER*SOS_LENFIC2 FICROUJ    ! Nom complet du fichier du modele de ROUJEAN.
      CHARACTER*SOS_LENFIC2 FICRH      ! Nom cpmplet du fichier du modele de RH.
      CHARACTER*SOS_LENFIC2 FICBREON   ! Nom complet du fichier du modele de BREON.
      CHARACTER*SOS_LENFIC2 FICNADAL   ! Nom complet du fichier du modele de NADAL.
      CHARACTER*SOS_LENFIC2 FICTMP     ! Nom complet du fichier temporaire pour la BPDF.
      CHARACTER*SOS_LENKEYWORD KEYWORD ! Mot clef du passage d'arguments.
      
      CHARACTER*LENLINE LIGNE_TEXTE    ! Ligne de caractères pour décodage de paramètres
      CHARACTER*1000 CVAL	       ! Chaine de caracteres pour la lecture de valeurs inutilisees
      
      LOGICAL EX		! Verifie existence d'un fichier.
      LOGICAL TRACE     	! = vrai,si ecriture dans le fichier trace

      INTEGER*4 NBMU_GAUSS      ! Nombre d'angles de Gauss utilises       
      INTEGER*4 LUM_NBMU        ! Nombre d'angles effectifs des simulations
      INTEGER*4 LDIRS		! longueur de la chaine de caractere DIRS.
      INTEGER*4 LDIRROUJEAN	! longueur de la chaine de caractere DIRROUJEAN
      INTEGER*4 LFICS		! longueur de la chaine de caractere FICS.     
      INTEGER*4 IER     	! code d'erreur =0 si pas d'erreur, =-1 sinon     
      INTEGER*4 J	        ! Indice de tableau                    
      INTEGER*4 OS_NB 		! Ordre des developpements des fonctions de phase 
				! en Polynomes de Legendre.
      INTEGER*4 OS_NS    	! Ordre maximal du developpement en polynomes de Legendre   
				! des elements de la matrice de Fresnel.
      INTEGER*4 OS_NM    	! Ordre maximal du developpement en series de Fourier  
				! de la fonction G des calculs de matrice de reflexion.

      INTEGER*4 IFIN            ! Nb de caractères significatifs d'une chaine. 

      INTEGER*4 NVAL 	        ! Integer pour la lecture de fichier ANGLES
  
      INTEGER*2 ISURF		! Type de surface     
      INTEGER*2 NBARGS		! Nombre d'arguments en entrée
      INTEGER*2 DEF_ISURF	! Indicateur d'affectation d'une valeur a ISURF 
      INTEGER*2 DEF_FICLOG	! Indicateur d'affectation d'une valeur a FICLOG 
      INTEGER*2 DEF_FICANGLES   ! Indicateur d'affectation d'une valeur a FICANGLES
      INTEGER*2 DEF_IND		! Indicateur d'affectation d'une valeur a IND 
      INTEGER*2 DEF_WIND	! Indicateur d'affectation d'une valeur a WIND 
      INTEGER*2 DEF_K0		! Indicateur d'affectation d'une valeur a K0 
      INTEGER*2 DEF_K1		! Indicateur d'affectation d'une valeur a K1
      INTEGER*2 DEF_K2		! Indicateur d'affectation d'une valeur a K2
      INTEGER*2 DEF_ALPHA	! Indicateur d'affectation d'une valeur a ALPHA
      INTEGER*2 DEF_BETA	! Indicateur d'affectation d'une valeur a BETA
      INTEGER*2 DEF_DIRS	! Indicateur d'affectation d'une valeur a DIRS
      INTEGER*2 DEF_DIRROUJEAN	! Indicateur d'affectation d'une valeur a DIRROUJEAN
      


      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
                                             !Cosinus des angles de Gauss (et utilisateur).
      DOUBLE PRECISION CHR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
                                             !Poids des angles de Gauss.				
          
      DOUBLE PRECISION WIND	! Vitesse du vent (m/s).
      
      DOUBLE PRECISION IND	! Indice de refraction du dioptre.	
      
      DOUBLE PRECISION K0,K1,K2	!Parametres du modele de Roujean
      
      DOUBLE PRECISION ALPHA	!Parametre alpha du modele de Nadal.
      DOUBLE PRECISION BETA	!Parametre beta du modele de Nadal.
   
      INTEGER*4 IARGC           !Type de la fonction IARGC
      
      INTEGER*4 NPARAM 		!Numero du parametre en entrée (type INTEGER*4 obligatoire
                                !pour le f77)


C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0
      
C* Lecture des parametres/arguments d'entree et decodage 
C------------------------------------------------------          

C --- Initialisation des indicateurs de definition des parametres (0 pour indefinis)
      DEF_ISURF=0
      DEF_FICLOG=0
      DEF_FICANGLES=0
      DEF_IND=0
      DEF_WIND=0
      DEF_K0=0
      DEF_K1=0
      DEF_K2=0
      DEF_ALPHA=0
      DEF_BETA=0
      DEF_DIRS=0
      DEF_DIRROUJEAN=0
       
C --- Nombre d'arguments en entree
      NBARGS = IARGC()

C --- Boucle de lecture des couples d'arguments       
      DO NPARAM=1,NBARGS, 2
         CALL GETARG(NPARAM,KEYWORD)
	  	  
	 IF (KEYWORD.EQ."-SURF.Type") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
	     READ(PROG_ARG,*,ERR=990) ISURF
	     DEF_ISURF=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Log") THEN
	     CALL GETARG (NPARAM+1,FICLOG)
	     DEF_FICLOG=1
	 ENDIF	  	 

	 IF (KEYWORD.EQ."-SURF.AngFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES)	!Nom du fichier des angles (chemin complet)
	     DEF_FICANGLES=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-SURF.Ind") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
	     READ(PROG_ARG,*,ERR=990) IND
	     DEF_IND=1
	 ENDIF
	  
	 IF (KEYWORD.EQ."-SURF.Glitter.Wind") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
             READ(PROG_ARG,*,ERR=990) WIND    
	     DEF_WIND=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Roujean.K0") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
             READ(PROG_ARG,*,ERR=990) K0    
	     DEF_K0=1
	 ENDIF
	  
	 IF (KEYWORD.EQ."-SURF.Roujean.K1") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
             READ(PROG_ARG,*,ERR=990) K1    
	     DEF_K1=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Roujean.K2") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
             READ(PROG_ARG,*,ERR=990) K2    
	     DEF_K2=1
	 ENDIF
	  	 
	 IF (KEYWORD.EQ."-SURF.Nadal.Alpha") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
             READ(PROG_ARG,*,ERR=990) ALPHA    
	     DEF_ALPHA=1
	 ENDIF	 

         IF (KEYWORD.EQ."-SURF.Nadal.Beta") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
             READ(PROG_ARG,*,ERR=990) BETA	
             DEF_BETA=1
         ENDIF	 

	 IF (KEYWORD.EQ."-SURF.Dir") THEN
	     CALL GETARG (NPARAM+1,DIRS)
	     DEF_DIRS=1
	 ENDIF	  

	 IF (KEYWORD.EQ."-SURF.DirBRDFtoBPDFcalculation") THEN
	     CALL GETARG (NPARAM+1,DIRROUJEAN)
	     DEF_DIRROUJEAN=1
	 ENDIF	  
	  	 	  
      ENDDO	!Fin de boucle de lecture des arguments



C --- Gestion de la coherence des arguments    

C     Argument FICANGLES obligatoire       
      IF (DEF_FICANGLES.EQ.0) GOTO 800
      
C     Argument ISURF obligatoire, controle de sa valeur  
      IF (DEF_ISURF.EQ.0) THEN
          GOTO 801
      ELSE    
         IF((ISURF.NE.1).AND.(ISURF.NE.3).AND.(ISURF.NE.4)
     &       .AND.(ISURF.NE.5).AND.(ISURF.NE.6)) GOTO 9901
      ENDIF 
       
C     Argument FICLOG obligatoire  
      IF (DEF_FICLOG.EQ.0) THEN
          GOTO 802
      ELSE  !Determination de la requete de suivi des traitements  
         IF ((FICLOG(1:1).EQ.'0').AND.(FICLOG(2:2).EQ.' ')) THEN
	    TRACE=.FALSE.
         ELSE
            TRACE=.TRUE.
         ENDIF
      ENDIF
       
      
C     Argument DIRS obligatoire        
      IF (DEF_DIRS.EQ.0) THEN
          GOTO 803
      ELSE    
         LDIRS = INDEX(DIRS,' ')
         LDIRS = LDIRS - 1
         IF(LDIRS.LE.0) LDIRS = SOS_LENDIR
      ENDIF

C     Argument DIRROUJEAN requis si ISURF >= 3      
      IF ( (DEF_DIRROUJEAN.EQ.0).AND.(ISURF.GE.3) ) THEN 
          GOTO 804  
      ELSE      
          LDIRROUJEAN = INDEX(DIRROUJEAN,' ')
          LDIRROUJEAN = LDIRROUJEAN - 1
          IF(LDIRROUJEAN.LE.0) LDIRROUJEAN = SOS_LENDIR
      ENDIF
 
C     Argument IND requis si ISURF = 1, 2, 4, 5 ou 6      
      IF ( (DEF_IND.EQ.0).AND.( (ISURF.EQ.1).OR.(ISURF.EQ.2)
     &      .OR.(ISURF.EQ.4).OR.(ISURF.EQ.5).OR.(ISURF.EQ.6)) ) THEN 
          GOTO 805  
      ENDIF      
      

            
C* Ouverture du fichier Trace 
C------------------------------------------------------
      IF (TRACE) THEN
      	 OPEN (99,FILE=FICLOG,ERR=911)
      ENDIF	 
	 

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
      
      READ(10,*,ERR=951) CVAL 
      READ(10,*,ERR=951) CVAL
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NB
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NS
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NM
      
      READ(10,*,ERR=951) CVAL
      
      DO J=1,LUM_NBMU
         READ(10,*,ERR=951) NVAL, RMU(J), CHR(J), NVAL
      ENDDO 

      DO J=1,LUM_NBMU
         RMU(-J) = -RMU(J)
	 CHR(-J) =  CHR(J)
      ENDDO 
      
      IF (TRACE) THEN
         WRITE (99,*,ERR=960) "ANGLES UTILISES"
	 WRITE (99,*,ERR=960) "****************"
	 WRITE (99,*,ERR=960) " "
	 WRITE(99,*,ERR=960) "Nombre d'angles au total : ", LUM_NBMU
	 WRITE(99,*,ERR=960) "Nombre d'angles de Gauss : ", NBMU_GAUSS
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
	 ENDIF
         WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,ERR=960) "Ordre max NB ", OS_NB
	 WRITE (99,*,ERR=960) "Ordre max NS ", OS_NS
	 WRITE (99,*,ERR=960) "Ordre max NM ", OS_NM
	 WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,err=960) "Liste des angles utilises :"
	 WRITE (99,*,err=960) "   Indice, mu, angle (degres), poids"
	 DO 1000 J=-LUM_NBMU,LUM_NBMU
	    IF (J.EQ.0) GOTO 1000
	    IF (J.LT.0) THEN 
	        WRITE(99,222,ERR=960) J, RMU(J),
     &	              -DACOS(-RMU(J))*180./SOS_PI, 
     &                 CHR(J)
            ELSE
	        WRITE(99,222,ERR=960) J, RMU(J),
     &	               DACOS(RMU(J))*180./SOS_PI, 
     &                 CHR(J)
	    ENDIF
1000     CONTINUE	
	 WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,ERR=960) " "  	
      ENDIF !Fichier Trace	 



C =====================================================================
C* Cas du GLITTER
C =====================================================================

      IF (ISURF.EQ.1) THEN

         IF (DEF_WIND.EQ.0) THEN
             GOTO 810
         ENDIF            

C* Constitution du nom du fichier GLITTER
C  --------------------------------------
	 CALL SOS_NOM_FICGLITTER(NBMU_GAUSS,FICANGLES_USER,IND,WIND,
     &	                         FICS)      
	 FICGLITTER = DIRS(1:LDIRS)//'/'//FICS	
	 LFICS=INDEX(FICGLITTER,' ')-1	  
	 WRITE(6,*) FICGLITTER(1:LFICS)
	 
C* Si le fichier GLITTER existe deja, le programme n'est pas lance
C  ---------------------------------------------------------------     
         INQUIRE(FILE=FICGLITTER(1:LFICS),EXIST=EX)
	  
         IF (EX) THEN  
            WRITE(6,*)'Fichier GLITTER deja calcule'
	  
	 ELSE   
      
C* Calcul du fichier Glitter
C  ------------------------- 
	    CALL SOS_GLITTER(LUM_NBMU,RMU,CHR,WIND,IND,
     &	                     OS_NB,OS_NS,OS_NM,
     &	                     FICGLITTER,TRACE,IER)
            IF (IER.NE.0) GOTO 993   ! cas d'erreur
	       
         ENDIF

      ENDIF ! Cas du Glitter
	 



C =====================================================================
C* Cas BRDF de ROUJEAN
C*   NB : Calculee egalement pour les modeles de BPDF.
C =====================================================================

      IF (ISURF.GE.3) THEN

         IF ( (DEF_K0.EQ.0).OR.(DEF_K1.EQ.0).OR.(DEF_K2.EQ.0) ) THEN
             GOTO 811
         ENDIF  
         

C* Controle de coherence par rapport au format F7.3
C  Valeurs de K0, K1,K2 forcees au format F7.3
C  -----------------------------------------	
         IF((K0.LT.0.).OR.(K0.GE.1000.)) GOTO 9910
	 IF((K1.LT.0.).OR.(K1.GE.1000.)) GOTO 9910
	 IF((K2.LT.0.).OR.(K2.GE.1000.)) GOTO 9910
	
	 K0 = DNINT(K0*1000.D+00) / 1000.D+00
	 K1 = DNINT(K1*1000.D+00) / 1000.D+00
	 K2 = DNINT(K2*1000.D+00) / 1000.D+00	 	 

	
C* Constitution du nom du fichier de ROUJEAN
C  -----------------------------------------
	 CALL SOS_NOM_FICROUJ(NBMU_GAUSS,FICANGLES_USER,K0,K1,K2,FICS) 	        
	 FICROUJ = DIRROUJEAN(1:LDIRROUJEAN)//'/'//FICS 
	 LFICS=INDEX(FICROUJ,' ')-1
	 WRITE(6,*) FICROUJ(1:LFICS)        
	 
C* Si le fichier existe deja, le programme n'est pas lance
C  --------------------------------------------------------	 
         INQUIRE(FILE=FICROUJ(1:LFICS),EXIST=EX)	 
	 
         IF (EX) THEN  
            WRITE(6,*)'Fichier du modele de ROUJEAN deja calcule'
	    
	 ELSE   
      
C* Calcul du fichier de ROUJEAN
C  -----------------------------
	    CALL SOS_ROUJEAN(LUM_NBMU,RMU,CHR,OS_NB,K0,K1,K2,
     &	                     FICROUJ,TRACE,IER)

            IF (IER.NE.0) GOTO 994   ! cas d'erreur
	       
         ENDIF
	 
      ENDIF ! Cas de la BRDF de Roujean



	 
C =====================================================================
C* Cas BPDF de Rondeaux - Herman
C =====================================================================

      IF (ISURF.EQ.4) THEN

C* Constitution du nom du fichier de Rondeaux - Herman
C  ---------------------------------------------------
	 CALL SOS_NOM_FICRH(NBMU_GAUSS,FICANGLES_USER,IND,K0,K1,K2,FICS) 	        
	 FICRH = DIRS(1:LDIRS)//'/'//FICS	 
	 LFICS=INDEX(FICRH,' ')-1
	 WRITE(6,*) FICRH(1:LFICS)	 
	 
C* Si le fichier existe deja, le programme n'est pas lance
C  -------------------------------------------------------	 
         INQUIRE(FILE=FICRH(1:LFICS),EXIST=EX)
	 	 
         IF (EX) THEN  
            WRITE(6,*)'Fichier du modele de RH deja calcule'
	  
	 ELSE   

C* Determination du nom du fichier de BPDF sans complement de BRDF
C* du modele de Roujean.     
            LFICS=INDEX(FICRH,' ') - 1
            IF (LFICS.LE.0) LFICS=SOS_LENDIR
            FICTMP=FICRH(1:LFICS)//'_BPDF'
      
            
C* Calcul du fichier de Rondeaux - Herman
C  --------------------------------------
	    CALL SOS_RONDEAUX_BREON(LUM_NBMU,RMU,CHR,1,IND,
     &	                            OS_NB,OS_NS,OS_NM,
     &	                            FICTMP,TRACE,IER)       
            IF (IER.NE.0) GOTO 995   ! cas d'erreur

		
C* Ajout de la BRDF de Roujean
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP,FICROUJ,LUM_NBMU,OS_NB,
     &	                             FICRH,IER)
	    IF (IER.NE.0) GOTO 998   ! cas d'erreur
	    
C* Destruction du fichier temporaire
C---------------------------------------   
	    OPEN(UNIT=3,FILE=FICTMP,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	    	                   
	        	       
         ENDIF ! Si fichier a creer.

      ENDIF ! Cas de la BPDF de Rondeaux - Herman



C =====================================================================
C* Cas BPDF de BREON
C =====================================================================

      IF (ISURF.EQ.5) THEN

C* Constitution du nom du fichier de BREON
C  ---------------------------------------
	 CALL SOS_NOM_FICBREON(NBMU_GAUSS,FICANGLES_USER,IND,K0,K1,K2,
     &	                       FICS) 	        
	 FICBREON = DIRS(1:LDIRS)//'/'//FICS
	 LFICS=INDEX(FICBREON,' ')-1
	 WRITE(6,*) FICBREON(1:LFICS)
	 
C* Si le fichier existe deja, le programme n'est pas lance
C  --------------------------------------------------------	 
         INQUIRE(FILE=FICBREON(1:LFICS),EXIST=EX)	 
	 
         IF (EX) THEN  
            WRITE(6,*)'Fichier du modele de BREON deja calcule'
	  
	 ELSE   
      
C* Determination du nom du fichier de BPDF sans complement de BRDF
C* du modele de Roujean.  
C  -------------------------------------------------------------------   
            LFICS=INDEX(FICBREON,' ') - 1
            IF (LFICS.LE.0) LFICS=SOS_LENDIR
            FICTMP=FICBREON(1:LFICS)//'_BPDF'
      
            
C* Calcul du fichier de Rondeaux - Herman
C  --------------------------------------
	    CALL SOS_RONDEAUX_BREON(LUM_NBMU,RMU,CHR,2,IND,
     &	                            OS_NB,OS_NS,OS_NM,
     &	                            FICTMP,TRACE,IER)       
            IF (IER.NE.0) GOTO 995   ! cas d'erreur

C* Ajout de la BRDF de Roujean
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP,FICROUJ,LUM_NBMU,OS_NB,
     &                               FICBREON,IER)
	    IF (IER.NE.0) GOTO 998   ! cas d'erreur
	    
C* Destruction du fichier temporaire
C---------------------------------------   
	    OPEN(UNIT=3,FILE=FICTMP,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	       
         ENDIF ! Si fichier a creer.

      ENDIF ! Cas de la BPDF de Breon 



C =====================================================================
C* Cas BPDF de NADAL
C =====================================================================

      IF (ISURF.EQ.6) THEN	

         IF ( (DEF_ALPHA.EQ.0).OR.(DEF_BETA.EQ.0) ) THEN
             GOTO 812
         ENDIF  
		
C* Controle de coherence par rapport aux formats F6.4 et F4.1
C  Valeurs de alpha et beta forcees a ces formats
C  -----------------------------------------	
         IF((ALPHA.LT.0.).OR.(ALPHA.GE.10.)) GOTO 9915
	 IF((BETA.LT.0.).OR.(BETA.GE.100.)) GOTO 9915

	 ALPHA = DNINT(ALPHA*10000.D+00) / 10000.D+00
	 BETA = DNINT(BETA*10.D+00) / 10.D+00	
		
C* Constitution du nom du fichier de NADAL
C  ---------------------------------------
	 CALL SOS_NOM_FICNADAL(NBMU_GAUSS,FICANGLES_USER,IND,
     &	                       ALPHA,BETA,K0,K1,K2,FICS) 	        
	 FICNADAL = DIRS(1:LDIRS)//'/'//FICS	 
	 LFICS=INDEX(FICNADAL,' ')-1
	 WRITE(6,*) FICNADAL(1:LFICS)
	 
C* Si le fichier existe deja, le programme n'est pas lance
C  --------------------------------------------------------	 
         INQUIRE(FILE=FICNADAL(1:LFICS),EXIST=EX)
	 	 
         IF (EX) THEN  
            WRITE(6,*)'Fichier du modele de NADAL deja calcule'
	  
	 ELSE   

C* Determination du nom du fichier de BPDF sans complement de BRDF
C* du modele de Nadal.     
            LFICS=INDEX(FICNADAL,' ') - 1
            IF (LFICS.LE.0) LFICS=SOS_LENDIR
            FICTMP=FICNADAL(1:LFICS)//'_BPDF'
            
C* Calcul du fichier de NADAL
C  ---------------------------
	    CALL SOS_NADAL(LUM_NBMU,RMU,CHR,IND,ALPHA,BETA,
     &	                   OS_NB,OS_NS,OS_NM, FICTMP,TRACE,IER)
            IF (IER.NE.0) GOTO 997   ! cas d'erreur

C* Ajout de la BRDF de Roujean
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP,FICROUJ,LUM_NBMU,OS_NB,
     &                               FICNADAL,IER)
	    IF (IER.NE.0) GOTO 998   ! cas d'erreur
	    
C* Destruction du fichier temporaire
C---------------------------------------   
	    OPEN(UNIT=3,FILE=FICTMP,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	       
         ENDIF ! Si fichier a creer.
	 
      ENDIF ! Cas de la BPDF de Nadal

     



C* Fermeture fichier
C-------------------
      IF(TRACE) CLOSE(99)      
      CLOSE(10) !Fichier des angles
        
C* Fin nominale 
C-------------------
      GOTO 9999
            
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  800 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.AngFile Value'
      CALL EXIT(1)
      GOTO 9999 
      
  801 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Type Value (format : I2)'
      CALL EXIT(1)
      GOTO 9999  
      
  802 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Log Value'
      CALL EXIT(1)
      GOTO 9999       
            
  803 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Dir Value'
      CALL EXIT(1)
      GOTO 9999        
                   
  804 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.DirBRDFtoBPDFcalculation Value'
      CALL EXIT(1)
      GOTO 9999        
                   
  805 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.ind Value'
      CALL EXIT(1)
      GOTO 9999        
             
  810 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Glitter.Wind Value (format : F4)'
      CALL EXIT(1)
      GOTO 9999        
           
  811 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Roujean.K0 Value (format : F7.3)',
     &           '     -SURF.Roujean.K1 Value (format : F7.3)',
     &           '     -SURF.Roujean.K2 Value (format : F7.3)'
      CALL EXIT(1)
      GOTO 9999        
           
  812 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Nadal.Alpha Value (format : F6.4)',
     &           '     -SURF.Nadal.Beta Value (format : F4.1)'
      CALL EXIT(1)
      GOTO 9999        
                    
  990 WRITE(6,*) '  ERROR on argument loading for SOS_SURFACE'
      WRITE(6,*) '  --> Check required formats.'
      CALL EXIT(1)
      GOTO 9999  
      
 9901 WRITE(6,*) '  ERROR on -SURF.Type value : out of available range.'
      CALL EXIT(1)
      GOTO 9999  
      
 9910 WRITE(6,*) '  ERROR on -SURF.Roujean.K* values : '
      WRITE(6,*) '  out of available range : ]0, 1000.]'
      CALL EXIT(1)
      GOTO 9999 
      
 9915 WRITE(6,*) '  ERROR on -SURF.Nadal.Alpha or .Beta values : '
      WRITE(6,*) '  out of available range : 0 < Alpha <= 10',
     &           ' & 0 < Beta <= 100. '
      CALL EXIT(1)
      GOTO 9999      
              
  911 WRITE(6,*) '  ERROR on logfile opening for SOS_SURFACE'
      CALL EXIT(1)
      GOTO 9999   
  
  950 WRITE(6,*) '  ERROR on angles file opening for SOS_SURFACE'
      CALL EXIT(1)
      GOTO 9999

  951 WRITE(6,*) '  ERROR on angles file reading for SOS_SURFACE'
      CALL EXIT(1)
      GOTO 9999
      
  960 WRITE(6,*) '  ERROR on logfile writing for SOS_SURFACE'
      IER=-1
      goto 9999
                     
  993 WRITE(6,*) '  ERROR on subroutine SOS_GLITTER'
      CALL EXIT(1)
      GOTO 9999
      
  994 WRITE(6,*) '  ERROR on subroutine SOS_ROUJEAN'
      CALL EXIT(1)
      GOTO 9999
      
 9951 WRITE(6,*) '  ERROR while deleting temp. file MAT_BPDF_FORMATOS'
      CALL EXIT(1)
      GOTO 9999      
      
  995 WRITE(6,*) '  ERROR on subroutine SOS_RONDEAUX_BREON'
      CALL EXIT(1)
      GOTO 9999
      
  997 WRITE(6,*) '  ERROR on subroutine SOS_NADAL'
      CALL EXIT(1)
      GOTO 9999
      
  998 WRITE(6,*) '  ERROR on subroutine SOS_BPDF_AJOUT_BRDF'
      CALL EXIT(1)
      GOTO 9999

C* Format
C---------   
  222 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X)

    
 9999 STOP
      END 	!FIN DU PROGRAMME SOS_SURFACE
      
      
      


C==============================================================================
C PROCEDURE: SOS_MAT_FRESNEL
C ==========
C      Cette procedure calcule le developpement de la matrice de Fresnel
C      dans la base des polynomes de Legendre.
C			
C      La matrice de Fresnel, exprimee pour un repere lie au plan de reflexion
C      s'ecrit :
C                  | R11(teta)  R12(teta)    0      |
C         F(teta)= | R12(teta)  R11(teta)    0      |
C                  |    0          0      R33(teta) |
C
C         avec :
C           _ teta : angle de reflexion par rapport a la normale au dioptre,
C           _ R11(teta) = 0.5 * (RL*RL + RR*RR),
C           _ R12(teta) = 0.5 * (RL*RL - RR*RR),
C           _ R33(teta) = RL*RR,
C
C         L'indice de refraction relatif du dioptre, entre le milieu de 
C         reflexion et le milieu de transmission, est IND.
C
C         On calcule RL et RR en fonction de IND et teta, par :
C           _ RL : le coefficient de reflexion parallele
C
C               RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C           _ RR : le coefficient de reflexion perpendiculaire
C
C               RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C
C
C      Le developpement donne les parametres ALPHA(k), BETA(k), GAMMA(k) et
C      ZETA(k) associes a la matrice de Fresnel.
C      Ils sont enregistres de k=0 a OS_NS dans le fichier "RES_FRESNEL".
C
C
C      NB: DEMARCHE D'IMPLEMENTATION :
C      *******************************
C          L'implementation du code est adaptee pour un calcul de la matrice 
C          de reflexion (routine SOS_MAT_REFLEXION) n'utilisant que des valeurs
C           positives des cos(Angle de Gauss). Pour y parvenir :
C          - On utilise non pas l'angle de diffusion vraie mais son angle
C            complementaire (2*angle de reflexion) : DEM 1
C            (routine SOS_MAT_FRESNEL).
C          - On utilise non pas le vrai coefficient de reflexion parallele
C            mais son oppose : DEM 2 (routine SOS_MAT_FRESNEL).
C
C          En aval, la routine calculant la matrice de reflexion introduit
C          un coefficient (-1)**N pour la sommation sur les termes de la 
C          series de Fourier pour assurer la coherence de l'ensemble : DEM 3.
C
C
C Description des parametres
C --------------------------
C       LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles.
C	RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double) : (E) Tableau des mu de Gauss 
C                               (angles complementaires des angles de diffusion).
C	CHR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double) : (E) Tableau des poids de Gauss.
C 	IND   (double) :   (E) Indice de refraction de l'eau / air.
C       OS_NS (I4) : (E) Ordre maximal a utiliser pour le developpement en
C                        polynomes de Legendre des elements de la matrice de Fresnel 
C                        et en series de Fourier des calculs de luminance.
C       TRACE (logical)  : (E) = vrai,si ecriture dans le fichier trace
C                               (No d'unite logique du fichier trace =99)
C       IER  (I4)  :  (S)  code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C Autre resultat fourni par la procedure
C -------------------------------------
C      Fichier temporaire RES_FRESNEL
C      Informations dans le fichier de trace (optionnel)
C         (No d'unite logique du fichier trace =99 , 
C          fichier cree/ouvert par le programme appelant)
C
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier temporaire RES_FRESNEL
C      Fichier Ascii, format d'un enregistrement 4E15.8, 
C      cree par SOS_MAT_FRESNEL
C      Fichier cree dans le repertoire de lancement du programme
C      Si un fichier de meme nom existe deja, il est ecrase
C      Nombre d'enregistrements : de k=0 a OS_NS
C      Chaque enregistrement contient les parametres ALPHA(k), BETA(k), 
C      GAMMA(k) et ZETA(k) associes a la matrice de Fresnel
C
C
C Common utilise:
C --------------
C     Aucun
C
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C
C==============================================================================
      SUBROUTINE SOS_MAT_FRESNEL(LUM_NBMU,RMU,CHR,IND,OS_NS,TRACE,IER)

      IMPLICIT NONE

C* Definition des variables                    
C*----------------------------------------------------------------------
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				! Tableau des mu de Gauss 
      				! (angles complementaires des angles de diffusion).
      DOUBLE PRECISION CHR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				! Tableau des poids de Gauss.

      DOUBLE PRECISION IND	! Indice de refraction de l'eau / air.
	
      DOUBLE PRECISION ALPHA(0:SOS_OS_NS_MAX)       ! Parametres alpha(k).
      DOUBLE PRECISION BETA(0:SOS_OS_NS_MAX)	! Parametres beta(k).
      DOUBLE PRECISION GAMMA(0:SOS_OS_NS_MAX)       ! Parametres gamma(k).
      DOUBLE PRECISION ZETA(0:SOS_OS_NS_MAX)	! Parametres zeta(k).

      DOUBLE PRECISION DELTA(0:SOS_OS_NS_MAX)	! Parametres delta(k).

      DOUBLE PRECISION A	!
      DOUBLE PRECISION B	!
      DOUBLE PRECISION C	!

      DOUBLE PRECISION RL	! Coefficient de reflexion parallele.
      DOUBLE PRECISION RR	! Coefficient de reflexion perpendiculaire.

      DOUBLE PRECISION R11(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
                                ! Element R11 de la matrice de Fresnel
      				! en fonction de l'angle de reflexion.
      DOUBLE PRECISION R12(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)     
      				! Element R12 de la matrice de Fresnel.
      DOUBLE PRECISION R33(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)    
      				! Element R33 de la matrice de Fresnel

      DOUBLE PRECISION IT	! Valeur de R11 recalculee a partir de sa
      				! decomposition en BETA(k).
      DOUBLE PRECISION QT	! Valeur de R12 recalculee a partir de sa
      				! decomposition en GAMMA(k).	
      DOUBLE PRECISION UT	! Valeur de R33 recalculee a partir de sa
      				! decomposition en DELTA(k).			

      DOUBLE PRECISION PL(-1:SOS_OS_NS_MAX+1)
                      !Polynomes de Legendre a l'ordre L calcule pour un angle
      DOUBLE PRECISION POL(0:SOS_OS_NS_MAX+1)
                      !Fonction de Legendre d'ordre S=2 et d'ordre L variable
		      !(calcule pour un angle)
		       			
      DOUBLE PRECISION D	!Variable pour le calcul de POL.
      DOUBLE PRECISION E	!Variable pour le calcul de POL.

      DOUBLE PRECISION X
      DOUBLE PRECISION XRMU	!Valeur utilisee du tableau des angles de Gauss

      DOUBLE PRECISION X2				
      DOUBLE PRECISION XX	
      DOUBLE PRECISION XXX

      DOUBLE PRECISION CO1     !Variable pour le calcul des ALPHA(k) et ZETA(k)
      DOUBLE PRECISION CO2     !Variable pour le calcul des ALPHA(k) et ZETA(k)
      DOUBLE PRECISION CO3     !Variable pour le calcul des ALPHA(k) et ZETA(k)

      DOUBLE PRECISION SOM1    !Variable pour le calcul des ALPHA(k) et ZETA(k)
      DOUBLE PRECISION SOM2    !Variable pour le calcul des ALPHA(k) et ZETA(k)
      DOUBLE PRECISION SOM3    !Variable pour le calcul des ALPHA(k) et ZETA(k)
      DOUBLE PRECISION SOM4    !Variable pour le calcul des ALPHA(k) et ZETA(k)

      LOGICAL TRACE        ! = vrai,si ecriture dans le fichier trace

      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations

      INTEGER*4 OS_NS   ! Ordre maximal du developpement en polynomes de Legendre   
			! des elements de la matrice de Fresnel.
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon  
      INTEGER*4 J	! Indice de tableau.
      INTEGER*4 K	! Indice de tableau.
      INTEGER*4	NN
      INTEGER*4	MM
      INTEGER*4 I	

C* Ouverture du fichier resultat
C---------------------------------
      OPEN(UNIT=3,FILE='RES_FRESNEL',ERR=991)

C* Initilisation
C---------------------------------	
      DO 101 K=0,OS_NS
         BETA(K)=0.D+00
         GAMMA(K)=0.D+00
         ALPHA(K)=0.D+00
         ZETA(K)=0.D+00
         DELTA(K)=0.D+00
  101 CONTINUE
  
C* Fichier trace 
C-----------------
      IF (TRACE) THEN
         WRITE(99,*,err=921)
	 WRITE(99,*,err=921) 
     &   ' Calcul des elements de la matrice de reflexion de Fresnel'
      END IF
			
C* Calcul des elements de la matrice de reflexion de Fresnel	
C--------------------------------------------------------------
      DO 104 J=-LUM_NBMU,LUM_NBMU
         IF(J.EQ.0) GOTO 104

C* Cosinus de l'angle complementaire de l'angle de diffusion :
C* angle de diffusion = PI - 2* angle de reflexion.
C-------------------------------------------------------------	 	
         C=RMU(J)				! DEM 1
	
C* Cosinus de l'angle de reflexion :	
C------------------------------------------------------
         C=DSQRT(.5*(1+C))


C* Coefficients de reflexion parallele (oppose a la valeur vraie)
C* et perpendiculaire
C-----------------------------------------------------------------	
         A=DSQRT(IND*IND-1.D+00+C*C)
	 B=IND*IND*C
         RL=-(B-A)/(B+A)			! DEM 2  	
         RR=(C-A)/(C+A)

C* Elements de la matrice de Fresnel
C-----------------------------------	  	   		
         R11(J)=.5*(RL*RL+RR*RR)
         R12(J)=.5*(RL*RL-RR*RR)
         R33(J)=RL*RR
	 
C* Fichier trace 
C-----------------
         IF (TRACE) THEN
	   WRITE(99,200,err=921) 
     &           ' Angle reflexion ',
     &           DACOS(C)*90./SOS_PI,
     &           ' R11 ',R11(J),' R12 ',R12(J),' R33 ',R33(J)
         END IF
	 
 104  CONTINUE


C* Calcul des parametres BETA(k) de la decomposition de la fonction
C* R11(teta) en polynomes de Legendre.
C------------------------------------------------------------------
      DO 4 J=-LUM_NBMU,LUM_NBMU
         IF(J.EQ.0) GO TO 4
         X=R11(J)*CHR(J)
         XRMU=RMU(J)

C* Calcul des polynomes de Legendre d'ordre K, pour la valeur de mu fixee.
C------------------------------------------------------
         PL(-1)=0.D+00
         PL(0)=1.D+00
         DO 5 K=0,OS_NS
            PL(K+1)=((2*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
            BETA(K)=BETA(K)+X*PL(K)
    5    CONTINUE
    4 CONTINUE

      DO 6 K=0,OS_NS
           BETA(K)=(2*K+1)*BETA(K)*.5
    6 CONTINUE
    
C* Fichier trace 
C-----------------
      IF (TRACE) THEN
           WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) ' Recombinaison de R11'
      
C* Controle de la recombinaison de R11  
C--------------------------------------
           DO 8 J=-LUM_NBMU,LUM_NBMU
              IT=0.D+00
              IF(J.EQ.0) GO TO 8
              XRMU=RMU(J)
              PL(-1)=0.D+00
              PL(0)=1.D+00
              DO 9 K=0,OS_NS
                PL(K+1)=((2*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
                IT=IT+BETA(K)*PL(K)
    9         CONTINUE
    
	      WRITE(99,203,err=921) ' RMU ',RMU(J),' R11 ',R11(J),
     &	                            ' R11 recomb  ',IT         
    8      CONTINUE
      END IF !Fichier Trace
      
C* Calcul des parametres GAMMA(k) et DELTA(k) de la decomposition des fonctions
C* R12 et R33 en polynomes ou fonctions de Legendre.
C-----------------------------------------------------------------------------
      DO 11 J=-LUM_NBMU,LUM_NBMU
         IF(J.EQ.0)  GO TO 11
         XXX=CHR(J)*R12(J)
         XX=CHR(J)*R33(J)
         POL(0)=0.D+00
         POL(1)=0.D+00
         XRMU=RMU(J)
         PL(-1)=0.D+00
         PL(0)=1.D+00
         POL(2)=3.*(1.-XRMU**2)/2./DSQRT(6.0D+00)
	
         DO 12 K=2,OS_NS
            D=(2.*K+1.)/DSQRT(1.D+00*(K+3.)*(K-1.))
            E=DSQRT(1.D+00*(K+2.)*(K-2.))/(2.*K+1.)
            POL(K+1)=D*(XRMU*POL(K)-E*POL(K-1))
            GAMMA(K)=GAMMA(K)+XXX*POL(K)
  12     CONTINUE

         DO 14 K=0,OS_NS
            PL(K+1)=((2.*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
            DELTA(K)=DELTA(K)+XX*PL(K)
  14     CONTINUE

   11 CONTINUE

      DO 15 K=0,OS_NS
         DELTA(K)=DELTA(K)*(2.*K+1.)*.5
         GAMMA(K)=GAMMA(K)*(2.*K+1.)*.5
   15 CONTINUE



C* Fichier trace 
C-----------------
      IF (TRACE) THEN
           WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) ' Recombinaison de R12'
      
C* Controle de la recombinaison de R12  
C--------------------------------------
           DO 81 J=-LUM_NBMU,LUM_NBMU
              QT=0.D+00
	      
	      IF(J.EQ.0)  GO TO 81
                 POL(0)=0.D+00
                 POL(1)=0.D+00
                 XRMU=RMU(J)
                 POL(2)=3.*(1.-XRMU**2)/2./DSQRT(6.0D+00)
   
                 DO 91 K=2,OS_NS
		    D=(2.*K+1.)/DSQRT(1.D+00*(K+3.)*(K-1.))
                    E=DSQRT(1.D+00*(K+2.)*(K-2.))/(2.*K+1.)
		    POL(K+1)=D*(XRMU*POL(K)-E*POL(K-1))
		    QT=QT+GAMMA(K)*POL(K)
  91             CONTINUE
    
	         WRITE(99,203,err=921) ' RMU ',RMU(J),' R12 ',R12(J),
     &	                            ' R12 recomb  ',QT 
  81       CONTINUE
  
  
           WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) ' Recombinaison de R33'
      
C* Controle de la recombinaison de R33  
C--------------------------------------
           DO 82 J=-LUM_NBMU,LUM_NBMU
	      UT=0.D+00
	      
	      IF(J.EQ.0)  GO TO 82
                 PL(-1)=0.D+00
                 PL(0)=1.D+00
		 XRMU=RMU(J)

                 DO 92 K=0,OS_NS
                    PL(K+1)=((2*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
                    UT=UT+DELTA(K)*PL(K)
  92             CONTINUE
    
	         WRITE(99,203,err=921) ' RMU ',RMU(J),' R33 ',R33(J),
     &                              ' R33 recomb  ',UT   
  82       CONTINUE
      END IF !Fichier Trace
      
      
      
C* Calcul des parametres ALPHA(k) et ZETA(k) en fonctions des BETA(k) et DELTA(k).
C-------------------------------------------------------------------------------
      DO 16 I=2,OS_NS
         CO1=4*(2*I+1.)/FLOAT(I)/(I-1.)/(I+1.)/(I+2.)
         CO2=I*(I-1.)/((I+1.)*(I+2.))
         CO3=CO2*DELTA(I)
         CO2=CO2*BETA(I)
         NN=I*.5
         MM=(I-1)*.5
         SOM1=0.D+00
         SOM2=0.D+00
         SOM3=0.D+00
         SOM4=0.D+00
         DO 17 J=1,NN
            X2=(I-1.)*(I-1.)-3.*(2*J-1.)*(I-J)
            SOM1=SOM1+X2*BETA(I-2*J)
            SOM2=SOM2+X2*DELTA(I-2*J)
  17     CONTINUE

         DO 18 J=0,MM
            X2=(I-1.)*(I-1.)-3.*J*(2*I-2*J-1.)
            SOM3=SOM3+X2*BETA(I-2*J-1)
            SOM4=SOM4+X2*DELTA(I-2*J-1)
  18     CONTINUE

         ZETA(I)=CO3-CO1*(SOM2-SOM3)
         ALPHA(I)=CO2-CO1*(SOM1-SOM4)
  16  CONTINUE


C* Enregistrement des resultats
C---------------------------------
      DO 20 K=0,OS_NS
         WRITE(3,207,err=992)ALPHA(K),BETA(K),GAMMA(K),ZETA(K)
   20 CONTINUE



C* Fichier trace 
C-----------------
      IF (TRACE) THEN
           WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) 'Ordre limite : ',OS_NS
	   WRITE(99,215,err=921)'K','ALPHA(K)','BETA(K)','GAMMA(K)',
     &                      'DELTA(I)','ZETA(K)'  
           DO 45 K=0,OS_NS
	     WRITE(99,205,err=921) K,ALPHA(K),BETA(K),GAMMA(K),DELTA(K),
     &	                           ZETA(K)
   45      CONTINUE
      END IF 
   
   
       
C* Fin nominale 
C-------------------
      IER=0
      goto 9999
      
C* Cas d'erreur : retour du status -1 au programme appelant
C------------------------------------------------------
  991 WRITE(6,*) '  ERROR on a file opening for SOS_MAT_FRESNEL'
      IER=-1
      goto 9999
  921 WRITE(6,*) '  ERROR on logfile writing for SOS_MAT_FRESNEL'
      IER=-1
      goto 9999
  992 WRITE(6,*) '  ERROR on a file reading/writing for SOS_MAT_FRESNEL'
      IER=-1 
      goto 9999

 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(3)
   
C* Format
C---------   
  207 FORMAT(4(E15.8))
  200 FORMAT(A17,F6.2,3(A5,E13.7))
  203 FORMAT(2(A5,E13.7),A12,E13.7)
  205 FORMAT(1X,I3,5(2X,E13.5))
  215 FORMAT(3X,A1,7X,A8,7X,A7,8X,A8,7X,A8,7X,A7)    
      RETURN
      END         !FIN DE LA PROCEDURE SOS_MAT_FRESNEL




C==============================================================================
C PROCEDURE: SOS_MAT_REFLEXION
C ==========
C      Cette procedure calcule la matrice de reflexion de surface 
C      pour le glitter.
C      Les angles d'incidence et de reflexion sont exprimes comme etant
C      inferieurs a 90 degres (mu de Gauss > 0).
C
C
C
C      NB: DEMARCHE D'IMPLEMENTATION :
C      *******************************
C          L'implementation du code est adaptee pour un calcul de la matrice 
C          de reflexion (routine SOS_MAT_REFLEXION) n'utilisant que des valeurs
C           positives des cos(Angle de Gauss). Pour y parvenir :
C          - On utilise non pas l'angle de diffusion vraie mais son angle
C            complementaire (2*angle de reflexion) : DEM 1
C            (routine SOS_MAT_FRESNEL).
C          - On utilise non pas le vrai coefficient de reflexion parallele
C            mais son oppose : DEM 2 (routine SOS_MAT_FRESNEL).
C
C          Cette routine calculant la matrice de reflexion introduit
C          un coefficient (-1)**K pour le calcul des termes de la 
C          series de Fourier pour assurer la coherence de l'ensemble : DEM 3.
C
C
C Description des parametres
C --------------------------
C     WIND	(double)    (E) Vitesse du vent (m/s).
C     LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles.
C     RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)    (E) Tableau des mu de Gauss 
C                             (angles complementaires des angles de diffusion).
C     OS_NB (I4) : (E) Ordre des developpements des fonctions de phase 
C		       en Polynomes de Legendre a utiliser.
C     OS_NS (I4) : (E) Ordre maximal a utiliser pour le developpement en
C                      polynomes de Legendre des elements de la matrice de Fresnel 
C                      et en series de Fourier des calculs de luminance.
C     OS_NM (I4) : (E) Ordre maximal du developpement en series de Fourier  
C                      de la fonction G des calculs de matrice de reflexion.
C     IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C Autre resultat fourni par la procedure
C -------------------------------------
C      Fichier temporaire RES_MAT_REFLEX
C
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier temporaire RES_GSF
C       Fichier binaire non formatte, cree par SOS_GSF
C       Fichier cree dans le repertoire de lancement du programme
C       Chaque enregistrement contient :
C           -numero d'angle Teta1 (valeur entiere)
C           -numero d'angle Teta2 (valeur entiere)
C           -ordre maximale utile de la decomposition : IL  (valeur entiere)
C           -Fonction G decomposee pour (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 valeurs double precision)
C
C   -->Contenu du fichier temporaire RES_FRESNEL
C      Fichier Ascii, format d'un enregistrement 4E15.8, 
C      cree par SOS_MAT_FRESNEL
C      Fichier cree dans le repertoire de lancement du programme
C      Nombre d'enregistrements : de k=0 a OS_NS
C      Chaque enregistrement contient les parametres ALPHA(k), BETA(k), 
C      GAMMA(k) et ZETA(k) associes a la matrice de Fresnel
C
C   -->Contenu du fichier temporaire RES_MAT_REFLEX
C       Fichier binaire non formatte, cree par SOS_MAT_REFLEXION
C       Fichier cree dans le repertoire de lancement du programme
C       Si un fichier de meme nom existe deja, il est ecrase
C       Enregistrement pour chaque couple (Teta1, Teta2) de
C       la decomposition en series de Fourier de la matrice de Fresnel:
C       M11,M12,M13,M21,M22,M23,M31,M32,M33 
C       (n valeurs reelles)
C
C
C Common utilise:
C --------------
C    Aucun
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C      Le fichier RES_GSF ou le fichier RES_FRESNEL n'existe pas
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C
C     
C  Routines utilises par la procedure:
C ------------------------------------
C  Ce programme fait appel aux routines:
C         -SOS_NOYAUX_FRESNEL
C
C=============================================================================
      SUBROUTINE SOS_MAT_REFLEXION(WIND,LUM_NBMU,RMU,OS_NB,OS_NS,OS_NM,
     &                             IER)

      IMPLICIT NONE

C* Definition des variables     
C*--------------------------------------------------------------------
      DOUBLE PRECISION WIND	! Vitesse du vent (m/s).
      DOUBLE PRECISION SIG 	! Carre du parametre Sigma de la theorie 
      				! de Cox et Munck.
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				! Tableau des mu de Gauss 
				!(angles complementaires des angles de diffusion)

      DOUBLE PRECISION ALPHA(0:SOS_OS_NS_MAX)   ! Parametres alpha(k).
      DOUBLE PRECISION BETA(0:SOS_OS_NS_MAX)	! Parametres beta(k).
      DOUBLE PRECISION GAMMA(0:SOS_OS_NS_MAX)   ! Parametres gamma(k).
      DOUBLE PRECISION ZETA(0:SOS_OS_NS_MAX)	! Parametres zeta(k).

      DOUBLE PRECISION G(0:SOS_OS_NM_MAX)	! Fonction G decomposee en series 
      					        ! de Fourier.

        ! Element Pij de la matrice de reflexion de Fresnel decomposee en series 
	! de Fourier pour un ordre IS = 0 a OS_NS.
	! Calcul pour la direction d'incidence Teta1 (RMU1)
	! et de reflexion Teta2 (RMU2).
      DOUBLE PRECISION BP(0:SOS_OS_NS_MAX,2)    !   Element P11.
      DOUBLE PRECISION GR(0:SOS_OS_NS_MAX,2)    !   Element P12 ou P21.
      DOUBLE PRECISION GT(0:SOS_OS_NS_MAX,2)    !   Element -P13 ou -P31.
      DOUBLE PRECISION ARR(0:SOS_OS_NS_MAX,2)   !   Element P22.
      DOUBLE PRECISION ART(0:SOS_OS_NS_MAX,2)   !   Element -P23 ou -P32.
      DOUBLE PRECISION ATT(0:SOS_OS_NS_MAX,2)   !   Element P33.

      DOUBLE PRECISION R111    	! Variable de calcul du terme M11(IS,1).
      DOUBLE PRECISION R121    	! Variable de calcul du terme M12(IS,1).
      DOUBLE PRECISION R122    	! Variable de calcul du terme M12(IS,2).
      DOUBLE PRECISION R131    	! Variable de calcul du terme M13(IS,1).
      DOUBLE PRECISION R132    	! Variable de calcul du terme M13(IS,2).
      DOUBLE PRECISION R211    	! Variable de calcul du terme M21(IS,1).
      DOUBLE PRECISION R212    	! Variable de calcul du terme M21(IS,2).
      DOUBLE PRECISION R221    	! Variable de calcul du terme M22(IS,1).
      DOUBLE PRECISION R222    	! Variable de calcul du terme M22(IS,2).
      DOUBLE PRECISION R231     ! Variable de calcul du terme M23(IS,1).
      DOUBLE PRECISION R232     ! Variable de calcul du terme M23(IS,2).
      DOUBLE PRECISION R311     ! Variable de calcul du terme M31(IS,1).
      DOUBLE PRECISION R312     ! Variable de calcul du terme M31(IS,2).
      DOUBLE PRECISION R321     ! Variable de calcul du terme M32(IS,1).
      DOUBLE PRECISION R322     ! Variable de calcul du terme M32(IS,2).
      DOUBLE PRECISION R331     ! Variable de calcul du terme M33(IS,1).
      DOUBLE PRECISION R332     ! Variable de calcul du terme M33(IS,2).

      DOUBLE PRECISION X,Y	! Coefficients multiplicatifs de la matrice de
      				! reflexion.

      	! Elements Mij de la matrice de reflexion en fonction de l'ordre IS
	! de la serie de Fourier, pour 2 directions.
	! L'indice de direction est :
	!     _ 1 pour une incidence sous Teta1 et une reflexion sous Teta2,
	!     _ 2 pour une incidence sous Teta2 et une reflexion sous Teta1.
      REAL M11(0:SOS_OS_NB_MAX,2)		
      REAL M12(0:SOS_OS_NB_MAX,2)	
      REAL M13(0:SOS_OS_NB_MAX,2)		
      REAL M21(0:SOS_OS_NB_MAX,2)	
      REAL M22(0:SOS_OS_NB_MAX,2)		
      REAL M23(0:SOS_OS_NB_MAX,2)		
      REAL M31(0:SOS_OS_NB_MAX,2)	
      REAL M32(0:SOS_OS_NB_MAX,2)		
      REAL M33(0:SOS_OS_NB_MAX,2)	

      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations
      
      INTEGER*4 OS_NB 	 ! Ordre des developpements des fonctions de phase 
			 ! en Polynomes de Legendre.
      INTEGER*4 OS_NS    ! Ordre maximal du developpement en polynomes de Legendre   
			 ! des elements de la matrice de Fresnel.
      INTEGER*4 OS_NM    ! Ordre maximal du developpement en series de Fourier  
			 ! de la fonction G des calculs de matrice de reflexion.
 
      INTEGER*4 I	! Indice de l'angle Teta1.
      INTEGER*4 J	! Indice de l'angle Teta2.

      INTEGER*4 L	! Ordre de la decomposition en polynomes de Legendre.
      INTEGER*4 LIM	! Ordre maximal de la decomposition en polynomes
      			! de Legendre.
      INTEGER*4 II,II2	! Indices enregistres avec la fonction G(IS).
      INTEGER*4 IJ	
      INTEGER*4 IS	! Ordre de la decomposition en series de Fourier.
      INTEGER*4 K
      INTEGER*4 I1	! Valeur K+IS
      INTEGER*4 I2	! Valeur ABS(K-IS)
      INTEGER*4 IM	! Valeur 1 ou -1.
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon  


C* Creation du fichier resultat de SOS_MAT_REFLEXION
C---------------------------------------------------------
      OPEN(UNIT=4,FILE='RES_MAT_REFLEX',FORM='UNFORMATTED',ERR=991)

C* Carre du parametre Sigma de la theorie de Cox et Munck.
C---------------------------------------------------------
      SIG = .003 + .00512*WIND

C* Lecture des parametres de la decomposition de la matrice de Fresnel
C* en polynomes ou fonctions de Legendre.
C----------------------------------------------------------------------
      OPEN(UNIT=2,FILE='RES_FRESNEL',STATUS='OLD',ERR=991)
      
      DO 951 L=0,OS_NS
         READ(2,1300,err=992) ALPHA(L),BETA(L),GAMMA(L),ZETA(L)
  951 CONTINUE
  
      CLOSE(2)
      
C* Ouverture du fichier RES_GSF
C---------------------------------------------------------
      OPEN(UNIT=3,FILE='RES_GSF',FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)

C* Boucle sur les directions d'incidence et de reflexion.
C--------------------------------------------------------
      DO 25 I=1,LUM_NBMU   !Teta1 : angle zenithal d'incidence  
                              !(0 < Teta1 < PI/2),
         DO 2 J=1,I           !Teta2 : angle zenithal de reflexion 
	                      !(0 < Teta2 < PI/2),

C* Lecture de la decomposition en series de Fourier de la fonction G
C* pour les angles Teta1 et Teta2,avec un ordre maximal utile L.
C-------------------------------------------------------------------
            READ(3,err=992)II,II2,LIM,(G(K),K=0,LIM)

C* Complement du tableau a 0.
C---------------------------------
            DO 246 IJ=LIM+1,OS_NM
  246          G(IJ)=0.
  
C------------------------------------------------------------------------
C* Decomposition en series de Fourier de la matrice de reflexion de Fresnel
C* pour les directions Teta1 et Teta2
C--------------------------------------------------------------------------

C*  Calcul des fonctions BP,GR,GT,ARR,ART et ATT en fonction de l'ordre IS
C*  pour les directions Teta1 et Teta2
C---------------------------------------------------------------------------
            CALL SOS_NOYAUX_FRESNEL(RMU(I),RMU(J),OS_NS,
     &                              ALPHA,BETA,GAMMA,ZETA,
     &                              BP,GR,GT,ARR,ART,ATT)

C* Calcul des elements de la matrice de reflexion pour chaque ordre IS.
C---------------------------------------------------------------------
            DO 24 IS=0,OS_NB
               X=G(IS)/(4*SIG)
	       
C*   Elements pour IS=0.
               R111     =X*BP(0,1)
               R121     =X*GR(0,1)
               R122     =X*GR(0,2)
               R131     =0.
               R132     =0.
               R231     =0.
               R232     =0.
               R211     =X*GR(0,2)
               R212     =X*GR(0,1)
               R221     =X*ARR(0,2)
               R222     =X*ARR(0,1)
               R311     =0.
               R312     =0.
               R321     =0.
               R322     =0.
               R331     =X*ATT(0,2)
               R332     =X*ATT(0,1)

C*   Elements pour la sommation de K=1 a OS_NS
               IM=1
               DO 3 K=1,OS_NS
                 IM=-IM						! DEM 3
                 I1=K+IS
                 I2=ABS(K-IS)
                 IF ((I1.GT.LIM).AND.(I2.GT.LIM)) GOTO 3	
                 X=IM*(G(I1)+G(I2))/(4*SIG)			! DEM 3
                 Y=IM*(G(I2)-G(I1))/(4*SIG)			! DEM 3
                 R111     =R111     +BP(K,1)*X
                 R121     =R121     +GR(K,1)*X
                 R122     =R122     +GR(K,2)*X
                 R131     =R131     +GT(K,1)*Y
                 R132     =R132     +GT(K,2)*Y
                 R211     =R211     +GR(K,2)*X
                 R212     =R212     +GR(K,1)*X
                 R221     =R221     +ARR(K,2)*X
                 R222     =R222     +ARR(K,1)*X
                 R231     =R231     +ART(K,2)*Y
                 R232     =R232     +ART(K,1)*Y
                 R311     =R311     +GT(K,2)*Y
                 R312     =R312     +GT(K,1)*Y
                 R321     =R321     +ART(K,1)*Y
                 R322     =R322     +ART(K,2)*Y
                 R331     =R331     +ATT(K,2)*X
                 R332     =R332     +ATT(K,1)*X
  3            CONTINUE

C* Affectation des resultats aux tableaux resultats Rij(IS)
 244           M11(IS,1)=R111
               M11(IS,2)=R111
               M12(IS,1)=R121
               M12(IS,2)=R122
               M13(IS,1)=R131
               M13(IS,2)=R132
               M21(IS,1)=R211
               M21(IS,2)=R212
               M22(IS,1)=R221
               M22(IS,2)=R222
               M23(IS,1)=R231
               M23(IS,2)=R232
               M31(IS,1)=-R311
               M31(IS,2)=-R312
               M32(IS,1)=-R321
               M32(IS,2)=-R322
               M33(IS,1)=-R331
               M33(IS,2)=-R332
 24         CONTINUE	!Fin de boucle sur IS

C* Enregistrement de l'ordre pour le couple (Teta1, Teta2) de
C* la decomposition en series de Fourier de la matrice de Fresnel.
C-----------------------------------------------------------------
            WRITE(4,err=992) M11,M12,M13,M21,M22,M23,M31,M32,M33

 2       CONTINUE	!Fin de boucle sur Teta2 (J)

 25   CONTINUE	!Fin de boucle sur Teta1 (I)
 
 
C* Fin nominale 
C-------------------
      IER=0
      goto 9999
      
C* Cas d'erreur : retour du status -1 au programme appelant
C----------------------------------------------------------    
  991 WRITE(6,*) '  ERROR on a file opening for SOS_MAT_REFLEXION'
      IER=-1
      goto 9999
  992 WRITE(6,*) '  ERROR on a file reading/writing for ',
     &           'SOS_MAT_REFLEXION '
      IER=-1 
      goto 9999

 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
 
C* Format
C---------   
 1300 FORMAT(4(E15.8))  
   
      RETURN
      END         !FIN DE LA PROCEDURE SOS_MAT_REFLEXION





C==============================================================================
C PROCEDURE: SOS_NOYAUX_FRESNEL
C ==========
C      Cette procedure calcule la decomposition en series de Fourier de
C      la matrice de Fresnel selon l'azimut, pour les directions Teta1 et Teta2.
C
C      On calcule les elements de matrice suivants :
C         BP(IS,mu)    GR(IS,mu)    GT(IS,mu)
C         ARR(IS,mu)   ART(IS,mu)   ATT(IS,mu)
C      pour les ordres IS de la decomposition en series de Fourier
C      et pour mu = cos(Teta1) ou cos(Teta2).
C
C
C Description des parametres
C --------------------------
C	RMU1   (double)	     (E) Cos(Teta1).
C	RMU2   (double)	     (E) Cos(Teta2).
C       OS_NS (I4)           (E) Ordre maximal a utiliser pour le developpement en
C                                polynomes de Legendre des elements de la matrice 
C                                de Fresnel  et en series de Fourier des calculs 
C                                de luminance.
C
C	ALPHA(0:SOS_OS_NS_MAX)  (double)	(E) Parametres alpha(k).
C	BETA(0:SOS_OS_NS_MAX)   (double)	(E) Parametres beta(k).
C	GAMMA(0:SOS_OS_NS_MAX)  (double)	(E) Parametres gamma(k).
C	ZETA(0:SOS_OS_NS_MAX)   (double)	(E) Parametres zeta(k).
C
C       Element Pij de la matrice de reflexion de Fresnel decomposee en series 
C       de Fourier pour un ordre IS = 0 a OS_NS.
C       Calcul pour la direction d'incidence Teta1 (RMU1)
C       et de reflexion Teta2 (RMU2).	
C	BP(0:SOS_OS_NS_MAX,2)   (double)	(S) Element P11.
C	GR(0:SOS_OS_NS_MAX,2)   (double)	(S) Element P12 ou P21.
C	GT(0:SOS_OS_NS_MAX,2)   (double)	(S) Element -P13 ou -P31.
C	ARR(0:SOS_OS_NS_MAX,2)  (double)	(S) Element P22.
C	ART(0:SOS_OS_NS_MAX,2)  (double)	(S) Element -P23 ou -P32.
C	ATT(0:SOS_OS_NS_MAX,2)  (double)	(S) Element P33.
C
C 
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Aucun
C
C
C=============================================================================
      SUBROUTINE SOS_NOYAUX_FRESNEL(RMU1,RMU2,OS_NS,
     &                              ALPHA,BETA,GAMMA,ZETA,
     &                              BP,GR,GT,ARR,ART,ATT)

      IMPLICIT NONE

C* Definition des variables                                  
C*----------------------------------------------------------------------
      DOUBLE PRECISION RMU1	! Cos(Teta1).
      DOUBLE PRECISION RMU2	! Cos(Teta2).

      DOUBLE PRECISION ALPHA(0:SOS_OS_NS_MAX)       ! Parametres alpha(k).
      DOUBLE PRECISION BETA(0:SOS_OS_NS_MAX)	! Parametres beta(k).
      DOUBLE PRECISION GAMMA(0:SOS_OS_NS_MAX)       ! Parametres gamma(k).
      DOUBLE PRECISION ZETA(0:SOS_OS_NS_MAX)	! Parametres zeta(k).

      DOUBLE PRECISION R(2)	! tableau contenant cos(Teta1) et cos(Teta2).

      	! Element Pij de la matrice de reflexion de Fresnel decomposee en series 
	! de Fourier pour un ordre IS = 0 a OS_NS.
	! Calcul pour la direction d'incidence Teta1 (RMU1)
	! et de reflexion Teta2 (RMU2).	
      DOUBLE PRECISION BP(0:SOS_OS_NS_MAX,2)    ! Element P11.
      DOUBLE PRECISION GR(0:SOS_OS_NS_MAX,2)    ! Element P12 ou P21.
      DOUBLE PRECISION GT(0:SOS_OS_NS_MAX,2)    ! Element -P13 ou -P31.
      DOUBLE PRECISION ARR(0:SOS_OS_NS_MAX,2)   ! Element P22.
      DOUBLE PRECISION ART(0:SOS_OS_NS_MAX,2)   ! Element -P23 ou -P32.
      DOUBLE PRECISION ATT(0:SOS_OS_NS_MAX,2)   ! Element P33.

      DOUBLE PRECISION PSL(-1:SOS_OS_NS_MAX,2)  
                                ! Fonction de Legendre PSL pour l'ordre S=IS
      				! en fonction de l'ordre L et de l'angle mu.
      DOUBLE PRECISION RSL(-1:SOS_OS_NS_MAX,2)  ! Fonction RSL.
      DOUBLE PRECISION TSL(-1:SOS_OS_NS_MAX,2)  ! Fonction TSL.

      DOUBLE PRECISION C	     ! Valeur ponctuelle du mu de Gauss.
      DOUBLE PRECISION A,B,D,E,F     ! Coefficients pour la recurrence 
      				     ! sur PSL, RSL et TSL.

      DOUBLE PRECISION X	! Variable intermediaire de calcul.
      DOUBLE PRECISION XX	! Variable intermediaire de calcul.
      DOUBLE PRECISION YY	! Variable intermediaire de calcul.

      DOUBLE PRECISION RAC3	! Racine de 3.
      DOUBLE PRECISION X26	! 2*racine(6)

      INTEGER*4 OS_NS   ! Ordre maximal du developpement en polynomes de Legendre   
			! des elements de la matrice de Fresnel.

      INTEGER*4 IS	! Ordre de la decomposition en series de Fourier.
      INTEGER*4 I
      INTEGER*4 J
      INTEGER*4 K
      INTEGER*4 L


C* Initialisation
C-------------------------
      RAC3=DSQRT(3.D+00)
      X26=2.*DSQRT(6.D+00)

C* Stockage de cos(Teta1) et cos(Teta2) dans le tableau R.
C------------------------------------------------------
      R(1)=RMU1
      R(2)=RMU2

C* Decomposition en series de Fourier de IS = 0 a OS_NS
C----------------------------------------------------------
      DO 1 IS=0,OS_NS

C-----------------------------------------------------
C* Initialisation des fonctions d'ordre IS pour les
C* ordres L = 0, 1 et 2.
C-----------------------------------------------------

C* Initialisation pour l'ordre IS = 0
C-------------------------------------
         IF (IS.NE.0) GOTO 700
         DO 25 J=1,2
            C=R(J)
            PSL(0,J)=1
            PSL(1,J)=C
            X=(3*C*C-1)*0.5
            PSL(2,J)=X
            RSL(1,J)=0
            X=3*(1-C*C)/X26
            RSL(2,J)=X
            TSL(1,J)=0.
            TSL(2,J)=0.
   25    CONTINUE
         GOTO 501

C* Initialisation pour l'ordre IS = 1
C-------------------------------------
  700    IF (IS.NE.1) GOTO 701
         DO 26 J=1,2
            C=R(J)
            X=1-C*C
            PSL(0,J)=0
            PSL(1,J)=DSQRT(X*0.5)
            PSL(2,J)=C*PSL(1,J)*RAC3
            TSL(1,J)=0.
            RSL(1,J)=0
            RSL(2,J)=-C*DSQRT(X)*0.5
            TSL(2,J)=-DSQRT(X)*0.5
   26    CONTINUE
         GOTO 501

C* Initialisation pour l'ordre IS >= 2
C----------------------------------------
  701    A=1
         DO 27 I=1,IS
            X=I
   27       A=A*DSQRT((I+IS)/X)*0.5
   
         B=A*DSQRT(IS/(IS+1.D+00))*DSQRT((IS-1.D+00)/(IS+2.)) 
         DO 28 J=1,2
             C=R(J)
             XX=1-C*C
             YY=IS*0.5
             PSL(IS-1,J)=0.
             RSL(IS-1,J)=0.
             TSL(IS-1,J)=0.
             X=XX**YY
             PSL(IS,J)=A*X
             YY=YY-1
             X=XX**YY
             RSL(IS,J) = B*(1+C*C)*X
             TSL(IS,J)= 2*B*C*X
   28    CONTINUE

C* Calcul des fonctions PSL(L,mu), RSL(L,mu) et TSL(L,mu)
C* par recurrence sur l'ordre L pour IS fixe.
C------------------------------------------------------
  501    K=2
         IF (IS.GT.2) K=IS

         DO 30 L=K,OS_NS-1
C* Coefficients pour la recurrence sur PSL
             A=(2*L+1.)/DSQRT((L+IS+1.D+00)*(L-IS+1.))
             B=DSQRT(DFLOAT((L+IS)*(L-IS)))/(2.*L+1.) 

C* Coefficients pour la recurrence sur RSL et TSL
             D=(L+1.)*(2*L+1.)/DSQRT((L+3.D+00)*(L-1.)*(L+IS+1.)
     &         *(L-IS+1.))
             E=DSQRT((L+2.D+00)*(L-2.)*(L+IS)*(L-IS))/(L*(2.*L+1.))
             F=2.*IS/(L*(L+1.))

C* Application de la recurrence pour chaque mu de Gauss
C------------------------------------------------------
             DO 31 J=1,2
                C=R(J)
                X=A*(C*PSL(L,J)-B*PSL(L-1,J))
                PSL(L+1,J)=X
                X=D*(C*RSL(L,J)-F*TSL(L,J)-E*RSL(L-1,J))
                RSL(L+1,J)=X
                X=D*(C*TSL(L,J)-F*RSL(L,J)-E*TSL(L-1,J))
                TSL(L+1,J)=X
   31        CONTINUE
   30    CONTINUE


C------------------------------------------------------
C* Calcul des elements de la matrice de reflexion
C------------------------------------------------------

C* Boucles sur les directions Teta1 et Teta2
C------------------------------------------------------
         DO 32 K=1,2
             J=3-K

C* Initialisation des tableaux
             BP(IS,K)=0.
             ARR(IS,K)=0.
             ATT(IS,K)=0
             GR(IS,K)=0.
             GT(IS,K)=0.
             ART(IS,K)=0

C* Boucle de sommation sur L = IS a OS_NS
C-----------------------------------------------
             DO 32 L=IS,OS_NS
                BP(IS,K)=BP(IS,K)+BETA(L)*PSL(L,J)*PSL(L,K)
	        GR(IS,K)=GR(IS,K)+GAMMA(L)*PSL(L,J)*RSL(L,K)
                GT(IS,K)=GT(IS,K)+GAMMA(L)*PSL(L,J)*TSL(L,K)
                ATT(IS,K)=ATT(IS,K)+ALPHA(L)*TSL(L,J)*TSL(L,K)+
     &                    ZETA(L)*RSL(L,J)*RSL(L,K)
                ARR(IS,K)=ARR(IS,K)+ZETA(L)*TSL(L,J)*TSL(L,K)+
     &		          ALPHA(L)*RSL(L,J)*RSL(L,K)
                ART(IS,K)=ART(IS,K)+ALPHA(L)*RSL(L,K)*TSL(L,J)+
     &                    ZETA(L)*RSL(L,J)*TSL(L,K)

   32    CONTINUE  ! boucle sur K et L

   1  CONTINUE
      RETURN
      END         !FIN DE LA PROCEDURE SOS_NOYAUX_FRESNEL






C==============================================================================
C PROCEDURE:  SOS_MISE_FORMAT
C ==========
C      Cette procedure reordonne les elements de la matrice Mij en une matrice
C      Pij utilisable par le code des OS.
C
C      Les elements Mij, enregistres dans le fichier "RES_MAT_REFLEX", donnent
C      pour chaque couple d'angles (Teta1, Teta2) la decomposition en series 
C      de Fourier associee. 
C
C      Les elements Pij, enregistres dans le fichier de sortie, donnent
C      pour chaque ordre de la decomposition en series de Fourier
C      la reflexion associee a l'ensemble des couples (Teta1, Teta2) .
C
C Description des parametres
C --------------------------
C      FICGLITTER  (CHARACTER*SOS_LENFIC2) : (E)  Nom complet du fichier resultat.
C      NBMU (I4)  : (E)  Nombre d'angles (positifs) effectivement utiles.
C      OS_NB (I4) : (E)  Ordre des developpements des fonctions de phase 
C		         en Polynomes de Legendre a utiliser.
C      IER  (I4)  : (S)  code d'erreur =0 si pas d'erreur, =-1 sinon
C
C
C Autre resultat fourni par la procedure
C -------------------------------------
C      Fichier resultat : FICGLITTER, fichier GLITTER
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier temporaire RES_MAT_REFLEX
C       Fichier binaire non formatte, cree par SOS_MAT_REFLEXION
C       Fichier cree dans le repertoire de lancement du programme
C       Enregistrement pour chaque couple (Teta1, Teta2) de
C       la decomposition en series de Fourier de la matrice de Fresnel:
C       M11,M12,M13,M21,M22,M23,M31,M32,M33 
C       (n valeurs reelles)
C
C   -->Contenu du fichier GLITTER
C       Fichier binaire non formatte, cree par SOS_MISE_FORMAT
C       Il contient la matrice de reflexion de surface 
C       (fichier GLITTER ou autres BRDF):
C       c'est a dire les elements Pij de la matrice de reflexion 
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
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C      Le fichier RES_MAT_REFLEX n'existe pas
C      Le fichier Glitter existe
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C
C=============================================================================
      SUBROUTINE SOS_MISE_FORMAT(FICGLITTER,NBMU,OS_NB,IER)

      IMPLICIT NONE

C* Definition des variables                    
C*--------------------------------------------------------------------
      CHARACTER*SOS_LENFIC2 FICGLITTER	! Nom du fichier resultat.

      	! Elements Mij de la matrice de reflexion en fonction de l'ordre IS
	! de la serie de Fourier, pour 2 directions.
	! L'indice de direction est :
	!     _ 1 pour une incidence sous Teta1 et une reflexion sous Teta2,
	!     _ 2 pour une incidence sous Teta2 et une reflexion sous Teta1.
      REAL M11(0:SOS_OS_NB_MAX,2)		
      REAL M12(0:SOS_OS_NB_MAX,2)		
      REAL M13(0:SOS_OS_NB_MAX,2)			
      REAL M21(0:SOS_OS_NB_MAX,2)	
      REAL M22(0:SOS_OS_NB_MAX,2)		
      REAL M23(0:SOS_OS_NB_MAX,2)		
      REAL M31(0:SOS_OS_NB_MAX,2)	
      REAL M32(0:SOS_OS_NB_MAX,2)	
      REAL M33(0:SOS_OS_NB_MAX,2)

      	! Elements Pij de la matrice de reflexion en fonction de la direction
	! d'incidence et de reflexion, pour un ordre IS du developpement
	! en series de Fourier
      REAL P11(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P12(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P13(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)	
      REAL P21(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P22(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P23(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P31(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P32(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P33(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)

      INTEGER*4 NBMU        ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB 	    ! Ordre des developpements des fonctions de phase 
			    ! en Polynomes de Legendre.
      INTEGER*4 IS	! Ordre de la decomposition en series de Fourier.
      INTEGER*4 I	! Indice de l'angle d'incidence Teta1.
      INTEGER*4 J	! Indice de l'angle de reflexion Teta2.
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon  
	
C* Ouverture fichier
C------------------------------------------------------------	
	
      OPEN(UNIT=1,FILE='RES_MAT_REFLEX',FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=2,FILE=FICGLITTER,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=991)

C* Boucle sur l'ordre de la decomposition en series de Fourier
C------------------------------------------------------------
      DO 1 IS=0,OS_NB

C* Boucle sur les directions d'incidence Teta1
C---------------------------------------------
         DO 2 I=1,NBMU

C* Boucle sur les directions de reflexion Teta2
C-----------------------------------------------

C*       Lecture du fichier RES_MAT_REFLEX enregistre pour la taille max
C*       des tableaux : SOS_OS_NBMU_MAX    
            DO 3 J=1,I
               READ(1,err=992) M11,M12,M13,M21,M22,M23,M31,M32,M33
               P11(I,J)=M11(IS,1)
               P11(J,I)=M11(IS,2)
               P12(I,J)=M12(IS,1)
               P12(J,I)=M12(IS,2)
               P13(I,J)=M13(IS,1)
               P13(J,I)=M13(IS,2)
               P21(I,J)=M21(IS,1)
               P21(J,I)=M21(IS,2)
               P22(I,J)=M22(IS,1)
               P22(J,I)=M22(IS,2)
               P23(I,J)=M23(IS,1)
               P23(J,I)=M23(IS,2)
               P31(I,J)=M31(IS,1)
               P31(J,I)=M31(IS,2)
               P32(I,J)=M32(IS,1)
               P32(J,I)=M32(IS,2)
               P33(I,J)=M33(IS,1)
               P33(J,I)=M33(IS,2)
 3          CONTINUE
 2       CONTINUE
 
C*       Enregistrement des tableaux de la matrice de reflexion 
C*       a l'ordre IS de la decomposition en series de Fourier.
C*       Bornes d'enregistrement des tableaux limitees au 
C*       dimensionnement utile -NBMU:NBMU
C-------------------------------------------------------------- 
         WRITE(2,err=992) ((P11(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P12(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P13(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P21(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P22(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P23(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P31(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P32(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P33(I,J),I=1,NBMU),J=1,NBMU)

C* Retour au debut du fichier avant passage a l'ordre IS suivant.
C---------------------------------------------------------------
         REWIND (1)
 1    CONTINUE


      
C* Fin nominale 
C-------------------
      IER=0
      goto 9999
      
C* Cas d'erreur : retour du status -1 au programme appelant
C------------------------------------------------------
  991 WRITE(6,*) '  ERROR on a file opening for SOS_MISE_FORMAT'
      IER=-1
      goto 9999
  992 WRITE(6,*) '  ERROR on a file reading/writing for SOS_MISE_FORMAT'
      IER=-1 
      goto 9999

 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(1)
      CLOSE(2)
      
      RETURN
      END         !FIN DE LA PROCEDURE SOS_MISE_FORMAT

      
      
      
C==============================================================================
C PROCEDURE: SOS_BPDF_AJOUT_BRDF
C ==========
C      Cette procedure lit deux fichiers de BRDF ou BPDF mis en forme pour les OS
C      et realise la sommation des termes pour creer un fichier final somme des deux
C      premiers.
C
C Description des parametres
C --------------------------
C       FICSURF1 (CHARACTER*SOS_LENFIC2) (E) : nom du 1er fichier de reflexion de surface, 
C                                              au format pour les OS.
C                                             (nom complet avec repertoire).
C       FICSURF2 (CHARACTER*SOS_LENFIC2) (E) : nom du 2nd fichier de reflexion de surface, 
C                                              au format pour les OS.
C                                             (nom complet avec repertoire).
C       NBMU (I4)  (E) : nombre d'angles (positifs) effectivement utiles.
C       OS_NB (I4) (E) : Ordre des developpements des fonctions de phase 
C		         en Polynomes de Legendre a utiliser.
C       FICS   (CHARACTER*SOS_LENFIC2)  (E)  : nom du fichier SURFACE resultat.
C                                             (nom complet avec repertoire).
C       IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier BRDF ou BPDF
C       Fichier binaire non formatte.
C       Il contient la matrice de reflexion de surface (fichier BRDF / BPDF):
C       c'est a dire les elements Rij de la matrice de reflexion 
C       pour OS_NB ordres de la decomposition en series de Fourier.
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
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C============================================================================== 
      SUBROUTINE SOS_BPDF_AJOUT_BRDF(FICSURF1,FICSURF2,NBMU,OS_NB,
     &                               FICS,IER)

      IMPLICIT NONE

C* Definition des variables   
C*--------------------------------------------------------------------------
      CHARACTER*SOS_LENFIC2 FICSURF1      ! nom du 1er fichier de reflexion 
      					  ! de surface, au format pour les OS.
      CHARACTER*SOS_LENFIC2 FICSURF2      ! nom du 2nd fichier de reflexion 
      					  ! de surface, au format pour les OS.
      CHARACTER*SOS_LENFIC2 FICS          ! nom du fichier resultat de reflexion 
      					  ! de surface, au format pour les OS.
					  
      	! Elements Pij de la matrice de reflexion en fonction de la direction
	! d'incidence et de reflexion, pour un ordre IS du developpement
	! en series de Fourier. Format d'entree des OS pour le 1er fichier.
      REAL P11_SURF1(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P12_SURF1(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P13_SURF1(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)	
      REAL P21_SURF1(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P22_SURF1(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P23_SURF1(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P31_SURF1(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P32_SURF1(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P33_SURF1(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
           
        ! Elements Pij de la matrice de reflexion en fonction de la direction
	! d'incidence et de reflexion, pour un ordre IS du developpement
	! en series de Fourier. Format d'entree des OS pour le 2nd fichier.
      REAL P11_SURF2(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P12_SURF2(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P13_SURF2(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)	
      REAL P21_SURF2(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P22_SURF2(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P23_SURF2(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P31_SURF2(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P32_SURF2(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P33_SURF2(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)   
          
        ! Elements Pij de la matrice de reflexion en fonction de la direction
	! d'incidence et de reflexion, pour un ordre IS du developpement
	! en series de Fourier. Format d'entree des OS pour le fichier resultat.
      REAL P11(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P12(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P13(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)	
      REAL P21(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P22(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P23(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P31(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P32(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL P33(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)

      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB 	! Ordre des developpements des fonctions de phase 
			! en Polynomes de Legendre.

      INTEGER*4 IS	! Ordre de la decomposition en series de Fourier.
      INTEGER*4 I	! Indice de l'angle d'incidence Teta1.
      INTEGER*4 J	! Indice de l'angle de reflexion Teta2.
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon	  					  					  
      	      
      

C* Ouverture des fichiers
C------------------------------------------------------------		    
      OPEN(UNIT=1,FILE=FICSURF1,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=2,FILE=FICSURF2,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=3,FILE=FICS,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=991)

	   
C* Boucle sur l'ordre de la decomposition en series de Fourier
C------------------------------------------------------------      
      DO 1 IS=0,OS_NB

C*      Lecture des fichiers qui ont ete enregistres avec des 
C*      bornes de tableaux limitees au dimensionnement 
C*      utile -NBMU:NBMU           
	   READ(1,err=992),
     &	       ((P11_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P12_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P13_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P21_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P22_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P23_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P31_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P32_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P33_SURF1(I,J),I=1,NBMU),J=1,NBMU)  
                 
	   READ(2,err=992),
     &	       ((P11_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P12_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P13_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P21_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P22_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P23_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P31_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P32_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P33_SURF2(I,J),I=1,NBMU),J=1,NBMU) 

           
C* Boucle sur les directions d'incidence Teta1
C---------------------------------------------
         DO 2 I=1,NBMU

C* Boucle sur les directions de reflexion Teta2
C-----------------------------------------------
            DO 3 J=1,NBMU
               P11(I,J) = P11_SURF1(I,J) + P11_SURF2(I,J)  
               P12(I,J) = P12_SURF1(I,J) + P12_SURF2(I,J)  
               P13(I,J) = P13_SURF1(I,J) + P13_SURF2(I,J)  
               P21(I,J) = P21_SURF1(I,J) + P21_SURF2(I,J)
               P22(I,J) = P22_SURF1(I,J) + P22_SURF2(I,J)
               P23(I,J) = P23_SURF1(I,J) + P23_SURF2(I,J)  
               P31(I,J) = P31_SURF1(I,J) + P31_SURF2(I,J)  
               P32(I,J) = P32_SURF1(I,J) + P32_SURF2(I,J)  
               P33(I,J) = P33_SURF1(I,J) + P33_SURF2(I,J)  
 3          CONTINUE
 2       CONTINUE
 
C*       Enregistrement des tableaux de la matrice de reflexion 
C*       a l'ordre IS de la decomposition en series de Fourier.
C*       Bornes d'enregistrement des tableaux limitees au 
C*       dimensionnement utile -NBMU:NBMU
C-------------------------------------------------------------- 
	 WRITE(3,err=992) ((P11(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P12(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P13(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P21(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P22(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P23(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P31(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P32(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P33(I,J),I=1,NBMU),J=1,NBMU)
     
 1    CONTINUE


C* Fin nominale 
C-------------------
      IER=0
      GOTO 9999
 
 
C* Cas d'erreur : retour du status -1 au programme appelant
C----------------------------------------------------------
  991 WRITE(6,*) '  ERROR on a file opening for SOS_BPDF_AJOUT_BRDF'
      IER=-1
      goto 9999
  992 WRITE(6,*) '  ERROR on a file reading/writing for ',
     &           'SOS_BPDF_AJOUT_BRDF'
      IER=-1 
      goto 9999

 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      
      RETURN            
      END      !FIN DE LA PROCEDURE SOS_BPDF_AJOUT_BRDF 
      
