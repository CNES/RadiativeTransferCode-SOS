C******************************************************************************
C* Copyright 2023, Centre National d'Etudes Spatiales (CNES)
C* 
C* This file is part of the SOS-ABS radiative transfer code.
C* 
C* SOS-ABS is free software: you can redistribute it and/or modify
C* it under the terms of the GNU General Public License as published by
C* the Free Software Foundation, either version 3 of the License, or
C* (at your option) any later version.
C* 
C* SOS-ABS is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C* GNU General Public License for more details.
C* 
C* You should have received a copy of the GNU General Public License
C* along with SOS-ABS. If not, see <http://www.gnu.org/licenses/>.
C******************************************************************************


C******************************************************************************
C* FICHIER: SOS_SURFACE.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module genere les matrices de réflexion pour différents types de
C*	 modèles.  
C*       
C*
C* AUTEUR: 
C*   Fondement algorithmique: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, B.Lafrance, X.Lenot
C*
C* DATE: creation 1.0 : 30/04/2002
C*
C* MOD:VERSION:1.1: Mise a niveau des codes du LOA (codes relatifs au Glitter).
C*                  Développements d'adaptation des codes aux autres modèles.
C* MOD:VERSION:2.0: Correction d'une erreur sur la polarisation de surface
C*                  pour la routine SOS_MAT_REFLEXION.
C* MOD:VERSION:2.1: 07/07/2008 
C*                  Modification du format de lecture du paramètre WIND,
C*                  pour correction d'une erreur de décodage si valeur entière.
C*			WIND passe de F4.1 en F4
C* MOD:VERSION:3.0: 17/11/2009 
C*                  * Modification du passage des arguments : 
C*                      - passage par couples (mot clef, valeur). 
C*                      - adaptation de la gestion des cas d'erreurs. 
C*                  * Conversion des messages d'erreurs en anglais. 
C* MOD:VERSION:3.1: 22/12/2009 
C*                  Modification du format de lecture des paramètres Alpha et Beta du modèle Nadal
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
C*     ==> Passage de CTE_OS_NBMU à CTE_OS_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'à LUM_NBMU égal au nombre effectif d'angles utilisés.
C*
C*   - Suppression des constantes CTE_OS_NB, CTE_OS_NS et CTE_OS_NM (de SOS.h)  
C*     au profit : d'un passage par argument des valeurs OS_NB, OS_NS, OS_NM (lus dans le fichier des angles)
C*
C*   - Modification de la taille des tableaux de decomposition en polynomes de Legendre :
C*     ==> Passage de CTE_OS_NB à CTE_OS_NB_MAX.
C*     ==> Passage de CTE_OS_NS à CTE_OS_NS_MAX.
C*     ==> Passage de CTE_OS_NM à CTE_OS_NM_MAX.
C*     ==> Boucle sur les ordre jusqu'à OS_NB, OS_NS ou OS_NM 
C*
C*   - Routines SOS_MISE_FORMAT & SOS_BPDF_AJOUT_BRDF : Modification de la  
C*     ligne d'enregistrement ou de lecture des tableaux de la  matrice de 
C*     réflexion a l'ordre IS de la decomposition en séries de Fourier :
C*     ==> limitation aux bornes utiles -LUM_NBMU:LUM_NBMU
C*       Rem : pas de modification des bornes d'enregistrement des fichiers 
C*             temporaire RES_FRESNEL, SOS_MAT_REFLEX car il ne s'agit pas 
C*             de produits finaux.
C*
C*   - Adaptation de l'appel des routines de denomination des fichiers
C*     (SOS_NOM_FICGLITTER, SOS_NOM_FICROUJ, SOS_NOM_FICRH, SOS_NOM_FICBREON, SOS_NOM_FICNADAL)
C*     aux nouveaux paramètres NBMU_GAUSS et FICANGLES_USER
C*
C*   - Adaptation de la routine SOS_BPDF_AJOUT_BRDF au passage du paramètre
C*     LUM_NBMU (nb d'angles utiles) et de OS_NB
C*
C*   - Adaptation de l'appel des routines SOS_GLITTER, SOS_ROUJEAN, SOS_RONDEAUX_BREON, 
C*     SOS_NADAL et SOS_BPDF_AJOUT_BRDF aux passages des paramètres
C*     LUM_NBMU (nb d'angles utiles), RMU (table des cosinus des angles) 
C*     et CHR (table des poids des angles), ainsi que OS_NB, OS_NS et OS_NM
C*
C* MOD:VERSION:4.1: 24/06/2010 
C*     Les etiquettes de FORMAT en ,X, passent en ,1X, pour eviter des erreurs de compilation.
C*
C* MOD:VERSION:4.2: 06/06/2013
C*
C*   -  Les paramètres d'entrée sont fournis sous la forme de liste de paramètres
C*      et non plus sous la forme de couples avec mot cle associe => tous les
C*	paramètres deviennent donc obligatoires
C*
C*   -  Modification de la valeur de Pi, declaree en constante dans le fichier
C*      (INCTE_PI remplace CTE_PI)
C*  
C* MOD:VERSION:4.3: 24/02/2016 
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*
C*    - Changement de la manière de lire les fichiers formatés (Angles, ...)
C*      pour être compatible avec le compilateur gfortran.
C* 
C*  
C* MOD:VERSION:5.0: 12/08/2020
C*   - Ajout des routines de calcul de la matrice BPDF de Maignan (avec BRDF de Roujean):
C*     --> Appel de la routine SOS_NOM_FICMAIGNAN (définie dans le fichier SOS_NOM_FIC_SURFACE.F)
C*     --> Introduction du nouveau paramètre d'entrée: COEF_C_MAIGNAN qui intervient dans le modèle de Maignan.
C*
C*   - Appel de la nouvelle routine SOS_SURFACE_BPDF au lieu de SOS_RONDEAUX_BREON et SOS_NADAL (qui sont supprimées)
C*     pour gérer le calcul des matrices de BRDF/BPDF de surface (modèles de Rondeaux, Bréon, Nadal et Maignan).
C*   
C*   - Introduction d'un répertoire DIR_TMP pour la localisation des fichiers temporaires.
C*
C*   -  Passage du nom des fichiers temporaires (RES_GSF, RES_FRESNEL, RES_MAT_REFLEX, FICTMP_BPDF),
C*      avec arborescence complète, aux sous-routines SOS_GLITTER, SOS_ROUJEAN, SOS_SURFACE_BPDF, 
C*      SOS_MAT_FRESNEL et SOS_MAT_REFLEXION.
C*
C*   - Modification des seuils pour les tests de valeurs limites sur les coefficients KO, K1 et K2 
C*     du modèle de Roujean : passe de 1000 à 10 
C*     (quelques cas avec K2 supérieurs à 1 selon l'article de ROUJEAN ET AL (A Bidirectional Reflectance Model 
C*      of the Earth's Surface for the Correction of Remote Sensing Data, JOURNAL OF GEOPHYSICAL RESEARCH,1992).
C*
C* MOD:VERSION:5.1: 31/08/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Conversion des traces en Anglais.
C*
C* MOD:VERSION:5.2: 26/09/2023
C*    - Suppression d'une virgule après la commande READ(...) en deux endroits.
C*      Impact nul sur les résultats de simulations, mais levant un warning à la compilation.
C*    - Modification de syntaxes pour compatibilité avec Fortran 2018 pour compilateur gcc/10.2.0:
C*         - Boucle de complétude à zéro dans la routine SOS_MAT_REFLEXION 
C*         - Boucles d'initialisation pour l'ordre IS >= 2 dans la routine SOS_NOYAUX_FRESNEL
C*         - Boucles sur Teta1 et Teta2 dans la routine SOS_NOYAUX_FRESNEL
C*
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Définition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    CTE_LENDIR : Longueur des noms de répertoires.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_LENKEYWORD : Taille maximale des Keywords pour le passage des arguments
C    INCTE_PI : Valeur de PI
C    CTE_OS_NBMU_MAX : Nombre maximal d'angles positifs.
C    CTE_OS_NB_MAX :  Valeur maximale pour l'ordre limite du développement en 
C                     polynomes de Legendre.
C    CTE_OS_NS_MAX :  Valeur maximale pour l'ordre limite du développement en polynomes 
C                     de Legendre des éléments de la matrice de Fresnel 
C                     et en séries de Fourier des luminances.
C    CTE_OS_NM_MAX :  Valeur maximale pour l'ordre limite  du développement 
C                     en séries de Fourier de la fonction G.
C
C Constantes spécifiques :
C    LENLINE  : Longueur maximale d'une ligne de données lue dans un fichier
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define LENLINE 1000

C==============================================================================
C PROGRAMME: SOS_SURFACE
C ==========  
C      Ce programme calcule les matrice de réflexion (BRDF et/ou BPDF) 
C      de différents modèles:
C	- Réflexion du glitter,
C	- Modèle de BRDF de Roujean pour simuler la réflexion en intensite 
C        de couverts vegetaux,
C	- Modèle simplifie de BPDF de Rondeaux et Herman (RH) pour simuler la 
C	  réflexion en polarisation des couverts vegetaux,
C	- Modèle simplifie de BPDF de Breon, pour simuler la réflexion
C        en polarisation des sols nus,
C	- Modèle de BPDF de Nadal dont l'intérêt est de permettre de simuler 
C	  la polarisation de différents types de surface, d'une manière 
C	  simplifiee mais représentative.
C        ==> PB D'APPLICABILITE AU FORMALISME DE SOS
C            MODELE NON UTILISABLE.
C	- Modèle de BPDF de Maignan représentatif des surfaces désertiques. 
C
C       Aux 4 derniers modèles de BPDF (Rondeaux, Breon, Nadal, Maignan), 
C       on ajoute la BRDF de Roujean.
C
C Donnees en entrée du programme
C ------------------------------
C     -->Les paramètres/arguments du programme
C
C	      ISURF  (I4) : (E) Type de surface ISURF 
C			       (format de decodage en entier I4)
C      		0 :surface lambertienne d'albédo RHO.
C			1 :surface lambertienne + glitter.
C			2 :surface lambertienne + réflexion de Fresnel.
C		       3 :surface lambertienne + BRDF de Roujean.
C			4 :surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
C			5 :surface lambertienne + BRDF de Roujean + BPDF de Breon.
C			6 :surface lambertienne + BRDF de Roujean + BPDF de Nadal.	
C			7 :surface lambertienne + BRDF de Roujean + BPDF de Maignan.		 	
C
C             FICANGLES (CHARACTER) : (E) nom complet du fichier des angles utilisés 
C                                      pour les calculs de luminance et BRDF
C                 		      (répertoire + nom fichier + extension)
C
C             IND (I4) : (E) Indice de réfraction de la surface / air.		 	
C
C             WIND (DOUBLE) : (E) ==> Associé à : ISURF = 1
C 	      		       Vitesse du vent en m/s.	
C		
C             K0_ROUJEAN (DOUBLE) : (E) ==> Associé à : ISURF >= 3
C 	       		     Paramètre k0 du modèle de BRDF de Roujean.       	
C		
C             K1_ROUJEAN (DOUBLE) : (E) ==> Associé à : ISURF >= 3
C 	      		    Paramètre k1 du modèle de BRDF de Roujean.       	
C		
C             K2_ROUJEAN (DOUBLE) : (E) ==> Associé à : ISURF >= 3
C 	      		    Paramètre k2 du modèle de BRDF de Roujean.       	
C
C             ALPHA_NADAL (DOUBLE) : (E) ==> Associé à : ISURF = 6
C 	      			Paramètre alpha du modèle de BPDF de Nadal.       	
C
C             BETA_NADAL (DOUBLE) : (E) ==> Associé à : ISURF = 6
C 	      		       Paramètre beta du modèle de BPDF de Nadal.       
C
C             COEF_C_MAIGNAN (DOUBLE) : (E) ==> Associé à : ISURF = 7
C 	      		       Paramètre C*exp(-NDVI) du modèle de BPDF de Maignan.  	
C
C             FICLOG (CHARACTER): (E) nom complet du fichier trace (répertoire + nom fichier + extension)
C                                     ou =0, si pas de  trace	 
C
C             DIRS (CHARACTER) : (E) Répertoire de stockage des fichiers de matrice de réflexion 
C                   		        de surface (fichier Glitter ou autres BRDF/BPDF).
C
C             DIRROUJEAN (CHARACTER) : (E) ==> Associé à : ISURF = 3 4 5 6 ou 7
C					répertoire de stockage des fichiers de BRDF de Roujean,
C                                  appliquées aux modèles de BPDF.
C
C             DIR_TMP  (CHARACTER)  : (E) répertoire temporaire de génération des fichiers
C
C      	IER (I4) : (S) code d'erreur =0 si pas d'erreur, =-1 sinon
C
C                                 
C
C    -->Fichiers
C       Le fichier Trace.
C       Le fichier résultat contenant l'expression, pour le code des OS, de la 
C       decomposition en séries de Fourier de la matrice de réflexion de surface.
C       Des fichiers intermédiaires (non conservés) : 
C           commentaires et description dans les routines appelées.
C       Le fichier des angles et paramètres internes pour les limites des développements.
C         
C       
C
C Resultats fournis par le programme
C ----------------------------------
C fichier BRDF / BPDF de la matrice de réflexion de surface
C pour une décomposition en séries de Fourier
C
C Variable d'environnement
C ------------------------
C   Aucune
C
C
C Description des fichiers utilisés
C ---------------------------------
C
C   -->Contenu du fichier BRDF / BPDF résultat
C       Fichier binaire non formatte, cree par SOS_MISE_FORMAT
C       Il contient la matrice de réflexion de surface 
C       (fichier GLITTER ou autres BRDF / BPDF):
C       c'est a dire les éléments Pij de la matrice de réflexion 
C              P11(LUM_NBMU,LUM_NBMU) (réel),
C              P12(LUM_NBMU,LUM_NBMU) (réel),
C              P13(LUM_NBMU,LUM_NBMU) (réel),
C              P21(LUM_NBMU,LUM_NBMU) (réel),
C              P22(LUM_NBMU,LUM_NBMU) (réel),
C              P23(LUM_NBMU,LUM_NBMU) (réel),
C              P31(LUM_NBMU,LUM_NBMU) (réel),
C              P32(LUM_NBMU,LUM_NBMU) (réel),
C              P33(LUM_NBMU,LUM_NBMU) (réel)
C
C        --> Enregistrement sur le domaine limité aux angles utiles
C                                  --------------------------------
C
C
C    -->Contenu du fichier des angles:
C       Fichier ascii formatté, crée par SOS_ANGLES.
C       Il contient l'information sur les angles à utiliser 
C       et sur les ordres des développements adaptes : 
C		  NB_TOTAL_ANGLES :  Nb d'angles total à utiliser
C               NB_GAUSS_ANGLES :  Nb d'angles de Gauss
C               ANGLES_USERFILE :  Nom du fichier des angles utilisateur
C	         SOLAR ZENITH ANGLE : Angle zenithal solaire (degrés)
C		  INTERNAL_IMUS :      Numéro d'indice de table correspondant
C               INTERNAL_OS_NB :   Ordre des développements des fonctions de phase 
C				       en Polynomes de Legendre à utiliser.
C               INTERNAL_OS_NS :   Ordre maximal à utiliser pour le développement en
C                                  polynomes de Legendre des éléments de la matrice de Fresnel 
C                                  et en séries de Fourier des calculs de luminance.
C               INTERNAL_OS_NM :   Ordre maximal du développement en séries de Fourier  
C                                  de la fonction G des calculs de matrice de réflexion.
C               INDEX COS_ANGLE WEIGHT OUTPUT: lignes donnant le numéro d'angle, son cosinus 
C                                              et son poids, s'il est donné en sortie "utilisateur"
C                                              au format I4,1X,2D21.14,1X,I4
C         
C
C Common utilisé:
C --------------
C  Aucun
C
C
C Cas d'erreur :
C ------------        
C     -Erreur dans un sous-programme
C     -Erreur à l'ouverture du fichier Trace
C     -Paramètres non valables
C
C     Affichage d'un message à l'écran, arrêt du programme et 
C     retour du status 1 au shell
C   
C
C  Sous programmes utilisés:
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
      SUBROUTINE SOS_SURFACE(ISURF, FICANGLES, IND, WIND,  
     &                       K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN,   
     &                       ALPHA_NADAL, BETA_NADAL,   
     &                       COEF_C_MAIGNAN, FICLOG,                
     &                       DIRS, DIRROUJEAN, DIR_TMP, IER)
      
      IMPLICIT NONE
      
C* Définition des variables   
C*-------------------------------------------------------------------------- 
         
      CHARACTER*CTE_LENDIR DIRS	     ! Chemin d'accès au fichier SURFACE résultat .
      CHARACTER*CTE_LENDIR DIRROUJEAN   ! Répertoire de stockage des fichiers de BRDF 
                                        ! de Roujean, appliquées aux modèles de BPDF.  
      CHARACTER*CTE_LENDIR DIR_TMP      ! Chemin d'accès aux fichiers SURFACE en cours de génération                                             
      CHARACTER*CTE_LENFIC1 FICS        ! Nom du fichier SURFACE (sans arborescence).
      CHARACTER*CTE_LENFIC2 FICLOG      ! Nom du fichier trace SURFACE.
      CHARACTER*CTE_LENFIC2 FICANGLES   ! Chemin et nom du fichier des angles effectifs
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! Chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*CTE_LENFIC2 FICSURF_TMP     ! Chemin et nom du fichier SURFACE en cours de génération
      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Nom du fichier RES_MAT_REFLEX
                                               ! avec arborescence complète.
      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Nom du fichier RES_GSF avec arborescence complète.
      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Nom du fichier RES_FRESNEL avec arborescence complète.
      CHARACTER*CTE_LENFIC2 FICGLITTER ! Nom complet du fichier Glitter.
      CHARACTER*CTE_LENFIC2 FICROUJ    ! Nom complet du fichier du modèle de ROUJEAN.
      CHARACTER*CTE_LENFIC2 FICRH      ! Nom cpmplet du fichier du modèle de RH.
      CHARACTER*CTE_LENFIC2 FICBREON   ! Nom complet du fichier du modèle de BREON.
      CHARACTER*CTE_LENFIC2 FICNADAL   ! Nom complet du fichier du modèle de NADAL.
      CHARACTER*CTE_LENFIC2 FICMAIGNAN ! Nom complet du fichier du modèle de NADAL.
      CHARACTER*CTE_LENFIC2 FICTMP_BPDF ! Nom complet du fichier temporaire pour la BPDF sans complément de la BRDF. 
      CHARACTER*CTE_LENCOM  COMMAND     ! Chaine de commande pour appel système.
  
      CHARACTER*LENLINE LIGNE_TEXTE    ! Ligne de caractères pour décodage de paramètres   
      CHARACTER*1000 CVAL	       ! Chaine de caractères pour la lecture de valeurs inutilisees
      
      LOGICAL EX		! Verifie existence d'un fichier.
      LOGICAL TRACE     	! = vrai,si écriture dans le fichier trace

      INTEGER*4 NBMU_GAUSS  ! Nombre d'angles de Gauss utilisés       
      INTEGER*4 LUM_NBMU    ! Nombre d'angles effectifs des simulations
      INTEGER*4 LDIRS	! Longueur de la chaine de caractère DIRS.
      INTEGER*4 LDIRTMP     ! Longueur de la chaine de caractères DIRTMP 
      INTEGER*4 LDIRROUJEAN	! Longueur de la chaine de caractère DIRROUJEAN
      INTEGER*4 LFICS	! Longueur de la chaine de caractère FICS.  
      INTEGER*4 LFICTMP	! Longueur de la chaine de caractères FICSURF_TMP          
      INTEGER*4 LCOM        ! Longueur de la chaine de caractères COMMAND      
      INTEGER*4 IER     	! Code d'erreur =0 si pas d'erreur, =-1 sinon     
      INTEGER*4 J	       ! Indice de tableau                    
      INTEGER*4 OS_NB 	! Ordre des développements des fonctions de phase 
				! en Polynomes de Legendre.
      INTEGER*4 OS_NS    	! Ordre maximal du développement en polynomes de Legendre   
				! des éléments de la matrice de Fresnel.
      INTEGER*4 OS_NM    	! Ordre maximal du développement en séries de Fourier  
				! de la fonction G des calculs de matrice de réflexion.
  
      INTEGER*2 ISURF	! Type de surface          

      INTEGER*4 IFIN        ! Nb de caractères significatifs d'une chaine. 
 

      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                             !Cosinus des angles de Gauss (et utilisateur).
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                             !Poids des angles de Gauss.				
          
      DOUBLE PRECISION WIND	! Vitesse du vent (m/s).
      
      DOUBLE PRECISION IND	! Indice de réfraction du dioptre.	
      
      DOUBLE PRECISION K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN	!Paramètres du modèle de Roujean
      
      DOUBLE PRECISION ALPHA_NADAL	!Paramètre alpha du modèle de Nadal.
      DOUBLE PRECISION BETA_NADAL	!Paramètre beta du modèle de Nadal.
      DOUBLE PRECISION COEF_C_MAIGNAN !Paramètre C*exp(-NDVI) du modèle de Maignan.
   
      
C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0
     
C --- Dénomination du fichier de BRDF/BPDF temporaire 
C     pour sa création dans l'espace de travail avant déplacement 
C     dans l'espace d'archivage.
      LDIRTMP = INDEX(DIR_TMP,'/TMP') +3
      FICSURF_TMP=DIR_TMP(1:LDIRTMP)//'/'//"SURFACE_FILE_BUILDING_TMP"
      LFICTMP = INDEX(FICSURF_TMP,' ') - 1

C     Fichier temporaire RES_GSF
      FIC_RES_GSF = DIR_TMP(1:LDIRTMP)//'/'//"RES_GSF"
      FIC_RES_FRESNEL = DIR_TMP(1:LDIRTMP)//'/'//"RES_FRESNEL"
      FIC_RES_MAT_REFLEX = DIR_TMP(1:LDIRTMP)//'/'//"RES_MAT_REFLEX"
      FICTMP_BPDF = DIR_TMP(1:LDIRTMP)//'/'//"FICTMP_BPDF"

C     Suppression des fichiers temporaires potentiellement existants
      INQUIRE(FILE=FICSURF_TMP,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF 

      INQUIRE(FILE=FIC_RES_GSF,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FIC_RES_GSF,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF

      INQUIRE(FILE=FIC_RES_FRESNEL,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FIC_RES_FRESNEL,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF

      INQUIRE(FILE=FIC_RES_MAT_REFLEX,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF

      INQUIRE(FILE=FICTMP_BPDF,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF


C* Ouverture du fichier Trace 
C------------------------------------------------------
      IF (FICLOG.EQ.'NO_LOG_FILE') THEN
	  TRACE=.FALSE.
      ELSE
         TRACE=.TRUE.
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
         READ(10,*,ERR=951) CVAL, RMU(J), CHR(J), CVAL
      ENDDO 

      DO J=1,LUM_NBMU
         RMU(-J) = -RMU(J)
	 CHR(-J) =  CHR(J)
      ENDDO 

	
      
      IF (TRACE) THEN
        WRITE (99,*,ERR=960) "USED ANGLES"
	 WRITE (99,*,ERR=960) "***********"
	 WRITE (99,*,ERR=960) " "
	 WRITE(99,*,ERR=960) "Total number of angles : ", LUM_NBMU
	 WRITE(99,*,ERR=960) "Number of Gauss angles : ", NBMU_GAUSS
	 WRITE (99,*,ERR=960) " "
        IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	     WRITE(99,*,ERR=960)  
     &	       "A file of user's angles is used :",
     &         FICANGLES_USER
     	     WRITE(99,*,ERR=960) "Number of user's angles :",
     &         (LUM_NBMU-NBMU_GAUSS)
	 ELSE
	     WRITE(99,*,ERR=960)  
     &	       "No file of user's angles"
	 ENDIF
        WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,ERR=960) "Max order NB ", OS_NB
	 WRITE (99,*,ERR=960) "Max order NS ", OS_NS
	 WRITE (99,*,ERR=960) "Max order NM ", OS_NM
	 WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,err=960) "List of used angles :"
	 WRITE (99,*,err=960) "   Index, mu, angle (degrees), weigth"
	 DO 1000 J=-LUM_NBMU,LUM_NBMU
	    IF (J.EQ.0) GOTO 1000
	    IF (J.LT.0) THEN 
	        WRITE(99,222,ERR=960) J, RMU(J),
     &	                  -DACOS(-RMU(J))*180./INCTE_PI, CHR(J)
            ELSE
	        WRITE(99,222,ERR=960) J, RMU(J),
     &	                   DACOS(RMU(J))*180./INCTE_PI, CHR(J)
	    ENDIF
1000    CONTINUE	
	 WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,ERR=960) " "  	
      ENDIF !Fichier Trace	 


c* Taille du chemin complet au répertoire de sortie
C* -----------------------------------------------------------------     

      LDIRS = INDEX(DIRS,' ')
      LDIRS = LDIRS - 1
      IF(LDIRS.LE.0) LDIRS = CTE_LENDIR

      LDIRROUJEAN = INDEX(DIRROUJEAN,' ')
      LDIRROUJEAN = LDIRROUJEAN - 1
      IF(LDIRROUJEAN.LE.0) LDIRROUJEAN = CTE_LENDIR


C =====================================================================
C* Cas du GLITTER
C =====================================================================

      IF (ISURF.EQ.1) THEN

	
C* Constitution du nom du fichier GLITTER
C  --------------------------------------
	 CALL SOS_NOM_FICGLITTER(NBMU_GAUSS,FICANGLES_USER,IND,WIND,
     &	                         FICS)      
	 FICGLITTER = DIRS(1:LDIRS)//'/'//FICS	
	 LFICS=INDEX(FICGLITTER,' ')-1	
	 WRITE(6,*) "File : ",FICGLITTER(1:LFICS)
	 
C* Si le fichier GLITTER existe déjà, le programme n'est pas lancé
C  ---------------------------------------------------------------     
         INQUIRE(FILE=FICGLITTER(1:LFICS),EXIST=EX)
	  
         IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE   
      
C* Calcul du fichier Glitter
C  ------------------------- 
           WRITE(6,*)'==> Starting the file computation'
	    CALL SOS_GLITTER(LUM_NBMU,RMU,CHR,WIND,IND,
     &	                     OS_NB,OS_NS,OS_NM,
     &                      FIC_RES_GSF,FIC_RES_FRESNEL,
     &                      FIC_RES_MAT_REFLEX,
     &	                     FICSURF_TMP,TRACE,IER)
            IF (IER.NE.0) GOTO 993   ! cas d'erreur

C* Déplacement de FICSURF_TMP vers FICGLITTER
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICGLITTER(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICGLITTER(1:LFICS),EXIST=EX)
            IF (EX) THEN !Si le fichier de BRDF/BPDF existe on détruit le fichier temporaire (généré par autre traitement)
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de BRDF/BPDF
               CALL SYSTEM(COMMAND)  
            ENDIF    
         ENDIF

	  
      ENDIF ! Cas du Glitter
	 



C =====================================================================
C* Cas BRDF de ROUJEAN
C*   NB : Calculée également pour les modèles de BPDF.
C =====================================================================

      IF (ISURF.GE.3) THEN


C* Controle de cohérence par rapport au format F7.3
C  Valeurs de K0, K1,K2 forcées au format F7.3
C  -----------------------------------------	
        IF((K0_ROUJEAN.LT.0.).OR.(K0_ROUJEAN.GE.10.)) GOTO 9910
	 IF((K1_ROUJEAN.LT.0.).OR.(K1_ROUJEAN.GE.10.)) GOTO 9910
	 IF((K2_ROUJEAN.LT.0.).OR.(K2_ROUJEAN.GE.10.)) GOTO 9910
	
	 K0_ROUJEAN = DNINT(K0_ROUJEAN*1000.D+00) / 1000.D+00
	 K1_ROUJEAN = DNINT(K1_ROUJEAN*1000.D+00) / 1000.D+00
	 K2_ROUJEAN = DNINT(K2_ROUJEAN*1000.D+00) / 1000.D+00	 	 

	
C* Constitution du nom du fichier de ROUJEAN
C  -----------------------------------------
	 CALL SOS_NOM_FICROUJ(NBMU_GAUSS,FICANGLES_USER,
     &                       K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS)     
	 FICROUJ = DIRROUJEAN(1:LDIRROUJEAN)//'/'//FICS 
	 LFICS=INDEX(FICROUJ,' ')-1      
	 WRITE(6,*) "File : ",FICROUJ(1:LFICS)
	 
C* Si le fichier existe déjà, le programme n'est pas lancé
C  --------------------------------------------------------	 
         INQUIRE(FILE=FICROUJ(1:LFICS),EXIST=EX)	 
	 
         IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	    
	 ELSE   
      
C* Calcul du fichier de ROUJEAN
C  -----------------------------
            WRITE(6,*)'==> Starting the file computation'
	     CALL SOS_ROUJEAN(LUM_NBMU,RMU,CHR,OS_NB,
     &	                      K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &	                      FIC_RES_MAT_REFLEX,FICSURF_TMP,TRACE,IER)

            IF (IER.NE.0) GOTO 994   ! cas d'erreur

C* Déplacement de FICSURF_TMP vers FICROUJ
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICROUJ(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICROUJ(1:LFICS),EXIST=EX)
            IF (EX) THEN !Si le fichier de BRDF/BPDF existe on détruit le fichier temporaire (généré par autre traitement)
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de BRDF/BPDF
               CALL SYSTEM(COMMAND)  
            ENDIF
	       
         ENDIF
	 
      ENDIF ! Cas de la BRDF de Roujean



	 
C =====================================================================
C* Cas BPDF de Rondeaux - Herman
C =====================================================================

      IF (ISURF.EQ.4) THEN

C* Constitution du nom du fichier de Rondeaux - Herman
C  ---------------------------------------------------
	 CALL SOS_NOM_FICRH(NBMU_GAUSS,FICANGLES_USER,IND,
     &	                    K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS)
	 FICRH = DIRS(1:LDIRS)//'/'//FICS	 
	 LFICS=INDEX(FICRH,' ')-1
	 WRITE(6,*) "File : ",FICRH(1:LFICS)
	 
C* Si le fichier existe déjà, le programme n'est pas lancé
C  -------------------------------------------------------	 
        INQUIRE(FILE=FICRH(1:LFICS),EXIST=EX)
	 	 
        IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE   
  
            
C* Calcul du fichier de Rondeaux - Herman
C  --------------------------------------
           WRITE(6,*)'==> Starting the file computation'
	    CALL SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                           ISURF,ALPHA_NADAL,BETA_NADAL,  
     &                           COEF_C_MAIGNAN,
     &	                          OS_NB,OS_NS,OS_NM,
     &                           FIC_RES_GSF,FIC_RES_FRESNEL,
     &                           FIC_RES_MAT_REFLEX,
     &	                          FICTMP_BPDF,TRACE,IER)     
            IF (IER.NE.0) GOTO 995   ! cas d'erreur


		
C* Ajout de la BRDF de Roujean
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP_BPDF,FICROUJ,LUM_NBMU,OS_NB,
     &	                              FICSURF_TMP,IER)
	     IF (IER.NE.0) GOTO 998   ! cas d'erreur

C* Déplacement de FICSURF_TMP vers FICRH
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICRH(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICRH(1:LFICS),EXIST=EX)
            IF (EX) THEN !Si le fichier de BRDF/BPDF existe on détruit le fichier temporaire (généré par autre traitement)
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de BRDF/BPDF
               CALL SYSTEM(COMMAND)  
            ENDIF


C* Destruction du fichier temporaire
C---------------------------------------   
	     OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	    	                   
	        	       
        ENDIF ! Si fichier à créer.

      ENDIF ! Cas de la BPDF de Rondeaux - Herman



C =====================================================================
C* Cas BPDF de BREON
C =====================================================================

      IF (ISURF.EQ.5) THEN

C* Constitution du nom du fichier de BREON
C  ---------------------------------------
	 CALL SOS_NOM_FICBREON(NBMU_GAUSS,FICANGLES_USER,IND,
     &	                       K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS)   
	 FICBREON = DIRS(1:LDIRS)//'/'//FICS
	 LFICS=INDEX(FICBREON,' ')-1
	 WRITE(6,*) FICBREON(1:LFICS)
	 
C* Si le fichier existe déjà, le programme n'est pas lancé
C  --------------------------------------------------------	 
        INQUIRE(FILE=FICBREON(1:LFICS),EXIST=EX)	 
	 
        IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE        
            
C* Calcul du fichier de Rondeaux - Herman
C  --------------------------------------
            WRITE(6,*)'==> Starting the file computation'
	     CALL SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                            ISURF,ALPHA_NADAL,BETA_NADAL,  
     &                            COEF_C_MAIGNAN,
     &	                           OS_NB,OS_NS,OS_NM,
     &                            FIC_RES_GSF,FIC_RES_FRESNEL,
     &                            FIC_RES_MAT_REFLEX,
     &	                           FICTMP_BPDF,TRACE,IER)  
            IF (IER.NE.0) GOTO 995   ! cas d'erreur

C* Ajout de la BRDF de Roujean
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP_BPDF,FICROUJ,LUM_NBMU,OS_NB,
     &                               FICSURF_TMP,IER)
	     IF (IER.NE.0) GOTO 998   ! cas d'erreur

C* Déplacement de FICSURF_TMP vers FICBREON
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICBREON(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICBREON(1:LFICS),EXIST=EX)
            IF (EX) THEN !Si le fichier de BRDF/BPDF existe on détruit le fichier temporaire (généré par autre traitement)
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de BRDF/BPDF
               CALL SYSTEM(COMMAND)  
            ENDIF
	    
C* Destruction du fichier temporaire
C---------------------------------------   
	     OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	       
        ENDIF ! Si fichier à créer.

      ENDIF ! Cas de la BPDF de Breon 



C =====================================================================
C* Cas BPDF de NADAL
C =====================================================================

      IF (ISURF.EQ.6) THEN	

		
C* Controle de cohérence par rapport aux formats F6.4 et F4.1
C  Valeurs de alpha et beta forcées à ces formats
C  -----------------------------------------	
        IF((ALPHA_NADAL.LT.0.).OR.(ALPHA_NADAL.GE.10.)) GOTO 9915
	 IF((BETA_NADAL.LT.0.).OR.(BETA_NADAL.GE.100.)) GOTO 9915

	 ALPHA_NADAL = DNINT(ALPHA_NADAL*10000.D+00) / 10000.D+00
	 BETA_NADAL = DNINT(BETA_NADAL*10.D+00) / 10.D+00	
		
C* Constitution du nom du fichier de NADAL
C  ---------------------------------------
	 CALL SOS_NOM_FICNADAL(NBMU_GAUSS,FICANGLES_USER,IND,
     &	                       ALPHA_NADAL,BETA_NADAL,
     &	                       K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS) 	        
	 FICNADAL = DIRS(1:LDIRS)//'/'//FICS	 
	 LFICS=INDEX(FICNADAL,' ')-1
	 WRITE(6,*) "File : ",FICNADAL(1:LFICS)
	 
C* Si le fichier existe déjà, le programme n'est pas lancé
C  --------------------------------------------------------	 
        INQUIRE(FILE=FICNADAL(1:LFICS),EXIST=EX)
	 	 
        IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE   
            
C* Calcul du fichier de NADAL
C  ---------------------------
            WRITE(6,*)'==> Starting the file computation'
	     CALL SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                            ISURF,ALPHA_NADAL,BETA_NADAL,  
     &                            COEF_C_MAIGNAN,
     &	                           OS_NB,OS_NS,OS_NM,
     &                            FIC_RES_GSF,FIC_RES_FRESNEL,
     &                            FIC_RES_MAT_REFLEX,
     &	                           FICTMP_BPDF,TRACE,IER)  
            IF (IER.NE.0) GOTO 997   ! cas d'erreur

C* Ajout de la BRDF de Roujean
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP_BPDF,FICROUJ,LUM_NBMU,OS_NB,
     &                               FICSURF_TMP,IER)
	     IF (IER.NE.0) GOTO 998   ! cas d'erreur

C* Déplacement de FICSURF_TMP vers FICNADAL
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICNADAL(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICNADAL(1:LFICS),EXIST=EX)
            IF (EX) THEN !Si le fichier de BRDF/BPDF existe on détruit le fichier temporaire (généré par autre traitement)
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de BRDF/BPDF
               CALL SYSTEM(COMMAND)  
            ENDIF  

C* Destruction du fichier temporaire
C---------------------------------------   
	     OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	       
        ENDIF ! Si fichier à créer.
	 
      ENDIF ! Cas de la BPDF de Nadal




C =====================================================================
C* Cas BPDF de MAIGNAN
C =====================================================================

      IF (ISURF.EQ.7) THEN

C* Constitution du nom du fichier de BREON
C  ---------------------------------------
	 CALL SOS_NOM_FICMAIGNAN(NBMU_GAUSS,FICANGLES_USER,IND,
     &	                         COEF_C_MAIGNAN,
     &	                         K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS)	        
	 FICMAIGNAN = DIRS(1:LDIRS)//'/'//FICS
	 LFICS=INDEX(FICMAIGNAN,' ')-1  
	 WRITE(6,*) "File : ",FICMAIGNAN(1:LFICS)
	 
C* Si le fichier existe déjà, le programme n'est pas lancé
C  --------------------------------------------------------	 
        INQUIRE(FILE=FICMAIGNAN(1:LFICS),EXIST=EX)	 
	 
        IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE        
            
C* Calcul du fichier de Rondeaux - Herman
C  --------------------------------------
            WRITE(6,*)'==> Starting the file computation'

	     CALL SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                            ISURF,ALPHA_NADAL,BETA_NADAL,  
     &                            COEF_C_MAIGNAN,
     &	                           OS_NB,OS_NS,OS_NM,
     &                            FIC_RES_GSF,FIC_RES_FRESNEL,
     &                            FIC_RES_MAT_REFLEX,
     &	                           FICTMP_BPDF,TRACE,IER)    
            IF (IER.NE.0) GOTO 996   ! cas d'erreur


C* Ajout de la BRDF de Roujean
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP_BPDF,FICROUJ,LUM_NBMU,OS_NB,
     &                               FICSURF_TMP,IER)
	     IF (IER.NE.0) GOTO 998   ! cas d'erreur

C* Déplacement de FICSURF_TMP vers FICMAIGNAN
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICMAIGNAN(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICMAIGNAN(1:LFICS),EXIST=EX)
            IF (EX) THEN !Si le fichier de BRDF/BPDF existe on détruit le fichier temporaire (généré par autre traitement)
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de BRDF/BPDF
               CALL SYSTEM(COMMAND)  
            ENDIF
	    
C* Destruction du fichier temporaire
C---------------------------------------   
	     OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	       
        ENDIF ! Si fichier à créer.

      ENDIF ! Cas de la BPDF de Maignan 




C* Fermeture fichier
C-------------------

      IF (TRACE) CLOSE(99)      
      CLOSE(10) !Fichier des angles
       
C* Fin nominale 
C-------------------
      GOTO 9999
            
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  800 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.AngFile Value'
      GOTO 9998 
      
  801 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Type Value (format : I2)'
      GOTO 9998  
      
  802 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Log Value'
      GOTO 9998       
            
  803 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Dir Value'
      GOTO 9998        
                   
  804 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.DirBRDFtoBPDFcalculation Value'
      GOTO 9998        
                   
  805 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.ind Value'
      GOTO 9998        
             
  810 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Glitter.Wind Value (format : F4)'
      GOTO 9998        
           
  811 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Roujean.K0 Value (format : F7.3)',
     &           '     -SURF.Roujean.K1 Value (format : F7.3)',
     &           '     -SURF.Roujean.K2 Value (format : F7.3)'
      GOTO 9998        
           
  812 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Nadal.Alpha Value (format : F6.4)',
     &           '     -SURF.Nadal.Beta Value (format : F4.1)'
      GOTO 9998        
                    
  990 WRITE(6,*) '  ERROR on argument loading for SOS_SURFACE'
      WRITE(6,*) '  --> Check required formats.'
      GOTO 9998  
      
 9901 WRITE(6,*) '  ERROR on -SURF.Type value : out of available range.'
      CALL EXIT(1)
      GOTO 9998  
      
 9910 WRITE(6,*) '  ERROR on -SURF.Roujean.K* values : '
      WRITE(6,*) '  out of available range : ]0, 1000.]'
      GOTO 9998 
      
 9915 WRITE(6,*) '  ERROR on -SURF.Nadal.Alpha or .Beta values : '
      WRITE(6,*) '  out of available range : 0 < Alpha <= 10',
     &           ' & 0 < Beta <= 100. '
      GOTO 9998      
              
  911 WRITE(6,*) '  ERROR on logfile opening for SOS_SURFACE'
      GOTO 9998   
  
  950 WRITE(6,*) '  ERROR on angles file opening for SOS_SURFACE'
      GOTO 9998

  951 WRITE(6,*) '  ERROR on angles file reading for SOS_SURFACE'
      GOTO 9998
      
  960 WRITE(6,*) '  ERROR on logfile writing for SOS_SURFACE'
      GOTO 9998
                     
  993 WRITE(6,*) '  ERROR on subroutine SOS_GLITTER'
      GOTO 9998
      
  994 WRITE(6,*) '  ERROR on subroutine SOS_ROUJEAN'
      GOTO 9998
      
 9951 WRITE(6,*) '  ERROR while deleting temp. file MAT_BPDF_FORMATOS'
      GOTO 9998      
      
  995 WRITE(6,*) '  ERROR on subroutine SOS_RONDEAUX_BREON'
      GOTO 9998
      
  996 WRITE(6,*) '  ERROR on subroutine SOS_MAIGNAN'
      GOTO 9998
      
  997 WRITE(6,*) '  ERROR on subroutine SOS_NADAL'
      GOTO 9998
      
  998 WRITE(6,*) '  ERROR on subroutine SOS_BPDF_AJOUT_BRDF'
      GOTO 9998


 1020 WRITE(6,*) '  ERROR : mv command'
      WRITE(6,*) '  => The command to move the temporary SURFACE file',
     &           ' to the final depository is maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998

9998  IER=-1
      GOTO 9999

C* Format
C---------   
  222 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X)

    
 9999 RETURN
 	 

      END 	!FIN DU PROGRAMME SOS_SURFACE
      
      
      


C==============================================================================
C PROCEDURE: SOS_MAT_FRESNEL
C ==========
C      Cette procédure calcule le développement de la matrice de Fresnel
C      dans la base des polynomes de Legendre.
C			
C      La matrice de Fresnel, exprimee pour un repere lie au plan de réflexion
C      s'ecrit :
C                  | R11(teta)  R12(teta)    0      |
C         F(teta)= | R12(teta)  R11(teta)    0      |
C                  |    0          0      R33(teta) |
C
C         avec :
C           _ teta : angle de réflexion par rapport a la normale au dioptre,
C           _ R11(teta) = 0.5 * (RL*RL + RR*RR),
C           _ R12(teta) = 0.5 * (RL*RL - RR*RR),
C           _ R33(teta) = RL*RR,
C
C         L'indice de réfraction relatif du dioptre, entre le milieu de 
C         réflexion et le milieu de transmission, est IND.
C
C         On calcule RL et RR en fonction de IND et teta, par :
C           _ RL : le coefficient de réflexion parallele
C
C               RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C           _ RR : le coefficient de réflexion perpendiculaire
C
C               RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C
C
C      Le développement donne les paramètres ALPHA(k), BETA(k), GAMMA(k) et
C      ZETA(k) associes à la matrice de Fresnel.
C      Ils sont enregistres de k=0 à OS_NS dans le fichier "RES_FRESNEL".
C
C
C      NB: DEMARCHE D'IMPLEMENTATION :
C      *******************************
C          L'implementation du code est adaptee pour un calcul de la matrice 
C          de réflexion (routine SOS_MAT_REFLEXION) n'utilisant que des valeurs
C           positives des cos(Angle de Gauss). Pour y parvenir :
C          - On utilise non pas l'angle de diffusion vraie mais son angle
C            complémentaire (2*angle de réflexion) : DEM 1
C            (routine SOS_MAT_FRESNEL).
C          - On utilise non pas le vrai coefficient de réflexion parallele
C            mais son oppose : DEM 2 (routine SOS_MAT_FRESNEL).
C
C          En aval, la routine calculant la matrice de réflexion introduit
C          un coefficient (-1)**N pour la sommation sur les termes de la 
C          séries de Fourier pour assurer la cohérence de l'ensemble : DEM 3.
C
C
C Description des paramètres
C --------------------------
C      LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles.
C	RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double) : (E) Tableau des mu de Gauss 
C                               (angles complémentaires des angles de diffusion).
C	CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double) : (E) Tableau des poids de Gauss.
C 	IND   (double) :   (E) Indice de réfraction de l'eau / air.
C      OS_NS (I4) : (E) Ordre maximal à utiliser pour le développement en
C                        polynomes de Legendre des éléments de la matrice de Fresnel 
C                        et en séries de Fourier des calculs de luminance.
C      FIC_RES_FRESNEL (CHARACTER)  : (S) Nom du fichier temporaire RES_FRESNEL
C                                         avec arborescence complète.
C      TRACE (logical)  : (E) = vrai,si écriture dans le fichier trace
C                               (No d'unite logique du fichier trace =99)
C      IER  (I4)  :  (S)  code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C Autre résultat fourni par la procédure
C -------------------------------------
C      Fichier temporaire RES_FRESNEL
C      Informations dans le fichier de trace (optionnel)
C         (No d'unite logique du fichier trace =99 , 
C          fichier cree/ouvert par le programme appelant)
C
C 
C Description des fichiers utilisés
C ---------------------------------
C   -->Contenu du fichier temporaire RES_FRESNEL
C      Fichier Ascii, format d'un enregistrement 4E15.8, 
C      cree par SOS_MAT_FRESNEL
C      Fichier cree dans le répertoire de lancement du programme
C      Si un fichier de meme nom existe déjà, il est ecrase
C      Nombre d'enregistrements : de k=0 a OS_NS
C      Chaque enregistrement contient les paramètres ALPHA(k), BETA(k), 
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
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C
C==============================================================================
      SUBROUTINE SOS_MAT_FRESNEL(LUM_NBMU,RMU,CHR,IND,OS_NS,
     &                           FIC_RES_FRESNEL,TRACE,IER)

      IMPLICIT NONE

C* Définition des variables                    
C*----------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Tableau des mu de Gauss 
      				! (angles complémentaires des angles de diffusion).
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Tableau des poids de Gauss.

      DOUBLE PRECISION IND	! Indice de réfraction de l'eau / air.
	
      DOUBLE PRECISION ALPHA(0:CTE_OS_NS_MAX)    ! Paramètres alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NS_MAX)	! Paramètres beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NS_MAX)    ! Paramètres gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NS_MAX)	! Paramètres zeta(k).

      DOUBLE PRECISION DELTA(0:CTE_OS_NS_MAX)	! Paramètres delta(k).

      DOUBLE PRECISION A	!
      DOUBLE PRECISION B	!
      DOUBLE PRECISION C	!

      DOUBLE PRECISION RL	! Coefficient de réflexion parallele.
      DOUBLE PRECISION RR	! Coefficient de réflexion perpendiculaire.

      DOUBLE PRECISION R11(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
                                ! Element R11 de la matrice de Fresnel
      				! en fonction de l'angle de réflexion.
      DOUBLE PRECISION R12(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)     
      				! Element R12 de la matrice de Fresnel.
      DOUBLE PRECISION R33(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)    
      				! Element R33 de la matrice de Fresnel

      DOUBLE PRECISION IT	! Valeur de R11 recalculée à partir de sa
      				! decomposition en BETA(k).
      DOUBLE PRECISION QT	! Valeur de R12 recalculée à partir de sa
      				! decomposition en GAMMA(k).	
      DOUBLE PRECISION UT	! Valeur de R33 recalculée à partir de sa
      				! decomposition en DELTA(k).			

      DOUBLE PRECISION PL(-1:CTE_OS_NS_MAX+1)
                      !Polynomes de Legendre a l'ordre L calcule pour un angle
      DOUBLE PRECISION POL(0:CTE_OS_NS_MAX+1)
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

      LOGICAL TRACE        ! = vrai,si écriture dans le fichier trace

      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations

      INTEGER*4 OS_NS   ! Ordre maximal du développement en polynomes de Legendre   
			   ! des éléments de la matrice de Fresnel.
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon  
      INTEGER*4 J	! Indice de tableau.
      INTEGER*4 K	! Indice de tableau.
      INTEGER*4	NN
      INTEGER*4	MM
      INTEGER*4 I	

      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Nom du fichier RES_FRESNEL
                                            ! avec arborescence complète


C* Ouverture du fichier résultat
C---------------------------------
      OPEN(UNIT=3,FILE=FIC_RES_FRESNEL,ERR=991)

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
     &	    'Calculation of coefficients of the Fresnel matrix ' //
     &	    'expansion into Legendre functions' 
      END IF
			
C* Calcul des éléments de la matrice de réflexion de Fresnel	
C--------------------------------------------------------------
      DO 104 J=-LUM_NBMU,LUM_NBMU
         IF(J.EQ.0) GOTO 104

C* Cosinus de l'angle complémentaire de l'angle de diffusion :
C* angle de diffusion = PI - 2* angle de réflexion.
C-------------------------------------------------------------	 	
         C=RMU(J)				! DEM 1
	
C* Cosinus de l'angle de réflexion :	
C------------------------------------------------------
         C=DSQRT(.5*(1+C))


C* Coefficients de réflexion parallele (oppose a la valeur vraie)
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
     &           ' Refracted angle:',DACOS(C)*90./INCTE_PI,
     &           ' R11 ',R11(J),' R12 ',R12(J),' R33 ',R33(J)
         END IF
	 
 104  CONTINUE


C* Calcul des paramètres BETA(k) de la decomposition de la fonction
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
	   WRITE(99,*,err=921) ' Recomposition of the R11 function'
      
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
     &	                            ' R11 recomp  ',IT         
    8      CONTINUE
      END IF !Fichier Trace
      
C* Calcul des paramètres GAMMA(k) et DELTA(k) de la decomposition des fonctions
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
	   WRITE(99,*,err=921) ' Recomposition of the R12 function'
      
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
     &	                            ' R12 recomp  ',QT 
  81       CONTINUE
  
  
           WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) ' Recomposition of the R33 function'
      
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
     &                              ' R33 recomp  ',UT   
  82       CONTINUE
      END IF !Fichier Trace
      
      
      
C* Calcul des paramètres ALPHA(k) et ZETA(k) en fonctions des BETA(k) et DELTA(k).
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


C* Enregistrement des résultats
C---------------------------------
      DO 20 K=0,OS_NS
         WRITE(3,207,err=992)ALPHA(K),BETA(K),GAMMA(K),ZETA(K)
   20 CONTINUE



C* Fichier trace 
C-----------------
      IF (TRACE) THEN
           WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) 'Limit order : ',OS_NS
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
C      Cette procédure calcule la matrice de réflexion de surface 
C      pour le glitter.
C      Les angles d'incidence et de réflexion sont exprimes comme etant
C      inferieurs à 90 degrés (mu de Gauss > 0).
C
C
C
C      NB: DEMARCHE D'IMPLEMENTATION :
C      *******************************
C          L'implementation du code est adaptee pour un calcul de la matrice 
C          de réflexion (routine SOS_MAT_REFLEXION) n'utilisant que des valeurs
C           positives des cos(Angle de Gauss). Pour y parvenir :
C          - On utilise non pas l'angle de diffusion vraie mais son angle
C            complémentaire (2*angle de réflexion) : DEM 1
C            (routine SOS_MAT_FRESNEL).
C          - On utilise non pas le vrai coefficient de réflexion parallele
C            mais son oppose : DEM 2 (routine SOS_MAT_FRESNEL).
C
C          Cette routine calculant la matrice de réflexion introduit
C          un coefficient (-1)**K pour le calcul des termes de la 
C          séries de Fourier pour assurer la cohérence de l'ensemble : DEM 3.
C
C
C Description des paramètres
C --------------------------
C     COEF  (double): (E) Coefficient multiplicateur de la matrice de réflexion.
C     LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles.
C     RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)    (E) Tableau des mu de Gauss 
C                             (angles complémentaires des angles de diffusion).
C     OS_NB (I4) : (E) Ordre des développements des fonctions de phase 
C		         en Polynomes de Legendre à utiliser.
C     OS_NS (I4) : (E) Ordre maximal à utiliser pour le développement en
C                      polynomes de Legendre des éléments de la matrice de Fresnel 
C                      et en séries de Fourier des calculs de luminance.
C     OS_NM (I4) : (E) Ordre maximal du développement en séries de Fourier  
C                      de la fonction G des calculs de matrice de réflexion.  
C     FIC_RES_GSF (CHARACTER)  : (E) Nom du fichier temporaire RES_GSF
C                                    avec arborescence complète.
C     FIC_RES_FRESNEL (CHARACTER)  : (E) Nom du fichier temporaire RES_FRESNEL
C                                        avec arborescence complète.
C     FIC_RES_MAT_REFLEX (CHARACTER)  : (S) Nom du fichier temporaire RES_MAT_REFLEX
C                                           avec arborescence complète.
C     IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C
C 
C Description des fichiers utilisés
C ---------------------------------
C   -->Contenu du fichier temporaire RES_GSF
C       Fichier binaire non formatte, cree par SOS_GSF
C       Fichier cree dans le répertoire de lancement du programme
C       Chaque enregistrement contient :
C           -numéro d'angle Teta1 (valeur entiere)
C           -numéro d'angle Teta2 (valeur entiere)
C           -ordre maximale utile de la decomposition : IL  (valeur entiere)
C           -Fonction G decomposee pour (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 valeurs double précision)
C
C   -->Contenu du fichier temporaire RES_FRESNEL
C      Fichier Ascii, format d'un enregistrement 4E15.8, 
C      cree par SOS_MAT_FRESNEL
C      Fichier cree dans le répertoire de lancement du programme
C      Nombre d'enregistrements : de k=0 a OS_NS
C      Chaque enregistrement contient les paramètres ALPHA(k), BETA(k), 
C      GAMMA(k) et ZETA(k) associes a la matrice de Fresnel
C
C   -->Contenu du fichier temporaire RES_MAT_REFLEX
C       Fichier binaire non formatte, cree par SOS_MAT_REFLEXION
C       Fichier cree dans le répertoire de lancement du programme
C       Si un fichier de meme nom existe déjà, il est ecrase
C       Enregistrement pour chaque couple (Teta1, Teta2) de
C       la decomposition en séries de Fourier de la matrice de Fresnel:
C       M11,M12,M13,M21,M22,M23,M31,M32,M33 
C       (n valeurs réelles)
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
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C
C     
C  Routines utilisées par la procédure:
C -------------------------------------
C  Ce programme fait appel aux routines:
C         -SOS_NOYAUX_FRESNEL
C
C=============================================================================
      SUBROUTINE SOS_MAT_REFLEXION(COEF,LUM_NBMU,RMU,OS_NB,OS_NS,OS_NM,
     &                             FIC_RES_FRESNEL, FIC_RES_GSF,
     &                             FIC_RES_MAT_REFLEX, IER)

      IMPLICIT NONE

C* Définition des variables     
C*--------------------------------------------------------------------
      DOUBLE PRECISION COEF !Coefficient multiplicateur de la matrice de réflexion.
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Tableau des mu de Gauss 
				!(angles complémentaires des angles de diffusion)

      DOUBLE PRECISION ALPHA(0:CTE_OS_NS_MAX)   ! Paramètres alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NS_MAX)    ! Paramètres beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NS_MAX)   ! Paramètres gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NS_MAX)    ! Paramètres zeta(k).

      DOUBLE PRECISION G(0:CTE_OS_NM_MAX)	! Fonction G decomposee en séries 
      					        ! de Fourier.

       ! Element Pij de la matrice de réflexion de Fresnel decomposee en séries 
	! de Fourier pour un ordre IS = 0 a OS_NS.
	! Calcul pour la direction d'incidence Teta1 (RMU1)
	! et de réflexion Teta2 (RMU2).
      DOUBLE PRECISION BP(0:CTE_OS_NS_MAX,2)    !   Element P11.
      DOUBLE PRECISION GR(0:CTE_OS_NS_MAX,2)    !   Element P12 ou P21.
      DOUBLE PRECISION GT(0:CTE_OS_NS_MAX,2)    !   Element -P13 ou -P31.
      DOUBLE PRECISION ARR(0:CTE_OS_NS_MAX,2)   !   Element P22.
      DOUBLE PRECISION ART(0:CTE_OS_NS_MAX,2)   !   Element -P23 ou -P32.
      DOUBLE PRECISION ATT(0:CTE_OS_NS_MAX,2)   !   Element P33.

      DOUBLE PRECISION R111    ! Variable de calcul du terme M11(IS,1).
      DOUBLE PRECISION R121    ! Variable de calcul du terme M12(IS,1).
      DOUBLE PRECISION R122    ! Variable de calcul du terme M12(IS,2).
      DOUBLE PRECISION R131    ! Variable de calcul du terme M13(IS,1).
      DOUBLE PRECISION R132    ! Variable de calcul du terme M13(IS,2).
      DOUBLE PRECISION R211    ! Variable de calcul du terme M21(IS,1).
      DOUBLE PRECISION R212    ! Variable de calcul du terme M21(IS,2).
      DOUBLE PRECISION R221    ! Variable de calcul du terme M22(IS,1).
      DOUBLE PRECISION R222    ! Variable de calcul du terme M22(IS,2).
      DOUBLE PRECISION R231    ! Variable de calcul du terme M23(IS,1).
      DOUBLE PRECISION R232    ! Variable de calcul du terme M23(IS,2).
      DOUBLE PRECISION R311    ! Variable de calcul du terme M31(IS,1).
      DOUBLE PRECISION R312    ! Variable de calcul du terme M31(IS,2).
      DOUBLE PRECISION R321    ! Variable de calcul du terme M32(IS,1).
      DOUBLE PRECISION R322    ! Variable de calcul du terme M32(IS,2).
      DOUBLE PRECISION R331    ! Variable de calcul du terme M33(IS,1).
      DOUBLE PRECISION R332    ! Variable de calcul du terme M33(IS,2).

      DOUBLE PRECISION X,Y	! Coefficients multiplicatifs de la matrice de
      				! réflexion.

      	! Elements Mij de la matrice de réflexion en fonction de l'ordre IS
	! de la série de Fourier, pour 2 directions.
	! L'indice de direction est :
	!     _ 1 pour une incidence sous Teta1 et une réflexion sous Teta2,
	!     _ 2 pour une incidence sous Teta2 et une réflexion sous Teta1.
      REAL M11(0:CTE_OS_NB_MAX,2)		
      REAL M12(0:CTE_OS_NB_MAX,2)	
      REAL M13(0:CTE_OS_NB_MAX,2)		
      REAL M21(0:CTE_OS_NB_MAX,2)	
      REAL M22(0:CTE_OS_NB_MAX,2)		
      REAL M23(0:CTE_OS_NB_MAX,2)		
      REAL M31(0:CTE_OS_NB_MAX,2)	
      REAL M32(0:CTE_OS_NB_MAX,2)		
      REAL M33(0:CTE_OS_NB_MAX,2)	

      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations
      
      INTEGER*4 OS_NB ! Ordre des développements des fonctions de phase 
			 ! en Polynomes de Legendre.
      INTEGER*4 OS_NS ! Ordre maximal du développement en polynomes de Legendre   
			 ! des éléments de la matrice de Fresnel.
      INTEGER*4 OS_NM ! Ordre maximal du développement en séries de Fourier  
			 ! de la fonction G des calculs de matrice de réflexion.
 
      INTEGER*4 I	! Indice de l'angle Teta1.
      INTEGER*4 J	! Indice de l'angle Teta2.

      INTEGER*4 L	! Ordre de la decomposition en polynomes de Legendre.
      INTEGER*4 LIM	! Ordre maximal de la decomposition en polynomes
      			! de Legendre.
      INTEGER*4 II,II2	! Indices enregistres avec la fonction G(IS).
      INTEGER*4 IJ	
      INTEGER*4 IS	! Ordre de la decomposition en séries de Fourier.
      INTEGER*4 K
      INTEGER*4 I1	! Valeur K+IS
      INTEGER*4 I2	! Valeur ABS(K-IS)
      INTEGER*4 IM	! Valeur 1 ou -1.
      INTEGER*4 IER  ! code d'erreur =0 si pas d'erreur, =-1 sinon  
  
      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX 
                     ! Nom du fichier RES_MAT_REFLEX
                     ! avec arborescence complète

      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Nom du fichier RES_GSF
                                        ! avec arborescence complète

      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Nom du fichier RES_FRESNEL
                                            ! avec arborescence complète

      
C* Creation du fichier résultat de SOS_MAT_REFLEXION
C---------------------------------------------------------
      OPEN(UNIT=4,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED',ERR=991)


C* Lecture des paramètres de la decomposition de la matrice de Fresnel
C* en polynomes ou fonctions de Legendre.
C----------------------------------------------------------------------
      OPEN(UNIT=2,FILE=FIC_RES_FRESNEL,STATUS='OLD',ERR=991)
      
      DO 951 L=0,OS_NS
         READ(2,1300,err=992) ALPHA(L),BETA(L),GAMMA(L),ZETA(L)
  951 CONTINUE
  
      CLOSE(2)
      
C* Ouverture du fichier RES_GSF
C---------------------------------------------------------
      OPEN(UNIT=3,FILE=FIC_RES_GSF,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)

C* Boucle sur les directions d'incidence et de réflexion.
C--------------------------------------------------------
      DO 25 I=1,LUM_NBMU   !Teta1 : angle zenithal d'incidence  
                           !(0 < Teta1 < PI/2),
         DO 2 J=1,I        !Teta2 : angle zenithal de réflexion 
	                    !(0 < Teta2 < PI/2),

C* Lecture de la decomposition en séries de Fourier de la fonction G
C* pour les angles Teta1 et Teta2,avec un ordre maximal utile L.
C-------------------------------------------------------------------
            READ(3,err=992)II,II2,LIM,(G(K),K=0,LIM)

C* Complément du tableau à 0.
C---------------------------------
            DO IJ=LIM+1,OS_NM
               G(IJ)=0.
            ENDDO
	      
C------------------------------------------------------------------------
C* Decomposition en séries de Fourier de la matrice de réflexion de Fresnel
C* pour les directions Teta1 et Teta2
C--------------------------------------------------------------------------

C*  Calcul des fonctions BP,GR,GT,ARR,ART et ATT en fonction de l'ordre IS
C*  pour les directions Teta1 et Teta2
C---------------------------------------------------------------------------
            CALL SOS_NOYAUX_FRESNEL(RMU(I),RMU(J),OS_NS,
     &                              ALPHA,BETA,GAMMA,ZETA,
     &                              BP,GR,GT,ARR,ART,ATT)

C* Calcul des éléments de la matrice de réflexion pour chaque ordre IS.
C---------------------------------------------------------------------
            DO 24 IS=0,OS_NB
               X=COEF*G(IS)/4.
	       
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
                 IM=-IM					! DEM 3
                 I1=K+IS
                 I2=ABS(K-IS)
                 IF ((I1.GT.LIM).AND.(I2.GT.LIM)) GOTO 3	
                 X=COEF*IM*(G(I1)+G(I2))/4.		! DEM 3
                 Y=COEF*IM*(G(I2)-G(I1))/4.		! DEM 3
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

C* Affectation des résultats aux tableaux résultats Rij(IS)
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
C* la decomposition en séries de Fourier de la matrice de Fresnel.
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
C      Cette procédure calcule la decomposition en séries de Fourier de
C      la matrice de Fresnel selon l'azimut, pour les directions Teta1 et Teta2.
C
C      On calcule les éléments de matrice suivants :
C         BP(IS,mu)    GR(IS,mu)    GT(IS,mu)
C         ARR(IS,mu)   ART(IS,mu)   ATT(IS,mu)
C      pour les ordres IS de la decomposition en séries de Fourier
C      et pour mu = cos(Teta1) ou cos(Teta2).
C
C
C Description des paramètres
C --------------------------
C	RMU1   (double)	 (E) Cos(Teta1).
C	RMU2   (double)	 (E) Cos(Teta2).
C       OS_NS (I4)           (E) Ordre maximal à utiliser pour le développement en
C                                polynomes de Legendre des éléments de la matrice 
C                                de Fresnel  et en séries de Fourier des calculs 
C                                de luminance.
C
C	ALPHA(0:CTE_OS_NS_MAX)  (double)	(E) Paramètres alpha(k).
C	BETA(0:CTE_OS_NS_MAX)   (double)	(E) Paramètres beta(k).
C	GAMMA(0:CTE_OS_NS_MAX)  (double)	(E) Paramètres gamma(k).
C	ZETA(0:CTE_OS_NS_MAX)   (double)	(E) Paramètres zeta(k).
C
C       Element Pij de la matrice de réflexion de Fresnel decomposee en séries 
C       de Fourier pour un ordre IS = 0 a OS_NS.
C       Calcul pour la direction d'incidence Teta1 (RMU1)
C       et de réflexion Teta2 (RMU2).	
C	BP(0:CTE_OS_NS_MAX,2)   (double)	(S) Element P11.
C	GR(0:CTE_OS_NS_MAX,2)   (double)	(S) Element P12 ou P21.
C	GT(0:CTE_OS_NS_MAX,2)   (double)	(S) Element -P13 ou -P31.
C	ARR(0:CTE_OS_NS_MAX,2)  (double)	(S) Element P22.
C	ART(0:CTE_OS_NS_MAX,2)  (double)	(S) Element -P23 ou -P32.
C	ATT(0:CTE_OS_NS_MAX,2)  (double)	(S) Element P33.
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

C* Définition des variables                                  
C*----------------------------------------------------------------------
      DOUBLE PRECISION RMU1	! Cos(Teta1).
      DOUBLE PRECISION RMU2	! Cos(Teta2).

      DOUBLE PRECISION ALPHA(0:CTE_OS_NS_MAX)    ! Paramètres alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NS_MAX)	! Paramètres beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NS_MAX)    ! Paramètres gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NS_MAX)	! Paramètres zeta(k).

      DOUBLE PRECISION R(2)	! tableau contenant cos(Teta1) et cos(Teta2).

      	! Element Pij de la matrice de réflexion de Fresnel decomposee en séries 
	! de Fourier pour un ordre IS = 0 a OS_NS.
	! Calcul pour la direction d'incidence Teta1 (RMU1)
	! et de réflexion Teta2 (RMU2).	
      DOUBLE PRECISION BP(0:CTE_OS_NS_MAX,2)    ! Element P11.
      DOUBLE PRECISION GR(0:CTE_OS_NS_MAX,2)    ! Element P12 ou P21.
      DOUBLE PRECISION GT(0:CTE_OS_NS_MAX,2)    ! Element -P13 ou -P31.
      DOUBLE PRECISION ARR(0:CTE_OS_NS_MAX,2)   ! Element P22.
      DOUBLE PRECISION ART(0:CTE_OS_NS_MAX,2)   ! Element -P23 ou -P32.
      DOUBLE PRECISION ATT(0:CTE_OS_NS_MAX,2)   ! Element P33.

      DOUBLE PRECISION PSL(-1:CTE_OS_NS_MAX,2)  
                            ! Fonction de Legendre PSL pour l'ordre S=IS
      				! en fonction de l'ordre L et de l'angle mu.
      DOUBLE PRECISION RSL(-1:CTE_OS_NS_MAX,2)  ! Fonction RSL.
      DOUBLE PRECISION TSL(-1:CTE_OS_NS_MAX,2)  ! Fonction TSL.

      DOUBLE PRECISION C	     ! Valeur ponctuelle du mu de Gauss.
      DOUBLE PRECISION A,B,D,E,F ! Coefficients pour la recurrence 
      				     ! sur PSL, RSL et TSL.

      DOUBLE PRECISION X	! Variable intermédiaire de calcul.
      DOUBLE PRECISION XX	! Variable intermédiaire de calcul.
      DOUBLE PRECISION YY	! Variable intermédiaire de calcul.

      DOUBLE PRECISION RAC3	! Racine de 3.
      DOUBLE PRECISION X26	! 2*racine(6)

      INTEGER*4 OS_NS   ! Ordre maximal du développement en polynomes de Legendre   
			   ! des éléments de la matrice de Fresnel.

      INTEGER*4 IS	! Ordre de la decomposition en séries de Fourier.
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

C* Decomposition en séries de Fourier de IS = 0 a OS_NS
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
         DO I=1,IS
            X=I
            A=A*DSQRT((I+IS)/X)*0.5
         ENDDO
	 
         B=A*DSQRT(IS/(IS+1.D+00))*DSQRT((IS-1.D+00)/(IS+2.)) 
         DO J=1,2
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
         ENDDO

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
C* Calcul des éléments de la matrice de réflexion
C------------------------------------------------------

C* Boucles sur les directions Teta1 et Teta2
C------------------------------------------------------
         DO K=1,2
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
             DO L=IS,OS_NS
                BP(IS,K)=BP(IS,K)+BETA(L)*PSL(L,J)*PSL(L,K)
	         GR(IS,K)=GR(IS,K)+GAMMA(L)*PSL(L,J)*RSL(L,K)
                GT(IS,K)=GT(IS,K)+GAMMA(L)*PSL(L,J)*TSL(L,K)
                ATT(IS,K)=ATT(IS,K)+ALPHA(L)*TSL(L,J)*TSL(L,K)+
     &                    ZETA(L)*RSL(L,J)*RSL(L,K)
                ARR(IS,K)=ARR(IS,K)+ZETA(L)*TSL(L,J)*TSL(L,K)+
     &		          ALPHA(L)*RSL(L,J)*RSL(L,K)
                ART(IS,K)=ART(IS,K)+ALPHA(L)*RSL(L,K)*TSL(L,J)+
     &                    ZETA(L)*RSL(L,J)*TSL(L,K)
             ENDDO  ! boucle sur  L
	     
         ENDDO  ! boucle sur K 

   1  CONTINUE
      RETURN
      END         !FIN DE LA PROCEDURE SOS_NOYAUX_FRESNEL






C==============================================================================
C PROCEDURE:  SOS_MISE_FORMAT
C ==========
C      Cette procédure reordonne les éléments de la matrice Mij en une matrice
C      Pij utilisable par le code des OS.
C
C      Les éléments Mij, enregistres dans le fichier "RES_MAT_REFLEX", donnent
C      pour chaque couple d'angles (Teta1, Teta2) la decomposition en séries 
C      de Fourier associee. 
C
C      Les éléments Pij, enregistres dans le fichier de sortie, donnent
C      pour chaque ordre de la decomposition en séries de Fourier
C      la réflexion associee a l'ensemble des couples (Teta1, Teta2) .
C
C Description des paramètres
C --------------------------
C      FIC_RES_MAT_REFLEX (CHARACTER)  : (E) Nom du fichier temporaire RES_MAT_REFLEX
C                                            avec arborescence complète.
C      FICGLITTER  (CHARACTER*CTE_LENFIC2) : (S)  Nom complet du fichier résultat.
C      NBMU (I4)  : (E)  Nombre d'angles (positifs) effectivement utiles.
C      OS_NB (I4) : (E)  Ordre des développements des fonctions de phase 
C		         en Polynomes de Legendre à utiliser.
C      IER  (I4)  : (S)  code d'erreur =0 si pas d'erreur, =-1 sinon
C
C
C Autre résultat fourni par la procédure
C -------------------------------------
C      Fichier résultat : FICGLITTER, fichier GLITTER
C 
C Description des fichiers utilisés
C ---------------------------------
C   -->Contenu du fichier temporaire RES_MAT_REFLEX
C       Fichier binaire non formatte, cree par SOS_MAT_REFLEXION
C       Fichier cree dans le répertoire de lancement du programme
C       Enregistrement pour chaque couple (Teta1, Teta2) de
C       la decomposition en séries de Fourier de la matrice de Fresnel:
C       M11,M12,M13,M21,M22,M23,M31,M32,M33 
C       (n valeurs réelles)
C
C   -->Contenu du fichier GLITTER
C       Fichier binaire non formatte, cree par SOS_MISE_FORMAT
C       Il contient la matrice de réflexion de surface 
C       (fichier GLITTER ou autres BRDF):
C       c'est a dire les éléments Pij de la matrice de réflexion 
C              P11(LUM_NBMU,LUM_NBMU) (réel), 
C              P12(LUM_NBMU,LUM_NBMU) (réel),
C              P13(LUM_NBMU,LUM_NBMU) (réel),
C              P21(LUM_NBMU,LUM_NBMU) (réel),
C              P22(LUM_NBMU,LUM_NBMU) (réel),
C              P23(LUM_NBMU,LUM_NBMU) (réel),
C              P31(LUM_NBMU,LUM_NBMU) (réel),
C              P32(LUM_NBMU,LUM_NBMU) (réel),
C              P33(LUM_NBMU,LUM_NBMU) (réel)
C
C        --> Enregistrements sur le domaine limité aux angles utiles
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
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C
C=============================================================================
      SUBROUTINE SOS_MISE_FORMAT(FIC_RES_MAT_REFLEX, FICGLITTER,
     &                           NBMU,OS_NB,IER)

      IMPLICIT NONE

C* Définition des variables                    
C*--------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICGLITTER	! Nom du fichier résultat.

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Nom du fichier RES_MAT_REFLEX
                                               ! avec arborescence complète

      	! Elements Mij de la matrice de réflexion en fonction de l'ordre IS
	! de la série de Fourier, pour 2 directions.
	! L'indice de direction est :
	!     _ 1 pour une incidence sous Teta1 et une réflexion sous Teta2,
	!     _ 2 pour une incidence sous Teta2 et une réflexion sous Teta1.
      REAL M11(0:CTE_OS_NB_MAX,2)		
      REAL M12(0:CTE_OS_NB_MAX,2)		
      REAL M13(0:CTE_OS_NB_MAX,2)			
      REAL M21(0:CTE_OS_NB_MAX,2)	
      REAL M22(0:CTE_OS_NB_MAX,2)		
      REAL M23(0:CTE_OS_NB_MAX,2)		
      REAL M31(0:CTE_OS_NB_MAX,2)	
      REAL M32(0:CTE_OS_NB_MAX,2)	
      REAL M33(0:CTE_OS_NB_MAX,2)

      	! Elements Pij de la matrice de réflexion en fonction de la direction
	! d'incidence et de réflexion, pour un ordre IS du développement
	! en séries de Fourier
      REAL P11(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P12(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P13(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)	
      REAL P21(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P22(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P23(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P31(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P32(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P33(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)

      INTEGER*4 NBMU     ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB    ! Ordre des développements des fonctions de phase 
			    ! en Polynomes de Legendre.
      INTEGER*4 IS	! Ordre de la decomposition en séries de Fourier.
      INTEGER*4 I	! Indice de l'angle d'incidence Teta1.
      INTEGER*4 J	! Indice de l'angle de réflexion Teta2.
      INTEGER*4 IER  ! code d'erreur =0 si pas d'erreur, =-1 sinon  
	
C* Ouverture fichier
C------------------------------------------------------------	
	
      OPEN(UNIT=1,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=2,FILE=FICGLITTER,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=991)

C* Boucle sur l'ordre de la decomposition en séries de Fourier
C------------------------------------------------------------
      DO 1 IS=0,OS_NB

C* Boucle sur les directions d'incidence Teta1
C---------------------------------------------
         DO 2 I=1,NBMU

C* Boucle sur les directions de réflexion Teta2
C-----------------------------------------------

C*       Lecture du fichier RES_MAT_REFLEX enregistre pour la taille max
C*       des tableaux : CTE_OS_NBMU_MAX    
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
 
C*       Enregistrement des tableaux de la matrice de réflexion 
C*       a l'ordre IS de la decomposition en séries de Fourier.
C*       Bornes d'enregistrement des tableaux limitées au 
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
C      Cette procédure lit deux fichiers de BRDF ou BPDF mis en forme pour les OS
C      et realise la sommation des termes pour creer un fichier final somme des deux
C      premiers.
C
C Description des paramètres
C --------------------------
C       FICSURF1 (CHARACTER*CTE_LENFIC2) (E) : nom du 1er fichier de réflexion de surface, 
C                                              au format pour les OS.
C                                             (nom complet avec répertoire).
C       FICSURF2 (CHARACTER*CTE_LENFIC2) (E) : nom du 2nd fichier de réflexion de surface, 
C                                              au format pour les OS.
C                                             (nom complet avec répertoire).
C       NBMU (I4)  (E) : nombre d'angles (positifs) effectivement utiles.
C       OS_NB (I4) (E) : Ordre des développements des fonctions de phase 
C		         en Polynomes de Legendre à utiliser.
C       FICS   (CHARACTER*CTE_LENFIC2)  (E)  : nom du fichier SURFACE résultat.
C                                             (nom complet avec répertoire).
C       IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C Description des fichiers utilisés
C ---------------------------------
C   -->Contenu du fichier BRDF ou BPDF
C       Fichier binaire non formatte.
C       Il contient la matrice de réflexion de surface (fichier BRDF / BPDF):
C       c'est a dire les éléments Rij de la matrice de réflexion 
C       pour OS_NB ordres de la decomposition en séries de Fourier.
C              P11(LUM_NBMU,LUM_NBMU) (réel),
C              P12(LUM_NBMU,LUM_NBMU) (réel),
C              P13(LUM_NBMU,LUM_NBMU) (réel),
C              P21(LUM_NBMU,LUM_NBMU) (réel),
C              P22(LUM_NBMU,LUM_NBMU) (réel),
C              P23(LUM_NBMU,LUM_NBMU) (réel),
C              P31(LUM_NBMU,LUM_NBMU) (réel),
C              P32(LUM_NBMU,LUM_NBMU) (réel),
C              P33(LUM_NBMU,LUM_NBMU) (réel)
C
C        --> Enregistrements sur le domaine limité aux angles utiles
C                                   --------------------------------
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C============================================================================== 
      SUBROUTINE SOS_BPDF_AJOUT_BRDF(FICSURF1,FICSURF2,NBMU,OS_NB,
     &                               FICS,IER)

      IMPLICIT NONE

C* Définition des variables   
C*--------------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICSURF1 ! nom du 1er fichier de réflexion 
      					  ! de surface, au format pour les OS.
      CHARACTER*CTE_LENFIC2 FICSURF2 ! nom du 2nd fichier de réflexion 
      					  ! de surface, au format pour les OS.
      CHARACTER*CTE_LENFIC2 FICS     ! nom du fichier résultat de réflexion 
      					  ! de surface, au format pour les OS.
					  
      	! Elements Pij de la matrice de réflexion en fonction de la direction
	! d'incidence et de réflexion, pour un ordre IS du développement
	! en séries de Fourier. Format d'entrée des OS pour le 1er fichier.
      REAL P11_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P12_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P13_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)	
      REAL P21_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P22_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P23_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P31_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P32_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P33_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
           
        ! Elements Pij de la matrice de réflexion en fonction de la direction
	! d'incidence et de réflexion, pour un ordre IS du développement
	! en séries de Fourier. Format d'entrée des OS pour le 2nd fichier.
      REAL P11_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P12_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P13_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)	
      REAL P21_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P22_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P23_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P31_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P32_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P33_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)   
          
        ! Elements Pij de la matrice de réflexion en fonction de la direction
	! d'incidence et de réflexion, pour un ordre IS du développement
	! en séries de Fourier. Format d'entrée des OS pour le fichier résultat.
      REAL P11(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P12(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P13(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)	
      REAL P21(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P22(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P23(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P31(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P32(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P33(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)

      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB   ! Ordre des développements des fonctions de phase 
			   ! en Polynomes de Legendre.

      INTEGER*4 IS	! Ordre de la decomposition en séries de Fourier.
      INTEGER*4 I	! Indice de l'angle d'incidence Teta1.
      INTEGER*4 J	! Indice de l'angle de réflexion Teta2.
      INTEGER*4 IER  ! code d'erreur =0 si pas d'erreur, =-1 sinon	  					  					  
      	      
      

C* Ouverture des fichiers
C------------------------------------------------------------		    
      OPEN(UNIT=1,FILE=FICSURF1,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=2,FILE=FICSURF2,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=3,FILE=FICS,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=991)

	   
C* Boucle sur l'ordre de la decomposition en séries de Fourier
C------------------------------------------------------------      
      DO 1 IS=0,OS_NB

C*      Lecture des fichiers qui ont ete enregistres avec des 
C*      bornes de tableaux limitées au dimensionnement 
C*      utile -NBMU:NBMU           
	   READ(1,err=992)
     &	       ((P11_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P12_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P13_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P21_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P22_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P23_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P31_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P32_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P33_SURF1(I,J),I=1,NBMU),J=1,NBMU)  
                 
	   READ(2,err=992)
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

C* Boucle sur les directions de réflexion Teta2
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
 
C*       Enregistrement des tableaux de la matrice de réflexion 
C*       a l'ordre IS de la decomposition en séries de Fourier.
C*       Bornes d'enregistrement des tableaux limitées au 
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
      
