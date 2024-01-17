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
C* FICHIER: SOS_PREPA_OS.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Module principal du code des OS
C*           - lit les paramètres de granulométries aérosols,
C*           - lit les paramètres angulaires 
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, B.Lafrance, X.Lenot
C*
C* DATE: 24/06/2013
C* MOD:VERSION:1.0: Creation du fichier par separation de la fonction SOS.F en 3
C*		      sous fonctions : SOS.F, SOS_PREPA_OS.F et SOS_FIN_OS.F
C*
C* MOD:VERSION:1.1: 24/02/2016 
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*
C*    - Changement de la manière de lire les fichiers formatés (Angles, Aérosols)
C*      pour être compatible avec le compilateur gfortran.
C*
C* MOD:VERSION:1.2: 17/07/2020 
C*    - Ajustements pour introduire le modèle de Maignan.
C*
C* MOD:VERSION:1.3: 07/07/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    CTE_LENDIR : Longueur des noms de répertoires.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_OS_NBMU_MAX : Nombre maximal d'angles positifs pour la résolution de l'ETR
C    CTE_OS_NB_MAX :  Valeur maximale pour l'ordre limite du développement en 
C                     polynomes de Legendre.
C
C Constantes spécifiques :
C    LENLINE  : Longueur maximale d'une ligne de données lue dans un fichier
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define LENLINE 1000

C==============================================================================
C PROGRAMME: SOS_PREPA_OS
C ==========      
C      Programme qui initie les variables pour le calcul des OS : lecture des
C      fichiers d'angles et d'aérosols. Identification du mode de calcul
C      de la réflectance de surface et du fichier associé.
C
C Données en entrée du programme 
C ------------------------------
C    -->Les paramètres/arguments du programme
C
C            FICANGLES (CHARACTER) : (E) nom complet du fichier des angles utilisés 
C                                      pour les calculs de luminance et BRDF
C                 		           (répertoire + nom fichier + extension)
C
C	      FICGRANU (CHARACTER) : (E) nom complet du fichier résultat de GRANU
C                                      (répertoire + nom fichier + extension)
C
C            ISURF  (I4)    : (E) Type de surface ISURF 
C			            (format de décodage en entier I4)
C      		0 :surface lambertienne d'albédo RHO.
C			1 :surface lambertienne + glitter.
C			2 :surface lambertienne + réflexion de Fresnel.
C		       3 :surface lambertienne + BRDF de Roujean.
C			4 :surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
C			5 :surface lambertienne + BRDF de Roujean + BPDF de Breon.
C			6 :surface lambertienne + BRDF de Roujean + BPDF de Nadal.		 	
C
C            IND_SURF (I4) : (E) Indice de réfraction de la surface / air.		 	
C
C            WIND (DOUBLE) : (E) ==> Associé à : ISURF = 1
C 	      		       Vitesse du vent en m/s.	
C		
C            K0_ROUJEAN (DOUBLE) : (E) ==> Associé à : ISURF >= 3
C 	       		      Paramètre k0 du modèle de BRDF de Roujean.       	
C		
C            K1_ROUJEAN (DOUBLE) : (E) ==> Associé à : ISURF >= 3
C 	      			     Paramètre k1 du modèle de BRDF de Roujean.       	
C		
C            K2_ROUJEAN (DOUBLE) : (E) ==> Associé à : ISURF >= 3
C 	      			     Paramètre k2 du modèle de BRDF de Roujean.       	
C
C            ALPHA_NADAL (DOUBLE) : (E) ==> Associé à : ISURF = 6
C 	      			     Paramètre alpha du modèle de BPDF de Nadal.       	
C
C            BETA_NADAL (DOUBLE) : (E) ==> Associé à : ISURF = 6
C 	      			     Paramètre beta du modèle de BPDF de Nadal.       	
C
C            COEF_C_MAIGNAN (DOUBLE) : (E) ==> Associé à : ISURF = 7
C 	      		            Paramètre C*exp(-NDVI) du modèle de BPDF de Maignan.  	
C
C
C            REPSURF (CHARACTER) : (E) ==> Associé à : ISURF= 1 3 4 5 ou 6 et FICSURF= "DEFAULT"
C             			     Répertoire de stockage des fichiers de matrice de réflexion 
C                   		            de surface (fichier Glitter ou autres BRDF/BPDF).
C
C            FICSURF (CHARACTER) : (E) ==> Associé à : ISURF= 1 3 4 5 ou 6 
C             			    nom complet du fichier contenant la matrice de réflexion 
C                   		           de surface  (fichier Glitter ou autres BRDF/BPDF).
C                   		           (répertoire + nom fichier + extension)
C                    		    => DEFAULT si on utilise un fichier de la base des BRDF / BPDF.	 
C
C            TRACE (LOGICAL) : (E) Si vrai => écriture des fichiers trace
C
C	      IDLOG (I4) : (E) Identifiant du fichier log
C
C            LUM_NBMU (I4) : (S) Nombre d'angles effectifs pour les calculs de luminances      			               
C
C            N0 (I4) : (S) Numéro de mu de Gauss pour l'angle zénithal solaire.      			
C			     (angle complémentaire : mus = -rmu(N0))
C			     (N0 > 0 => pour en tenir compte).
C
C            PIZ (DOUBLE) : (S) Albédo de simple diffusion.
C
C            PIZTR (DOUBLE) : (S) Albédo de simple diffusion.après troncature
C
C            A (DOUBLE) : (S) Coefficient de troncature
C
C            RMU (DOUBLE) : (S) Tableau 1D de 2*CTE_OS_NBMU_MAX + 1 éléments
C			           contenant le cosinus des angles de Gauss et des
C			           angles fournis par l'utilisateur
C
C            GA (DOUBLE) : (S) Tableau 1D de 2*CTE_OS_NBMU_MAX + 1 éléments
C			          contenant le poids associé à chaque angle de Gauss et chaque
C			          angle fourni par l'utilisateur
C 
C            OS_NB (I4) : (S) Ordre des développements des fonctions de phase 
C			         en Polynomes de Legendre.      			               
C
C            IMAT_SURF (I4) :(S) Indice précisant si les simulations
C      			     utilisent un modèle de réflexion BRDF ou BPDF :
C			            1 : pour calcul de BRDF ou BPDF (il faut utiliser N0 > 0).
C			            sinon : pas de simulation de BRDF ou BPDF. 	
C
C            IFRESNEL (I4) : (S) Indice précisant si les simulations
C      			     introduisent la réflexion de Fresnel sur mer plate :
C			      1 : pour un calcul de la réflexion de Fresnel.
C			      Sinon : pas de prise en compte.    
C
C            IGLI (I4) : (S) Indice précisant si les simulations introduisent le Glitter
C			        1 : pour un calcul de Glitter
C			        Sinon : pas de prise en compte.    
C 
C            IROUJEAN (I4) : (S) Indice précisant si les simulations
C      			     introduisent le modèle de réflexion de Roujean :
C			            1 : pour un calcul dde BRDF Roujean (il faut utiliser N0 > 0).
C			            Sinon : pas de prise en compte.    
C
C            IRONDEAUX (I4) : (S) Indice précisant si les simulations
C      			      introduisent le modèle de réflexion de Rondeaux - Herman :
C			             1 : pour un calcul dde BRDF Rondeaux - Herman (il faut utiliser N0 > 0).
C			             Sinon : pas de prise en compte.    
C
C            IBREON (I4) : (S) Indice précisant si les simulations
C      			   introduisent le modèle de réflexion de Breon :
C			          1 : pour un calcul dde BRDF Breon (il faut utiliser N0 > 0).
C			          Sinon : pas de prise en compte.    
C
C            IMAIGNAN (I4) : (S) Indice précisant si les simulations
C      			     introduisent le modèle de réflexion de Maignan :
C			            1 : pour un calcul dde BRDF Maignan (il faut utiliser N0 > 0).
C			            Sinon : pas de prise en compte.    
C
C            INADAL (I4) : (S) Indice précisant si les simulations
C      			   introduisent le modèle de réflexion de Nadal :
C			          1 : pour un calcul dde BRDF Nadal (il faut utiliser N0 > 0).
C			          Sinon : pas de prise en compte.    
C
C            ALPHA (DOUBLE) : (S) Tableau 1D de CTE_OS_NB_MAX + 1 éléments    
C				      contenant les coefficients Alpha de la décomposition en fonctions 
C				      de Legendre des éléments Pij de la fonction de phase  			               
C
C            BETA (DOUBLE) : (S) Tableau 1D de CTE_OS_NB_MAX + 1 éléments    
C				     contenant les coefficients Beta de la décomposition en fonctions 
C				     de Legendre des éléments Pij de la fonction de phase  			               
C
C            GAMMA (DOUBLE) : (S) Tableau 1D de CTE_OS_NB_MAX + 1 éléments    
C				      contenant les coefficients Gamma de la décomposition en fonctions 
C				      de Legendre des éléments Pij de la fonction de phase  			               
C
C            ZETA (DOUBLE) : (S) Tableau 1D de CTE_OS_NB_MAX + 1 éléments    
C				     contenant les coefficients Zeta de la décomposition en fonctions 
C				     de Legendre des éléments Pij de la fonction de phase  			               
C
C            IND_ANGOUT_FIN (DOUBLE) : (S) Indices des tableaux d'angles pour lesquels il doit y avoir  
C                                      une sortie dans le fichier résultat des angles utilisateur
C      				    (1 : pour une sortie requise)
C                                      Indices à partir de 0 pour passage à Python.
C				 				
C            IER (I4) : (S) code d'erreur =0 si pas d'erreur, =-1 sinon
C 
C
C    -->Fichiers
C       Le fichier des angles et paramètres internes pour les limites des développements.
C       Le fichier Matrice de réflexion de surface
C       Le fichier Aérosols
C
C
C Résultats fournis par le programme
C ----------------------------------
C   Contenu du fichier des angles
C   Contenu du fichier aérosols       
C   Type de modelisation de la réflexion de surface
C   Nom du fichier de surface
C   Fichier de trace (optionnel)
C              (No d'unite logique du fichier trace =IDLOG)
C
C
C
C 
C Description des fichiers utilisés
C ---------------------------------
C
C   -->Contenu du fichier BRDF/BPDF
C       Fichier binaire non formatté, créé par SOS_SURFACE
C       Il contient la matrice de réflexion de surface 
C       (fichier Glitter ou autres BRDF/BPDF):
C       c'est à dire les éléments Pij de la matrice de réflexion :
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
C        --> Enregistrements sur le domaine limite aux angles utiles
C                                   --------------------------------
C
C
C    -->Contenu du fichier des angles:
C       Fichier ascii formatté, créé par SOS_ANGLES.
C       Il contient l'information sur les angles à utiliser 
C       et sur les ordres des développements adaptes : 
C		NB_TOTAL_ANGLES :    nb d'angles total à utiliser
C             NB_GAUSS_ANGLES :    nb d'angles de Gauss
C             ANGLES_USERFILE :    Nom du fichier des angles utilisateur
C	       SOLAR ZENITH ANGLE : Angle zénithal solaire (degrés)
C		INTERNAL_IMUS :      Numéro d'indice de table correspondant
C             INTERNAL_OS_NB :   Ordre des développements des fonctions de phase 
C				     en Polynomes de Legendre à utiliser.
C             INTERNAL_OS_NS :   Ordre maximal à utiliser pour le développement en
C                                polynomes de Legendre des éléments de la matrice de Fresnel 
C                                et en séries de Fourier des calculs de luminance.
C                                (non utilisé dans ce programme)
C             INTERNAL_OS_NM :   Ordre maximal du développement en séries de Fourier  
C                                de la fonction G des calculs de matrice de réflexion.
C                                (non utilisé dans ce programme)
C             INDEX COS_ANGLE WEIGHT OUTPUT: lignes donnant le numéro d'angle, son cosinus 
C                                            et son poids, s'il est donné en sortie "utilisateur"
C                                            au format I4,1X,2D21.14,1X,I4
C         
C   -->Contenu du fichier de MIE
C       Fichier binaire non formatté, créé par SOS_MIE
C       1 enregistrement d'entête avec
C                Indice de réfraction (partie réelle)
C                Indice de réfraction (partie imaginaire)
C                Paramètre de taille limite
C                Nombre d'angles de Gauss positifs
C
C      1 enregistrement pour chaque valeur du paramètre de taille avec:
C        - Paramètre de taille (réel),
C        - Coefficient d'efficacite en extinction (réel),
C        - Coefficient d'efficacite en diffusion (réel),
C        - Facteur d'anisotropie (double précision),
C        - Tableau de la fonction de phase Mie(mu) (réel),
C        - Tableau de la fonction de phase polarisee Qmie(mu) (réel),
C        - Tableau de la fonction de phase polarisee Umie(mu) (réel).
C
C    --> Les tableaux ont une taille utile limitée à MIE_NBMU : 
C        ils sont créés et lus avec cette limitation.
C
C
C Common utilise:
C --------------
C  Aucun
C
C
C Cas d'erreur :
C ------------        
C     -Décodage des paramètres incorrect
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
C   - SOS_NOM_FICGLITTER
C   - SOS_NOM_FICROUJ
C   - SOS_NOM_FICRH
C   - SOS_NOM_FICBREON
C   - SOS_NOM_FICNADAL
C   - SOS_NOM_FICMAIGNAN
C
C============================================================================
       SUBROUTINE SOS_PREPA_OS(FICANGLES, FICGRANU, 
     &	       ISURF, IND_SURF,WIND, 
     &		K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN, 
     &		ALPHA_NADAL, BETA_NADAL, COEF_C_MAIGNAN,
     &		REPSURF,FICSURF,TRACE,IDLOG,
     &		LUM_NBMU,N0,PIZ,PIZTR,A,RMU,GA,OS_NB,IMAT_SURF,
     &		IGLI,IFRESNEL,IROUJEAN,IRONDEAUX,IBREON,INADAL,IMAIGNAN,
     &		ALPHA,BETA,GAMMA,ZETA,IND_ANGOUT_FIN,IER)      
      
      IMPLICIT NONE

C* Définition des variables                         
C*-----------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                   !Cosinus des angles de Gauss (et utilisateur).
      DOUBLE PRECISION GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                   !Poids des angles de Gauss (et utilisateur).				
 
 
      DOUBLE PRECISION TETAS	! Angle zénithal solaire en degrés (0 < tetas < 90).
      				       ! (angle complémentaire : mus= -cos(tetas))
      
      DOUBLE PRECISION IND_SURF   ! Indice de réfraction de la surface / air.
      
      DOUBLE PRECISION WIND       ! Vitesse du vent en m/s.
      
      DOUBLE PRECISION K0_ROUJEAN ! Paramètre k0 du modèle de BRDF de Roujean.
      DOUBLE PRECISION K1_ROUJEAN ! Paramètre k1 du modèle de BRDF de Roujean.
      DOUBLE PRECISION K2_ROUJEAN ! Paramètre k2 du modèle de BRDF de Roujean.
      
      DOUBLE PRECISION ALPHA_NADAL ! Paramètre alpha du modèle de BPDF de Nadal.
      DOUBLE PRECISION BETA_NADAL  ! Paramètre beta du modèle de BPDF de Nadal.

      DOUBLE PRECISION COEF_C_MAIGNAN !Paramètre C*exp(-NDVI) du modèle de Maignan.
            
      DOUBLE PRECISION A	! Coefficient de troncature.
      DOUBLE PRECISION PIZ	! Albédo de simple diffusion.
      DOUBLE PRECISION PIZTR	! Albédo de simple diffusion après troncature.

      DOUBLE PRECISION ALPHA(0:CTE_OS_NB_MAX)	! Paramètres alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NB_MAX)	! Paramètres beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NB_MAX)	! Paramètres gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Paramètres zeta(k).
								               
     
      INTEGER*4 N0	! Numéro de mu de Gauss pour l'angle
      			! zénithal solaire.
			! (angle complémentaire : mus = -rmu(N0))
			! (N0 > 0 => pour en tenir compte).

      INTEGER*2 ISURF	! Type de surface :
			       ! 0 :surface lambertienne d'albédo RHO.
			       ! 1 :surface lambertienne + glitter.
			       ! 2 :surface lambertienne + réflexion de Fresnel.
			       ! 3 :surface lambertienne + BRDF de Roujean.
			       ! 4 :surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
			       ! 5 :surface lambertienne + BRDF de Roujean + BPDF de Breon.
			       ! 6 :surface lambertienne + BRDF de Roujean + BPDF de Nadal.
				
      INTEGER*4 IGLI	! Indice précisant si les simulations
      			! introduisent le Glitter :
			!    1 => pour calcul du Glitter (il faut alors utiliser N0 > 0).
			!    Sinon => pas de simulation du Glitter.            
	
      INTEGER*4 IFRESNEL ! Indice précisant si les simulations
      			    ! introduisent la réflexion de Fresnel sur mer plate :
			    !    1 => pour un calcul de la réflexion de Fresnel.
			    !    Sinon => pas de prise en compte.            
	
      INTEGER*4 IROUJEAN ! Indice précisant si les simulations
      			    ! utilisent le modèle de réflexion de Roujean :
			    !    1 => pour calcul de BRDF Roujean (il faut utiliser N0 > 0).
			    !    Sinon => pas de simulation de BRDF Roujean. 
			 
      INTEGER*4 IRONDEAUX ! Indice précisant si les simulations
      			     ! utilisent le modèle de réflexion de Rondeaux - Herman :
			     !    1 => pour calcul de BPDF RH (il faut utiliser N0 > 0).
			     !    Sinon => pas de simulation de BPDF RH. 			 
			
      INTEGER*4 IBREON   ! Indice précisant si les simulations
      			    ! utilisent le modèle de réflexion de Breon :
			    !    1 => pour calcul de BPDF Breon (il faut utiliser N0 > 0).
			    !    Sinon => pas de simulation de BPDF Breon. 

      INTEGER*4 INADAL   ! Indice précisant si les simulations
      		           ! utilisent le modèle de réflexion de Nadal :
		           !    1 => pour calcul de BPDF Nadal (il faut utiliser N0 > 0).
			    !    Sinon => pas de simulation de BPDF Nadal. 	

      INTEGER*4 IMAIGNAN ! Indice précisant si les simulations
      		           ! utilisent le modèle de réflexion de Maignan :
		           !    1 => pour calcul de BPDF Nadal (il faut utiliser N0 > 0).
			    !    Sinon => pas de simulation de BPDF Maignan. 
			 
	
      INTEGER*4 IMAT_SURF ! Indice précisant si les simulations
      			     ! utilisent un modèle de réflexion BRDF ou BPDF :
			     !    1 => pour calcul de BRDF ou BPDF (il faut utiliser N0 > 0).
			     !    Sinon => pas de simulation de BRDF ou BPDF. 	
			  			  			 											

      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon

      INTEGER*4 K       ! Indice de tableau
      INTEGER*4 I,J     ! Indice de tableau
      
      INTEGER*4 LDIRSURF ! longueur de la chaine de caractère REPSURF.
      INTEGER*4 LFICSURF ! longueur de la chaine de caractère FICSURF.
      
      INTEGER*4 NBMU_GAUSS   ! Nombre d'angles de Gauss utilisés       
      INTEGER*4 LUM_NBMU     ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB 	 ! Ordre des développements des fonctions de phase 
			        ! en Polynomes de Legendre.      

      INTEGER*4 IND_ANGOUT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                            ! Indices des tableaux d'angles pour lesquels il doit y avoir 
      				! une sortie dans le fichier résultat des angles utilisateur 
				! (1 : pour une sortie requise)	
			
      INTEGER*4 IND_ANGOUT_FIN(0:CTE_OS_NBMU_MAX)  
                            ! Indices des tableaux d'angles pour lesquels il doit y avoir une sortie dans le fichier résultat des angles utilisateur 
                            ! ==> 1 : pour une sortie requise

      INTEGER*4 IFIN        ! Nb de caractères significatifs d'une chaine. 

			               
      
      LOGICAL TRACE   ! = vrai,si écriture dans le fichier trace
      LOGICAL EX	! Vérifie l'existence d'un fichier

      CHARACTER*CTE_LENDIR REPSURF	! Répertoire de stockage des matrice de réflexion de
      				       ! surface (fichier Glitter ou autres BRDF/BPDF).
      CHARACTER*CTE_LENFIC2 FICSURF ! Fichier contenant la matrice de réflexion de
      				        ! surface (fichier Glitter ou autres BRDF/BPDF).
      CHARACTER*CTE_LENFIC2 FICANGLES           ! chemin et nom du fichier des angles effectifs    
      CHARACTER*CTE_LENFIC2 FICANGLES_USER      ! chemin et nom du fichier des angles "utilisateurs"
      
      CHARACTER*CTE_LENFIC2 FICGRANU	! nom du fichier d'entrée AEROSOLS

      CHARACTER*CTE_LENFIC1 FICS ! nom recomposé du fichier SURFACE à partir 
      				     ! des arguments de simulation de la BRDF / BPDF.

      CHARACTER*LENLINE LIGNE_TEXTE    ! Ligne de caractères pour décodage de paramètres      
      CHARACTER*1000 CVAL	           ! Chaine de caractères pour la lecture de valeurs inutilisées
             
      INTEGER*4 IDLOG	     ! Numéro identifiant du fichier Trace	    
     	  
      	        		 

C* Initialisation
C-----------------
      IER=0
      IGLI=0
      IFRESNEL=0
      IROUJEAN=0
      IRONDEAUX=0
      IBREON=0
      INADAL=0
      IMAIGNAN=0
      IMAT_SURF=0
      IF (ISURF.EQ.1) IGLI=1
      IF (ISURF.EQ.2) IFRESNEL=1
      IF (ISURF.GE.3) IROUJEAN=1
      IF (ISURF.EQ.4) IRONDEAUX=1
      IF (ISURF.EQ.5) IBREON=1
      IF (ISURF.EQ.6) INADAL=1
      IF (ISURF.EQ.7) IMAIGNAN=1
      IF ((IGLI.EQ.1).OR.(IROUJEAN.EQ.1).OR.(IRONDEAUX.EQ.1).
     &     OR.(IBREON.EQ.1).OR.(INADAL.EQ.1).
     &     OR.(IMAIGNAN.EQ.1)) IMAT_SURF=1


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
         READ(10,*,ERR=951) CVAL, RMU(J), GA(J), IND_ANGOUT(J)
      ENDDO 

      
      DO J=1,LUM_NBMU
         RMU(-J) = -RMU(J)
	  GA(-J) =  GA(J)
	  IND_ANGOUT(-J) = IND_ANGOUT(J)
	  IND_ANGOUT_FIN(J-1) = IND_ANGOUT(J)	!Passage à indexation à partir de 0 pour Python
      ENDDO 
      
      !Controle d'applicabilite
      IF (ABS(N0).GT.LUM_NBMU) GOTO 1000
      IF ((TETAS.LT.0.D-00).OR.(TETAS.GE.90.D-00)) GOTO 1001
      IF ((ISURF.NE.0).AND.(N0.LT.0)) GOTO 1012
       
      IF (TRACE) THEN
        WRITE (IDLOG,*,ERR=960) "USED ANGLES"
	 WRITE (IDLOG,*,ERR=960) "***********"
	 WRITE (IDLOG,*,ERR=960) " "
	 WRITE (IDLOG,*,ERR=960) "Total number of angles : ",LUM_NBMU
	 WRITE (IDLOG,*,ERR=960) "Number of Gauss angles : ",NBMU_GAUSS
	 WRITE (IDLOG,*,ERR=960) "Solar incidence : "
	 WRITE (IDLOG,*,ERR=960) " --> Solar zenith angle: ", TETAS
	 WRITE (IDLOG,*,ERR=960) " --> Index in the angles table: ", N0 
	 WRITE (IDLOG,*,ERR=960) " "
        IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	     WRITE(IDLOG,*,ERR=960)  
     &	       "A file of user's angles is used :", FICANGLES_USER
     	     WRITE(IDLOG,*,ERR=960) 
     &       "Number of user's angles + solar angle SZA :",
     &         (LUM_NBMU-NBMU_GAUSS)
	 ELSE
	     WRITE(IDLOG,*,ERR=960)  
     &	       "No file of user's angles"
	 ENDIF !Fin test Pas de fichier d'angles utilisateur
        WRITE (IDLOG,*,ERR=960) " " 
	 WRITE (IDLOG,*,ERR=960) "Max order NB ", OS_NB
	 WRITE (IDLOG,*,ERR=960) " " 
	 WRITE (IDLOG,*,err=960) "List of used angles :"
	 WRITE (IDLOG,*,err=960)  
     &  "Index        Mu             angle (deg)    Weigth"// 
     &  "         user's angle"
	 DO 1200 J=-LUM_NBMU,LUM_NBMU
	    IF (J.EQ.0) GOTO 1200
	    IF (J.LT.0) THEN 
	        WRITE(IDLOG,223,ERR=960) J, RMU(J),
     &	                     -DACOS(-RMU(J))*180./INCTE_PI, GA(J),
     &                       IND_ANGOUT(-J)
            ELSE
	        WRITE(IDLOG,223,ERR=960) J, RMU(J),
     &	                      DACOS(RMU(J))*180./INCTE_PI, GA(J), 
     &      		        IND_ANGOUT(J) 
	    ENDIF
1200    CONTINUE	
	 WRITE (IDLOG,*,ERR=960) " " 
	 WRITE (IDLOG,*,ERR=960) " "  	
      ENDIF !Fichier Trace	 

	     	
      
C* Détermination du fichier de BRDF ou BPDF de surface
C  ---------------------------------------------------

      IF (IMAT_SURF.EQ.1) THEN

C* Restitution du nom du fichier de Réflexion de surface
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
     &                             FICS)
	     ENDIF
	     
	     IF (ISURF.EQ.5) THEN
	         CALL SOS_NOM_FICBREON(NBMU_GAUSS,FICANGLES_USER,
     &	                               IND_SURF,
     &                                K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                                FICS)
	     ENDIF
	     
	     IF (ISURF.EQ.6) THEN
	         CALL SOS_NOM_FICNADAL(NBMU_GAUSS,FICANGLES_USER,
     &                                IND_SURF,ALPHA_NADAL,BETA_NADAL,
     &	                               K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                                FICS)
	     ENDIF

	     IF (ISURF.EQ.7) THEN
	         CALL SOS_NOM_FICMAIGNAN(NBMU_GAUSS,FICANGLES_USER,
     &                                  IND_SURF,COEF_C_MAIGNAN,
     &                                  K0_ROUJEAN,K1_ROUJEAN,
     &	                                 K2_ROUJEAN,FICS) 	
	     ENDIF

	     
	     LDIRSURF = INDEX(REPSURF,' ')
            LDIRSURF = LDIRSURF - 1
            IF(LDIRSURF.LE.0) LDIRSURF = CTE_LENDIR
            FICSURF = REPSURF(1:LDIRSURF)//'/'//FICS
	 
	 ENDIF        
		
C* Controle de l'existence du fichier de Matrice de réflexion
C  ----------------------------------------------------------	
	 
        LFICSURF=INDEX(FICSURF,' ')-1 
        INQUIRE(FILE=FICSURF(1:LFICSURF),EXIST=EX)
        IF (EX.EQV..FALSE.) GOTO 1020   !Cas d'erreur  
	       
      ENDIF ! Fin du test pour l'utilisation du fichier de Matrice de réflexion


      
C* Lecture du fichier AEROSOLS
C  ---------------------------
      OPEN(UNIT=1,FILE=FICGRANU,STATUS='OLD',ERR=913)
	 
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
         READ(1,*,ERR=923) ALPHA(K),BETA(K),GAMMA(K),ZETA(K)
      ENDDO      


C* Calcul de l'albédo de simple diffusion non tronquée 
C* à partir de la valeur tronquée
C  ---------------------------------------------------
      PIZ=PIZTR/(1+0.5*A*(PIZTR-1))
            
     
      
      
C* Fermeture fichier
C-------------------
      CLOSE(1)
      CLOSE(10)

      
C* Format 
C----------
  223 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X,I2)  


C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------    

  950 WRITE(6,*) '  SOS_PREPA_OS : ERROR_950'
      WRITE(6,*) '  ERROR on angles file opening for SOS'
      GOTO 9998

  951 WRITE(6,*) '  SOS_PREPA_OS : ERROR_951'
      WRITE(6,*) '  ERROR on angles file reading for SOS'
      GOTO 9998

  960 WRITE(6,*) '  SOS_PREPA_OS : ERROR_960'      
      WRITE(6,*) '  ERROR on logfile writing for SOS'
      GOTO 9998
             
  913 WRITE(6,*) '  SOS_PREPA_OS : ERROR_913'      
      WRITE(6,*) '  ERROR on AEROSOLS file opening for SOS'
      GOTO 9998    

  922 WRITE(6,*) '  SOS_PREPA_OS : ERROR_922'      
      WRITE(6,*) '  ERROR on PROFILE file reading for SOS'
      GOTO 9998
 
  923 WRITE(6,*) '  SOS_PREPA_OS : ERROR_923'        
      WRITE(6,*) '  ERROR on AEROSOLS file reading for SOS'
      GOTO 9998 

 1000 WRITE(6,*) '  SOS_PREPA_OS : ERROR_1000'                 
      WRITE(6,*) '  ERROR on angles file value for IMUS' 
      WRITE(6,*) '  Out of available range.'
      WRITE(6,*) '  --> Defined value : ',N0
      WRITE(6,*) '      Have to be define between 1 and ',LUM_NBMU
      GOTO 9998

 1001 WRITE(6,*) '  SOS_PREPA_OS : ERROR_1001'                 
      WRITE(6,*) '  ERROR on TETAS value : out of available range'
      GOTO 9998      

 1012 WRITE(6,*) '  SOS_PREPA_OS : ERROR_1012'                 
      WRITE(6,*) '  ERROR on angles file value for IMUS' 
      WRITE(6,*) '  Out of available range.'
      WRITE(6,*) '  --> ISURF > 0 requires IMUS > 0 '
      WRITE(6,*) '      IMUS corresponding to the number '
      WRITE(6,*) '      of the solar angle in the angles file.'
      GOTO 9998     
                                                                                                       
 1020 WRITE(6,*) '  SOS_PREPA_OS : ERROR_1020'                 
      WRITE(6,*) '  ERROR on FICSURF value : file does not exist.'
      WRITE(6,*) '  ISURF: ',ISURF
      WRITE(6,*) '  IMAT_SURF: ',IMAT_SURF     
      GOTO 9998 
        
 9998 IER=-1
    
 9999 RETURN
 
      END	!FIN DU PROGRAMME SOS




