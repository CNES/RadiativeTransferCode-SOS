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
C* FICHIER: SOS_ROUJEAN.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Genere les fichiers de BRDF de Roujean.
C*
C* AUTEUR: 
C*    CS GROUP France, B.Lafrance.
C* DATE: 30/04/2002
C*
C* MOD:VERSION:1.0: Creation du programme
C* MOD:VERSION:1.1: Amelioration des commentaires
C*
C* MOD:VERSION:2.0: 10/01/2010 
C*
C*   - Suppression de la constante CTE_OS_FICANGLE (de SOS.h) et de la lecture du fichier correspondant 
C*     des angles de Gauss au profit : d'un passage par argument des tables d'angles (RMU) et poids (CHR).
C*
C*   - Modification de la taille des tableaux de luminance :
C*     ==> Passage de CTE_OS_NBMU à CTE_OS_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'à LUM_NBMU (nouveau argument d'entrée) 
C*         égal au nombre effectif d'angles utilises.
C*
C*   - Suppresion des constantes CTE_OS_NB  (de SOS.h)  
C*     au profit : d'un passage par argument de la valeur  OS_NB.
C*
C*   - Modification de la taille des tableaux de decomposition en polynomes de Legendre :
C*     ==> Passage de CTE_OS_NB à CTE_OS_NB_MAX.
C*     ==> Boucle sur les ordre jusqu'à OS_NB  (nouveaux argument d'entrée)
C*
C*   - Routine SOS_FSF_ROUJEAN : Modification de l'appel de la routine pour inclure l'argument LUM_NBMU   
C*     (nb d'angles utiles) et pour limiter les boucles aux bornes utiles LUM_NBMU.
C*     et pour inclure l'argument OS_NB.    
C*
C*   - Routine SOS_MISE_FORMAT_RJ : Modification de la ligne d'enregistrement des  
C*     tableaux de la  matrice de reflexion a l'ordre IS de la decomposition 
C*     en series de Fourier :
C*     ==> limitation aux bornes utiles -LUM_NBMU:LUM_NBMU
C*       Rem : pas de modification des bornes d'enregistrement du fichier 
C*             temporaire RES_MAT_REFLEX car il ne s'agit pas d'un produit final.
C*     Utilisation de OS_NB.
C*
C* MOD:VERSION:2.1: 06/06/2013
C*   -  Modification de la valeur de Pi, declaree en constante dans le fichier
C*      (INCTE_PI remplace CTE_PI)
C*
C* MOD:VERSION:2.2: 24/02/2016 
C*    - Ajustements mineurs pour le respect strict des 72 colonnes
C*
C* MOD:VERSION:2.3: 18/05/2020
C*
C*   - Ajout d'un test dans la routine SOS_CALC_F_ROUJEAN 
C*     pour le calcul du cosinus de l'angle Zeta de sorte à éviter des 
C*     arrondis numériques donnant une valeur 1+Epsilon ou -1-Epsilon.
C*
C*   - Introduction d'une nouvelle constante de seuil : CTE_SEUIL_NUM
C*     (définie dans inc/SOS.h)
C* 
C*   - Modification des valeurs stockées de la fonction de BRDF F 
C*     (par la routine SOS_CALC_F_ROUJEAN) : 
C*     ==> On stocke F * cos(Angle Réflexion) * cos(Angle Incidence)
C*         sans introduire de seuil sur les angles pour le calcul des cosinus
C*         de pondération. En effet, les calculs de flux incident en surface
C*         dans le coeur de code des OS (SOS_OS) nécessitent d'avoir les angles
C*         vrais.
C*
C* MOD:VERSION:2.4: 31/08/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Conversion des traces en Anglais.
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    CTE_LENDIR : Longueur des noms de repertoires.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_OS_NBMU_MAX : Nombre maximal d'angles positifs.
C    INCTE_PI : valeur de PI
C    CTE_OS_NB_MAX :  Valeur maximale pour l'ordre limite du développement en 
C                     polynomes de Legendre.
C    CTE_PH_NU : Nombre d'angles azimutaux (2**NQ)
C    CTE_TETAS_LIM_ROUJEAN : Angle zenithal solaire limite pour les calculs
C                            de BRDF du modèle de Roujean.
C    CTE_TETAV_LIM_ROUJEAN : Angle zenithal de visee limite pour les calculs
C                            de BRDF du modèle de Roujean.
C    CTE_SEUIL_SF_ROUJEAN : Valeur seuil pour l'arret de la decomposition 
C                           en series de Fourier de la BRDF de Roujean.
C    CTE_SEUIL_NUM : Valeur de seuil de précision numérique
C
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROCEDURE: SOS_ROUJEAN
C ==========      
C      Cette procédure genere les fichiers de BRDF du modèle de Roujean.
C
C
C Donnees en entrée du programme 
C ------------------------------
C    -->Les paramètres/arguments du programme
C       LUM_NBMU (I4) (E) : Nombre d'angles (positifs) effectivement utiles
C	 RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) (E) : Cosinus des angles.
C	 CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) (E) : Poids des angles.
C       OS_NB (I4) : (E) Ordre des développements des fonctions de phase 
C		           en Polynomes de Legendre a utiliser.
C 	 K0  (DOUBLE)  (E) : Coefficient du modèle de Roujean.         	
C 	 K1  (DOUBLE)  (E) : Coefficient du modèle de Roujean.         	
C 	 K2  (DOUBLE)  (E) : Coefficient du modèle de Roujean.    
C       FIC_RES_MAT_REFLEX (CHARACTER*CTE_LENFIC2)  :  
C          (E) Nom du fichier temporaire RES_MAT_REFLEX avec arborescence complète.     	
C       FICROUJ  (CHARACTER*CTE_LENFIC2) (E) : nom complet du fichier Roujean 
C                                              (arborescence complete)
C       TRACE  (LOGICAL) (E) : vrai,si écriture dans le fichier trace.
C       IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon
C
C    -->Fichier
C       Le fichier d'angles de Gauss 
C
C Resultat fourni par le programme
C --------------------------------
C   Fichier BRDF de Roujean
C   Fichier de trace (optionnel), 
C               si un fichier de meme nom existe deja, il est ecrase
C              (No d'unite logique du fichier trace =99)
C
C 
C Variable d'environnement
C ------------------------
C   Aucune
C
C
C Description des fichiers utilises
C ---------------------------------
C
C   -->Contenu du fichier de BRDF de Roujean
C       Fichier binaire non formatte, cree par SOS_MISE_FORMAT_RJ
C       Il contient la matrice de reflexion de surface :
C       c'est a dire les elements Pij de la matrice de reflexion 
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
C   -->Contenu du fichier temporaire RES_MAT_REFLEX
C       Fichier binaire non formatte,
C       Fichier cree dans le repertoire de lancement du programme
C       Si un fichier de meme nom existe deja, il est ecrase
C       Enregistrement pour chaque couple (Teta1, Teta2) de la 
C       decomposition en series de Fourier de la BRDF de Roujean :
C       M11 (n valeurs réelles)
C
C    Le fichier temporaire est ecrase s'il existe puis 
C    detruit en fin de programme
C
C
C Common utilise:
C --------------
C  Aucun
C
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C      Erreur si le fichier BRDF de Roujean existe deja
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, arret du programme et 
C     retour du status 1 au shell
C        
C   
C
C  Sous programmes utilises:
C --------------------------
C  Cette procédure fait appel aux routines:
C      - SOS_FSF_ROUJEAN
C      - SOS_MISE_FORMAT_RJ
C
C
C=============================================================================

      SUBROUTINE SOS_ROUJEAN(LUM_NBMU,RMU,CHR,OS_NB,K0,K1,K2,
     &                       FIC_RES_MAT_REFLEX,FICROUJ,TRACE,IER)

      IMPLICIT NONE


C* Definition des variables   
C*--------------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICROUJ  ! (E) Nom complet du fichier ROUJEAN resultat.

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Nom du fichier RES_MAT_REFLEX
                                               ! avec arborescence complète

      LOGICAL TRACE         ! (E) vrai si écriture dans le fichier trace 

      DOUBLE PRECISION K0	! (E) Paramètre du modèle de Roujean.
      DOUBLE PRECISION K1	! (E) Paramètre du modèle de Roujean.
      DOUBLE PRECISION K2	! (E) Paramètre du modèle de Roujean.
      	
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Tableau des mu de Gauss
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Poids de Gauss.
      
      
      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB    ! Ordre des développements des fonctions de phase 
			    ! en Polynomes de Legendre.
      INTEGER*4 I	   ! Indice de tableau.
      INTEGER*4 IFIN    ! nb caracteres significatifs d'une chaine
      INTEGER*4 IER     ! (S) code d'erreur =0 si pas d'erreur, =-1 sinon

C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0

C* Ouverture du fichier Trace 
C------------------------------------------------------
      IF (TRACE) THEN
         WRITE(99,*,err=921)
	  WRITE(99,*,err=921) '*****************************************'
         WRITE(99,*,err=921) 'Calculation of the ROUJEAN file :' 
	  WRITE(99,*,err=921) FICROUJ 
         WRITE(99,*,err=921) ' '
         WRITE(99,*,err=921)'Simulation parameters : '
         WRITE(99,*,err=921)' k0 =  ',K0
	  WRITE(99,*,err=921)' k1 =  ',K1
	  WRITE(99,*,err=921)' k2 =  ',K2
	  WRITE(99,*,err=921)
	  WRITE(99,*,err=921)' Limit of viewing zenith angle = ',
     &	                     CTE_TETAV_LIM_ROUJEAN
	  WRITE(99,*,err=921)' Limit of solar zenith angle = ',
     &	                     CTE_TETAS_LIM_ROUJEAN
         WRITE(99,*,err=921)
	  WRITE(99,*,err=921)	 	 
      END IF
        


C* Decomposition en series de Fourier de la fonction F(Teta1,Teta2,phi)
C* et enregistrement de la BRDF developpee en series de Fourier.
C* (BRDF convertie en luminance normalisée et adaptee au code des OS) 
C* => Creation du fichier temporaire RES_MAT_REFLEX.
C---------------------------------------------------------------------
      CALL SOS_FSF_ROUJEAN(LUM_NBMU,RMU,OS_NB,K0,K1,K2,
     &                     FIC_RES_MAT_REFLEX,TRACE,IER)
      IF (IER.NE.0) GOTO 995   ! cas d'erreur      
      
      
C* Rearrangement des resultats pour etre compatible avec le code des OS.
C-----------------------------------------------------------------------
      CALL SOS_MISE_FORMAT_RJ(FIC_RES_MAT_REFLEX,FICROUJ,
     &                        LUM_NBMU,OS_NB,IER)
      IF (IER.NE.0) GOTO 996   ! cas d'erreur

C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  912 WRITE(6,*) 'Erreur a l''ouverture du fichier angles de Gauss'
      IER=-1
      GOTO 9999
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Surface'
      IER=-1
      GOTO 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1
      GOTO 9999      
  995 WRITE(6,*) 'Erreur dans la routine SOS_FSF_ROUJEAN'
      IER=-1
      GOTO 9999
  996 WRITE(6,*) 'Erreur dans la routine SOS_MISE_FORMAT_RJ'
      IER=-1
      GOTO 9999  
      
C* Format
C---------   
 2    FORMAT(2(D21.14))    
              
 9999 CONTINUE

C* Destruction des fichiers temporaires
C---------------------------------------   
      OPEN(UNIT=3,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED')
      CLOSE(3,STATUS='DELETE') 
                                                      
      RETURN
      END    !FIN DE LA PROCEDURE SOS_ROUJEAN 





C==============================================================================
C PROCEDURE: SOS_FSF_ROUJEAN
C ==========
C      Cette procédure realise la decomposition en series de Fourier
C      (sur phi) de la fonction F(Teta1,Teta2,phi) :
C
C	 F(Teta1,Teta2,phi) = K0 + K1*F1(Teta1,Teta2,phi) + K2*F2(Teta1,Teta2,phi)
C
C	 avec Teta1 : angle zenithal d'incidence : tetas  (0 < Teta1 < PI/2),
C             Teta2 : angle zenithal de reflexion : tetav (0 < Teta2 < PI/2),
C             phi   : angle azimutal relatif entre les directions
C                     d'incidence et de reflexion selon la convention 
C                     du modèle de Roujean (0 <= phi <= PI),
C
C	      K0, K1, K2 : coefficients du modèle de Roujean,
C			
C             F1 = [ (PI-phi)*cos(phi)+sin(phi) ]*tgs*tgv / 2*PI
C                 -[ tgs + tgv + SQRT{tgs**2+tgv**2-2*tgs*tgv*cos(phi)} ] / PI
C
C		  en posant tgs = tan(tetas) et tgv = tan(tetav),
C
C            F2 = [ (PI/2-zeta)*cos(zeta)+sin(zeta) ] * (4/3*PI) 
C                 / [ cos(tetas)+cos(tetav) ]   - 1/3
C
C                 avec zeta l'angle de phase defini par :
C                 cos(zeta) = cos(tetas)*cos(tetav)+sin(tetas)*sin(tetav)*cos(phi). 
C
C
C      La BRDF du modèle de Roujean est convertie de reflectance en une 
C      grandeur utile pour les OS : luminance normalisée * cos(Teta2),
C      c'est a dire : BRDF * cos(Teta1) * cos(Teta2).			
C
C
C      La procédure enregistre les resultats dans le fichier "RES_MAT_REFLEX" 
C      sous la forme:
C         Fonction F decomposee pour (Teta1,Teta2) : (E(IS),IS=0,CTE_OS_NB)
C
C
C   REM : phiOS = PI - phiRoujean
C   *****************************
C
C Description des paramètres
C --------------------------
C      LUM_NBMU (I4) (E) : Nombre d'angles (positifs) effectivement utiles.
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)	 (E) : Tableau des mu de Gauss.
C      OS_NB (I4) : (E) Ordre des développements des fonctions de phase 
C		        en Polynomes de Legendre a utiliser.
C      K0	(double)  (E) : Coefficient du modèle de Roujean.
C      K1	(double)  (E) : Coefficient du modèle de Roujean.
C      K2	(double)  (E) : Coefficient du modèle de Roujean.
C      FIC_RES_MAT_REFLEX (CHARACTER)  : (S) Nom du fichier temporaire RES_MAT_REFLEX
C                                            avec arborescence complète.
C      TRACE    (logical) (E) : = vrai,si écriture dans le fichier trace (unite 99)
C      IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C Autre resultat fourni par la procédure
C -------------------------------------
C      Fichier temporaire RES_MAT_REFLEX
C
C 
C Description des fichiers utilises
C ---------------------------------
C
C   -->Contenu du fichier temporaire RES_MAT_REFLEX
C       Fichier binaire non formatte,
C       Fichier cree dans le repertoire de lancement du programme
C       Si un fichier de meme nom existe deja, il est ecrase
C       Enregistrement pour chaque couple (Teta1, Teta2) de la 
C       decomposition en series de Fourier de la BRDF de Roujean :
C       M11 (n valeurs réelles)
C
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C     -Lecture/écriture dans un fichier
C     -BRDF < 0
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C     
C  Routines utilises par la procédure:
C ------------------------------------
C  Ce programme fait appel aux routines:
C        - SOS_CALC_F_ROUJEAN
C
C==============================================================================
      SUBROUTINE SOS_FSF_ROUJEAN(LUM_NBMU,RMU,OS_NB,K0,K1,K2,
     &                           FIC_RES_MAT_REFLEX,TRACE,IER)
 	
      IMPLICIT NONE

C* Definition des variables    
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Tableau des mu de Gauss

      DOUBLE PRECISION K0	! Coefficient du modèle de Roujean.
      DOUBLE PRECISION K1	! Coefficient du modèle de Roujean.
      DOUBLE PRECISION K2	! Coefficient du modèle de Roujean.

      DOUBLE PRECISION E(0:CTE_OS_NB_MAX) ! Resultat de la decomposition en 
      				              ! series de Fourier de la fonction F.
	
      DOUBLE PRECISION PHIOS	! Azimut relatif pour la convention OS (radians).
      DOUBLE PRECISION PHIRJ	! Azimut relatif pour la convention Roujean (radians).
      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).

      DOUBLE PRECISION F	! Valeur de la fonction F(phiRJ).

      DOUBLE PRECISION FMIN	! Valeur minimale retenue pour la fonction 
      				! F(phiRJ) pour Theta1 et Theta2 fixes.
      DOUBLE PRECISION FMAX	! Valeur maximale de la fonction F(phiRJ).

      DOUBLE PRECISION U(0:CTE_PH_NU)  ! Fonction F(phiRJ) entre 0 et PI
      			                  ! pour Theta1 et Theta2 fixes.

      DOUBLE PRECISION Q	! Pas en phi.

      DOUBLE PRECISION T1	! Recombinaison de la fonction F(phi) a 
      				! partir de E(IS).

      DOUBLE PRECISION Y
      DOUBLE PRECISION B1	! Valeur test pour l'arret de la decomposition
      				! en series de Fourier.
      DOUBLE PRECISION B1_PREC  ! Valeur de B1 pour l'ordre IS-1.

      DOUBLE PRECISION PHIFMIN	!Azimut relatif donnant Fmin (selon la convention OS).
      DOUBLE PRECISION PHIFMAX	!Azimut relatif donnant Fmax (selon la convention OS).
     
	
      DOUBLE PRECISION FRECOMB	! Valeur de la fonction F recombinee à partir de sa
      				       ! decomposition en series de Fourier.
      DOUBLE PRECISION ECART	! Ecart relatif entre la valeur exacte de F  
      				       ! et celle recombinee.			
      DOUBLE PRECISION ECARTMAX ! Ecart relatif maximal entre la valeur exacte de F  
      				    ! et celle recombinee.
      
      DOUBLE PRECISION PHIECARTMAX ! Angle azimutal pour la configuration donnant ECARTMAX.
      DOUBLE PRECISION TVECARTMAX  ! Angle zenithal de visee pour ECARTMAX.	
      DOUBLE PRECISION TSECARTMAX  ! Angle zenithal solaire pour ECARTMAX.	 
      DOUBLE PRECISION F_ECARTMAX  ! Valeur exacte de F donnant ECARTMAX.
      DOUBLE PRECISION FRECOMB_ECARTMAX  ! Valeur recombinee de F donnant ECARTMAX.
      DOUBLE PRECISION SOMECART	! Sommation des ecarts relatifs.

      DOUBLE PRECISION ECARTMAX_C1C2 ! Ecart relatif maximal entre la valeur exacte de F  
      				         ! et celle recombinee 
				         !=> pour le couple (Teta1, Teta2) en cours.
      DOUBLE PRECISION PHIECARTMAX_C1C2 ! Angle azimutal pour ECARTMAX_C1C2.
      DOUBLE PRECISION F_ECARTMAX_C1C2  ! Valeur exacte de F donnant ECARTMAX_C1C2.
      DOUBLE PRECISION FRECOMB_ECARTMAX_C1C2  ! Valeur recombinee de F donnant ECARTMAX_C1C2.      
            
      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations 
      INTEGER*4 OS_NB ! Ordre des développements des fonctions de phase 
			 ! en Polynomes de Legendre.
 	    
      INTEGER*4 IL	! Ordre maximale utile pour la fonction E.
      INTEGER*4 I1	! Indice de numéro d'angle Teta1.
      INTEGER*4 I2	! Indice de numéro d'angle Teta2.
      INTEGER*4 IS	! Ordre IS de la decomposition en series de Fourier.
      INTEGER*4 IS2	! Ordre IS de la decomposition en series de Fourier.
      INTEGER*4 IER  ! code d'erreur =0 si pas d'erreur, =-1 sinon  
      INTEGER*4 I
      INTEGER*4 COMPT	! Nombre de configurations pour la sommation des ecarts relatifs.

      	! Elément M11 de la matrice de reflexion en fonction de l'ordre IS
	! de la serie de Fourier, pour une une incidence sous Teta1 et une 
	! reflexion sous Teta2
      REAL M11(0:CTE_OS_NB_MAX)		
      
      LOGICAL TRACE

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX 
                     ! Nom du fichier RES_MAT_REFLEX
                     ! avec arborescence complète
     
	
C* Initialisation
C---------------------------------------------------------	
      ECARTMAX = -1.D+300
      SOMECART = 0.
      COMPT = 0
      

C* Creation du fichier resultat de SOS_FSF_ROUJEAN
C---------------------------------------------------------
      OPEN(UNIT=4,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED',ERR=991)

C* Boucle sur l'angle Teta1 > 0.
C---------------------------------
      DO  3 I1=1 ,LUM_NBMU

C* Cosinus et sinus de Teta1 : angle d'incidence
C---------------------------------
         C1=RMU(I1)
         S1=DSQRT(1-C1*C1)

C*   Boucle sur l'angle teta2 > 0.
C---------------------------------
         DO 4 I2=1 ,LUM_NBMU	

C*   Cosinus et sinus de Teta2 : angle de reflexion
C---------------------------------
           C2=RMU(I2)
           S2=DSQRT(1-C2*C2)

C*   Determination de F pour phi = 0 a PI par pas PI/CTE_PH_NU
C-------------------------------------------------------------         
           Q=INCTE_PI/CTE_PH_NU
           DO 6 I=0,CTE_PH_NU
              PHIOS=Q*I
	       PHIRJ = INCTE_PI - PHIOS	!Adaptation convention Roujean
              CALL SOS_CALC_F_ROUJEAN(K0,K1,K2,C1,S1,C2,S2,PHIRJ,F)     
              U(I)=F
	      
	      IF(F.LT.0.D+00) GOTO 993	!BRDF < 0 => cas d'erreur.
 6         CONTINUE

		
C*   Determination des valeurs min et max de F
C*   pour I1 et I2 fixes mais PHI variant => pour fichier Trace
C--------------------------------------------------------------
	   IF (TRACE) THEN
	       FMAX=-1.D+300
	       FMIN=1.D+300
	       DO I=0,CTE_PH_NU
	          IF(U(I).GT.FMAX) THEN
                   FMAX = U(I)
		     PHIFMAX = Q*I		!Azimuth convention OS
                 ENDIF   
	          IF(U(I).LT.FMIN) THEN
		     FMIN = U(I)
		     PHIFMIN = Q*I		!Azimuth convention OS
	          ENDIF   	 
	       ENDDO

	       WRITE(99,1000,err=921) ACOS(C1)*180./INCTE_PI, 
     &	                              ACOS(C2)*180./INCTE_PI,
     &                               PHIFMIN*180./INCTE_PI,FMIN,
     &				         PHIFMAX*180./INCTE_PI,FMAX	          	 
           ENDIF !Trace 
	      


	 
C*   Realisation de la decomposition en serie de Fourier.
C*   Boucle sur l'ordre de la serie.                     
C--------------------------------------------------------------------

C--------- Initialisations 
C--------------------------------       
          DO IS=0,OS_NB
	      E(IS) = 0.D+00
	   ENDDO
	   
	   B1_PREC = 1.D+300		


C--------- Boucle sur l'ordre IS 
C-------------------------------- 	   
	   DO 13 IS=0,OS_NB
              

C--------- Determination de E(IS) 
C--------------------------------      
	      Y = 0.D+00      
	      DO I=0,CTE_PH_NU
	         PHIOS = I*Q		 		!Azimuth convention OS
	         Y=Y+U(I)*DCOS(IS*PHIOS)
	      ENDDO
	      E(IS) = Y*Q / INCTE_PI
	      
C--------- Determination de l'ecart maximal entre F exact et recombine : B1
C---------------------------------------------------------------------      
	      B1=0.D+00			!initialisation
	      DO I=0,CTE_PH_NU
	         PHIOS=Q*I
	         PHIRJ = INCTE_PI - PHIOS	!Adaptation convention Roujean
                CALL SOS_CALC_F_ROUJEAN(K0,K1,K2,C1,S1,C2,S2,PHIRJ,F)     

C-------------- Calcul de F recombine à partir de la fonction E(IS) pour PHI fixe. 
                 T1 = E(0)
		 DO IS2 = 1, IS		 
		     T1 = T1+2.*E(IS2)*DCOS(IS2*PHIOS)	!Azimuth convention OS		     
		 ENDDO !Fin boucle sur IS2
		 
		 DO IS2 = 0, IS		 
		     IF(IS2.EQ.0) THEN
	                T1=E(0)			
	             ELSE 
	                T1=T1+2.*E(IS2)*DCOS(IS2*PHIOS)	!Azimuth convention OS
	             ENDIF		     
		 ENDDO !Fin boucle sur IS2
       
                 B1=MAX( DABS((T1-F)/F) , B1)
		 		 
	      ENDDO !Fin boucle sur I (phi = 0 a PI)
	           	     

C--------- Test d'arret de la decomposition :
C---------    * Arret si l'ecart relatif maximal entre F exacte et recombinee est
C---------      inferieur au seuil CTE_SEUIL_SF_ROUJEAN
C---------      ou s'il augmente par rapport a l'ordre precedent.
C---------    * Poursuite de la boucle sur IS sinon.
C-------------------------------------------------------------------------------- 		
              IF (B1.GT.CTE_SEUIL_SF_ROUJEAN) GOTO 14
              IL=IS	      
              GOTO 17	!Sortie de boucle sur IS
	      
 14           IF(B1.LT.B1_PREC) GOTO 141
              IL=IS-1
	      GOTO 17	!Sortie de boucle sur IS
	      
 141	      B1_PREC=B1	        
 13        CONTINUE	!Fin de boucle sur IS     
	   IL=OS_NB


	     



C------------------------------------------------------------------------
C* Stockage de la decomposition en series de Fourier de la matrice de 
C* reflexion de Roujean pour les directions Teta1 et Teta2
C--------------------------------------------------------------------------

C* Conversion de type pour chaque ordre IS : double précision -> réel
C---------------------------------------------------------------------
17          DO IS=0,OS_NB
                M11(IS)=E(IS)  		               
            ENDDO

C* Enregistrement pour le couple (Teta1, Teta2) de la decomposition
C* en series de Fourier de BRDF (convertie en luminance normalisée
C* et adaptee au code des OS).
C-----------------------------------------------------------------
            WRITE(4,err=992) M11



C* Fichier Trace : Calculs de paramètres d'analyse de qualite
C*                 de la decomposition en series de Fourier
C*                 et enregistrements.
C-----------------------------------------------------------------
           IF (TRACE) THEN

	      ECARTMAX_C1C2=0.
	      	      
	      Q=INCTE_PI/CTE_PH_NU
	      DO I=0,CTE_PH_NU
	         PHIOS=I*Q
		  PHIRJ = INCTE_PI - PHIOS		!Adaptation convention Roujean
	         CALL SOS_CALC_F_ROUJEAN(K0,K1,K2,C1,S1,C2,S2,PHIRJ,F)
		 
                FRECOMB = E(0)
	         DO IS=1,IL
	            FRECOMB = FRECOMB + 2.*DCOS(IS*PHIOS)*E(IS)
	         ENDDO
            
	         IF(F.NE.0.D+00) THEN
	            ECART = DABS((FRECOMB-F)/F)   

C*       Analyse pour tous les azimuts du couple courant (Teta1, Teta2)	
		     IF(ECART.GT.ECARTMAX_C1C2) THEN
		        ECARTMAX_C1C2 = ECART
		        PHIECARTMAX_C1C2 = PHIOS*180./INCTE_PI	!Azimuth convention OS
	               F_ECARTMAX_C1C2 = F
	               FRECOMB_ECARTMAX_C1C2 = FRECOMB
		     ENDIF 
		    
C*       Analyse pour tous les azimuts et tous les couples (Teta1, Teta2)		    
	            IF(ECART.GT.ECARTMAX) THEN
	               ECARTMAX = ECART
	               PHIECARTMAX = PHIOS*180./INCTE_PI	!Azimuth convention OS
	               TVECARTMAX = ACOS(C2)*180./INCTE_PI
	               TSECARTMAX = ACOS(C1)*180./INCTE_PI
	               F_ECARTMAX = F
	               FRECOMB_ECARTMAX = FRECOMB
	            ENDIF   
	    
	            SOMECART=SOMECART + ECART
	            COMPT = COMPT + 1
		 ENDIF
	      ENDDO !Fin boucle sur I (phi)
         	
	      	
	      IF(F_ECARTMAX_C1C2.NE.0.D+00) THEN	 
	         WRITE(99,2000,err=921) IL,PHIECARTMAX_C1C2,
     &		 			    F_ECARTMAX_C1C2,
     &		                         FRECOMB_ECARTMAX_C1C2,
     &					    ECARTMAX_C1C2*100.
	      ELSE
	      	 WRITE(99,2500,err=921) IL,PHIECARTMAX_C1C2,
     &		 			   F_ECARTMAX_C1C2,
     &		                        FRECOMB_ECARTMAX_C1C2
	      ENDIF	 
	      WRITE(99,*,err=921)   	    
                 
	      
           END IF !Fin Trace
      
      	
		                            
 4       CONTINUE	!Fin de boucle sur Teta2 (I2)
 
 3    CONTINUE	!Fin de boucle sur Teta1 (I1)
 
 


C* Fichier Trace : Enregistrements des paramètres d'analyse de qualite
C*                 de la decomposition en series de Fourier.
C----------------------------------------------------------------- 
      IF (TRACE) THEN
         WRITE(99,*,err=921)
	  WRITE(99,*,err=921)
         WRITE(99,*,err=921)  
     &	  "Analysis of the Fourier expansion of the ROUJEAN's BRDF: "
	
	 IF(COMPT.GT.1) THEN
	    WRITE(99,*,err=921)
	    WRITE(99,*,err=921)
	    WRITE(99,*,err=921)  
     &	    'Mean relative difference from the recombination'//
     &	    ' (% excluding null values):', SOMECART*100/COMPT
	    WRITE(99,*,err=921)
	    WRITE(99,*,err=921)  
     &	    'Max relative difference from the recombination (%):',
     &	    ECARTMAX*100
           WRITE(99,*,err=921)   
     &	    'Related delta phi angle (SOS convention) :',PHIECARTMAX
           WRITE(99,*,err=921)    
     &	    'Related viewing angle :',TVECARTMAX
	    WRITE(99,*,err=921)    
     &	    'Related solar angle :',TSECARTMAX
	    WRITE(99,*,err=921)    
     &	    'Related value of BRDF :',F_ECARTMAX
	    WRITE(99,*,err=921)   
     &	    'Related value of the recombined BRDF :',FRECOMB_ECARTMAX 
        ELSE
	    WRITE(99,*,err=921)
	    WRITE(99,*,err=921) 'WARNING : BRDF forced to zero !'
	    WRITE(99,*,err=921) 'Unsuitable simulation parameters'
	 ENDIF 
	 
	 WRITE(99,*,err=921)
	 WRITE(99,*,err=921)
	 
      END IF 
      
      

       
      
      
C* Fin nominale 
C-------------------
      IER=0
      GOTO 9999
      
C* Cas d'erreur : retour du status -1 au programme appelant
C------------------------------------------------------
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Surface'
      IER=-1
      GOTO 9999
  991 WRITE(6,*) 'SOS_FSF_ROUJEAN : Erreur a l''ouverture d''un fichier'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'SOS_FSF_ROUJEAN',
     &		 'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      GOTO 9999
  993 WRITE(6,*) 'SOS_FSF_ROUJEAN : ATTENTION : PRESENCE DE BRDF < 0 !',
     &		 'Paramètres de simulation inadaptes.'
      IER=-1 
      GOTO 9999
      
 9999 CONTINUE


C* Fermeture fichier :'RES_MAT_REFLEX'
C-------------------------------------
      CLOSE(4)
 
 
C* Format
C--------

1000  FORMAT(6hTeta1:,1X,F5.2,2X,6hTeta2:,1X,F5.2,2X,
     &15h; (PHIos,Fmin):,1X,F6.2,2X,E11.4,2X,
     &15h; (PHIos,Fmax):,1X,F6.2,2X,E11.4)    
 
2000  FORMAT(5X,12hLimit order:,I3,3h ; ,14h Error max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4,2X,
     &12h, Ecart (%):,1X,F7.3) 
     
2500  FORMAT(5X,12hLimit order:,I3,3h ; ,14h Error max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4)                               
      RETURN
      END         !FIN DE LA PROCEDURE SOS_FSF_ROUJEAN





C==============================================================================
C PROCEDURE: SOS_CALC_F_ROUJEAN
C ==========
C      Cette procédure calcule la BRDF du modèle de Roujean:
C	 F(Teta1,Teta2,phi) = K0 + K1*F1(Teta1,Teta2,phi) + K2*F2(Teta1,Teta2,phi)
C
C	 avec Teta1 : angle zenithal d'incidence : tetas  (0 < Teta1 < PI/2),
C             Teta2 : angle zenithal de reflexion : tetav(0 < Teta2 < PI/2),
C             phi   : angle azimutal relatif entre les directions
C                     d'incidence et de reflexion selon la convention 
C                     du modèle de Roujean  (0 <= phi <= PI),
C
C	      K0, K1, K2 : coefficients du modèle de Roujean,
C			
C            F1 = [ (PI-phi)*cos(phi)+sin(phi) ]*tgs*tgv / 2*PI
C                 -[ tgs + tgv + SQRT{tgs**2+tgv**2-2*tgs*tgv*cos(phi)} ] / PI
C
C		  en posant tgs = tan(tetas) et tgv = tan(tetav),
C
C            F2 = [ (PI/2-zeta)*cos(zeta)+sin(zeta) ] * (4/3*PI) 
C                 / [ cos(tetas)+cos(tetav) ]   - 1/3
C
C                 avec zeta l'angle de phase defini par :
C                 cos(zeta) = cos(tetas)*cos(tetav)+sin(tetas)*sin(tetav)*cos(phi). 
C
C
C      La BRDF du modèle de Roujean est convertie de reflectance en une 
C      grandeur utile pour les OS : luminance normalisée * cos(Teta2),
C      c'est a dire : BRDF * cos(Teta1) * cos(Teta2).			
C
C Description des paramètres
C --------------------------
C       K0	(double)   (E) Coefficient du modèle de Roujean.
C       K1	(double)   (E) Coefficient du modèle de Roujean.
C       K2	(double)   (E) Coefficient du modèle de Roujean.
C	 C1	(double)   (E) Cos(Teta1).
C	 S1	(double)   (E) Sin(Teta1).
C	 C2	(double)   (E) Cos(Teta2).
C	 S2	(double)   (E) Sin(Teta2).
C       PHI     (double)   (E) Azimut relatif pour la convention
C                              du modèle de Roujean (radians).
C 	 F	(double)   (S) Valeur de F(Teta1,Teta2,phi).
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
C==============================================================================
      SUBROUTINE SOS_CALC_F_ROUJEAN(K0,K1,K2,C1,S1,C2,S2,PHI,F)

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------------------------------------
      DOUBLE PRECISION K0	! (E) Paramètre du modèle de Roujean.
      DOUBLE PRECISION K1	! (E) Paramètre du modèle de Roujean.
      DOUBLE PRECISION K2	! (E) Paramètre du modèle de Roujean.

      DOUBLE PRECISION C1	! (E)  Cos(Teta1).
      DOUBLE PRECISION S1	! (E)  Sin(Teta1).
      DOUBLE PRECISION C2	! (E)  Cos(Teta2).
      DOUBLE PRECISION S2	! (E)  Sin(Teta2).
      DOUBLE PRECISION PHI	! (E) Azimut relatif pour la convention 
      				!     Roujean (radians).
      
      DOUBLE PRECISION xPHI	! Valeur de PHI réellement utilisee :
      				! => valeur entre 0 et PI.
      
      DOUBLE PRECISION xC1	! Valeur de C1 réellement utilisee.
      DOUBLE PRECISION xS1	! Valeur de S1 réellement utilisee.
      DOUBLE PRECISION xC2	! Valeur de C2 réellement utilisee.
      DOUBLE PRECISION xS2	! Valeur de S2 réellement utilisee.
      
      DOUBLE PRECISION F	! (S) Valeur de la BRDF de Roujean 
      				! pour la geometrie d'observation
				! (convertie en luminance normalisée * cos(teta2).	
      
      DOUBLE PRECISION F1	! Fonction de la BRDF de ROUJEAN.
      DOUBLE PRECISION F2	! Fonction de la BRDF de ROUJEAN.
      
      DOUBLE PRECISION COSPHI	! Cosinus de l'azimut relatif.
      DOUBLE PRECISION TANTS	! Tangente de l'angle d'incidence.
      DOUBLE PRECISION TANTV	! Tangente de l'angle de reflexion.
            
      DOUBLE PRECISION ZETA	! Angle de phase.
      DOUBLE PRECISION COSZETA	! Cosinus de l'angle de phase entre les
      				       ! directions d'incidence et de reflexion.

      DOUBLE PRECISION PI	! Valeur de PI
      

C* Traitements
C---------------------------------    
      PI=4.D+00*DATAN(1.D+00)

C*   Controle des angles                 
C---------------------------------------------------------      

C*   Conversion de PHI entre 0 et PI
      xPHI = PHI
      IF(xPHI.LT.0.D+00) xPHI = - xPHI
      IF(xPHI.GT.PI) xPHI = 2.*PI - xPHI

      
C*   Controle des angles zenitaux de visee et d'eclairement
      xC1=C1
      xS1=S1
      xC2=C2
      xS2=S2  
      
      IF (DACOS(C1)*180./PI.GT.CTE_TETAS_LIM_ROUJEAN) THEN
          xC1 = DCOS(CTE_TETAS_LIM_ROUJEAN*PI/180.D+00)
 	   xS1 = DSIN(CTE_TETAS_LIM_ROUJEAN*PI/180.D+00)
      ENDIF
      IF (DACOS(C2)*180./PI.GT.CTE_TETAV_LIM_ROUJEAN) THEN
          xC2 = DCOS(CTE_TETAV_LIM_ROUJEAN*PI/180.D+00)
 	   xS2 = DSIN(CTE_TETAV_LIM_ROUJEAN*PI/180.D+00)
      ENDIF	
	
	
C*   Calcul du cos(phi) et des tangentes de Tetas et Tetav                
C---------------------------------------------------------     	
      COSPHI = DCOS(xPHI)
      TANTS =  xS1 / xC1
      TANTV =  xS2 / xC2

C*   Calcul de F1                 
C---------------------------------------------------------      
      F1 = 0.5 * ( (PI-xPHI)*COSPHI + DSIN(xPHI) ) * TANTS * TANTV 
      F1 = F1 - TANTS - TANTV
      F1 = F1 - DSQRT(TANTS*TANTS + TANTV*TANTV -2.*TANTV*TANTS*COSPHI)
      F1 = F1 / PI
      
      
C*   Calcul de F2                 
C---------------------------------------------------------            
      COSZETA = xC1*xC2+xS1*xS2*COSPHI  

      !Gestion du cas d'une petite erreur d'arrondie numérique donnant un dépassement de 1 ou -1
      IF ( DABS(DABS(COSZETA)-1.D+00).LE.CTE_SEUIL_NUM ) THEN
         IF ( (COSZETA.GE.(1.-CTE_SEUIL_NUM))
     &       .AND.(COSZETA.LE.(1.+CTE_SEUIL_NUM)) ) THEN
            COSZETA = 1.D+00
         ELSE
            COSZETA = -1.D+00
         ENDIF
      ENDIF

      ZETA = DACOS(COSZETA)

      F2 = 4.*( (PI/2. - ZETA)*COSZETA + DSIN(ZETA) ) 
     & / ( 3.*PI*(xC1+xC2))
      F2 = F2 - (1.D+00/3.D+00)



C*   Calcul de F               
C---------------------------------------------------------            
      F = K0 + K1*F1 + K2*F2

C*   Conversion reflectance -> luminance normalisée (*C1)
C*   et adaptation aux OS -> *C2.  
C*
C*   NB : 
C*    - L'adaptation au code des OS par multiplication de 
C*      cos(tetav) doit utiliser la valeur correspondante au 
C*      fichier des angles de Gauss car c'est cette valeur
C*      qui sera prise dans les calculs du code des OS
C*      (on prend C2 et non pas xC2).
C*    - De même, la conversion de réflectance en luminance normalisée  
C*      est faite pour la valeur de cos(TetaIncident) vraie, et non pas
C*      seuillée par CTE_TETAS_LIM_ROUJEAN, car les calculs de flux
C*      dans SOS_OS nécessitent d'utiliser l'angle d'incidence vrai : 
C*      on prend C1 et non pas xC1.
                   
C--------------------------------------------------------- 	
      F = F*C2*C1	     

      RETURN
      END         !FIN DE LA PROCEDURE SOS_CALC_F_ROUJEAN










C==============================================================================
C PROCEDURE:  SOS_MISE_FORMAT_RJ
C ==========
C      Cette procédure reordonne les elements de la matrice Mij en une matrice
C      Pij utilisable par le code des OS.
C
C      Les elements M11, enregistres dans le fichier "RES_MAT_REFLEX", donnent
C      pour chaque couple d'angles (Teta1, Teta2) la decomposition en series 
C      de Fourier associee. Les autres elements (de polarisation) ne sont pas 
C      precalcules car ils sont poses nuls.
C
C      Les elements Pij, enregistres dans le fichier de sortie, donnent
C      pour chaque ordre de la decomposition en series de Fourier
C      la reflexion associee a l'ensemble des couples (Teta1, Teta2) .
C
C Description des paramètres
C --------------------------
C     FIC_RES_MAT_REFLEX (CHARACTER)  : (E) Nom du fichier temporaire RES_MAT_REFLEX
C                                           avec arborescence complète.

C      FIC  (CHARACTER*CTE_LENFIC2)  (E)  Nom complet du fichier resultat.
C      NBMU (I4)  : (E)  Nombre d'angles (positifs) effectivement utiles.
C      OS_NB (I4) : (E)  Ordre des développements des fonctions de phase 
C		         en Polynomes de Legendre a utiliser.
C      IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon
C
C
C Autre resultat fourni par la procédure
C -------------------------------------
C      Fichier resultat : FIC, fichier de BRDF du modèle de Roujean.
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier temporaire RES_MAT_REFLEX
C       Fichier binaire non formatte, cree par SOS_FSF_ROUJEAN
C       Fichier cree dans le repertoire de lancement du programme
C       Enregistrement pour chaque couple (Teta1, Teta2) de
C       la decomposition en series de Fourier de : M11 (n valeurs réelles)
C
C   -->Contenu du fichier BRDF ROUJEAN
C       Fichier binaire non formatte, cree par SOS_MISE_FORMAT_RJ
C       Il contient la matrice de reflexion de surface 
C       c'est a dire les elements Pij de la matrice de reflexion 
C              P11(CTE_OS_NBMU,CTE_OS_NBMU) (réel),
C              P12(CTE_OS_NBMU,CTE_OS_NBMU) (réel),
C              P13(CTE_OS_NBMU,CTE_OS_NBMU) (réel),
C              P21(CTE_OS_NBMU,CTE_OS_NBMU) (réel),
C              P22(CTE_OS_NBMU,CTE_OS_NBMU) (réel),
C              P23(CTE_OS_NBMU,CTE_OS_NBMU) (réel),
C              P31(CTE_OS_NBMU,CTE_OS_NBMU) (réel),
C              P32(CTE_OS_NBMU,CTE_OS_NBMU) (réel),
C              P33(CTE_OS_NBMU,CTE_OS_NBMU) (réel)
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C      Le fichier RES_MAT_REFLEX n'existe pas
C      Le fichier de BRDF Roujean existe
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C
C=============================================================================
      SUBROUTINE SOS_MISE_FORMAT_RJ(FIC_RES_MAT_REFLEX,FIC,
     &                              NBMU,OS_NB,IER)

      IMPLICIT NONE

C* Definition des variables                    
C*--------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FIC	! Nom du fichier resultat.

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX 
                     ! Nom du fichier RES_MAT_REFLEX
                     ! avec arborescence complète

      	! Elément M11 de la matrice de reflexion en fonction de l'ordre IS
	! de la serie de Fourier, pour une incidence sous Teta1 et une 
	! reflexion sous Teta2.
      REAL M11(0:CTE_OS_NB_MAX)		
      

      	! Eléments Pij de la matrice de reflexion en fonction de la direction
	! d'incidence et de reflexion, pour un ordre IS du développement
	! en series de Fourier
      REAL P11(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P12(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P13(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)	
      REAL P21(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P22(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P23(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P31(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P32(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P33(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)

      INTEGER*4 IS	! Ordre de la decomposition en series de Fourier.
      INTEGER*4 I	! Indice de l'angle d'incidence Teta1.
      INTEGER*4 J	! Indice de l'angle de reflexion Teta2.
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon  
      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB ! Ordre des développements des fonctions de phase 
			 ! en Polynomes de Legendre.
       
C* Ouverture fichier
C------------------------------------------------------------
      OPEN(UNIT=1,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=2,FILE=FIC,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=991)

C* Boucle sur l'ordre de la decomposition en series de Fourier
C------------------------------------------------------------
      DO 1 IS=0,OS_NB

C* Boucle sur les directions d'incidence Teta1
C---------------------------------------------
         DO 2 I=1,NBMU

C* Boucle sur les directions de reflexion Teta2
C-----------------------------------------------
            DO 3 J=1,NBMU
               READ(1,err=992) M11
               P11(I,J)=M11(IS)
               P12(I,J)=0. 
               P13(I,J)=0.  
               P21(I,J)=0.
               P22(I,J)=0. 
               P23(I,J)=0.
               P31(I,J)=0.
               P32(J,I)=0.
               P33(J,I)=0.
 3          CONTINUE
 2       CONTINUE

C*       Enregistrement des tableaux de la matrice de reflexion 
C*       a l'ordre IS de la decomposition en series de Fourier.
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
      GOTO 9999
      
C* Cas d'erreur : retour du status -1 au programme appelant
C------------------------------------------------------
  991 WRITE(6,*) 'SOS_MISE_FORMAT_RJ : '//
     &           'Erreur a l''ouverture d''un fichier'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'SOS_MISE_FORMAT_RJ : '//
     &           'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      GOTO 9999

 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(1)
      CLOSE(2)
      
      RETURN
      END         !FIN DE LA PROCEDURE SOS_MISE_FORMAT_RJ





