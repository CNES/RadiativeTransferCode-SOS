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
C* FICHIER: SOS_SURFACE_BPDF.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Génère les fichiers de BPDF des modèles de Rondeaux-Herman, Bréon et Maignan.
C*
C* AUTEUR: 
C*   CS GROUP France, B.Lafrance.
C* DATE: 12/08/2020
C*
C* MOD:VERSION:1.0: Création du programme.
C*    - Récupération des routines existantes : SOS_GSF_RONDEAUX_BREON, SOS_F21SF_NADAL, SOS_CALC_F21_NADAL_SUR_FRESNEL.
C*    - Développement de routines nouvelles: SOS_GSF_MAIGNAN et SOS_CALCG_MAIGNAN.
C*
C* MOD:VERSION:1.1: 31/08/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Conversion des traces en Anglais.
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisées :
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_OS_NBMU_MAX : Nombre maximal d'angles positifs.
C    INCTE_PI : Valeur de PI.
C    CTE_OS_NB_MAX :  Valeur maximale pour l'ordre limite du développement en 
C                     polynomes de Legendre.
C    CTE_OS_NS_MAX :  Valeur maximale pour l'ordre limite du développement en polynomes 
C                     de Legendre des éléments de la matrice de Fresnel 
C                     et en séries de Fourier des luminances.
C    CTE_OS_NM_MAX :  Valeur maximale pour l'ordre limite  du développement 
C                     en séries de Fourier de la fonction G.
C    CTE_PH_TEST : Facteur pour la comparaison de GMIN et GMAX
C    CTE_PH_NU : Nombre d'angle azimutal (2**NQ)
C    CTE_PH_NQ : Puissance de 2 de CTE_PH_NU (NU=2**NQ)
C
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROCEDURE: SOS_SURFACE_BPDF
C ==========      
C      Cette procédure génère les fichiers de réflexion de surface des modèles de BPDF

C          Modèle de Rondeaux :
C           => Rp(Teta1, Teta2, phi) =  G(Teta1, Teta2, phi) * Fp(Teta1, Teta2, phi) 
C                avec G(Teta1, Teta2, phi) = 1 / [ 4 * (cosTeta1+cosTeta2) ]
C
C          Modèle de Bréon :
C           => Rp(Teta1, Teta2, phi) =  G(Teta1, Teta2, phi) * Fp(Teta1, Teta2, phi) 
C                avec G(Teta1, Teta2, phi) = 1 / [ 4*cosTeta1*cosTeta2 ]
C
C          Modèle de Maignan :
C           => Rp(Teta1, Teta2, phi) =  G(Teta1, Teta2, phi) * Fp(Teta1, Teta2, phi) 
C                avec G(Teta1, Teta2, phi) = C*exp[-tan(TetaN)] / [ 4 * (cosTeta1+cosTeta2) ]
C
C	     avec Teta1 : angle zenithal d'incidence  (0 < Teta1 < PI/2),
C                Teta2 : angle zenithal de réflexion (0 < Teta2 < PI/2),
C                phi   : angle azimutal relatif entre les directions
C                        d'incidence et de réflexion,
C                Fp(Teta1, Teta2, phi) : matrice de réflexion de Fresnel,
C                TetaN : angle d'incidence par rapport à la normale de la facette 
C                        qui pour une incidence TetaN donne une réflexion sous Teta2 pour l'angle azimutal phi,
C                C     : coefficient constant du modèle de Maignan (intègre le terme exp(-NDVI) du modèle).
C
C
C Données en entrée du programme 
C ------------------------------
C    -->Les paramètres/arguments du programme
C       LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C	 RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (E) Cosinus des angles.
C	 CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (E) Poids des angles.
C       IND   (DOUBLE)  : (E) Indice de réfraction de l'eau / air.
C	 ISURF  (I2) : (E) Type de surface ISURF 
C			       (format de decodage en entier I2)
C			4 :surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
C			5 :surface lambertienne + BRDF de Roujean + BPDF de Breon.
C			6 :surface lambertienne + BRDF de Roujean + BPDF de Nadal.	
C			7 :surface lambertienne + BRDF de Roujean + BPDF de Maignan.	    
C 	 ALPHA   (DOUBLE)  (E) : Paramètre alpha du modèle de Nadal ==> Associé à ISURF = 6         	
C 	 BETA    (DOUBLE)  (E) : Paramètre beta du modèle de Nadal  ==> Associé à ISURF = 6
C       COEF_C  (DOUBLE)  (E) : Paramètre C*exp(-NDVI) du modèle de BPDF de Maignan
C                               ==> Associé à : ISURF = 7   		        	  
C       OS_NB (I4) : (E) Ordre des développements des fonctions de phase 
C		           en Polynomes de Legendre à utiliser.
C       OS_NS (I4) : (E) Ordre maximal à utiliser pour le développement en
C                        polynomes de Legendre des éléments de la matrice de Fresnel 
C                        et en séries de Fourier des calculs de luminance.
C       OS_NM (I4) : (E) Ordre maximal du développement en séries de Fourier  
C                        de la fonction G des calculs de matrice de réflexion.
C       FIC_RES_FRESNEL (CHARACTER*CTE_LENFIC2)  : 
C                 (E) Nom du fichier temporaire RES_FRESNEL avec arborescence complète.
C       FIC_RES_GSF (CHARACTER*CTE_LENFIC2)  :  
C                 (E) Nom du fichier temporaire RES_GSF avec arborescence complète.
C       FIC_RES_FRESNEL (CHARACTER*CTE_LENFIC2)  :  
C                 (E) Nom du fichier temporaire RES_FRESNEL avec arborescence complète.
C       FIC_RES_MAT_REFLEX (CHARACTER*CTE_LENFIC2)  :  
C                 (E) Nom du fichier temporaire RES_MAT_REFLEX avec arborescence complète.
C       FICBPDF  (CHARACTER*CTE_LENFIC2) : (S) nom complet du fichier résultat.
C       TRACE  (LOGICAL)   ! = vrai,si écriture dans le fichier trace.
C       IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon
C
C
C Résultat fourni par le programme
C --------------------------------
C   Fichier BRDF / BPDF
C   Fichier de trace (optionnel), 
C               si un fichier de même nom existe déjà, il est écrasé
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
C   -->Contenu du fichier Résultat
C       Fichier binaire non formatté, créé par SOS_MISE_FORMAT
C       Il contient la matrice de réflexion de surface ,
C       c'est à dire les éléments Pij de la matrice de réflexion 
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
C   -->Contenu du fichier temporaire RES_GSF
C       Fichier binaire non formatte, créé par SOS_GSF_RONDEAUX_BREON, SOS_F21SF_NADAL ou SOS_GSF_MAIGNAN
C       Fichier créé dans le répertoire de lancement du programme
C       Si un fichier de même nom existe déjà, il est écrasé
C       Chaque enregistrement contient :
C           -numéro d'angle Teta1 (valeur entière)
C           -numéro d'angle Teta2 (valeur entière)
C           -ordre maximale utile de la décomposition : IL  (valeur entière)
C           -Fonction G décomposée pour (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 valeurs double précision)
C
C   -->Contenu du fichier temporaire RES_FRESNEL
C      Fichier Ascii, format d'un enregistrement 4E15.8, 
C      créé par SOS_MAT_FRESNEL
C      Fichier créé dans le répertoire de lancement du programme
C      Si un fichier de même nom existe déjà, il est écrasé
C      Nombre d'enregistrements : de k=0 a OS_NS
C      Chaque enregistrement contient les paramètres ALPHA(k), BETA(k), 
C      GAMMA(k) et ZETA(k) associés à la matrice de Fresnel
C
C   -->Contenu du fichier temporaire RES_MAT_REFLEX
C       Fichier binaire non formatte, créé par SOS_MAT_REFLEXION
C       Fichier créé dans le répertoire de lancement du programme
C       Si un fichier de même nom existe déjà, il est écrasé
C       Enregistrement pour chaque couple (Teta1, Teta2) de
C       la décomposition en séries de Fourier de la matrice de Fresnel:
C       M11,M12,M13,M21,M22,M23,M31,M32,M33 
C       (n valeurs réelles)
C
C    Les 3 fichiers temporaires sont écrasés s'ils existent, puis 
C    détruits en fin de programme.
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
C      Erreur si le fichier Glitter existe déjà
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, arrêt du programme et 
C     retour du status 1 au shell
C        
C
C  Sous programmes utilises:
C --------------------------
C  Cette procédure fait appel aux routines:
C      - SOS_GSF_RONDEAUX_BREON
C      - SOS_F21SF_NADAL et SOS_CALC_F21_NADAL_SUR_FRESNEL 
C      - SOS_GSF_MAIGNAN et SOS_CALCG_MAIGNAN
C      - SOS_MAT_FRESNEL 
C      - SOS_MAT_REFLEXION 
C      - SOS_MISE_FORMAT
C
C=============================================================================

      SUBROUTINE SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                            ISURF,ALPHA,BETA,COEF_C,
     &                            OS_NB,OS_NS,OS_NM,
     &                            FIC_RES_GSF,FIC_RES_FRESNEL,
     &                            FIC_RES_MAT_REFLEX,
     &                            FICBPDF,TRACE,IER)

      IMPLICIT NONE


C* Definition des variables   
C*--------------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICBPDF	! Nom complet du fichier BPDF/BRDF résultat.

      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Nom du fichier RES_GSF
                                        ! avec arborescence complète

      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Nom du fichier RES_FRESNEL
                                            ! avec arborescence complète

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Nom du fichier RES_MAT_REFLEX
                                               ! avec arborescence complète

      DOUBLE PRECISION IND	! Indice de réfraction de l'eau / air

      DOUBLE PRECISION ALPHA	!Paramètre alpha du modèle de Nadal.
      DOUBLE PRECISION BETA	       !Paramètre beta du modèle de Nadal.
      DOUBLE PRECISION COEF_C      !Paramètre C*exp(-NDVI) du modèle de Maignan.
	
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Tableau des mu de Gauss
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Poids de Gauss.
      
      LOGICAL TRACE     ! = vrai,si écriture dans le fichier trace

      INTEGER*2 ISURF   ! Type de surface     

      INTEGER*4 OS_NB   ! Ordre des développements des fonctions de phase 
			   ! en Polynomes de Legendre.
      INTEGER*4 OS_NS   ! Ordre maximal du développement en polynomes de Legendre   
			   ! des éléments de la matrice de Fresnel.
      INTEGER*4 OS_NM   ! Ordre maximal du développement en séries de Fourier  
			   ! de la fonction G des calculs de matrice de réflexion.
			
      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations
      INTEGER*4 I	    ! Indice de tableau.
      INTEGER*4 IFIN     ! nb caracteres significatifs d'une chaine
      INTEGER*4 IER      ! code d'erreur =0 si pas d'erreur, =-1 sinon




C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0

     
      
C* Ouverture du fichier Trace 
C------------------------------------------------------
      IF (TRACE) THEN
         WRITE(99,*,err=921)	 
         WRITE(99,*,err=921) 'Calculation of the MAIGNAN BPDF file'  
	  WRITE(99,*,err=921) FICBPDF  
         WRITE(99,*,err=921) ' '
         WRITE(99,*,err=921)'Simulation parameters ' 
         WRITE(99,*,err=921)
     &        ' Refractive index surface /air ',IND
         WRITE(99,*,err=921)' OS_NB ',OS_NB
	  WRITE(99,*,err=921)' OS_NS ',OS_NS
	  WRITE(99,*,err=921)' OS_NM ',OS_NM
	  WRITE(99,*,err=921) ' ' 	 	 
      END IF
        

C* Decomposition en séries de Fourier de la fonction G(Teta1,Teta2,phi).
C---------------------------------------------------------------------
      IF (ISURF.EQ.4) THEN	!Modèle BPDf de Rondeaux
          CALL SOS_GSF_RONDEAUX_BREON(LUM_NBMU,RMU,1,OS_NM,
     &                                FIC_RES_GSF,TRACE,IER)
      ENDIF

      IF (ISURF.EQ.5) THEN	!Modèle BPDf de Bréon
          CALL SOS_GSF_RONDEAUX_BREON(LUM_NBMU,RMU,2,OS_NM,
     &                                FIC_RES_GSF,TRACE,IER)
      ENDIF

C			4 :surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
C			5 :surface lambertienne + BRDF de Roujean + BPDF de Breon.



      IF (ISURF.EQ.6) THEN
          CALL SOS_F21SF_NADAL(LUM_NBMU,RMU,IND,ALPHA,BETA,OS_NB,
     &                         FIC_RES_GSF,TRACE,IER)
      ENDIF

      IF (ISURF.EQ.7) THEN
          CALL SOS_GSF_MAIGNAN(LUM_NBMU,RMU,COEF_C,OS_NM,
     &                         FIC_RES_GSF,IER)
      ENDIF
     
      IF (IER.NE.0) GOTO 995   ! cas d'erreur
      
      
C* Développement en fonctions de Legendre de la matrice de Fresnel.
C------------------------------------------------------------------
      CALL SOS_MAT_FRESNEL(LUM_NBMU,RMU,CHR,IND,OS_NS,
     &                     FIC_RES_FRESNEL,TRACE,IER)
      IF (IER.NE.0) GOTO 996   ! cas d'erreur
      
C* Calcul de la matrice de réflexion developpee en séries de Fourier.
C------------------------------------------------------------------
      CALL SOS_MAT_REFLEXION(1.D+00,LUM_NBMU,RMU,OS_NB,OS_NS,OS_NM,
     &                       FIC_RES_FRESNEL,FIC_RES_GSF,
     &                       FIC_RES_MAT_REFLEX,IER)
      IF (IER.NE.0) GOTO 997   ! cas d'erreur
      
C* Rearrangement des résultats pour etre compatible avec le code des OS.
C-----------------------------------------------------------------------
      CALL SOS_MISE_FORMAT(FIC_RES_MAT_REFLEX, FICBPDF,
     &                     LUM_NBMU,OS_NB,IER)
      IF (IER.NE.0) GOTO 998   ! cas d'erreur


C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  912 WRITE(6,*) 'Erreur à l''ouverture du fichier angles de Gauss'
      IER=-1
      GOTO 9999
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Glitter'
      IER=-1
      GOTO 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1
      GOTO 9999      
  995 IF(ISURF.EQ.6) THEN
         WRITE(6,*) 'Erreur dans la routine SOS_F21SF_NADAL'
      ENDIF
      IF(ISURF.EQ.7) THEN
         WRITE(6,*) 'Erreur dans la routine SOS_GSF_MAIGNAN'
      ENDIF
      IER=-1
      GOTO 9999
  996 WRITE(6,*) 'Erreur dans la routine SOS_MAT_FRESNEL'
      IER=-1
      GOTO 9999
  997 WRITE(6,*) 'Erreur dans la routine SOS_MAT_REFLEXION'
      IER=-1
      GOTO 9999
  998 WRITE(6,*) 'Erreur dans la routine SOS_MISE_FORMAT'
      IER=-1
      GOTO 9999         
           
C* Format
C---------   
 2    FORMAT(2(D21.14))    
              
 9999 CONTINUE

C* Destruction des fichiers temporaires
C---------------------------------------   
      OPEN(UNIT=3,FILE=FIC_RES_GSF,FORM='UNFORMATTED') 
      CLOSE(3,STATUS='DELETE')     
      OPEN(UNIT=3,FILE=FIC_RES_FRESNEL) 
      CLOSE(3,STATUS='DELETE')       
      OPEN(UNIT=3,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED')
      CLOSE(3,STATUS='DELETE') 
                                                      
      RETURN
      END    !FIN DU PROGRAMME SOS_SURFACE_BPDF




C==============================================================================
C PROCEDURE: SOS_GSF_RONDEAUX_BREON
C ==========
C      Cette procédure realise la decomposition en series de Fourier
C      (sur phi) de la fonction G(Teta1,Teta2,phi) :
C
C         Modèle de Rondeaux :G(Teta1, Teta2, phi) = 1. / [ 4 * (cosTeta1+cosTeta2) ]
C
C         Modèle de Breon :   G(Teta1, Teta2, phi) = 1. / [ 4*cosTeta1*cosTeta2 ]
C
C	    avec  Teta1 : angle zenithal d'incidence  (0 < Teta1 < PI/2),
C                Teta2 : angle zenithal de reflexion (0 < Teta2 < PI/2),
C                phi   : angle azimutal relatif entre les directions
C                        d'incidence et de reflexion,
C			
C
C      Enregistre les resultats dans le fichier "RES_GSF" sous la forme:
C	   Numéro d'angle Teta1 : I1,
C         Numéro d'angle Teta2 : I2,
C         Ordre maximal utile de la decomposition: IL
C            (IL = 0 pour ces fonctions independantes  de l'azimut),
C         Fonction G decomposee pour (Teta1,Teta2) : (E(J),J=0,IL)
C
C Description des paramètres
C --------------------------
C      LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles.
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)	 (E) : Tableau des mu de Gauss.
C      IRONDEAUX (integer*4) (E) : Indice de modèle 1 <-> Rondeaux, Breon sinon.
C      OS_NM (I4) : (E) Ordre maximal du développement en series de Fourier  
C                       de la fonction G des calculs de matrice de reflexion.
C      FIC_RES_GSF (CHARACTER*CTE_LENFIC2)  :  
C                 (E) Nom du fichier temporaire RES_GSF avec arborescence complète.
C      TRACE (logical) (E) : Option de Trace.
C      IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon.
C
C Autre resultat fourni par la procédure
C -------------------------------------
C      Fichier temporaire RES_GSF
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier temporaire RES_GSF
C       Fichier binaire non formatte, cree par SOS_GSF_RONDEAUX_BREON
C       Fichier cree dans le repertoire de lancement du programme
C       Si un fichier de meme nom existe deja, il est ecrase
C       Chaque enregistrement contient :
C           -numéro d'angle Teta1 (valeur entiere)
C           -numéro d'angle Teta2 (valeur entiere)
C           -ordre maximale utile de la decompositon : IL  (valeur entiere)
C           -Fonction G decomposee pour (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 valeurs double précision)
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
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C==============================================================================
      SUBROUTINE SOS_GSF_RONDEAUX_BREON(LUM_NBMU,RMU,IRONDEAUX,
     &                                  OS_NM,FIC_RES_GSF,TRACE,IER)
 	
      IMPLICIT NONE

C* Definition des variables    
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Tableau des mu de Gauss

      DOUBLE PRECISION E(0:CTE_OS_NM_MAX)	! Resultat de la decomposition en 
      				        ! series de Fourier de la fonction G.

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).

      LOGICAL TRACE     ! = vrai,si écriture dans le fichier trace.

      INTEGER*4 OS_NM   ! Ordre maximal du développement en series de Fourier  
			   ! de la fonction G des calculs de matrice de reflexion.      
      INTEGER*4 LUM_NBMU  ! Nombre d'angles effectifs des simulations
      INTEGER*4 IRONDEAUX ! 1 si simulation du modèle de Rondeaux.
      		            ! Sinon simulation du modèle de Breon.
      INTEGER*4 IL	! Ordre maximale utile pour la fonction E.
      INTEGER*4 I1	! Indice de numéro d'angle Teta1.
      INTEGER*4 I2	! Indice de numéro d'angle Teta2.
      INTEGER*4 IER  ! code d'erreur =0 si pas d'erreur, =-1 sinon  
      INTEGER*4 J	! Indice de tableau.
      
      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Nom du fichier RES_GSF
                                        ! avec arborescence complète


C* Ouverture du fichier resultat
C--------------------------------
      OPEN(UNIT=30,FILE=FIC_RES_GSF,FORM='UNFORMATTED',ERR=991)

C* Fichier Trace
C--------------------------------
      IF (TRACE) THEN
         IF(IRONDEAUX.EQ.1) THEN
	     WRITE(99,*,err=921) 'Values of the function G Rondeaux'     
	  ELSE
	     WRITE(99,*,err=921) 'Values of the function G Breon'       
         ENDIF
	    	 
      ENDIF  
	   
C* Boucle sur l'angle Teta1 > 0.
C---------------------------------
      DO  3 I1=1 ,LUM_NBMU

C* Cosinus et sinus de Teta1
C---------------------------------
         C1=RMU(I1)
         S1=DSQRT(1-C1*C1)

C*   Boucle sur l'angle teta2 > 0.
C---------------------------------
         DO 4 I2=1 ,I1

C*   Cosinus et sinus de Teta2
C---------------------------------
           C2=RMU(I2)
           S2=DSQRT(1-C2*C2)
          

C*   Fonction G a l'ordre IS = 0 pour le modèle de Rondeaux.
C*   => Independant de l'azimut : Terme unique IS = 0.
C---------------------------------
           IF(IRONDEAUX.EQ.1) THEN
	      E(0) = 1. / (1./C1 + 1./C2)
	   ELSE 
	      E(0) = 1.
	   ENDIF   
	   IL=0

C* Ecriture des resultats dans le fichier Trace
C--------------------------------
	   IF (TRACE) THEN
	       WRITE(99,1000,err=921) ACOS(C1)*180./INCTE_PI, 
     &	                              ACOS(C2)*180./INCTE_PI,
     &                               E(0)         	 
           ENDIF  
	      
                 

C*   Enregistrement des resultats
C---------------------------------
 17        WRITE(30,err=992)I1,I2,IL,(E(J),J=0,IL)

 4       CONTINUE	!Fin de boucle sur Teta2 (I2)
 
 3    CONTINUE	!Fin de boucle sur Teta1 (I1)
      
C* Fin nominale 
C-------------------
      IER=0
      GOTO 9999
      

C* Cas d'erreur : retour du status -1 au programme appelant
C------------------------------------------------------
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Surface'
      IER=-1
      GOTO 9999
  991 WRITE(6,*) 'SOS_GSF_RONDEAUX_BREON: ',
     &           'Erreur a l''ouverture d''un fichier'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'SOS_GSF_RONDEAUX_BREON: ',
     &		 'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      GOTO 9999

 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(30)

C* Format
C--------
1000  FORMAT(6hTeta1:,1X,F5.2,2X,6hTeta2:,1X,F5.2,2X,
     &2hG:,1X,E11.4)    
           
      RETURN
      END         !FIN DE LA PROCEDURE SOS_GSF_RONDEAUX_BREON



C==============================================================================
C PROCEDURE: SOS_F21SF_NADAL
C ==========
C      Cette procédure realise la decomposition en series de Fourier
C      (sur phi) de la fonction F(Teta1,Teta2,phi) :
C
C	 F(Teta1,Teta2,phi) = F12nadal(Teta1,Teta2,phi) 
C                            / F12fresnel(Teta1,Teta2,phi)
C	 avec Teta1 : angle zenithal d'incidence : tetas  (0 < Teta1 < PI/2),
C             Teta2 : angle zenithal de reflexion : tetav (0 < Teta2 < PI/2),
C             phi   : angle azimutal relatif entre les directions
C                     d'incidence et de reflexion selon la convention 
C                     des OS (0 <= phi <= PI),
C			
C             F12fresnel = 0.5 * (RR*RR - RL*RL) 
C             selon la convention de signe du modèle de Nadal.
C
C                avec RL le coefficient de reflexion parallele,
C                calcules en fonction de l'indice de réfraction de la surface
C                et de l'angle de reflexion teta :
C                RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                         / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C                RR le coefficient de reflexion perpendiculaire :
C                RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                         / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C                
C
C             F12nadal = alpha * 
C                  [ 1 -exp(-beta*F12fresnel / (cosTeta1 + cosTeta2) ) ]
C
C                 avec alpha et beta les paramètres du modèle de Nadal.
C
C      La BPDF du modèle de Nadal est convertie de reflectance polarisee 
C      en une grandeur utile pour les OS : luminance normalisée * cos(Teta2),
C      c'est a dire : BPDF * cos(Teta1) * cos(Teta2).			
C
C
C      La procédure enregistre les resultats dans le fichier "RES_MAT_REFLEX" 
C      sous la forme:
C         Fonction F decomposee pour (Teta1,Teta2) : (E(IS),IS=0,OS_NB)
C
C Description des paramètres
C --------------------------
C      LUM_NBMU (I4) (E) : Nombre d'angles (positifs) effectivement utiles.
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)	 (E) Tableau des mu de Gauss.
C      IND     (double)  (E) : Indice de réfraction de la surface.      	
C      ALPHA   (double)  (E) : Paramètre alpha du modèle de Nadal.         	
C      BETA    (double)  (E) : Paramètre beta du modèle de Nadal.   
C      OS_NB (I4) (E) : Ordre des développements des fonctions de phase 
C		          en Polynomes de Legendre a utiliser.     	
C      TRACE    (logical) (E) = vrai,si écriture dans le fichier trace
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
C       decomposition en series de Fourier de fonction F21 du modèle 
C       de Nadal.
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
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C     
C  Routines utilises par la procédure:
C ------------------------------------
C  Ce programme fait appel aux routines:
C        - SOS_CALC_F21_NADAL_SUR_FRESNEL
C
C==============================================================================
      SUBROUTINE SOS_F21SF_NADAL(LUM_NBMU,RMU,IND,ALPHA,BETA,
     &                           OS_NB,FIC_RES_GSF,TRACE,IER)
 	
      IMPLICIT NONE

C* Definition des variables    
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Tableau des mu de Gauss

      DOUBLE PRECISION IND	 ! Indice de réfraction surface / air
      DOUBLE PRECISION ALPHA !Paramètre alpha du modèle de Nadal.
      DOUBLE PRECISION BETA	 !Paramètre beta du modèle de Nadal.
      
      DOUBLE PRECISION E(0:CTE_OS_NB_MAX) ! Resultat de la decomposition en 
      				              ! series de Fourier de la fonction F.
	
      DOUBLE PRECISION PHI	! Azimut relatif pour la convention OS (radians).

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).

      DOUBLE PRECISION F	! Valeur de la fonction F(phi).
      DOUBLE PRECISION FMIN	! Valeur minimale retenue pour la fonction 
      				! F(phi) pour Theta1 et Theta2 fixes.
      DOUBLE PRECISION FMAX	! Valeur maximale de la fonction F(phi).

      DOUBLE PRECISION U(0:CTE_PH_NU)  ! Fonction F(phi) entre 0 et PI
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
      				    ! et celle recombinee 
				    !=> pour tous les couples (Teta1, Teta2).
      
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
      
      LOGICAL TRACE

      INTEGER*4 OS_NB ! Ordre des développements des fonctions de phase 
			 ! en Polynomes de Legendre.
				
      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations
      	    
      INTEGER*4 J	! Indice de tableau.
      INTEGER*4 IL	! Ordre maximale utile pour la fonction E.
      INTEGER*4 I1	! Indice de numéro d'angle Teta1.
      INTEGER*4 I2	! Indice de numéro d'angle Teta2.
      INTEGER*4 IS	! Ordre IS de la decomposition en series de Fourier.
      INTEGER*4 IS2	! Ordre IS de la decomposition en series de Fourier.
      INTEGER*4 IER  ! code d'erreur =0 si pas d'erreur, =-1 sinon  
      INTEGER*4 I
      INTEGER*4 COMPT	! Nombre de configurations pour la sommation des ecarts relatifs.
     
      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Nom du fichier RES_GSF
                                        ! avec arborescence complète      
	
C* Initialisation
C---------------------------------------------------------	
      ECARTMAX = -1.D+300
      SOMECART = 0.
      COMPT = 0
      
     
C* Creation du fichier resultat de SOS_MAT_REFLEXION
C---------------------------------------------------------
      OPEN(UNIT=30,FILE=FIC_RES_GSF,FORM='UNFORMATTED',ERR=991)

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
	
              PHI=Q*I
              CALL SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,
     &	                                           C1,S1,C2,S2,PHI,F)   
              U(I)=F
	      
 6         CONTINUE

		
C*   Determination des valeurs min et max de F
C*   pour I1 et I2 fixes mais PHI variant.
C----------------------------------------------
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
           ENDIF  !trace
	      
	 
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
	         PHI = I*Q		 		!Azimuth convention OS
	         Y=Y+U(I)*DCOS(IS*PHI)
	      ENDDO
	      E(IS) = Y*Q / INCTE_PI
	      
C--------- Determination de l'ecart maximal entre F exact et recombine : B1
C---------------------------------------------------------------------      
	      B1=0.D+00			!initialisation
	      DO I=0,CTE_PH_NU
	         PHI = I*Q		 		!Azimuth convention OS
	         CALL SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,
     &	                                          C1,S1,C2,S2,PHI,F)

C-------------- Calcul de F recombine à partir de la fonction E(IS) pour PHI fixe. 
                T1=E(0)	
		  DO IS2 = 1, IS		  
	              T1=T1+2.*E(IS2)*DCOS(IS2*PHI)	!Azimuth convention OS          		     
		  ENDDO !Fin boucle sur IS2
       
                B1=MAX( DABS((T1-F)/F) , B1)
		 		 
	      ENDDO !Fin boucle sur I (phi = 0 a PI)
	           	     

C--------- Test d'arret de la decomposition :
C---------    * Arret si l'ecart relatif maximal entre F exacte et recombinee est
C---------      inferieur au seuil CTE_SEUIL_SF_NADAL
C---------      ou s'il augmente par rapport a l'ordre precedent.
C---------    * Poursuite de la boucle sur IS sinon.
C-------------------------------------------------------------------------------- 		
              IF (B1.GT.CTE_SEUIL_SF_NADAL) GOTO 14
              IL=IS	      
              GOTO 17	!Sortie de boucle sur IS
	      
 14           IF(B1.LT.B1_PREC) GOTO 141
              IL=IS-1
	      GOTO 17	!Sortie de boucle sur IS
	      
 141	      B1_PREC=B1	        
 13        CONTINUE	!Fin de boucle sur IS     
	   IL=OS_NB



C* Enregistrement pour le couple (Teta1, Teta2) de la decomposition
C* en series de Fourier de BPDF (convertie en luminance normalisée
C* et adaptee au code des OS).
C-----------------------------------------------------------------
 17        WRITE(30,err=992)I1,I2,IL,(E(J),J=0,IL)


C* Fichier Trace : Calculs de paramètres d'analyse de qualite
C*                 de la decomposition en series de Fourier
C*                 et enregistrements.
C-----------------------------------------------------------------
           IF (TRACE) THEN

	      ECARTMAX_C1C2=0.
	            	      
	      Q=INCTE_PI/CTE_PH_NU
	      DO I=0,CTE_PH_NU
	         PHI=I*Q
	         CALL SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,
     &	                                             C1,S1,C2,S2,PHI,F)   
		 
                FRECOMB = E(0)
	         DO IS=1,IL
	            FRECOMB = FRECOMB + 2.*DCOS(IS*PHI)*E(IS)
	         ENDDO
            
	         IF(F.NE.0.D+00) THEN
	            ECART = DABS((FRECOMB-F)/F)  
		    
C*       Analyse pour tous les azimuts du couple courant (Teta1, Teta2)	
		     IF(ECART.GT.ECARTMAX_C1C2) THEN
		        ECARTMAX_C1C2 = ECART
		        PHIECARTMAX_C1C2 = PHI*180./INCTE_PI	!Azimuth convention OS
	               F_ECARTMAX_C1C2 = F
	               FRECOMB_ECARTMAX_C1C2 = FRECOMB
		     ENDIF 
		    
C*       Analyse pour tous les azimuts et tous les couples (Teta1, Teta2)		    
	            IF(ECART.GT.ECARTMAX) THEN
	               ECARTMAX = ECART
	               PHIECARTMAX = PHI*180./INCTE_PI	!Azimuth convention OS
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
     &					    F_ECARTMAX_C1C2,
     &		                         FRECOMB_ECARTMAX_C1C2,
     &					    ECARTMAX_C1C2*100.
	      ELSE
	      	  WRITE(99,2500,err=921) IL,PHIECARTMAX_C1C2,
     &				 	    F_ECARTMAX_C1C2,
     &		                         FRECOMB_ECARTMAX_C1C2
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
     &	  "Analysis of the Fourier expansion of the NADAL's BPDF" 
	
	  IF(COMPT.GT.1) THEN
	     WRITE(99,*,err=921)
	     WRITE(99,*,err=921)
	     WRITE(99,*,err=921) 
     &	          'Mean relative difference from the recombination'//
     &	          ' (% excluding null values):', SOMECART*100/COMPT
	     WRITE(99,*,err=921)
	     WRITE(99,*,err=921)   
     &	          'Max relative difference from the recombination (%):',
     &	          ECARTMAX*100
            WRITE(99,*,err=921)   
     &	          'Related delta phi angle (SOS convention) :', 
     &	          PHIECARTMAX
            WRITE(99,*,err=921)    
     &	          'Related viewing angle :',TVECARTMAX
	     WRITE(99,*,err=921)    
     &	          'Related solar angle :',TSECARTMAX
	     WRITE(99,*,err=921)   
     &	          'Related value of BRDF :',F_ECARTMAX
	     WRITE(99,*,err=921)
     &	          'Related value of the recombined BRDF :',
     &	          FRECOMB_ECARTMAX
         ELSE
	     WRITE(99,*,err=921)
	     WRITE(99,*,err=921) 'WARNING : nulle BPDF !'
	     WRITE(99,*,err=921) 'Simulation parameters to be checked'
	 ENDIF 
	 
	 WRITE(99,*,err=921)
	 WRITE(99,*,err=921)
	 
      END IF 
      
     
      
C* Fin nominale 
C-------------------
      IER=0
      goto 9999
      
C* Cas d'erreur : retour du status -1 au programme appelant
C------------------------------------------------------
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Surface'
      IER=-1
      goto 9999
  991 WRITE(6,*) 'SOS_F21SF_NADAL : Erreur a l''ouverture d''un fichier'
      IER=-1
      goto 9999
  992 WRITE(6,*) 'SOS_F21SF_NADAL',
     &		 'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      goto 9999
      
 9999 CONTINUE


C* Fermeture fichier :'RES_MAT_REFLEX'
C-------------------------------------
      CLOSE(30)
 
 
C* Format
C--------

1000  FORMAT(6hTeta1:,1X,F5.2,2X,6hTeta2:,1X,F5.2,2X,
     &15h; (PHIos,Fmin):,1X,F6.2,2X,E11.4,2X,
     &15h; (PHIos,Fmax):,1X,F6.2,2X,E11.4)    
      
2000  FORMAT(5X,12hLimit Order:,I3,3h ; ,14h Error max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4,2X,
     &12h, Ecart (%):,1X,F7.3)         
     
2500  FORMAT(5X,12hLimit Order:,I3,3h ; ,14h Error max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4)         
                           
      RETURN
      END         !FIN DE LA PROCEDURE SOS_F21SF_NADAL



C==============================================================================
C PROCEDURE: SOS_CALC_F21_NADAL_SUR_FRESNEL
C		
C      Cette procédure calcule la fonction G(Teta1, Teta2, phi) pour le
C      modèle de Nadal :
C
C      G(Teta1, Teta2, phi) 
C        =  F21Nadal(Teta1, Teta2, phi) / F21Fresnel(Teta1, Teta2, phi)
C
C        avec F21Nadal(Teta1, Teta2, phi)
C             = Alpha * [ 1 -exp(-Beta*F21Fresnel / ( cosTeta1+cosTeta2) ) ]
C
C        et F21Fresnel(Teta1, Teta2, phi) = 0.5 * (RR*RR - RL*RL) 
C        selon la convention de signe du modèle de Nadal.
C
C            RL est le coefficient de reflexion parallele
C
C               RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C            RR est le coefficient de reflexion perpendiculaire
C
C               RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C            calcules en fonction de l'indice de réfraction de la surface
C            de l'angle de reflexion teta.
C
C
C
C      La BPDF du modèle de Nadal est convertie de reflectance polarisee 
C      en une grandeur utile pour les OS : luminance normalisée * cos(Teta2),
C      c'est a dire : BPDF * cos(Teta1) * cos(Teta2).			
C
C Description des paramètres
C --------------------------
C       IND	(double)   (E) Indice de réfraction surface / air.
C       ALPHA	(double)   (E) Paramètre alpha du modèle de Nadal.
C       BETA	(double)   (E) Paramètre beta du modèle de Nadal.	
C	 C1	(double)   (E) Cos(Teta1).
C	 S1	(double)   (E) Sin(Teta1).
C	 C2	(double)   (E) Cos(Teta2).
C	 S2	(double)   (E) Sin(Teta2).
C       PHI   (double)   (E) Azimut relatif pour la convention OS (radians).
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
      SUBROUTINE SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,
     &	                                        C1,S1,C2,S2,PHI,F)   

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------------------------------------
      DOUBLE PRECISION IND	! (E)  Indice de réfraction surface / air
      DOUBLE PRECISION ALPHA! (E)  Paramètre alpha du modèle de Nadal.
      DOUBLE PRECISION BETA	! (E)  Paramètre beta du modèle de Nadal.
      DOUBLE PRECISION C1	! (E)  Cos(Teta1).
      DOUBLE PRECISION S1	! (E)  Sin(Teta1).
      DOUBLE PRECISION C2	! (E)  Cos(Teta2).
      DOUBLE PRECISION S2	! (E)  Sin(Teta2).
      DOUBLE PRECISION PHI	! (E)  Azimut relatif pour la convention 
      				!      OS (radians).
            
      DOUBLE PRECISION F	! (S) Valeur de la fonction F21nadal/F21fresnel
      				! pour la geometrie d'observation
				! (convertie en luminance normalisée * cos(teta2).	
      
      DOUBLE PRECISION F21FRESNEL	! Valeur de la fonction F21 de Fresnel.
      DOUBLE PRECISION F21NADAL	! Valeur de la fonction F21 de Nadal.
      
      DOUBLE PRECISION COEF

      DOUBLE PRECISION COSDIF	! Cosinus de l'angle de diffusion.      
      DOUBLE PRECISION COSW	! Cosinus de l'angle de reflexion W
      				! (angle d'incidence par rapport a la 
				!  normale du dioptre)
      
      DOUBLE PRECISION V	! Sinus carre de W.
      DOUBLE PRECISION X	! 

      DOUBLE PRECISION IND2	! Carre de l'indice de réfraction surface / air      
      DOUBLE PRECISION RL	! Coefficient de reflexion parallele.
      DOUBLE PRECISION RR	! Coefficient de reflexion perpendiculaire.

      DOUBLE PRECISION PI	! Valeur de PI
      

C* Traitements
C---------------------------------    
      PI=4.D+00*DATAN(1.D+00)

      COEF = 4.*C1*PI
	
C*   Calcul du cosinus de l'angle de diffusion                
C---------------------------------------------------------     	
      COSDIF=-C1*C2+DSQRT(1.-C1*C1)*DSQRT(1.-C2*C2)*DCOS(PHI)

C*   Calcul du cosinus de l'angle de reflexion                
C---------------------------------------------------------   
      COSW=DSQRT(.5*(1-COSDIF))
      
C*   Calcul du sinus au carre de l'angle de reflexion                
C--------------------------------------------------------- 
      V=.5*(1+COSDIF)  
     
C*   Calcul des coefficients de reflexion parallele et
C*   perpendiculaire selon la theorie de Fresnel                
C---------------------------------------------------------         
      IND2=IND*IND         
      X=DSQRT(IND2-V)     
      RL=(IND2*COSW-X)/(IND2*COSW+X)
      RR=(COSW-X)/(COSW+X)
    

C*   Calcul de F21 de Fresnel (convention de signe de Nadal)                
C--------------------------------------------------------- 
      F21FRESNEL = 0.5D+00 * (RR*RR-RL*RL)      
      
C*   Calcul de F21 de Nadal                 
C--------------------------------------------------------- 
      F21NADAL = -BETA * F21FRESNEL / (C1 + C2)     
      F21NADAL = ALPHA * ( 1.D+00 - DEXP(F21NADAL) )

C*   Calcul de F               
C---------------------------------------------------------            
      IF(F21FRESNEL.LT.1.0D-10) THEN  !Développement limité au voisinage 
      					   !de 0 pour F21Fresnel
         F = ALPHA*BETA / (C1 + C2)
      ELSE	 
         F = F21NADAL / F21FRESNEL
      ENDIF


C*   Conversion reflectance -> luminance normalisée (*C1)
C*   et adaptation aux OS -> *C2.  
C--------------------------------------------------------- 
      F = F*COEF*C2*C1		

      RETURN
      
      END         !FIN DE LA PROCEDURE SOS_CALC_F21_NADAL_SUR_FRESNEL







C==============================================================================
C procédure: SOS_GSF_MAIGNAN
C ==========
C
C      Cette procédure realise la décomposition en séries de Fourier
C      (sur phi) de la fonction G(Teta1,Teta2,phi) :
C
C	 G(Teta1,Teta2,phi) = SIG * exp [-tan(TetaN)] / [ 4 * (cosTeta1+cosTeta2) ]
C
C	 avec  Teta1 : angle zénithal d'incidence  (0 < Teta1 < PI/2),
C             Teta2 : angle zénithal de réflexion (0 < Teta2 < PI/2),
C             phi   : angle azimutal relatif entre les directions
C                     d'incidence et de réflexion,
C
C	       TetaN : angle zénithal de la normale à la facette de vague,
C                     calcule à partir de Teta1, Teta2 et phi.
C			
C             SIG   : carre du paramètre Sigma de la théorie de Cox et Munck.
C
C
C      Enregistre les résultats dans le fichier "RES_GSF" sous la forme:
C	   Numéro d'angle Teta1 : I1,
C         Numéro d'angle Teta2 : I2,
C         Ordre maximal utile de la décomposition: IL,
C         Fonction G décomposée pour (Teta1,Teta2) : (E(J),J=0,IL)
C
C Description des paramètres
C --------------------------
C      LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles.
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double) : (E) Tableau des mu de Gauss.
C      OS_NM (I4) : (E) Ordre maximal du développement en séries de Fourier  
C                       de la fonction G des calculs de matrice de réflexion.
C      FIC_RES_GSF  (CHARACTER)  : (E) Nom du fichier temporaire RES_GSF
C                                      avec arborescence complète.
C
C      IER  (I4) : (S) code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C Autre résultat fourni par la procédure
C -------------------------------------
C      Fichier temporaire RES_GSF
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier temporaire RES_GSF
C       Fichier binaire non formatte, créé par SOS_GSF
C       Fichier créé dans le répertoire de lancement du programme
C       Si un fichier de même nom existe déjà, il est écrasé
C       Chaque enregistrement contient :
C           -numéro d'angle Teta1 (valeur entière)
C           -numéro d'angle Teta2 (valeur entière)
C           -ordre maximale utile de la decompositon : IL  (valeur entière)
C           -Fonction G décomposée pour (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 valeurs double précision)
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
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C     
C  Routines utilises par la procédure:
C ------------------------------------
C  Ce programme fait appel a la routine:
C        - SOS_CALCG_MAIGNAN
C
C==============================================================================
      SUBROUTINE SOS_GSF_MAIGNAN(LUM_NBMU,RMU,COEF_C,OS_NM,
     &                           FIC_RES_GSF,IER)
      IMPLICIT NONE

C* Definition des variables    
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Tableau des mu de Gauss

      DOUBLE PRECISION COEF_C	! Coefficient C du modèle de BPDF de Maignan

      DOUBLE PRECISION E(0:CTE_OS_NM_MAX) ! Resultat de la décomposition en 
      				              ! séries de Fourier de la fonction G.

      DOUBLE PRECISION PHI	! Azimut relatif (radians).

      DOUBLE PRECISION G	! Valeur de la fonction G(phi).

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).
      DOUBLE PRECISION S12	! Sin(Teta1) * Sin(Teta2).

      DOUBLE PRECISION GMIN	! Valeur minimale retenue pour la fonction 
      				! G(phi) pour Theta1 et Theta2 fixes.
      DOUBLE PRECISION GMAX	! Valeur maximale de la fonction G(phi) = G(0).

      DOUBLE PRECISION PHIB	! Angle azimutal phi maximal pour la 
      				! description de la fonction G(phi).
			
      DOUBLE PRECISION PHI1	! Angle intermediaire pour le calcul de PHIB.
      DOUBLE PRECISION PHI2	! Angle intermediaire pour le calcul de PHIB.
						
      DOUBLE PRECISION U(0:CTE_PH_NU)  ! Fonction G(phi) entre 0 et PHIB
      			                  ! pour Theta1 et Theta2 fixes.

      DOUBLE PRECISION Q	! Pas en phi.

      DOUBLE PRECISION T2	! Recombinaison de la fonction G(phi) a 
      				! partir de E(IS).
      DOUBLE PRECISION T1	! Sommation des E(IS).

      DOUBLE PRECISION X,Y,Z
      DOUBLE PRECISION XT,B1

      INTEGER*4 OS_NM  ! Ordre maximal du développement en séries de Fourier  
			  ! de la fonction G des calculs de matrice de réflexion.
			 
      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations
      
      INTEGER*4 IL	! Ordre maximale utile pour la fonction E.
      INTEGER*4 I1	! Indice de numéro d'angle Teta1.
      INTEGER*4 I2	! Indice de numéro d'angle Teta2.
      INTEGER*4 IS	! Ordre IS de la décomposition en séries de Fourier.
      INTEGER*4 IER  ! code d'erreur =0 si pas d'erreur, =-1 sinon  

      INTEGER*4 I
      INTEGER*4 IA,IP,K,J
      
      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Nom du fichier RES_GSF
                                        ! avec arborescence complète



C* Ouverture du fichier résultat
C--------------------------------
      OPEN(UNIT=30,FILE=FIC_RES_GSF,FORM='UNFORMATTED',ERR=991)

C* Boucle sur l'angle Teta1 > 0.
C---------------------------------
      DO  3 I1=1 ,LUM_NBMU

C* Cosinus et sinus de Teta1
C---------------------------------
         C1=RMU(I1)
         S1=DSQRT(1-C1*C1)

C*   Boucle sur l'angle teta2 > 0.
C---------------------------------
         DO 4 I2=1 ,I1

C*   Cosinus et sinus de Teta2
C---------------------------------
           C2=RMU(I2)
           S2=DSQRT(1-C2*C2)
           S12=S1*S2

C*   Determination de G pour phi = 0.
C-------------------------------------
           PHI=0.0
           CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)
           U(0)=G
           GMAX=G

C*   Determination de G pour phi = INCTE_PI.
C-----------------------------------------
           PHI=INCTE_PI
           CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)
           U(CTE_PH_NU)=G
           GMIN=G
		

C*   Evaluation de l'importance de l'ecart entre Gmin et Gmax.
C*      1er cas: L'ecart est modere.
C*               Pour Teta1 et Teta2 fixes, on definit la fonction
C*               G(phi) entre 0 et PI par pas de PI/CTE_PH_NU.
C*      2nd cas: L'ecart est important.
C*               La fonction G(phi) devient insignifiante devant Gmax
C*		   à partir d'un angle PHIB a determiner.
C*               Pour Teta1 et Teta2 fixes, on definit la fonction
C*               G(phi) entre 0 et PHIB par pas de PHIB/CTE_PH_NU.
C----------------------------------------------------------------------
           X=CTE_PH_TEST*GMIN
           IF(X.GE.GMAX) THEN
C* 1er cas: L'ecart entre Gmin et Gmax est modere.
C--------------------------------------------------
              PHIB=INCTE_PI

              Q=INCTE_PI/CTE_PH_NU
              DO 6 I=1,CTE_PH_NU
                 PHI=Q*I
                 CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)
                 U(I)=G
 6            CONTINUE

           ELSE
C* 2eme cas: L'ecart entre Gmin et Gmax est important.
C*          On recherche l'angle PHIB au dela duquel G(phi) devient
C*          insignifiante devant Gmax.
C------------------------------------------------------------------
C*      Initialisation
              PHI1=0
              PHI2=INCTE_PI
	
C*      Estimation de PHIB
C* (boucle sur l'etiquette 9)
C---------------------------	
 9            PHIB=.5*(PHI1+PHI2)

C*      Calcul de G(PHIB)
C---------------------------
              CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHIB,COEF_C,G)
	
C*      Controle de l'importance de l'ecart entre G(PHIB) et Gmax.
C*      Cas 2A: L'ecart reste important.
C*              PHIB doit etre re-estime.
C*              On resserre PHI1 ou PHI2, selon le signe de l'ecart.
C*              On estime a nouveau PHIB comme la moyenne de PHI1 et PHI2
C*              et on calcule la nouvelle valeur G(PHIB).
C*      Cas 2B: L'ecart est modere.
C*              PHIB est une estimation satisfaisante de la valeur de
C*              l'angle phi jusqu'à laquelle G(phi) reste significative
C*              devant Gmax. La fonction G(phi) est decrite entre 0 et 
C*              PHIB par pas de PHIB/CTE_PH_NU.	
C-------------------------------------------------------------------
              X=CTE_PH_TEST*G
              IF (DABS(X-GMAX).LT..01*GMAX) GOTO 8
	
C*      Cas 2A: L'ecart entre G(PHIB) et GMAX reste important.
C*              On affine la valeur de PHIB.
C-------------------------------------------------------------
              IF (X.LE.GMAX) GOTO 88
              PHI1=PHIB
              GOTO 9
 88           PHI2=PHIB
              GOTO 9
	
C*      Cas 2B: L'ecart entre G(PHIB) et GMAX est modere.
C*              L'estimation de PHIB est satisfaisante.
C*              On definit G(phi) entre 0 et PHIB par pas de 
C*              PHIB/CTE_PH_NU.
C*              Gmin devient G(PHIB).	
C--------------------------------------------------------------------
 8            Q=PHIB/CTE_PH_NU

              DO 10 I=1,CTE_PH_NU
                 PHI=Q*I
                 CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)
                 U(I)=G
 10           CONTINUE
              GMIN=U(CTE_PH_NU)
              T1=0.
              T2=0.
           ENDIF


C*   Realisation de la décomposition en serie de Fourier.
C*   Boucle sur l'ordre de la serie.                     	
C--------------------------------------------------------------------
           DO 13 IS=0,OS_NM
              Z=.5*(GMAX+GMIN*DCOS(IS*PHIB))

              IA=1
              DO 14 I=1,CTE_PH_NQ
                 IA=2*IA
                 IP=CTE_PH_NU/IA
                 Y=0
                 DO 15 J=1,IA,2
                    K=IP*J
                    Y=Y+U(K)*DCOS(IS*K*Q)
 15              CONTINUE
                 Y=2*Y/IA
	
                 XT=DABS(Z-Y)/Z
                 IF (XT.LT..0001) GOTO 16
                 Z=.5*(Y+Z)
 14           CONTINUE

 16           E(IS)=PHIB*Z/INCTE_PI
              IF (IS.NE.0) GOTO 133
              T1=E(0)
              T2=E(0)
              GOTO 13
	
 133          T1=T1+2*E(IS)
              T2=T2+2*DCOS(IS*PHIB*.5)*E(IS)

C*   Test d'arret de la décomposition - poursuite de la boucle sur IS sinon.
C---------------------------------------------------------------------------
              B1=DABS(T1-GMAX)/GMAX
              IF (B1.GT..001) GOTO 13
              IL=IS
              GOTO 17	!Sortie de boucle
 13        CONTINUE	!Fin de boucle sur IS
           IL=OS_NM

C*   Enregistrement des résultats
C---------------------------------
 17        WRITE(30,err=992)I1,I2,IL,(E(J),J=0,IL)

 4       CONTINUE	!Fin de boucle sur Teta2 (I2)
 
 3    CONTINUE	!Fin de boucle sur Teta1 (I1)
      
C* Fin nominale 
C-------------------
      IER=0
      GOTO 9999
      
C* Cas d'erreur : retour du status -1 au programme appelant
C------------------------------------------------------
  991 WRITE(6,*) 'SOS_GSF Erreur a l''ouverture d''un fichier'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'SOS_GSF',
     &		 'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      GOTO 9999

 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(30)
      
      RETURN
      END         !FIN DE LA PROCEDURE SOS_GSF_MAIGNAN





C==============================================================================
C PROCEDURE: SOS_CALCG_MAIGNAN
C ==========
C       La fonction de pondération de la matrice de Fresnel pour le modèle
C       de BPDF de Maignan est :
C	 G(Teta1,Teta2,phi) = C * exp [-tan(TetaI)] / [ 4 * (cosTeta1+cosTeta2) ]
C
C	 avec Teta1 : angle zénithal d'incidence  (0 < Teta1 < PI/2),
C            Teta2 : angle zénithal de réflexion (0 < Teta2 < PI/2),
C            phi   : angle azimutal relatif entre les directions
C                    d'incidence et de réflexion,
C
C	      TetaI : angle d'incidence par rapport à la normale à la facette de surface,
C                    calculé à partir de Teta1, Teta2 et phi.
C
C	      C : coefficient du site (indépendant de la géométrie d'observation).	
C
C      Cette procédure calcule : G*4*cosTeta1*cosTeta2 
C      soit C * exp [-tan(TetaI)] / [ 1/cosTeta1 + 1/cosTeta2) ]
C
C Description des paramètres
C --------------------------
C	C1	(double)   (E) Cos(Teta1).
C	C2	(double)   (E) Cos(Teta2).
C	S12	(double)   (E) Sin(Teta1) * Sin(Teta2).
C	PHI	(double)   (E) Azimut relatif (radians).
C	COEF_C	(double)   (E) Coefficient C du modèle de Maignan.
C 	G	(double)   (S) Valeur de G(TetaN).
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
      SUBROUTINE SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------------------------------------
      DOUBLE PRECISION PHI	! Azimut relatif (radians).

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S12	! Sin(Teta1) * Sin(Teta2).      
      DOUBLE PRECISION COS_2I  ! Cosinus de deux fois l'angle d'incidence.
      DOUBLE PRECISION TAN2_I  ! Carré de la tangente de l'angle d'incidence.
      DOUBLE PRECISION COEF_C ! Coefficient C du modèle de BPDF de Maignan
      DOUBLE PRECISION G	! Valeur de G(Teta1,Teta2,phi).	

C* Calcul
C---------------------------------
      COS_2I = C1*C2-S12*DCOS(PHI)
      TAN2_I = (1-COS_2I) / (1+COS_2I)

      !Gestion du cas d'erreur d'arrondis numériques (tangente carrée négative impossible).
      IF (TAN2_I.LT.0.D+00) TAN2_I=0.D+00

      G= COEF_C *DEXP(-DSQRT(TAN2_I)) / (1./C1 + 1./C2)

      RETURN
      END         !FIN DE LA PROCEDURE SOS_CALCG_MAIGNAN








