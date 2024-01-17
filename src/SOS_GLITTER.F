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
C* FICHIER: SOS_GLITTER.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Genere les fichiers "glitter" pour le calcul du glitter.
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, B.Lafrance
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0: Mise a niveau du code des OS du LOA
C* MOD:VERSION:1.1: Mise à jour des commentaires
C*                  
C* MOD:VERSION:2.0: 10/01/2010 
C*
C*   - Suppression de la constante SOS_OS_FICANGLE (de SOS.h) et de la lecture du fichier correspondant 
C*     des angles de Gauss au profit : d'un passage par argument des tables d'angles (RMU) et poids (CHR).
C*
C*   - Modification de la taille des tableaux de luminance :
C*     ==> Passage de CTE_OS_NBMU à CTE_OS_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'à LUM_NBMU (nouveau argument d'entrée) 
C*         égal au nombre effectif d'angles utilises.
C*
C*   - Suppression des constantes CTE_OS_NB, CTE_OS_NS et CTE_OS_NM (de SOS.h)  
C*     au profit : d'un passage par argument des valeurs OS_NB, OS_NS, OS_NM.
C*
C*   - Modification de la taille des tableaux de décomposition en polynomes de Legendre :
C*     ==> Passage de CTE_OS_NB à CTE_OS_NB_MAX.
C*     ==> Passage de CTE_OS_NS à CTE_OS_NS_MAX.
C*     ==> Passage de CTE_OS_NM à CTE_OS_NM_MAX.
C*     ==> Boucle sur les ordre jusqu'à OS_NB, OS_NS ou OS_NM (nouveaux argument d'entrée)
C*
C*   - Routine SOS_GSF : Modification de l'appel de la routine pour inclure l'argument LUM_NBMU (nb d'angles utiles)  
C*     et pour limiter les boucles aux bornes utiles LUM_NBMU, et pour inclure l'argument OS_NM.  
C*
C*   - Adaptation de l'appel des routines SOS_MAT_FRESNEL, SOS_MAT_REFLEXION  
C*     et SOS_MISE_FORMAT, aux passages des paramètres LUM_NBMU (nb d'angles utiles),
C*     RMU (table des cosinus des angles) et CHR (table des poids des angles),
C*     ainsi que OS_NB, OS_NS et OS_NM.
C*
C* MOD:VERSION:6.0: 06/06/2013
C*
C*   - Modification de la valeur de Pi, déclarée en constante dans le fichier
C*     (INCTE_PI remplace CTE_PI). 
C*
C* MOD:VERSION:7.0: 16/06/2020
C*
C*   - Passage du nom des fichiers temporaires (RES_GSF, RES_FRESNEL, RES_MAT_REFLEX) 
C*     avec arborescence complète.
C*   - Modification en conséquence des arguments fournis à la routine SOS_GSF.
C*   - Appel de la routine SOS_GSF avec le paramètre SIG (de Cox&Munk) et de la routine 
C*     SOS_MAT_REFLEXION avec le paramètre (1./SIG) au lieu de la vitesse de vent.
C*
C* MOD:VERSION:7.1: 31/08/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Conversion des traces en Anglais.
C*
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
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
C PROCEDURE: SOS_GLITTER
C ==========      
C      Cette procédure genere les fichiers de reflexion de surface
C      pour le calcul du glitter.
C
C
C Donnees en entrée du programme 
C ------------------------------
C    -->Les paramètres/arguments du programme
C       LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C	 RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (E) Cosinus des angles.
C	 CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (E) Poids des angles.
C 	 WIND  (DOUBLE)  : (E) Vitesse du vent en m/s. 
C       IND   (DOUBLE)  : (E) Indice de réfraction de l'eau / air.
C       OS_NB (I4) : (E) Ordre des développements des fonctions de phase 
C		           en Polynomes de Legendre à utiliser.
C       OS_NS (I4) : (E) Ordre maximal à utiliser pour le développement en
C                        polynomes de Legendre des éléments de la matrice de Fresnel 
C                        et en séries de Fourier des calculs de luminance.
C       OS_NM (I4) : (E) Ordre maximal du développement en séries de Fourier  
C                        de la fonction G des calculs de matrice de reflexion.
C       FIC_RES_FRESNEL (CHARACTER*CTE_LENFIC2)  : 
C                 (E) Nom du fichier temporaire RES_FRESNEL avec arborescence complète.
C       FIC_RES_GSF (CHARACTER*CTE_LENFIC2)  :  
C                 (E) Nom du fichier temporaire RES_GSF avec arborescence complète.
C       FIC_RES_FRESNEL (CHARACTER*CTE_LENFIC2)  :  
C                 (E) Nom du fichier temporaire RES_FRESNEL avec arborescence complète.
C       FIC_RES_MAT_REFLEX (CHARACTER*CTE_LENFIC2)  :  
C                 (E) Nom du fichier temporaire RES_MAT_REFLEX avec arborescence complète.
C       FICGLITTER  (CHARACTER*CTE_LENFIC2) : (S) nom complet du fichier Glitter résultat.
C       TRACE  (LOGICAL)   ! = vrai,si écriture dans le fichier trace.
C       IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon
C
C
C Resultat fourni par le programme
C --------------------------------
C   Fichier Glitter
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
C   -->Contenu du fichier Glitter
C       Fichier binaire non formatte, cree par SOS_MISE_FORMAT
C       Il contient la matrice de reflexion de surface ,
C       c'est a dire les éléments Pij de la matrice de reflexion 
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
C       Fichier binaire non formatte, cree par SOS_GSF
C       Fichier cree dans le repertoire de lancement du programme
C       Si un fichier de meme nom existe deja, il est ecrase
C       Chaque enregistrement contient :
C           -numéro d'angle Teta1 (valeur entiere)
C           -numéro d'angle Teta2 (valeur entiere)
C           -ordre maximale utile de la décomposition : IL  (valeur entiere)
C           -Fonction G décomposée pour (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 valeurs double précision)
C
C   -->Contenu du fichier temporaire RES_FRESNEL
C      Fichier Ascii, format d'un enregistrement 4E15.8, 
C      cree par SOS_MAT_FRESNEL
C      Fichier cree dans le repertoire de lancement du programme
C      Si un fichier de mSOS_RONDEAUX_BREONeme nom existe deja, il est ecrase
C      Nombre d'enregistrements : de k=0 a OS_NS
C      Chaque enregistrement contient les paramètres ALPHA(k), BETA(k), 
C      GAMMA(k) et ZETA(k) associes a la matrice de Fresnel
C
C   -->Contenu du fichier temporaire RES_MAT_REFLEX
C       Fichier binaire non formatte, cree par SOS_MAT_REFLEXION
C       Fichier cree dans le repertoire de lancement du programme
C       Si un fichier de meme nom existe deja, il est ecrase
C       Enregistrement pour chaque couple (Teta1, Teta2) de
C       la décomposition en séries de Fourier de la matrice de Fresnel:
C       M11,M12,M13,M21,M22,M23,M31,M32,M33 
C       (n valeurs réelles)
C
C    Les 3 fichiers temporaires sont ecrases s'ils existent, puis 
C    detruits en fin de programme.
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
C      Erreur si le fichier Glitter existe deja
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, arret du programme et 
C     retour du status 1 au shell
C        
C
C  Sous programmes utilises:
C --------------------------
C  Cette procédure fait appel aux routines:
C      - SOS_GSF
C      - SOS_MAT_FRESNEL 
C      - SOS_MAT_REFLEXION 
C      - SOS_MISE_FORMAT
C
C=============================================================================

      SUBROUTINE SOS_GLITTER(LUM_NBMU,RMU,CHR,WIND,IND,
     &                       OS_NB,OS_NS,OS_NM,
     &                       FIC_RES_GSF,FIC_RES_FRESNEL,
     &                       FIC_RES_MAT_REFLEX,
     &                       FICGLITTER,TRACE,IER)

      IMPLICIT NONE


C* Definition des variables   
C*--------------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICGLITTER	! Nom complet du fichier GLITTER résultat.

      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Nom du fichier RES_GSF
                                        ! avec arborescence complète

      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Nom du fichier RES_FRESNEL
                                            ! avec arborescence complète

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Nom du fichier RES_MAT_REFLEX
                                               ! avec arborescence complète
            
      DOUBLE PRECISION WIND	! Vitesse du vent (m/s).

      DOUBLE PRECISION SIG 	! Carré du paramètre Sigma de la théorie de Cox et Munck.
      
      DOUBLE PRECISION IND	! Indice de réfraction de l'eau / air
	
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Tableau des mu de Gauss
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Poids de Gauss.
      
      LOGICAL TRACE     ! = vrai,si écriture dans le fichier trace

      INTEGER*4 OS_NB   ! Ordre des développements des fonctions de phase 
			   ! en Polynomes de Legendre.
      INTEGER*4 OS_NS   ! Ordre maximal du développement en polynomes de Legendre   
			   ! des éléments de la matrice de Fresnel.
      INTEGER*4 OS_NM   ! Ordre maximal du développement en séries de Fourier  
			   ! de la fonction G des calculs de matrice de reflexion.
			
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
         WRITE(99,*,err=921) 'Calculation of Sun Glint file'  
         WRITE(99,*,err=921) ' '
         WRITE(99,*,err=921)'Simulation parameters '
         WRITE(99,*,err=921)' Wind velocity  (m/s) ',WIND
         WRITE(99,*,err=921)
     &        ' Refractive index of sea water ',IND
         WRITE(99,*,err=921)' OS_NB ',OS_NB
	  WRITE(99,*,err=921)' OS_NS ',OS_NS
	  WRITE(99,*,err=921)' OS_NM ',OS_NM
	  WRITE(99,*,err=921) ' ' 	 	 
      END IF
        
C* Carre du paramètre Sigma de la théorie de Cox et Munk.
C--------------------------------------------------------
      SIG = .003 + .00512*WIND

C* Decomposition en séries de Fourier de la fonction G(Teta1,Teta2,phi).
C---------------------------------------------------------------------
      CALL SOS_GSF(LUM_NBMU,RMU,SIG,OS_NM,FIC_RES_GSF,IER)
      IF (IER.NE.0) GOTO 995   ! cas d'erreur
      
C* Développement en fonctions de Legendre de la matrice de Fresnel.
C------------------------------------------------------------------
      CALL SOS_MAT_FRESNEL(LUM_NBMU,RMU,CHR,IND,OS_NS,
     &                     FIC_RES_FRESNEL,TRACE,IER)
      IF (IER.NE.0) GOTO 996   ! cas d'erreur
      
C* Calcul de la matrice de reflexion developpee en séries de Fourier.
C------------------------------------------------------------------
      CALL SOS_MAT_REFLEXION((1./SIG),LUM_NBMU,RMU,OS_NB,OS_NS,OS_NM,
     &                       FIC_RES_FRESNEL,FIC_RES_GSF,
     &                       FIC_RES_MAT_REFLEX,IER)
      IF (IER.NE.0) GOTO 997   ! cas d'erreur
      
C* Rearrangement des résultats pour etre compatible avec le code des OS.
C-----------------------------------------------------------------------
      CALL SOS_MISE_FORMAT(FIC_RES_MAT_REFLEX, FICGLITTER,
     &                     LUM_NBMU,OS_NB,IER)
      IF (IER.NE.0) GOTO 998   ! cas d'erreur


C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  912 WRITE(6,*) 'Erreur a l''ouverture du fichier angles de Gauss'
      IER=-1
      GOTO 9999
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Glitter'
      IER=-1
      GOTO 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1
      GOTO 9999      
  995 WRITE(6,*) 'Erreur dans la routine SOS_GSF'
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
      END    !FIN DU PROGRAMME SOS_GLITTER  





C==============================================================================
C PROCEDURE: SOS_GSF
C ==========
C      Cette procédure realise la décomposition en séries de Fourier
C      (sur phi) de la fonction G(Teta1,Teta2,phi) :
C
C	 G(Teta1,Teta2,phi) = exp [-tan(TetaN)**2 / SIG] / cos(TetaN)**4
C
C	 avec  Teta1 : angle zénithal d'incidence  (0 < Teta1 < PI/2),
C             Teta2 : angle zénithal de reflexion (0 < Teta2 < PI/2),
C             phi   : angle azimutal relatif entre les directions
C                     d'incidence et de reflexion,
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
C      SIG	(double) : (E)  Carré du paramètre Sigma de la théorie de Cox et Munck.
C      OS_NM (I4) : (E) Ordre maximal du développement en séries de Fourier  
C                       de la fonction G des calculs de matrice de reflexion.
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
C       Fichier binaire non formatte, cree par SOS_GSF
C       Fichier cree dans le repertoire de lancement du programme
C       Si un fichier de meme nom existe deja, il est ecrase
C       Chaque enregistrement contient :
C           -numéro d'angle Teta1 (valeur entiere)
C           -numéro d'angle Teta2 (valeur entiere)
C           -ordre maximale utile de la decompositon : IL  (valeur entiere)
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
C        - SOS_CALCG
C
C==============================================================================
      SUBROUTINE SOS_GSF(LUM_NBMU,RMU,SIG,OS_NM,FIC_RES_GSF,IER)
 	
      IMPLICIT NONE

C* Definition des variables    
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Tableau des mu de Gauss

      DOUBLE PRECISION E(0:CTE_OS_NM_MAX) ! Resultat de la décomposition en 
      				              ! séries de Fourier de la fonction G.
	
      DOUBLE PRECISION SIG 	! Carre du paramètre Sigma de la théorie de Cox et Munck.

      DOUBLE PRECISION PHI	! Azimut relatif (radians).

      DOUBLE PRECISION G	! Valeur de la fonction G(phi).

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).
      DOUBLE PRECISION C12	! Cos(Teta1) * Cos(Teta2).
      DOUBLE PRECISION S12	! Sin(Teta1) * Sin(Teta2).
      DOUBLE PRECISION CS12	! 0.5 * [ Cos(Teta1) + Cos(Teta2) ]**2 .

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

      INTEGER*4 OS_NM ! Ordre maximal du développement en séries de Fourier  
			 ! de la fonction G des calculs de matrice de reflexion.
			 
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
           C12=C1*C2
           S12=S1*S2
           CS12=(C1+C2)
           CS12=.5*CS12*CS12

C*   Determination de G pour phi = 0.
C-------------------------------------
           PHI=0.0
           CALL SOS_CALCG(CS12,C12,S12,SIG,PHI,G)
           U(0)=G
           GMAX=G

C*   Determination de G pour phi = INCTE_PI.
C-----------------------------------------
           PHI=INCTE_PI
           CALL SOS_CALCG(CS12,C12,S12,SIG,PHI,G)
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
                 CALL SOS_CALCG(CS12,C12,S12,SIG,PHI,G)
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
              CALL SOS_CALCG(CS12,C12,S12,SIG,PHIB,G)
	
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
                 CALL SOS_CALCG(CS12,C12,S12,SIG,PHI,G)
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
      END         !FIN DE LA PROCEDURE SOS_GSF





C==============================================================================
C PROCEDURE: SOS_CALCG
C ==========
C      Cette procédure calcule :
C	 G(Teta1,Teta2,phi) = exp [-tan(TetaN)**2 / SIG] / cos(TetaN)**4
C
C	 avec  Teta1 : angle zénithal d'incidence  (0 < Teta1 < PI/2),
C             Teta2 : angle zénithal de reflexion (0 < Teta2 < PI/2),
C             phi   : angle azimutal relatif entre les directions
C                     d'incidence et de reflexion,
C
C	       TetaN : angle zénithal de la normale à la facette de vague,
C                     calcule à partir de Teta1, Teta2 et phi.
C			
C             SIG   : carre du paramètre Sigma de la théorie de Cox et Munck.
C
C			
C
C Description des paramètres
C --------------------------
C	PHI	(double)   (E) Azimut relatif (radians).
C	C12	(double)   (E) Cos(Teta1) * Cos(Teta2).
C	S12	(double)   (E) Sin(Teta1) * Sin(Teta2).
C	CS12	(double)   (E) 0.5 * [ Cos(Teta1) + Cos(Teta2) ]**2 .
C 	SIG 	(double)   (E) Carre du paramètre Sigma de la théorie
C    			        de Cox et Munck.
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
      SUBROUTINE SOS_CALCG(CS12,C12,S12,SIG,PHI,G)

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------------------------------------
      DOUBLE PRECISION PHI	! Azimut relatif (radians).

      DOUBLE PRECISION C12	! Cos(Teta1) * Cos(Teta2).
      DOUBLE PRECISION S12	! Sin(Teta1) * Sin(Teta2).
      DOUBLE PRECISION CS12	! 0.5 * [ Cos(Teta1) + Cos(Teta2) ]**2 .

      DOUBLE PRECISION SIG 	! Carre du paramètre Sigma de la théorie
      				! de Cox et Munck.

      DOUBLE PRECISION G	! Valeur de G(TetaN).	
      
      DOUBLE PRECISION COSTETAD	! Cosinus de l'angle de diffusion entre les
      				       ! directions d'incidence et de reflexion.
      DOUBLE PRECISION X	! 1 / cos(TetaN)**2.


C* Calcul
C---------------------------------
      COSTETAD=-C12+S12*DCOS(PHI)
      X=(1-COSTETAD)/CS12
      G=X*X*DEXP (-(X-1)/SIG)
      
      RETURN
      END         !FIN DE LA PROCEDURE SOS_CALCG








