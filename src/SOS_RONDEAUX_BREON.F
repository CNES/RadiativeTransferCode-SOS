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
C* FICHIER: SOS_RONDEAUX_BREON.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Genere les fichiers de BPDF du modele de Rondeaux et de Breon.
C*
C* AUTEUR: 
C*   CS SI, B.Lafrance.
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0: Creation du programme.
C*             1.1: Mise a jour des commentaires 
C*                  
C* MOD:VERSION:2.0: 10/01/2010 
C*
C*   - Suppression de la constante SOS_OS_FICANGLE (de SOS.h) et de la lecture du fichier correspondant 
C*     des angles de Gauss au profit : d'un passage par argument des tables d'angles (RMU) et poids (CHR).
C*
C*   - Modification de la taille des tableaux de luminance :
C*     ==> Passage de SOS_OS_NBMU a SOS_OS_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'a LUM_NBMU (nouveau argument d'entree) 
C*         egal au nombre effectif d'angles utilises.
C*
C*   - Suppresion des constantes SOS_OS_NB, SOS_OS_NS et SOS_OS_NM (de SOS.h)  
C*     au profit : d'un passage par argument des valeurs OS_NB, OS_NS, OS_NM.
C*
C*   - Modification de la taille des tableaux de decomposition en polynomes de Legendre :
C*     ==> Passage de SOS_OS_NM a SOS_OS_NM_MAX.
C*     ==> Boucle sur les ordre jusqu'a OS_NM (nouvel argument d'entree)
C*
C*   - Routine SOS_GSF_RONDEAUX_BREON : Modification de l'appel de la routine pour inclure l'argument 
C*     LUM_NBMU (nb d'angles utiles)  et pour limiter les boucles aux bornes utiles LUM_NBMU,
C*     et pour inclure l'argument OS_NM. 
C*
C*   - Adaptation de l'appel des routines SOS_MAT_FRESNEL, SOS_MAT_REFLEXION , 
C*     et SOS_MISE_FORMAT, aux passages des parametres LUM_NBMU (nb d'angles utiles),
C*     RMU (table des cosinus des angles) et CHR (table des poids des angles),
C*     ainsi que OS_NB, OS_NS et OS_NM.
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    SOS_LENDIR : Longueur des noms de repertoires.
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_OS_NBMU_MAX : Nombre maximal d'angles positifs.
C    SOS_PI : Valeur de PI
C    SOS_OS_NB_MAX :  Valeur maximale pour l'ordre limite du developpement en 
C                     polynomes de Legendre.
C    SOS_OS_NS_MAX :  Valeur maximale pour l'ordre limite du developpement en polynomes 
C                     de Legendre des elements de la matrice de Fresnel 
C                     et en series de Fourier des luminances.
C    SOS_OS_NM_MAX :  Valeur maximale pour l'ordre limite  du developpement 
C                     en series de Fourier de la fonction G.
C
C----------------------------------------------------------------------------
#include "SOS.h"

C==============================================================================
C PROCEDURE: SOS_RONDEAUX_BREON
C ==========      
C      Cette procedure genere les fichiers de reflexion de surface pour le 
C      modele de Rondeaux ou Breon.
C
C      La matrice de reflectance polarisee s'ecrit :
C          
C          Modele de Rondeaux :
C           => Rp(Teta1, Teta2, phi) =  G(Teta1, Teta2, phi) * Fp(Teta1, Teta2, phi) 
C                avec G(Teta1, Teta2, phi) = 1 / [ 4 * (cosTeta1+cosTeta2) ]
C
C          Modele de Breon :
C           => Rp(Teta1, Teta2, phi) =  G(Teta1, Teta2, phi) * Fp(Teta1, Teta2, phi) 
C                avec G(Teta1, Teta2, phi) = 1 / [ 4*cosTeta1*cosTeta2 ]
C
C	     avec Teta1 : angle zenithal d'incidence  (0 < Teta1 < PI/2),
C                 Teta2 : angle zenithal de reflexion (0 < Teta2 < PI/2),
C                 phi   : angle azimutal relatif entre les directions
C                         d'incidence et de reflexion,
C                 Fp(Teta1, Teta2, phi) : matrice de reflexion de Fresnel.
C
C
C       Le programme :
C          - decompose la fonction G en series de Fourier sur l'azimut 
C             => Stockage dans RES_GSF.
C          - decompose la matrice de Fresnel Fp en series de Fourier sur l'azimut 
C             => Stockage dans RES_FRESNEL.
C          - calcule la matrice G*Fp decomposee en series de Fourier sur l'azimut 
C            a partir de RES_GSF et RES_FRESNEL.
C             => Stockage dans RES_MAT_REFLEX.
C          - rearrange le contenu du fichier RES_MAT_FRESNEL pour avoir des 
C            enregistrements par ordre de la serie (compatible avec l'usage des OS).
C             => Generation du fichier de BRDF/BPDF resultat.
C
C
C        NB : Les routines utilisees sont les memes que celles pour le calcul du glitter
C             (decomposition de la matrice de Fresnel, calcul de la matrice de reflexion
C              de surface)
C             => Ceci impose d'introduire une vitesse de vent pour la routine SOS_MAT_REFLEXION.
C                Afin de compenser l'application d'une constante dans le calcul de la matrice
C                de reflexion on multiplie par l'inverse de cette constante les fonctions G
C                issues de la decomposition en series de Fourier pour les modeles de Rondeaux
C                ou Breon. La routine SOS_GSF_RONDEAUX_BREON requiet donc la meme vitesse
C                de vent que la routine SOS_MAT_REFLEXION.
C
C Donnees en entree du programme 
C ------------------------------
C    -->Les parametres/arguments du programme
C       LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C	RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (DOUBLE) : (E) Cosinus des angles.
C	CHR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (DOUBLE) : (E) Poids des angles.
C	IMODELE  (INTEGER*4)  : 1 <-> modele de Rondeaux, Breon sinon.      	
C       IND   (DOUBLE)  : Indice de refraction de l'eau / air.
C       OS_NB (I4) : (E) Ordre des developpements des fonctions de phase 
C		         en Polynomes de Legendre a utiliser.
C       OS_NS (I4) : (E) Ordre maximal a utiliser pour le developpement en
C                        polynomes de Legendre des elements de la matrice de Fresnel 
C                      et en series de Fourier des calculs de luminance.
C       OS_NM (I4) : (E) Ordre maximal du developpement en series de Fourier  
C                        de la fonction G des calculs de matrice de reflexion.
C       FICBPDF  (CHARACTER*SOS_LENFIC2) : nom complet du fichier Resultat.
C       TRACE  (LOGICAL) : vrai,si ecriture dans le fichier trace.
C       IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon
C
C
C Resultat fourni par le programme
C --------------------------------
C   Fichier BRDF / BPDF
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
C   -->Contenu du fichier Resultat
C       Fichier binaire non formatte, cree par SOS_MISE_FORMAT
C       Il contient la matrice de reflexion de surface 
C       (fichier GLITTER ou autres BRDF/BPDF):
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
C   -->Contenu du fichier temporaire RES_GSF
C       Fichier binaire non formatte, cree par SOS_GSF_RONDEAUX_BREON
C       Fichier cree dans le repertoire de lancement du programme
C       Si un fichier de meme nom existe deja, il est ecrase
C       Chaque enregistrement contient :
C           -numero d'angle Teta1 (valeur entiere)
C           -numero d'angle Teta2 (valeur entiere)
C           -ordre maximale utile de la decompositon : IL  (valeur entiere)
C           -Fonction G decomposee pour (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 valeurs double precision)
C
C   -->Contenu du fichier temporaire RES_FRESNEL
C      Fichier Ascii, format d'un enregistrement 4E15.8, 
C      cree par SOS_MAT_FRESNEL
C      Fichier cree dans le repertoire de lancement du programme
C      Si un fichier de meme nom existe deja, il est ecrase
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
C    Les 3 fichiers temporaires sont ecrases s'ils existent puis 
C    detruits en fin de programme
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
C      Erreur si le fichier Resultat existe deja
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, arret du programme et 
C     retour du status 1 au shell
C        
C   
C
C  Sous programmes utilises:
C --------------------------
C  Cette procedure fait appel aux routines:
C      - SOS_GSF_RONDEAUX_BREON
C      - SOS_MAT_FRESNEL 
C      - SOS_MAT_REFLEXION 
C      - SOS_MISE_FORMAT
C
C
C=============================================================================

      SUBROUTINE SOS_RONDEAUX_BREON(LUM_NBMU,RMU,CHR,IMODELE,IND,
     &                              OS_NB,OS_NS,OS_NM, 
     &                              FICBPDF,TRACE,IER)

      IMPLICIT NONE

C* Definition des variables   
C*--------------------------------------------------------------------------
      CHARACTER*SOS_LENFIC2 FICBPDF	! Nom complet du fichier resultat.
            
      DOUBLE PRECISION WIND	! Vitesse de vent (artefact pour conserver les
      				!                  routines pre existantes).
      DOUBLE PRECISION IND	! Indice de refraction surface / air
	
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   ! Tableau des mu de Gauss
      DOUBLE PRECISION CHR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   ! Poids de Gauss.
      
      LOGICAL TRACE     ! = vrai,si ecriture dans le fichier trace

      INTEGER*4 OS_NB 	! Ordre des developpements des fonctions de phase 
			! en Polynomes de Legendre.
      INTEGER*4 OS_NS   ! Ordre maximal du developpement en polynomes de Legendre   
			! des elements de la matrice de Fresnel.
      INTEGER*4 OS_NM   ! Ordre maximal du developpement en series de Fourier  
			! de la fonction G des calculs de matrice de reflexion.
			
      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations
      INTEGER*4 IMODELE  ! 1 si simulation du modele de Rondeaux.
      		         ! Sinon simulation du modele de Breon.
      INTEGER*4 I	 ! Indice de tableau.
      INTEGER*4 IFIN     ! nb caracteres significatifs d'une chaine
      INTEGER*4 IER      ! code d'erreur =0 si pas d'erreur, =-1 sinon


	

C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0

C* Initialisation de la vitesse du vent
C------------------------------------------------------
      WIND = 0.D+00
      
      
C* Ouverture du fichier Trace 
C------------------------------------------------------
      IF (TRACE) THEN
         WRITE(99,*,err=921)
	 WRITE(99,*,err=921) '*****************************************'
         IF(IMODELE.EQ.1) THEN 
	    WRITE(99,*,err=921) 'Calcul du fichier RONDEAUX'
	 ELSE
	    WRITE(99,*,err=921) 'Calcul du fichier BREON'
	 ENDIF 
	 WRITE(99,*,err=921) FICBPDF     
         WRITE(99,*,err=921) ' '
         WRITE(99,*,err=921)'PARAMETRES DE SIMULATION '
         WRITE(99,*,err=921)
     &        ' Indice de refraction surface / air ',IND
	 WRITE(99,*,err=921) ' ' 	 	 
      END IF
        


C* Decomposition en series de Fourier de la fonction G(Teta1,Teta2,phi).
C*    => Creation du fichier RES_GSF.
C---------------------------------------------------------------------
      CALL SOS_GSF_RONDEAUX_BREON(LUM_NBMU,RMU,IMODELE,WIND,OS_NM,
     &                            TRACE,IER)
      IF (IER.NE.0) GOTO 995   ! cas d'erreur
      
C* Developpement en fonctions de Legendre de la matrice de Fresnel.
C*    => Creation du fichier RES_FRESNEL.
C------------------------------------------------------------------
      CALL SOS_MAT_FRESNEL(LUM_NBMU,RMU,CHR,IND,OS_NS,TRACE,IER)
      IF (IER.NE.0) GOTO 996   ! cas d'erreur
      
C* Calcul de la matrice de reflexion developpee en series de Fourier.
C*    => Creation du fichier RES_MAT_FRESNEL 
C*       a partir de RES_GSF et RES_FRESNEL.
C------------------------------------------------------------------
      CALL SOS_MAT_REFLEXION(WIND,LUM_NBMU,RMU,OS_NB,OS_NS,OS_NM,IER)
      IF (IER.NE.0) GOTO 997   ! cas d'erreur
      
C* Rearrangement des resultats pour etre compatible avec le code des OS.
C*  => resultat dans un fichier temporaire n'incluant que la BPDF.
C-----------------------------------------------------------------------  
      CALL SOS_MISE_FORMAT(FICBPDF,LUM_NBMU,OS_NB,IER)
      IF (IER.NE.0) GOTO 998   ! cas d'erreur
       
	    
C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  912 WRITE(6,*) 'Erreur a l''ouverture du fichier angles de Gauss'
      IER=-1
      GOTO 9999
  921 WRITE(6,*) 'Erreur d''ecriture sur le fichier trace SURFACE'
      IER=-1
      GOTO 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1
      GOTO 9999      
  995 WRITE(6,*) 'Erreur dans la routine SOS_GSF_RONDEAUX_BREON'
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
      OPEN(UNIT=3,FILE='RES_GSF',FORM='UNFORMATTED') 
      CLOSE(3,STATUS='DELETE')     
      OPEN(UNIT=3,FILE='RES_FRESNEL') 
      CLOSE(3,STATUS='DELETE')       
      OPEN(UNIT=3,FILE='RES_MAT_REFLEX',FORM='UNFORMATTED')
      CLOSE(3,STATUS='DELETE') 
                                                            
      RETURN
      END    !FIN DE LA PROCEDURE SOS_RONDEAUX_BREON





C==============================================================================
C PROCEDURE: SOS_GSF_RONDEAUX_BREON
C ==========
C      Cette procedure realise la decomposition en series de Fourier
C      (sur phi) de la fonction G(Teta1,Teta2,phi) :
C
C         Modele de Rondeaux :G(Teta1, Teta2, phi) = SIG / [ 4 * (cosTeta1+cosTeta2) ]
C
C         Modele de Breon :   G(Teta1, Teta2, phi) = SIG / [ 4*cosTeta1*cosTeta2 ]
C
C	    avec Teta1 : angle zenithal d'incidence  (0 < Teta1 < PI/2),
C                Teta2 : angle zenithal de reflexion (0 < Teta2 < PI/2),
C                phi   : angle azimutal relatif entre les directions
C                        d'incidence et de reflexion,
C			
C                SIG   : carre du parametre Sigma de la theorie de Cox et Munck.
C
C
C      Enregistre les resultats dans le fichier "RES_GSF" sous la forme:
C	  Numero d'angle Teta1 : I1,
C         Numero d'angle Teta2 : I2,
C         Ordre maximal utile de la decomposition: IL
C            (IL = 0 pour ces fonctions independantes  de l'azimut),
C         Fonction G decomposee pour (Teta1,Teta2) : (E(J),J=0,IL)
C
C Description des parametres
C --------------------------
C      LUM_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles.
C      RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)	 (E) : Tableau des mu de Gauss.
C      IRONDEAUX (integer*4) (E) : Indice de modele 1 <-> Rondeaux, Breon sinon.
C      WIND	(double)  (E)  : Vitesse du vent (m/s).
C      OS_NM (I4) : (E) Ordre maximal du developpement en series de Fourier  
C                       de la fonction G des calculs de matrice de reflexion.
C      TRACE (logical) (E) : Option de Trace.
C      IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon.
C
C Autre resultat fourni par la procedure
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
C           -numero d'angle Teta1 (valeur entiere)
C           -numero d'angle Teta2 (valeur entiere)
C           -ordre maximale utile de la decompositon : IL  (valeur entiere)
C           -Fonction G decomposee pour (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 valeurs double precision)
C
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
      SUBROUTINE SOS_GSF_RONDEAUX_BREON(LUM_NBMU,RMU,IRONDEAUX,WIND,
     &                                  OS_NM,TRACE,IER)
 	
      IMPLICIT NONE

C* Definition des variables    
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! Tableau des mu de Gauss

      DOUBLE PRECISION WIND	! Vitesse du vent (m/s).

      DOUBLE PRECISION E(0:SOS_OS_NM_MAX)	! Resultat de la decomposition en 
      				        ! series de Fourier de la fonction G.
	
      DOUBLE PRECISION SIG 	! Carre du parametre Sigma de la theorie 
                                ! de Cox et Munck.

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).

      LOGICAL TRACE     ! = vrai,si ecriture dans le fichier trace.

      INTEGER*4 OS_NM   ! Ordre maximal du developpement en series de Fourier  
			! de la fonction G des calculs de matrice de reflexion.      
      INTEGER*4 LUM_NBMU  ! Nombre d'angles effectifs des simulations
      INTEGER*4 IRONDEAUX ! 1 si simulation du modele de Rondeaux.
      		          ! Sinon simulation du modele de Breon.
      INTEGER*4 IL	! Ordre maximale utile pour la fonction E.
      INTEGER*4 I1	! Indice de numero d'angle Teta1.
      INTEGER*4 I2	! Indice de numero d'angle Teta2.
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon  
      INTEGER*4 J	! Indice de tableau.
      

C* Carre du parametre Sigma de la theorie de Cox et Munck.
C--------------------------------------------------------
      SIG = .003 + .00512*WIND

C* Ouverture du fichier resultat
C--------------------------------
      OPEN(UNIT=30,FILE='RES_GSF',FORM='UNFORMATTED',ERR=991)

C* Fichier Trace
C--------------------------------
      IF (TRACE) THEN
          IF(IRONDEAUX.EQ.1) THEN
	     WRITE(99,*,err=921) 'Valeurs de la fonction G Rondeaux'       
	  ELSE
	     WRITE(99,*,err=921) 'Valeurs de la fonction G Breon'       
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
          

C*   Fonction G a l'ordre IS = 0 pour le modele de Rondeaux.
C*   => Independant de l'azimut : Terme unique IS = 0.
C---------------------------------
           IF(IRONDEAUX.EQ.1) THEN
	      E(0) = SIG / (1./C1 + 1./C2)
	   ELSE 
	      E(0) = SIG
	   ENDIF   
	   IL=0

C* Ecriture des resultats dans le fichier Trace
C--------------------------------
	   IF (TRACE) THEN
	       WRITE(99,1000,err=921) ACOS(C1)*180./SOS_PI, 
     &	                              ACOS(C2)*180./SOS_PI,
     &                                E(0)         	 
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
  921 WRITE(6,*) 'Erreur d''ecriture sur le fichier trace Surface'
      IER=-1
      GOTO 9999
  991 WRITE(6,*) 'SOS_GSF_RONDEAUX_BREON: ',
     &           'Erreur a l''ouverture d''un fichier'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'SOS_GSF_RONDEAUX_BREON: ',
     &		 'Erreur de lecture ou d''ecriture sur un fichier'
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








