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
C* FICHIER: SOS_NADAL.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Genere les fichiers de BPDF selon le modele de Nadal.
C*
C* AUTEUR: 
C*    CS SI, B.Lafrance.
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0: Creation du programme
C* MOD:VERSION:1.1: Mise a jour des commentaires.
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
C*   - Routine SOS_F21SF_NADAL : Modification de l'appel de la routine pour inclure l'argument LUM_NBMU 
C*     (nb d'angles utiles) et pour limiter les boucles aux bornes utiles LUM_NBMU 
C*     et pour inclure l'argument OS_NB. 
C*
C*   - Adaptation de l'appel des routines SOS_MAT_FRESNEL, SOS_MAT_REFLEXION , 
C*     et SOS_MISE_FORMAT, aux passages des parametres LUM_NBMU (nb d'angles utiles),
C*     RMU (table des cosinus des angles) et CHR (table des poids des angles) ,
C*     ainsi que OS_NB, OS_NS et OS_NM.
C*
C*
C* MOD:VERSION:2.1: 08/12/2015 
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*      (incluant au décodage des constante) : requis pour compilation gfortran
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
C    SOS_PH_NU : Nombre d'angles azimutaux (2**NQ)
C    SOS_SEUIL_SF_NADAL : Valeur seuil pour l'arret de la decomposition 
C                         en series de Fourier de la fonction 
C                        (F21 Nadal) / (F21 Fresnel).
C
C----------------------------------------------------------------------------
#include "SOS.h"

C==============================================================================
C PROCEDURE: SOS_NADAL
C ==========      
C      Cette procedure genere les fichiers de BPDF du modele semi-empirique 
C      de Nadal.
C
C
C      La matrice de reflectance polarisee du modele de Nadal s'ecrit :
C          
C           Rp(Teta1, Teta2, phi) =  F(Teta1, Teta2, phi) * Fp(Teta1, Teta2, phi) 
C               avec F(Teta1, Teta2, phi) = 1 / [ 4 * (cosTeta1+cosTeta2) ]
C
C	     avec Teta1 : angle zenithal d'incidence  (0 < Teta1 < PI/2),
C                 Teta2 : angle zenithal de reflexion (0 < Teta2 < PI/2),
C                 phi   : angle azimutal relatif entre les directions
C                         d'incidence et de reflexion,
C                 Fp(Teta1, Teta2, phi) : matrice de reflexion de Fresnel.
C
C    	          F(Teta1,Teta2,phi) = F12nadal(Teta1,Teta2,phi) 
C                                     / F12fresnel(Teta1,Teta2,phi)
C	           
C                    F12fresnel = 0.5 * (RR*RR - RL*RL) 
C                    selon la convention de signe du modele de Nadal.
C
C                       avec RL le coefficient de reflexion parallele,
C                       calcules en fonction de l'indice de refraction de la surface
C                       et de l'angle de reflexion teta :
C                       RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                                / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C                       RR le coefficient de reflexion perpendiculaire :
C                       RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                                / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C                
C
C                    F12nadal = alpha * 
C                         [ 1 -exp(-beta*F12fresnel / (cosTeta1 + cosTeta2) ) ]
C
C                         avec alpha et beta les parametres du modele de Nadal.
C
C
C       Le programme :
C          - decompose la fonction F en series de Fourier sur l'azimut 
C             => Stockage dans RES_GSF.
C          - decompose la matrice de Fresnel Fp en series de Fourier sur l'azimut 
C             => Stockage dans RES_FRESNEL.
C          - calcule la matrice F*Fp decomposee en series de Fourier sur l'azimut 
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
C                de reflexion on multiplie par l'inverse de cette constante les fonctions F
C                issues de la decomposition en series de Fourier pour les modeles de Rondeaux
C                ou Breon. La routine SOS_GSF_RONDEAUX_BREON requiet donc la meme vitesse
C                de vent que la routine SOS_MAT_REFLEXION.
C
C
C Donnees en entree du programme 
C ------------------------------
C    -->Les parametres/arguments du programme
C       LUM_NBMU (I4) (E) : Nombre d'angles (positifs) effectivement utiles
C	RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (DOUBLE) (E) : Cosinus des angles.
C	CHR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (DOUBLE) (E) : Poids des angles.
C 	IND     (DOUBLE)  (E) : Indice de refraction de la surface.      	
C 	ALPHA   (DOUBLE)  (E) : Parametre alpha du modele de Nadal.         	
C 	BETA    (DOUBLE)  (E) : Parametre beta du modele de Nadal.          	
C       FICBPDF (CHARACTER*SOS_LENFIC2) (E) : nom complet du fichier Resultat 
C                                    (arborescence complete)
C       OS_NB (I4) : (E) Ordre des developpements des fonctions de phase 
C		         en Polynomes de Legendre a utiliser.
C       OS_NS (I4) : (E) Ordre maximal a utiliser pour le developpement en
C                        polynomes de Legendre des elements de la matrice de Fresnel 
C                      et en series de Fourier des calculs de luminance.
C       OS_NM (I4) : (E) Ordre maximal du developpement en series de Fourier  
C                        de la fonction G des calculs de matrice de reflexion.
C       TRACE  (LOGICAL) (E) : vrai,si ecriture dans le fichier trace.
C       IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon
C
C
C Resultat fourni par le programme
C --------------------------------
C   Fichier BPDF de Nadal (sans intensite).
C   Fichier de trace (optionnel), 
C               si un fichier de meme nom existe deja, il est complete
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
C       Fichier binaire non formatte, cree par SOS_F21SF_NADAL
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
C      Nombre d'enregistrements : de k=0 a SOS_OS_NS
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
C  Sous programmes utilises:
C --------------------------
C  Cette procedure fait appel aux routines:
C      - SOS_F21SF_NADAL
C      - SOS_MAT_FRESNEL 
C      - SOS_MAT_REFLEXION 
C      - SOS_MISE_FORMAT
C
C=============================================================================

      SUBROUTINE SOS_NADAL(LUM_NBMU,RMU,CHR,IND,ALPHA,BETA,
     &                     OS_NB,OS_NS,OS_NM, FICBPDF,TRACE,IER)

      IMPLICIT NONE

C* Definition des variables   
C*--------------------------------------------------------------------------
      CHARACTER*SOS_LENFIC2 FICBPDF	! Nom complet du fichier resultat. 
            
      DOUBLE PRECISION IND	! Indice de refraction surface / air
      DOUBLE PRECISION ALPHA	! Parametre alpha du modele de Nadal.
      DOUBLE PRECISION BETA	! Parametre beta du modele de Nadal.
      	
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   ! Tableau des mu de Gauss
      DOUBLE PRECISION CHR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   ! Poids de Gauss.

      DOUBLE PRECISION WIND	! Vitesse de vent (artefact pour conserver les
      				!                  routines pre existantes).
								      
      LOGICAL TRACE     ! = vrai,si ecriture dans le fichier trace

      INTEGER*4 OS_NB 	! Ordre des developpements des fonctions de phase 
			! en Polynomes de Legendre.
      INTEGER*4 OS_NS   ! Ordre maximal du developpement en polynomes de Legendre   
			! des elements de la matrice de Fresnel.
      INTEGER*4 OS_NM   ! Ordre maximal du developpement en series de Fourier  
			! de la fonction G des calculs de matrice de reflexion.

      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations
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
         WRITE(99,*,err=921) 'Calcul du fichier BPDF de NADAL'	    
         WRITE(99,*,err=921) FICBPDF 
	 WRITE(99,*,err=921) ' '
         WRITE(99,*,err=921)'PARAMETRES DE SIMULATION '
         WRITE(99,*,err=921)
     &        ' Indice de refraction surface / air ',IND
         WRITE(99,*,err=921) ' Alpha =  ',ALPHA
         WRITE(99,*,err=921) ' Beta  =  ',BETA     
	 WRITE(99,*,err=921) ' ' 	 	 
      END IF
        

C* Decomposition en series de Fourier de la fonction F21 de Nadal.
C*    => Creation du fichier RES_GSF.
C---------------------------------------------------------------------
      CALL SOS_F21SF_NADAL(LUM_NBMU,RMU,IND,ALPHA,BETA,WIND,OS_NB,
     &                     TRACE,IER)
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
      goto 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  912 WRITE(6,*) 'Erreur a l''ouverture du fichier angles de Gauss'
      IER=-1
      goto 9999
  921 WRITE(6,*) 'Erreur d''ecriture sur le fichier trace SURFACE'
      IER=-1
      goto 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1
      goto 9999      
  995 WRITE(6,*) 'Erreur dans la routine SOS_F21SF_NADAL'
      IER=-1
      goto 9999
  996 WRITE(6,*) 'Erreur dans la routine SOS_MAT_FRESNEL'
      IER=-1
      GOTO 9999
  997 WRITE(6,*) 'Erreur dans la routine SOS_MAT_REFLEXION'
      IER=-1
      GOTO 9999
  998 WRITE(6,*) 'Erreur dans la routine SOS_MISE_FORMAT'
      IER=-1
      goto 9999              
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
      END    !FIN DE LA PROCEDURE SOS_NADAL






C==============================================================================
C PROCEDURE: SOS_F21SF_NADAL
C ==========
C      Cette procedure realise la decomposition en series de Fourier
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
C             selon la convention de signe du modele de Nadal.
C
C                avec RL le coefficient de reflexion parallele,
C                calcules en fonction de l'indice de refraction de la surface
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
C                 avec alpha et beta les parametres du modele de Nadal.
C
C      La BPDF du modele de Nadal est convertie de reflectance polarisee 
C      en une grandeur utile pour les OS : luminance normalisee * cos(Teta2),
C      c'est a dire : BPDF * cos(Teta1) * cos(Teta2).			
C
C
C      La procedure enregistre les resultats dans le fichier "RES_MAT_REFLEX" 
C      sous la forme:
C         Fonction F decomposee pour (Teta1,Teta2) : (E(IS),IS=0,OS_NB)
C
C Description des parametres
C --------------------------
C      LUM_NBMU (I4) (E) : Nombre d'angles (positifs) effectivement utiles.
C      RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)	 (E) Tableau des mu de Gauss.
C      IND     (double)  (E) : Indice de refraction de la surface.      	
C      ALPHA   (double)  (E) : Parametre alpha du modele de Nadal.         	
C      BETA    (double)  (E) : Parametre beta du modele de Nadal. 
C      WIND    (double)  (E) : Vitesse de vent (artefact pour conserver les 
C                                               routines existantes).    
C      OS_NB (I4) (E) : Ordre des developpements des fonctions de phase 
C		        en Polynomes de Legendre a utiliser.     	
C      TRACE    (logical) (E) = vrai,si ecriture dans le fichier trace
C      IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C Autre resultat fourni par la procedure
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
C       decomposition en series de Fourier de fonction F21 du modele 
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
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C     
C  Routines utilises par la procedure:
C ------------------------------------
C  Ce programme fait appel aux routines:
C        - SOS_CALC_F21_NADAL_SUR_FRESNEL
C
C==============================================================================
      SUBROUTINE SOS_F21SF_NADAL(LUM_NBMU,RMU,IND,ALPHA,BETA,WIND,
     &                           OS_NB,TRACE,IER)
 	
      IMPLICIT NONE

C* Definition des variables    
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! Tableau des mu de Gauss

      DOUBLE PRECISION IND	! Indice de refraction surface / air
      DOUBLE PRECISION ALPHA	!Parametre alpha du modele de Nadal.
      DOUBLE PRECISION BETA	!Parametre beta du modele de Nadal.
      DOUBLE PRECISION WIND	! Vitesse de vent (artefact pour conserver les
      				!                  routines pre existantes).
      DOUBLE PRECISION SIG 	! Carre du parametre Sigma de la theorie 
                                ! de Cox et Munck.	
      
      DOUBLE PRECISION E(0:SOS_OS_NB_MAX) ! Resultat de la decomposition en 
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

      DOUBLE PRECISION U(0:SOS_PH_NU)  ! Fonction F(phi) entre 0 et PI
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
     
	
      DOUBLE PRECISION FRECOMB	! Valeur de la fonction F recombinee a partir de sa
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

      INTEGER*4 OS_NB 	! Ordre des developpements des fonctions de phase 
			! en Polynomes de Legendre.
				
      INTEGER*4 LUM_NBMU ! Nombre d'angles effectifs des simulations
      	    
      INTEGER*4 J	! Indice de tableau.
      INTEGER*4 IL	! Ordre maximale utile pour la fonction E.
      INTEGER*4 I1	! Indice de numero d'angle Teta1.
      INTEGER*4 I2	! Indice de numero d'angle Teta2.
      INTEGER*4 IS	! Ordre IS de la decomposition en series de Fourier.
      INTEGER*4 IS2	! Ordre IS de la decomposition en series de Fourier.
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon  
      INTEGER*4 I
      INTEGER*4 COMPT	! Nombre de configurations pour la sommation des ecarts relatifs.

      
	
C* Initialisation
C---------------------------------------------------------	
      ECARTMAX = -1.D+300
      SOMECART = 0.
      COMPT = 0
      
C* Carre du parametre Sigma de la theorie de Cox et Munck.
C--------------------------------------------------------
      SIG = .003 + .00512*WIND
      
C* Creation du fichier resultat de SOS_MAT_REFLEXION
C---------------------------------------------------------
      OPEN(UNIT=30,FILE='RES_GSF',FORM='UNFORMATTED',ERR=991)

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

C*   Determination de F pour phi = 0 a PI par pas PI/SOS_PH_NU
C-------------------------------------------------------------
          
           Q=SOS_PI/SOS_PH_NU
           DO 6 I=0,SOS_PH_NU
	
              PHI=Q*I
              CALL SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,SIG,
     &	                                          C1,S1,C2,S2,PHI,F)   
              U(I)=F
	      
 6         CONTINUE

		
C*   Determination des valeurs min et max de F
C*   pour I1 et I2 fixes mais PHI variant.
C----------------------------------------------
	   IF (TRACE) THEN
	       FMAX=-1.D+300
	       FMIN=1.D+300
	       DO I=0,SOS_PH_NU
	          IF(U(I).GT.FMAX) THEN
                     FMAX = U(I)
		     PHIFMAX = Q*I		!Azimuth convention OS
                  ENDIF   
	          IF(U(I).LT.FMIN) THEN
		     FMIN = U(I)
		     PHIFMIN = Q*I		!Azimuth convention OS
	          ENDIF   	 
	       ENDDO

	       WRITE(99,1000,err=921) ACOS(C1)*180./SOS_PI, 
     &	                              ACOS(C2)*180./SOS_PI,
     &                                PHIFMIN*180./SOS_PI,
     &				      FMIN,
     &				      PHIFMAX*180./SOS_PI,
     &				      FMAX	          	 
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
	      DO I=0,SOS_PH_NU
	         PHI = I*Q		 		!Azimuth convention OS
	         Y=Y+U(I)*DCOS(IS*PHI)
	      ENDDO
	      E(IS) = Y*Q / SOS_PI
	      
C--------- Determination de l'ecart maximal entre F exact et recombine : B1
C---------------------------------------------------------------------      
	      B1=0.D+00			!initialisation
	      DO I=0,SOS_PH_NU
	         PHI = I*Q		 		!Azimuth convention OS
	         CALL SOS_CALC_F21_NADAL_SUR_FRESNEL(
     &	                   IND,ALPHA,BETA,SIG,C1,S1,C2,S2,PHI,F)

C-------------- Calcul de F recombine a partir de la fonction E(IS) pour PHI fixe. 
                 T1=E(0)	
		 DO IS2 = 1, IS		  
	              T1=T1+2.*E(IS2)*DCOS(IS2*PHI)	!Azimuth convention OS          		     
		 ENDDO !Fin boucle sur IS2
       
                 B1=MAX( DABS((T1-F)/F) , B1)
		 		 
	      ENDDO !Fin boucle sur I (phi = 0 a PI)
	           	     

C--------- Test d'arret de la decomposition :
C---------    * Arret si l'ecart relatif maximal entre F exacte et recombinee est
C---------      inferieur au seuil SOS_SEUIL_SF_NADAL
C---------      ou s'il augmente par rapport a l'ordre precedent.
C---------    * Poursuite de la boucle sur IS sinon.
C-------------------------------------------------------------------------------- 		
              IF (B1.GT.SOS_SEUIL_SF_NADAL) GOTO 14
              IL=IS	      
              GOTO 17	!Sortie de boucle sur IS
	      
 14           IF(B1.LT.B1_PREC) GOTO 141
              IL=IS-1
	      GOTO 17	!Sortie de boucle sur IS
	      
 141	      B1_PREC=B1	        
 13        CONTINUE	!Fin de boucle sur IS     
	   IL=OS_NB



C* Enregistrement pour le couple (Teta1, Teta2) de la decomposition
C* en series de Fourier de BPDF (convertie en luminance normalisee
C* et adaptee au code des OS).
C-----------------------------------------------------------------
 17        WRITE(30,err=992)I1,I2,IL,(E(J),J=0,IL)


C* Fichier Trace : Calculs de parametres d'analyse de qualite
C*                 de la decomposition en series de Fourier
C*                 et enregistrements.
C-----------------------------------------------------------------
           IF (TRACE) THEN

	      ECARTMAX_C1C2=0.
	            	      
	      Q=SOS_PI/SOS_PH_NU
	      DO I=0,SOS_PH_NU
	         PHI=I*Q
	         CALL SOS_CALC_F21_NADAL_SUR_FRESNEL(
     &	                    IND,ALPHA,BETA,SIG,C1,S1,C2,S2,PHI,F)   
		 
                 FRECOMB = E(0)
	         DO IS=1,IL
	            FRECOMB = FRECOMB + 2.*DCOS(IS*PHI)*E(IS)
	         ENDDO
            
	         IF(F.NE.0.D+00) THEN
	            ECART = DABS((FRECOMB-F)/F)  
		    
C*       Analyse pour tous les azimuts du couple courant (Teta1, Teta2)	
		    IF(ECART.GT.ECARTMAX_C1C2) THEN
		       ECARTMAX_C1C2 = ECART
		       PHIECARTMAX_C1C2 = PHI*180./SOS_PI	!Azimuth convention OS
	               F_ECARTMAX_C1C2 = F
	               FRECOMB_ECARTMAX_C1C2 = FRECOMB
		    ENDIF 
		    
C*       Analyse pour tous les azimuts et tous les couples (Teta1, Teta2)		    
	            IF(ECART.GT.ECARTMAX) THEN
	               ECARTMAX = ECART
	               PHIECARTMAX = PHI*180./SOS_PI	!Azimuth convention OS
	               TVECARTMAX = ACOS(C2)*180./SOS_PI
	               TSECARTMAX = ACOS(C1)*180./SOS_PI
	               F_ECARTMAX = F
	               FRECOMB_ECARTMAX = FRECOMB
	            ENDIF   
	    
	            SOMECART=SOMECART + ECART
	            COMPT = COMPT + 1
		 ENDIF
	      ENDDO !Fin boucle sur I (phi)
         	

	      IF(F_ECARTMAX_C1C2.NE.0.D+00) THEN	 
	         WRITE(99,2000,err=921) IL,PHIECARTMAX_C1C2,
     &				        F_ECARTMAX_C1C2,
     &		                        FRECOMB_ECARTMAX_C1C2,
     &					ECARTMAX_C1C2*100.
	      ELSE
	      	 WRITE(99,2500,err=921) IL,PHIECARTMAX_C1C2,
     &				        F_ECARTMAX_C1C2,
     &		                        FRECOMB_ECARTMAX_C1C2
	      ENDIF	 
	      WRITE(99,*,err=921)   	       
                 
	      
           END IF !Fin Trace
      
      	
		                            
 4       CONTINUE	!Fin de boucle sur Teta2 (I2)
 
 3    CONTINUE	!Fin de boucle sur Teta1 (I1)
 
 


C* Fichier Trace : Enregistrements des parametres d'analyse de qualite
C*                 de la decomposition en series de Fourier.
C----------------------------------------------------------------- 
      IF (TRACE) THEN
         WRITE(99,*,err=921)
	 WRITE(99,*,err=921)
         WRITE(99,*,err=921) 'Analyse de la decomposition en series'// 
     &	                     ' de Fourier de la BPDF NADAL :' 
	
	 IF(COMPT.GT.1) THEN
	  WRITE(99,*,err=921)
	  WRITE(99,*,err=921)
	  WRITE(99,*,err=921) 'Ecart relatif moy. pour la recombinaison'
     &	                      // ' (% hors valeurs nulles):',
     &	                       SOMECART*100/COMPT
	  WRITE(99,*,err=921)
	  WRITE(99,*,err=921) 'Ecart relatif max. pour la recomb. (%):',
     &	                       ECARTMAX*100
          WRITE(99,*,err=921) 'Angle delta phi associe'//
     &	                      ' (convention OS):',PHIECARTMAX
          WRITE(99,*,err=921) 'Angle de visee associe :',TVECARTMAX
	  WRITE(99,*,err=921) 'Angle solaire associe :',TSECARTMAX
	  WRITE(99,*,err=921) 'Valeur de la BPDF associee:',F_ECARTMAX
	  WRITE(99,*,err=921) 'Valeur de la BPDF recomb. associee:',
     &	                      FRECOMB_ECARTMAX 
         ELSE
	  WRITE(99,*,err=921)
	  WRITE(99,*,err=921) 'ATTENTION : BPDF nulle !'
	  WRITE(99,*,err=921) 'Parametres de simulation a controler.'
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
  921 WRITE(6,*) 'Erreur d''ecriture sur le fichier trace Surface'
      IER=-1
      goto 9999
  991 WRITE(6,*) 'SOS_F21SF_NADAL : Erreur a l''ouverture d''un fichier'
      IER=-1
      goto 9999
  992 WRITE(6,*) 'SOS_F21SF_NADAL',
     &		 'Erreur de lecture ou d''ecriture sur un fichier'
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
      
2000  FORMAT(5X,13hOrdre Limite:,I3,3h ; ,15h Erreur max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4,2X,
     &12h, Ecart (%):,1X,F7.3)         
     
2500  FORMAT(5X,13hOrdre Limite:,I3,3h ; ,15h Erreur max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4)         
                           
      RETURN
      END         !FIN DE LA PROCEDURE SOS_F21SF_NADAL







C==============================================================================
C PROCEDURE: SOS_CALC_F21_NADAL_SUR_FRESNEL
C		
C      Cette procedure calcule la fonction G(Teta1, Teta2, phi) pour le
C      modele de Nadal :
C
C      G(Teta1, Teta2, phi) 
C        =  F21Nadal(Teta1, Teta2, phi) / F21Fresnel(Teta1, Teta2, phi)
C
C        avec F21Nadal(Teta1, Teta2, phi)
C             = Alpha * [ 1 -exp(-Beta*F21Fresnel / ( cosTeta1+cosTeta2) ) ]
C
C        et F21Fresnel(Teta1, Teta2, phi) = 0.5 * (RR*RR - RL*RL) 
C        selon la convention de signe du modele de Nadal.
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
C            calcules en fonction de l'indice de refraction de la surface
C            de l'angle de reflexion teta.
C
C
C
C      La BPDF du modele de Nadal est convertie de reflectance polarisee 
C      en une grandeur utile pour les OS : luminance normalisee * cos(Teta2),
C      c'est a dire : BPDF * cos(Teta1) * cos(Teta2).			
C
C Description des parametres
C --------------------------
C       IND	(double)   (E) Indice de refraction surface / air.
C       ALPHA	(double)   (E) Parametre alpha du modele de Nadal.
C       BETA	(double)   (E) Parametre beta du modele de Nadal.
C       SIG	(double)   (E) Carre du parametre Sigma de Cox et Munck.	
C	C1	(double)   (E) Cos(Teta1).
C	S1	(double)   (E) Sin(Teta1).
C	C2	(double)   (E) Cos(Teta2).
C	S2	(double)   (E) Sin(Teta2).
C       PHI     (double)   (E) Azimut relatif pour la convention OS (radians).
C 	F	(double)   (S) Valeur de F(Teta1,Teta2,phi).
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
      SUBROUTINE SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,SIG,
     &	                                        C1,S1,C2,S2,PHI,F)   

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------------------------------------
      DOUBLE PRECISION IND	! (E)  Indice de refraction surface / air
      DOUBLE PRECISION ALPHA	! (E)  Parametre alpha du modele de Nadal.
      DOUBLE PRECISION BETA	! (E)  Parametre beta du modele de Nadal.
      DOUBLE PRECISION SIG 	! (E)  Carre du parametre Sigma de la theorie 
      				!      de Coxet Munck.	
      DOUBLE PRECISION C1	! (E)  Cos(Teta1).
      DOUBLE PRECISION S1	! (E)  Sin(Teta1).
      DOUBLE PRECISION C2	! (E)  Cos(Teta2).
      DOUBLE PRECISION S2	! (E)  Sin(Teta2).
      DOUBLE PRECISION PHI	! (E)  Azimut relatif pour la convention 
      				!      OS (radians).
            
      DOUBLE PRECISION F	! (S) Valeur de la fonction F21nadal/F21fresnel
      				! pour la geometrie d'observation
				! (convertie en luminance normalisee * cos(teta2).	
      
      DOUBLE PRECISION F21FRESNEL	! Valeur de la fonction F21 de Fresnel.
      DOUBLE PRECISION F21NADAL		! Valeur de la fonction F21 de Nadal.
      
      DOUBLE PRECISION COEF

      DOUBLE PRECISION COSDIF	! Cosinus de l'angle de diffusion.      
      DOUBLE PRECISION COSW	! Cosinus de l'angle de reflexion W
      				! (angle d'incidence par rapport a la 
				!  normale du dioptre)
      
      DOUBLE PRECISION V	! Sinus carre de W.
      DOUBLE PRECISION X	! 

      DOUBLE PRECISION IND2	! Carre de l'indice de refraction surface / air      
      DOUBLE PRECISION RL	! Coefficient de reflexion parallele.
      DOUBLE PRECISION RR	! Coefficient de reflexion perpendiculaire.

      DOUBLE PRECISION PI	! Valeur de PI
      

C* Traitements
C---------------------------------    
      PI=4.D+00*DATAN(1.D+00)

      COEF = 4.*C1*PI*SIG
	
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
      IF(F21FRESNEL.LT.1.0D-10) THEN	!Developpement limite au voisinage 
      					!de 0 pour F21Fresnel
         F = ALPHA*BETA / (C1 + C2)
      ELSE	 
         F = F21NADAL / F21FRESNEL
      ENDIF


C*   Conversion reflectance -> luminance normalisee (*C1)
C*   et adaptation aux OS -> *C2.  
C--------------------------------------------------------- 
      F = F*COEF*C2*C1		

      RETURN
      
      END         !FIN DE LA PROCEDURE SOS_CALC_F_ROUJEAN










