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


C****************************************************************************
C* FICHIER: SOS_OS.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module gere les procedures permettant les calculs de resolution
C*       de l'Equation de Transfert Radiatif, en polarisation, par ordres
C*       successifs de diffusion.
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et developpement : CS SI, B.Lafrance, C. Casteras.
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0: Mise a niveau du code des OS du LOA
C* MOD:VERSION:2.0: Correction d'une erreur sur la polarisation de surface
C* MOD:VERSION:2.1: 02/07/2008 Ajustement du format d'affichage du profil pour le 
C*                             fichier de Trace 
C*                             ==> Numero de niveau sur 4 caracteres au lieu de 3.
C* MOD:VERSION:3.0: 10/01/2010 
C*   - Modification de la taille des tableaux de champs de luminance :
C*     ==> Passage de SOS_OS_NBMU a SOS_OS_NBMU_MAX.
C*
C*   - Modification de la taille des tableaux des ordres de developpement :
C*     ==> Passage de SOS_OS_NB  a SOS_OS_NB_MAX.
C*
C*   - Ajout du parametre d'entree : NBMU (nb d'angles effectifs pour les simulations)
C*   - Ajout du parametre d'entree : OS_NB (ordre max des developpements des fonctions 
C		                    de phase en Polynomes de Legendre a utiliser).
C*
C*   - Modification de la ligne d'enregistrement des tableaux de luminances 
C*     a l'ordre IS de la decomposition en series de Fourier 
C*     ==> limitation aux bornes utiles -NBMU:NBMU
C*                  
C*   - Modification de la ligne de lecture des tableaux de la  matrice de 
C*     reflexion a l'ordre IS de la decomposition en series de Fourier :
C*     ==> limitation aux bornes utiles -NBMU:NBMU
C* 
C*   - Adaptation des routines SOS_MAT_FRESNEL_PLAN_REFL, SOS_NOYAUX, 
C*     SOS_INTEGR_EPOPT, SOS_FSOURCE_ORDRE1, SOS_FSOURCE_ORDREIG,
C*     SOS_FSOURCE_DIFF_FRESNEL1, SOS_PARAM_CONV, SOS_ARRET_DIFFUS, 
C*     SOS_ARRET_FOURIER et SOS_AJOUT_QUEUE au passage du parametre
C*     LUM_NBMU (nb d'angles utiles).
C*  
C*   - Adaptation de la routine SOS_NOYAUX au passage de OS_NB
C*     (pour boucler sur les ordres de developpement).
C*
C*   - Le calcul du flux diffus incident (normalise par l'eclairement solaire) EMOINS
C*     est passe en argument de sortie 
C*     (pour fournir la transmission diffuse si le sol est noir)
C*
C*   - Correction d'une erreur de dimensionnement de boucle pour l'annulation de la  
C*     polarisation de surface (IPOLAR=0 & IMAT_SURF=1) : demarrage de boucle 0 --> 1
C*
C* MOD:VERSION:3.1: 08/12/2015
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*      (incluant au décodage des constante) : requis pour compilation gfortran
C****************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    SOS_LENDIR : Longueur des noms de repertoires.
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_PI : valeur de PI                           
C    SOS_OS_NB_MAX :  Valeur maximale pour l'ordre limite du developpement en 
C                     polynomes de Legendre (dimensionnement des tables).
C    SOS_OS_NT : Nombre de couches du profil atmospherique.
C    SOS_OS_NBMU_MAX : Nombre maximal d'angles positifs pour la resolution de l'ETR
C    SOS_OS_IBOR : Ordre minimal du developpement en series de Fourier
C    SOS_PH_SEUIL_CV_SG : Seuil de convergence en serie geometrique.
C    SOS_PH_SEUIL_SUMDIF : Seuil d'arret des diffusions multiples.
C    SOS_PH_SEUIL_SF : Seuil pour l'arret de la decomposition en series 
C                      de Fourier.
C----------------------------------------------------------------------------
#include "SOS.h"



C*==========================================================================
C PROCEDURE: SOS_OS
C ==========
C      Cette procedure calcule le champ de rayonnement diffus (I, Q, U) dans
C      une atmosphere diffusante. 
C      Le champ resultat est decomposee en series de Fourier.
C
C Description des parametres
C ---------------------------
C     NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C
C     RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (DOUBLE) : (E) Cosinus des angles.
C     GA(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (DOUBLE)  : (E) Poids des angles.
C
C     OS_NB (I4) : (E) Ordre max des developpements des fonctions de phase 
C		       en Polynomes de Legendre a utiliser.
C
C     FICSURF (CHARACTER*SOS_LENFIC2)  (E) : 
C                           nom du fichier contenant la matrice de reflexion 
C                           de surface (fichier Glitter ou autres BRDF/BPDF).
C                           (repertoire + nom fichier + extension)
C     FICOS   (CHARACTER*SOS_LENFIC2)  (E) : 
C                           nom du fichier resultat des OS 
C                           (repertoire + nom fichier + extension)
C     N0  (I4)   (E) : Numero de mu de Gauss pour l'angle zenithal solaire.
C                     (angle complementaire : mus = -rmu(N0))
C                     (N0 > 0 => pour en tenir compte).
C     TETAS (double)   (E) : Angle zenithal solaire.
C                           (angle complementaire : mus= -cos(tetas))
C                           (N0 < 0 => pour en tenir compte).
C     RO (double)   (E) : Albedo de surface.
C     IMAT_SURF  (I4)   (E) : Option de BRDF / BPDF de surface
C      		        0 :surface lambertienne d'albedo RHO.
C			1 :surface lambertienne + glitter.
C			2 :surface lambertienne + reflexion de Fresnel.
C		        3 :surface lambertienne + BRDF de Roujean.
C			4 :surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
C			5 :surface lambertienne + BRDF de Roujean + BPDF de Breon.
C			6 :surface lambertienne + BRDF de Roujean + BPDF de Nadal.	
C     IFRESNEL  (I4)  (E) : Option de calcul de reflexion de Fresnel 
C                           sur mer plate
C                        = 0 : pas de reflexion sur mer plate.
C                        = 1 : calcul de la reflexion de Fresnel sur mer plate.
C     IND_SURF (double) (E) : Indice de refraction du dioptre surface / atmosphere.
C     H(0:SOS_OS_NT) (double)   (E) : Epaisseur optique a chaque niveau.
C     XDEL(0:SOS_OS_NT) (double)   (E) : Pourcentage d'aerosols par niveau.
C     YDEL(0:SOS_OS_NT) (double)   (E) : Pourcentage de molecules par niveau.
C     RON (double)   (E) : Facteur de depolarisation (moleculaire).
C     ALPHA(0:SOS_OS_NB_MAX) (double)   (E) : Parametres alpha(k).
C     BETA(0:SOS_OS_NB_MAX)  (double)   (E) : Parametres beta(k).
C     GAMMA(0:SOS_OS_NB_MAX) (double)   (E) : Parametres gamma(k).
C     ZETA(0:SOS_OS_NB_MAX)  (double)   (E) : Parametres zeta(k).
C     IOUT  (I4)   (E) : Type de simulation :
C                       = -1 =>  Sortie standard:
C                                Luminance emergeante au TOA
C                                + descendante au niveau du sol.
C                       = n =>  Luminance montante et descendante
C                               pour le niveau n (0 <= n <= SOS_OS_NT).
C     IGMAX  (I4)   (E) : Ordre maximal de diffusion.
C     IBORM  (I4)   (E) : Ordre maximal du developpement en series de Fourier
C                         en cours de traitement :
C                         _ 2 pour une atmosphere purement moleculaire,
C                         _ OS_NB en presence d'aerosols.
C     IPOLAR  (I4)   (E) : Option de simulation sans polarisation.
C      			   _ 0 : calculs sans polarisation.
C			   _ 1 : calculs normaux avec polarisation.
C     TRACE (logical)   (E) : = vrai,si ecriture dans le fichier trace
C                               (No d'unite logique du fichier trace =99)
C     EMOINS (double)   (S) : Flux diffus eclairant le sol, normalise par
C			      le flux solaire incident au TOA.
C     IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C
C Autre resultat fourni par la procedure
C -------------------------------------
C   Fichier binaire resultat des OS
C   Informations dans le fichier de trace (optionnel)
C   (No d'unite logique du fichier trace =99 , 
C    fichier cree/ouvert par le programe appelant)
C
C 
C Variable d'environnement
C ------------------------
C   Aucune
C
C
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier Vent
C       Fichier binaire non formatte, cree par SOS_SURFACE
C       Il contient la matrice de reflexion de surface 
C       (fichier Glitter ou autres BRDF/BPDF):
C       c'est a dire les elements Rij de la matrice de reflexion 
C              R11(NBMU,NBMU) (reel),
C              R12(NBMU,NBMU) (reel),
C              R13(NBMU,NBMU) (reel),
C              R21(NBMU,NBMU) (reel),
C              R22(NBMU,NBMU) (reel),
C              R23(NBMU,NBMU) (reel),
C              R31(NBMU,NBMU) (reel),
C              R32(NBMU,NBMU) (reel),
C              R33(NBMU,NBMU) (reel)
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
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C       Erreur si le fichier de surface n'existe pas
C     -Lecture/ecriture dans un fichier
C     -Incoherence du niveau de sortie dans le profil atmospherique
C     -Erreur dans une routine appelee
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C     Remarque : Si le fichier resultat des OS existe deja, il est ecrase
C                par le nouveau fichier cree par cette routine
C                (pas de message et pas d'erreur)
C
C      Cas d'un warning 
C     (message sans arret du programme et sans retour d'erreur)
C            Arret des calculs des OS pour une incidence au limbe.
C
C
C  Routines utilises par la procedure:
C ------------------------------------
C  Ce programme fait appel aux routines:
C      - SOS_MAT_FRESNEL_PLAN_REFL
C      - SOS_NOYAUX
C      - SOS_INTEGR_EPOPT
C      - SOS_FSOURCE_ORDRE1  
C      - SOS_FSOURCE_ORDREIG
C      - SOS_FSOURCE_DIFF_FRESNEL1 
C      - SOS_PARAM_CONV 
C      - SOS_ARRET_DIFFUS
C      - SOS_ARRET_FOURIER
C      - SOS_AJOUT_QUEUE
C
C===========================================================================
      SUBROUTINE SOS_OS(NBMU,RMU,GA,OS_NB,
     &                  FICSURF,FICOS,
     &                  N0,TETAS,RO,IMAT_SURF,IFRESNEL,IND_SURF,
     &                  H,XDEL,YDEL,RON,
     &                  ALPHA,BETA,GAMMA,ZETA,IOUT,IGMAX,IBORM,
     &                  IPOLAR,TRACE,EMOINS,IER)

      IMPLICIT NONE


C* Definition des variables   
C*------------------------------------------------------
				
      DOUBLE PRECISION TETAS	! Angle zenithal solaire en degres (0 < tetas < 90).
      				! (angle complementaire : mus= -cos(tetas))
      				! (N0 < 0 => pour en tenir compte).
	
      DOUBLE PRECISION RO	! Albedo de surface.
	
      DOUBLE PRECISION IND_SURF  ! Indice de refraction surface / atmosphere.
      				
      DOUBLE PRECISION H(0:SOS_OS_NT)	  !Epaisseur optique pour chaque niveau.
      DOUBLE PRECISION XDEL(0:SOS_OS_NT)  !Pourcentage d'aerosols par niveau.
      DOUBLE PRECISION YDEL(0:SOS_OS_NT)  !Pourcentage de molecules par niveau.

      DOUBLE PRECISION RON	! Facteur de depolarisation (moleculaire).

      DOUBLE PRECISION ALPHA(0:SOS_OS_NB_MAX)	! Parametres alpha(k).
      DOUBLE PRECISION BETA(0:SOS_OS_NB_MAX)	! Parametres beta(k).
      DOUBLE PRECISION GAMMA(0:SOS_OS_NB_MAX)	! Parametres gamma(k).
      DOUBLE PRECISION ZETA(0:SOS_OS_NB_MAX)	! Parametres zeta(k).

        				
      DOUBLE PRECISION I3OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				! Parametre de Stokes I pour un ordre S
      				! du developpement en series de Fourier,
				! en fonction de l'angle mu.
      DOUBLE PRECISION Q3OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				! Parametre de Stokes Q.
      DOUBLE PRECISION U3OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				! Parametre de Stokes U.

      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
                                               ! Tableau des mu de Gauss.
      DOUBLE PRECISION GA(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   
                                               ! Tableau des poids de Gauss.

      DOUBLE PRECISION RRMU	! Double inversion de mu.

      DOUBLE PRECISION TAB      ! Valeur du cosinus de l'angle zenithal 
                                ! solaire (mus).

C*  - Reflexion lambertienne						
      DOUBLE PRECISION XR(0:SOS_OS_NBMU_MAX) ! champ de luminance de reflexion 
                                         ! lambertienne au sol.

      DOUBLE PRECISION LSOL   ! Luminance de reflexion lambertienne au sol.	

C*  - Contribution de la reflexion directe de surface (BRDF/BPDF) 
C*    sur le signal au TOA et pour le niveau de sortie choisi
      DOUBLE PRECISION RII(SOS_OS_NBMU_MAX)	     ! Contribution pour I TOA. 	
      DOUBLE PRECISION RQQ(SOS_OS_NBMU_MAX)	     ! Contribution pour Q TOA.
      DOUBLE PRECISION RUU(SOS_OS_NBMU_MAX)	     ! Contribution pour U TOA.
      DOUBLE PRECISION RIIOUT(SOS_OS_NBMU_MAX)   ! Contribution pour I 
                                             ! au niveau de sortie. 
      DOUBLE PRECISION RQQOUT(SOS_OS_NBMU_MAX)   ! Contribution pour Q 
                                             ! au niveau de sortie
      DOUBLE PRECISION RUUOUT(SOS_OS_NBMU_MAX)   ! Contribution pour U 
                                             ! au niveau de sortie.

      DOUBLE PRECISION A	! Transmission dans la direction de visee
      				! entre le sol et le niveau de sortie.

      DOUBLE PRECISION RR	! Attenuation du faisceau solaire direct au sol
      				! divisee par mus.

      DOUBLE PRECISION F11(0:SOS_OS_NBMU_MAX)  ! Element de la matrice de Fresnel.
      DOUBLE PRECISION F12(0:SOS_OS_NBMU_MAX)  ! Element de la matrice de Fresnel.
      DOUBLE PRECISION F33(0:SOS_OS_NBMU_MAX)  ! Element de la matrice de Fresnel.

C*  - Pour la diffusion moleculaire :         	
      DOUBLE PRECISION BETA0	! Coefficients du developpement de la matrice 
      DOUBLE PRECISION BETA2	! de phase moleculaire en fonctions de Legendre
      DOUBLE PRECISION GAMMA2	
      DOUBLE PRECISION ALPHA2	

      DOUBLE PRECISION AAA	! Facteur correctif de la depolarisation 
                                ! des molecules.

      DOUBLE PRECISION XPL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Fonction de Legendre PSL à l'ordre L=2 et S=IS
      			! en fonction de l'angle mu.     		
      DOUBLE PRECISION XRL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Fonction RSL a l'ordre L=2 et S=IS en fonction de mu.
      DOUBLE PRECISION XTL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Fonction TSL a l'ordre L=2 et S=IS en fonction de mu.

C*  - Pour la diffusion aerosol :
	! Element Pij de la matrice de phase des aerosols pour un ordre
	! IS de la decomposition en series de Fourier.
	! Soit mup la direction d'incidence (-1 < mup < 1)
	! et mu celle de diffusion (-1 < mu < 1),
	! la diffusion est exprimee par Pij(mu,mup).	
      DOUBLE PRECISION 
     &    BP(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! P12 ou P21.
      DOUBLE PRECISION 
     &    GT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! -P13 ou -P31.
      DOUBLE PRECISION 
     &    ARR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! P22.
      DOUBLE PRECISION 
     &    ART(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! -P23 ou -P32.
      DOUBLE PRECISION 
     &    ATT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! P33.

      DOUBLE PRECISION CH(0:SOS_OS_NT)	! Attenuation du faisceau solaire 
                                        ! direct / 4 pour le melange aerosols
					! + molecules 

C* Variables associees a l'integration spatiale et verticale : 
      	! Champ de luminance d'ordre IG-1 ou IG en fonction du niveau du profil
	! et de la direction.       	
      DOUBLE PRECISION I1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Champ I.
      DOUBLE PRECISION Q1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Champ Q.
      DOUBLE PRECISION U1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Champ U.

      DOUBLE PRECISION 
     &    I1fresnel(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			! Champ I pour la diffusion primaire de la
			! lumiere du rayonnement solaire remontant 
			! apres une reflexion directe sur mer plate,  
			! selon les lois de Fresnel.
      DOUBLE PRECISION 
     &    Q1fresnel(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Champ Q.
      DOUBLE PRECISION 
     &    U1fresnel(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Champ U.
      
	! Vecteur fonction source en fonction du niveau du profil
	! et de la direction.
      DOUBLE PRECISION I2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour I.
      DOUBLE PRECISION Q2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour Q.
      DOUBLE PRECISION U2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour U.

C* Variable associee aux tests de sorties de boucles :
      DOUBLE PRECISION Z1	! Variable seuil pour les tests d'arret.
     	
C*  - Boucle sur l'ordre de la decomposition en series de Fourier :
	! Tableaux de sommation des ordres S de la decomposition
        ! en series de Fourier, en fonction de l'angle mu.
      DOUBLE PRECISION I4(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      DOUBLE PRECISION Q4(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      DOUBLE PRECISION U4(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)

        ! Tableaux de sommation des ordres S de la decomposition
        ! en series de Fourier, avec permutation du signe des ordres S,
	! en fonction de l'angle mu.
      DOUBLE PRECISION I5(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      DOUBLE PRECISION Q5(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      DOUBLE PRECISION U5(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	

      DOUBLE PRECISION SIGN	! Valeur 1 ou -1.
      DOUBLE PRECISION COEF	! Coefficient 1 ou 2.

C*  - Boucle sur les diffusions :
	! Fonctions pour le test de convergence en series geometrique,
	! champ emergeant au TOA (mu > 0) et descendant au sol (mu < 0).
	! Les fonctions sont etablies pour un ordre S du developpement
	! en series de Fourier, en fonction de l'angle mu.	
      DOUBLE PRECISION I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)		
                                ! Parametre de Stokes I pour le cumul des
      				! ordres de diffusion precedent IG en cours.
      DOUBLE PRECISION Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes Q.
      DOUBLE PRECISION U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes U.
            			
      DOUBLE PRECISION A1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
       				! Parametre de Stokes I pour l'ordre IG-2 
				! du champ emergeant au TOA (mu>0),
				! ou descendant au niveau du sol (mu<0). 
      DOUBLE PRECISION B1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG-2.
      DOUBLE PRECISION C1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes U pour l'ordre IG-2.

      DOUBLE PRECISION D1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG-1.	
      DOUBLE PRECISION E1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG-1.
      DOUBLE PRECISION F1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG-1.

      DOUBLE PRECISION G1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG.	
      DOUBLE PRECISION H1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG.
      DOUBLE PRECISION P1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG.

	! Memes fonctions pour le niveau de sortie des resultats:
      DOUBLE PRECISION D1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG-1.	
      DOUBLE PRECISION E1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG-1.
      DOUBLE PRECISION F1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG-1.

      DOUBLE PRECISION G1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG.	
      DOUBLE PRECISION H1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG.
      DOUBLE PRECISION P1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG.

      DOUBLE PRECISION EMOINS	! Flux diffus eclairant le sol, normalise par
      				! le flux solaire incident au TOA.
      DOUBLE PRECISION EPLUS	! Flux diffus emergeant au TOA de l'atmosphere.

      DOUBLE PRECISION II2 	! Valeur intermediaire pour le calcul de I2 
                                ! (mu > 0).
      DOUBLE PRECISION QQ2 	! Valeur intermediaire pour le calcul de Q2 
                                ! (mu > 0).
      DOUBLE PRECISION UU2 	! Valeur intermediaire pour le calcul de U2 
                                ! (mu > 0).
      DOUBLE PRECISION XI1	! Valeur de I du champ montant d'ordre IG-1.
      DOUBLE PRECISION XQ1	! Valeur de Q du champ montant d'ordre IG-1.
      DOUBLE PRECISION XU1	! Valeur de U du champ montant d'ordre IG-1.

      DOUBLE PRECISION Z	! Valeur ponctuelle du tableau des poids 
                                ! de Gauss.


C*  - Elements Rij de la matrice de reflexion de surface
      REAL*4 R11(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL*4 R12(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL*4 R13(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL*4 R21(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL*4 R22(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL*4 R23(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL*4 R31(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL*4 R32(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)
      REAL*4 R33(SOS_OS_NBMU_MAX,SOS_OS_NBMU_MAX)

								
      CHARACTER*SOS_LENFIC2 FICSURF  ! Fichier contenant la matrice de reflexion
      				     ! de surface (fichier Glitter ou autres BRDF/BPDF).

      CHARACTER*SOS_LENFIC2 FICOS    ! Nom du fichier resultat des OS.
      

      LOGICAL TRACE        ! = vrai,si ecriture dans le fichier trace

      INTEGER*4 OS_NB 		! Ordre max des developpements des fonctions de phase 
				! en Polynomes de Legendre.

      INTEGER*4 N0		! Numero de mu de Gauss pour l'angle
      				! zenithal solaire.
				! (angle complementaire : mus = -rmu(N0))
				! (N0 > 0 => pour en tenir compte).
               
      INTEGER*4 IMAT_SURF ! Option de BRDF / BPDF de surface :
      			  !   0 : surface lambertienne d'albedo RHO.
			  !   1 : surface lambertienne + Glitter.
			  !   2 : surface lambertienne + reflexion de Fresnel.
		          !   3 : surface lambertienne + BRDF de Roujean.
			  !   4 : surface lambertienne + BRDF de Roujean + BPDF de Rondeaux.
			  !   5 : surface lambertienne + BRDF de Roujean + BPDF de Breon.
			  !   6 : surface lambertienne + BRDF de Roujean + BPDF de Nadal.
			  
      INTEGER*4 IFRESNEL ! Indice precisant si les simulations
      			 ! introduisent la reflexion de Fresnel sur mer plate :
			 !    1 => pour un calcul de la reflexion de Fresnel.
			 !    Sinon => pas de prise en compte.            
			 
      INTEGER*4 IOUT	! Type de simulation :
			!   = -1 =>  Sortie standard:
			!            Luminance emergeante au TOA
			!          + descendante au niveau du sol.
			!   =  n =>  Luminance montante et descendante
			!            pour le niveau n (0 <= n <= SOS_OS_NT).
 	
      INTEGER*4 IGMAX	! Ordre maximal de diffusion.
      
      INTEGER*4 IBORM	! Ordre maximal du developpement en series de Fourier
      			!  en cours de traitement :
			!   _ 2 pour une atmosphere purement moleculaire,
			!   _ OS_NB en presence d'aerosols.

      INTEGER*4 IPOLAR	! Option de simulation sans polarisation.
      			! 0 : calculs sans polarisation.
			! 1 : calculs normaux avec polarisation.
			
      INTEGER*4 JK	! Numero de stockage de mus dans le tableau RMU
      			! (0 par defaut).

      INTEGER*4 NUP_OUT	   ! Niveau de sortie de I, Q et U pour le 
                           ! champ montant.
      INTEGER*4 NDOWN_OUT  ! Niveau de sortie de I, Q et U pour le 
                           ! champ descendant.

      INTEGER*4 IS	   ! Ordre en cours du developpement en series 
                           ! de Fourier associe a la matrice de phase
			   
      INTEGER*4 IG	   ! Ordre d'interaction en cours 
                           ! (diffusion + reflexion).
								
      INTEGER*4 IER        !code d'erreur =0 si pas d'erreur, =-1 sinon
      INTEGER*4 IFIN       !nb caracteres significatifs d'une chaine
      INTEGER*4 I,K,J      !indices de tableaux
      
      INTEGER*4 NBMU       ! nombre d'angles utiles	

C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0

C* Ouverture des fichiers 
C------------------------------------------------------
C     Fichier BRDF / BPDF de surface

      IF (IMAT_SURF.EQ.1) THEN
          OPEN(UNIT=12,FILE=FICSURF,FORM='UNFORMATTED',
     &         STATUS='OLD',ERR=910)      
      ENDIF  
     
C     Fichier resultat des OS
C     Si ce fichier existe deja, il est ecrase par le nouveau fichier
      OPEN(UNIT=10,FILE=FICOS,FORM='UNFORMATTED',ERR=913)


C* Prise en compte du facteur de depolarisation pour les parametres
C* de la fonction de phase moleculaire.
C---------------------------------------------------------------------
      AAA=RON/(2-RON)
      AAA=(1-AAA)/(1+2*AAA)

      BETA0=1.
      BETA2=0.5*AAA
      GAMMA2=-AAA*DSQRT(1.5D+00)
      ALPHA2=3.*AAA


C* Coupure de la polarisation atmospherique.
C---------------------------------------------------------------------
      IF (IPOLAR.EQ.0) THEN
          GAMMA2 = 0.D+00
          ALPHA2 = 0.D+00
	  
	  DO K=0,OS_NB
	     ALPHA(K) = 0.D+00
	     GAMMA(K) = 0.D+00
             ZETA(K)  = 0.D+00
   	  ENDDO
      
      ENDIF !Coupure de la polarisation
      

C* Definition de mus (cosinus de l'angle zenithal solaire)
C* selon que la simulation est faite pour un angle de Gauss
C* ou pour un angle quelconque.
C----------------------------------------------------------
      IF (N0.GT.0) THEN
         TAB=-RMU(N0)
      ELSE
         TAB=-DCOS(SOS_PI*TETAS/180.)
      ENDIF

C* Stockage de mus pour l'indice 0 du tableau des mu de GAUSS (RMU).
C------------------------------------------------------------------
      JK=0
      RMU(JK)=TAB

C* Fichier trace 
C------------------------------------------------------  
      IF (TRACE) THEN
        WRITE(99,*,err=921) ' '  
        WRITE(99,*,err=921) ' '	   
        WRITE(99,*,err=921) 'ROUTINE OS  '
	WRITE(99,*,err=921) '============='
        WRITE(99,*,err=921) ' '
        WRITE(99,*,err=921) 'PARAMETRES DE SIMULATIONS'
	WRITE(99,*,err=921) '------------------------- '
	WRITE(99,*,err=921)' '   
        WRITE(99,*,err=921)' Incidence solaire '
	WRITE(99,*,err=921)'    N0 : ',N0
	IF (N0.LE.0) THEN
	   WRITE(99,*,err=921)'    Angle zenithal solaire : ',TETAS 
	ELSE
	   WRITE(99,*,err=921)
     &       '    Angle zenithal solaire : ', 
     &            DACOS(RMU(N0))*180.D+00/SOS_PI
	ENDIF
	WRITE(99,*,err=921)'    Valeur de mus : ' ,TAB
	WRITE(99,*,err=921)
	WRITE(99,*,err=921)'    Transmission directe TOA -> surface : '
	WRITE(99,*,err=921)'    (avec AOT tronquee si troncature) ',
     &                          DEXP(H(SOS_OS_NT)/(TAB))
        WRITE(99,*,err=921) 
        IF (IPOLAR.EQ.1) THEN
	   WRITE(99,*,err=921)'    Simulation de la polarisation'
	ELSE
	   WRITE(99,*,err=921)'    Simulation SANS polarisation' 
	ENDIF
	WRITE(99,*,err=921)' ' 
	WRITE(99,*,err=921)' Surface' 
	WRITE(99,*,err=921)'    Albedo de surface : ',RO 
        IF (IMAT_SURF.EQ.1) THEN
	   WRITE(99,*,err=921)'    Calcul avec BRDF / BPDF de surface'
	   WRITE(99,*,err=921)
     &       '    Fichier de la matrice de reflexion : ', FICSURF
	ELSE
	   WRITE(99,*,err=921)'    Pas de BRDF / BPDF de surface' 
	ENDIF
	IF (IFRESNEL.EQ.1) THEN
	   WRITE(99,*,err=921)'    Calcul des reflexions sur mer plate'
	   WRITE(99,*,err=921)
     &	   '    Indice de refraction surface / atmosphere', IND_SURF  
	ELSE
	   WRITE(99,*,err=921)'    Pas de reflexions sur mer plate' 
	ENDIF				
	WRITE(99,*,err=921)' '
	WRITE(99,*,err=921)' Profil atmospherique ' 
	WRITE(99,*,err=921)'  (avec troncature si appliquee)' 
	WRITE(99,*,err=921)'    Nombre de couches : ',SOS_OS_NT 
	WRITE(99,*,err=921)
     &	     '      Couche  % aerosols    % molecules   Ep. opt.'
        DO 3 J=0,SOS_OS_NT
	   WRITE(99,'(8X,I4,5X,F6.4,8X,F6.4,8X,F6.4)
     &	        ',err=921) J,XDEL(J),YDEL(J),H(J)
   3    CONTINUE
	WRITE(99,*,err=921)' '	   
	WRITE(99,*,err=921)'    Coefficients de la matrice de phase :'
        WRITE(99,*,err=921) ' '  
	WRITE(99,35,err=921)'I','ALPHA(I)','BETA(I)','GAMMA(I)',
     &			    'ZETA(I)'  
        DO K=0,OS_NB
	   WRITE(99,25,err=921) K,ALPHA(K),BETA(K),GAMMA(K),ZETA(K)
   	ENDDO
	WRITE(99,*,err=921)' '	
	WRITE(99,*,err=921)
     &	     '    Coefficients de la matrice de phase moleculaire'
	WRITE(99,*,err=921)'         Facteur de depolarisation : ',RON
	WRITE(99,*,err=921)
     &       '         Facteur correctif de la depolarisation : ',AAA
	WRITE(99,*,err=921)'         Beta0 = ',BETA0
	WRITE(99,*,err=921)'         Beta2 = ',BETA2	
	WRITE(99,*,err=921)'         Gamma2 = ',GAMMA2
	WRITE(99,*,err=921)'         Alpha2 = ',ALPHA2	
	WRITE(99,*,err=921)' '	
	WRITE(99,*,err=921)'RESULTATS'
	WRITE(99,*,err=921)'---------'
      END IF


C* Arret des calculs pour une incidence au limbe.
C------------------------------------------------------
      IF(TAB.EQ.0.D+00) GOTO 995

C* Test coherence du niveau de sortie dans le profil atmospherique
C----------------------------------------------------------------
      IF ((IOUT.LT.-1).OR.(IOUT.GT.SOS_OS_NT)) GOTO 996

C* Definition des niveaux du profil pour le champ montant et descendant.
C-----------------------------------------------------------------------
      IF (IOUT.EQ.-1) THEN
         NDOWN_OUT = SOS_OS_NT		!Sol.
	 NUP_OUT = 0			!TOA.
      ELSE
         NDOWN_OUT = IOUT		!Niveau utilisateur.
	 NUP_OUT = IOUT			!Niveau utilisateur.
      ENDIF



C* Calcul des coefficients de reflexion de la matrice de Fresnel
C---------------------------------------------------------------
      IF (IFRESNEL.EQ.1) THEN
          CALL SOS_MAT_FRESNEL_PLAN_REFL(NBMU,RMU,IND_SURF,IPOLAR,
     &	                                 F11,F12,F33)
	  
	  IF (TRACE) THEN
	     WRITE(99,*,err=921) 
     &       ' Matrice de reflexion de Fresnel : '
             DO J= 0,NBMU
                WRITE(99,40,err=921) 
     &		      ' Angle reflexion ',DACOS(RMU(J))*180./SOS_PI,
     &	              ' F11 ',F11(J),' F12 ',F12(J),' F33 ',F33(J)   
             ENDDO
             WRITE(99,*,err=921) ' '	   
	  ENDIF
	  
      ENDIF !reflexion de Fresnel sur mer plate


C* Attenuation du faisceau solaire direct (divisee par 4).
C------------------------------------------------------
      DO 84 I=0,SOS_OS_NT
         CH(I)=DEXP(-H(I)/(-TAB))/4.
  84  CONTINUE


C* Initialisation des tableaux pour le test d'arret de la decomposition
C* en series de Fourier.
C-----------------------------------------------------------------------
      DO 22 J=-NBMU,NBMU
         Q4(J)=0.
         U4(J)=0.
         I4(J)=0.

         Q5(J)=0.
         U5(J)=0.
         I5(J)=0.
 22   CONTINUE


C--------------------------------------------------------------
C
C       ******  DECOMPOSITION EN SERIE DE FOURIER   ******
C
C--------------------------------------------------------------
      SIGN=-1.
	
C* Boucle sur les ordres de la decomposition en series de Fourier
C------------------------------------------------------
      DO 24 IS=SOS_OS_IBOR,IBORM
 	 
         SIGN=-SIGN
 
C* Fichier trace 
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(99,*,err=921) 
     &       ' Ordre de la decomposition en serie de Fourier : IS = ',IS
            WRITE(99,*,err=921) 
     &       ' -----------------------------------------------'
            WRITE(99,*,err=921) ' '	   
         ENDIF

C* Suppression du coefficient moleculaire BETA0 pour le calcul
C* de la matrice de phase de diffusion moleculaire pour les
C* ordres IS > 0.
C--------------------------------------------------------------
         IF(IS.GT.0) BETA0=0.

C* Initialisation des tableaux donnant les parametres I, Q, U 
C* pour l'ordre IS de la decomposition en series de Fourier.
C------------------------------------------------------
         DO 23 J=-NBMU,NBMU
            I3(J)=0.
            Q3(J)=0.
            U3(J)=0.
            I3OUT(J)=0.
            Q3OUT(J)=0.
            U3OUT(J)=0.
   23    CONTINUE

C* Lecture des elements de la matrice de reflexion pour 
C* l'ordre IS en cours.
C--------------------------------------------------------------
         IF (IMAT_SURF.EQ.1) THEN
C*      Lecture des fichiers qui ont ete enregistres avec des 
C*      bornes de tableaux limitees au dimensionnement 
C*      utile -NBMU:NBMU         
	     READ(12,err=920),
     &	       ((R11(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R12(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R13(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R21(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R22(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R23(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R31(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R32(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R33(I,J),I=1,NBMU),J=1,NBMU)  
     	     
C*    Coupure de la polarisation de surface.
             IF (IPOLAR.EQ.0) THEN 	  
	         DO K=1,NBMU
		    DO J=1,NBMU
	               R12(J,K) = 0.
		       R13(J,K) = 0.
		       R21(J,K) = 0.
		       R22(J,K) = 0.
		       R23(J,K) = 0.
		       R31(J,K) = 0.
		       R32(J,K) = 0.
		       R33(J,K) = 0.	     
   	            ENDDO
		 ENDDO        
             ENDIF !IPOLAR=0
      
	 ENDIF  !IMAT_SURF=1
	   

C* Appel de la routine de calcul des noyaux de diffusion pour 
C* l'ordre IS de la decomposition en series de Fourier.
C--------------------------------------------------------------
         CALL SOS_NOYAUX(IS,NBMU,RMU,OS_NB,ALPHA,BETA,GAMMA,ZETA,
     &                   XPL,XRL,XTL,BP,GR,GT,ARR,ART,ATT)



C--------------------------------------------------------------
C
C   ******  CALCULS DE DIFFUSION ET REFLEXION PRIMAIRE   ******
C
C--------------------------------------------------------------

C* Calcul de la fonction source de diffusion primaire : I2, Q2, U2
C----------------------------------------------------
         CALL SOS_FSOURCE_ORDRE1(IS,NBMU,JK,XDEL,YDEL,
     &                           BETA0,BETA2,GAMMA2,
     &                           XPL,XRL,XTL,BP,GR,GT,CH,
     &                           I2,Q2,U2)       


C* Condition limite au sol
C-------------------------
         DO 108 K=1,NBMU
            I1(SOS_OS_NT,K)=0.
            Q1(SOS_OS_NT,K)=0.
            U1(SOS_OS_NT,K)=0.
            XR(K)=0.

C* Reflexion lambertienne au sol. Stockage dans XR.
C------------------------------------------------
            IF (RO.EQ.0..OR.IS.NE.0) GOTO 1022
            I1(SOS_OS_NT,K)=-RO*TAB*DEXP(H(SOS_OS_NT)/TAB)
            XR(K)=I1(SOS_OS_NT,K)

C* Reflexion de surface avec BRDF / BPDF
C---------------------------------------
 1022       IF (IMAT_SURF.EQ.1) THEN
               RR = DEXP(H(SOS_OS_NT)/TAB) / RMU(K)

               I1(SOS_OS_NT,K)=  I1(SOS_OS_NT,K) + R11(N0,K)*RR
               Q1(SOS_OS_NT,K)=  R21(N0,K)*RR
               U1(SOS_OS_NT,K)=  R31(N0,K)*RR    	
            ENDIF
	    
  108    CONTINUE


C* Integration du champ de luminance montant (mu > 0) du sol au TOA.
C* et du champ de luminance descendant (mu < 0) du TOA au sol.
C--------------------------------------------------------------
         CALL SOS_INTEGR_EPOPT(NBMU,RMU,H,I2,Q2,U2,I1,Q1,U1)



C* Reflexion de Fresnel :
C* --------------------
C* Pour un eclairement solaire direct, la luminance montante (L1fresnel)
C* obtenue par reflexion de Fresnel est reduite a la direction mu = -mus.
C* Il n'est pas possible de la developper en series de Fourier.
C* Par contre, on developpe en series de Fourier la luminance de
C* premiere diffusion associee a L1fresnel.
C--------------------------------------------------------------
         IF (IFRESNEL.EQ.1) THEN

C*   Calcul de la fonction source pour la diffusion de L1fresnel 	 
	     CALL SOS_FSOURCE_DIFF_FRESNEL1(
     &   	       F11(0),F12(0),XDEL,YDEL,
     &                 BETA0,BETA2,GAMMA2,ALPHA2,BP,GR,GT,ARR,ART,
     &                 XPL,XRL,XTL,IS,NBMU,TAB,H,I2,Q2,U2)
     

C*   Integration sur le profil atmospherique. 
C*   Calcul du champ de luminance pour la diffusion primaire de L1fresnel.
	 
             DO 111 K=1,NBMU
                 I1fresnel(SOS_OS_NT,K)=0.	!Initialisation au sol.
                 Q1fresnel(SOS_OS_NT,K)=0.
                 U1fresnel(SOS_OS_NT,K)=0. 
111	     CONTINUE
     
 	     CALL SOS_INTEGR_EPOPT(NBMU,RMU,H,I2,Q2,U2,I1fresnel,
     &	                           Q1fresnel,U1fresnel)
     

C* Sommation des termes de diffusion primaire :
C*   - du faisceau solaire direct + reflexion lambertienne
C*   - du faisceau solaire reflechi sur mer plate, selon les lois de Fresnel. 
            DO 112 I=0,SOS_OS_NT
               DO 113 K=-NBMU,NBMU
                     I1(I,K)=I1(I,K)+I1fresnel(I,K)
                     Q1(I,K)=Q1(I,K)+Q1fresnel(I,K)
                     U1(I,K)=U1(I,K)+U1fresnel(I,K) 
113	       CONTINUE		     
112	    CONTINUE

	 ENDIF	!Fin des calculs pour IFRESNEL=1
	 

	 
C* Stockage dans RII, RQQ et RUU de la contribution de la reflexion
C* directe de surface (BRDF / BPDF) dans le signal pour le niveau 
C* de sortie souhaite (TOA : Nup_out=0, sol : Nup_out = NT).
C--------------------------------------------------------------
         DO 4325 K=1,NBMU
            RII(K)=0. 	
            RQQ(K)=0.
            RUU(K)=0.
            RIIOUT(K)=0. 	
            RQQOUT(K)=0.
            RUUOUT(K)=0.
 4325    CONTINUE
		
         IF (IMAT_SURF.EQ.1) THEN

            DO 4321 K=1 , NBMU
               A = -( H(SOS_OS_NT) - H(NUP_OUT) ) / RMU(K)
               A= DEXP(A)
               RIIOUT(K) = A*(I1(SOS_OS_NT,K)-XR(K))
               RQQOUT(K) = A*Q1(SOS_OS_NT,K)
               RUUOUT(K )= A*U1(SOS_OS_NT,K)
	
	       A = -H(SOS_OS_NT)  / RMU(K)
               A= DEXP(A)
	       RII(K) = A*(I1(SOS_OS_NT,K)-XR(K))
               RQQ(K) = A*Q1(SOS_OS_NT,K)
               RUU(K )= A*U1(SOS_OS_NT,K)
 4321       CONTINUE

         ENDIF


C------------------------------------------------------------
C* Stockage des resultats de calcul de l'interaction d'ordre 1 
C* (diffusion primaire et reflexion directe).
C-------------------------------------------------------------

C* Stockage du champ descendant d'ordre IG = 1.
C------------------------------------------------
         DO K=-NBMU,-1

C* Au sol pour le test de convergence :
            I3(K)=I1(SOS_OS_NT,K)
	    Q3(K)=Q1(SOS_OS_NT,K)
            U3(K)=U1(SOS_OS_NT,K)
            D1(K)=I1(SOS_OS_NT,K)
            E1(K)=Q1(SOS_OS_NT,K)
            F1(K)=U1(SOS_OS_NT,K)
	
C* Au niveau de sortie des resultats :
           I3OUT(K)=I1(NDOWN_OUT,K)
	   Q3OUT(K)=Q1(NDOWN_OUT,K)
           U3OUT(K)=U1(NDOWN_OUT,K)
         ENDDO


C* Stockage du champ montant d'ordre IG = 1.
C-------------------------------------------
         DO K=1,NBMU

C* Au TOA pour le test de convergence :
            I3(K)=I1(0,K)
	    Q3(K)=Q1(0,K)
            U3(K)=U1(0,K)
	    D1(K)=I1(0,K)
            E1(K)=Q1(0,K)
            F1(K)=U1(0,K)
	 
C* Au niveau de sortie des resultats : 	
            I3OUT(K)=I1(NUP_OUT,K)
	    Q3OUT(K)=Q1(NUP_OUT,K)
            U3OUT(K)=U1(NUP_OUT,K)
         ENDDO


C--------------------------------------------------------------
C
C   ******  BOUCLE SUR LES INTERACTIONS SUCCESSIVES   ******
C
C--------------------------------------------------------------

         IG=1

  503    IG=IG+1
  
C* Test d'arret des calculs si l'ordre maximal de diffusion est atteint
C-----------------------------------------------------------------------
         IF (IG.GT.IGMAX) GOTO 505


C* Calcul de la fonction source d'ordre IG : I2, Q2, U2
C----------------------------------------------------
         CALL SOS_FSOURCE_ORDREIG(IS,NBMU,XDEL,YDEL,
     &                            BETA0,BETA2,GAMMA2,ALPHA2,
     &                            XPL,XRL,XTL,I1,Q1,U1,
     &                            BP,GR,GT,ARR,ART,ATT,GA,
     &                            I2,Q2,U2)
     
	
C* Condition limite au sol
C-------------------------
         DO 1081 K=1,NBMU
            I1(SOS_OS_NT,K)=0.
            Q1(SOS_OS_NT,K)=0.
            U1(SOS_OS_NT,K)=0.
            XR(K)=0.
 1081    CONTINUE	    

	    
C* Integration de la reflexion sur un sol lambertien de 
C* la luminance descendante d'ordre IG-1.
C-------------------------------------------------------
         LSOL=0.
         DO 9876 J=1,NBMU
            LSOL=LSOL+GA(J)*I1(SOS_OS_NT,-J)*RMU(J)            
 9876    CONTINUE
         LSOL=2*LSOL*RO

C* Reflexion lambertienne au sol. Stockage dans XR.
C--------------------------------------------------
         IF (RO.EQ.0..OR.IS.NE.0) GOTO 2121
  
         DO 2122 J=1,NBMU
            I1(SOS_OS_NT,J) = LSOL  
            XR(J) = LSOL
 2122    CONTINUE

C* Reflexion de surface (BRDF / BPDF)
C------------------------------------		
 2121   IF (IMAT_SURF.EQ.1) THEN
 
            DO 960 K=1,NBMU
               II2=0.
               QQ2=0.
               UU2=0.
               RRMU=2/RMU(K)
	 
               DO 961 J=1,NBMU
                  Z= GA(J)
                  XI1 = I1(SOS_OS_NT,-J)
                  XQ1 = Q1(SOS_OS_NT,-J)
                  XU1 = U1(SOS_OS_NT,-J)

                  II2 = II2 +Z*( XI1*R11(J,K)+XQ1*R12(J,K)+XU1*R13(J,K)) 
                  QQ2 = QQ2 +Z*( XI1*R21(J,K)+XQ1*R22(J,K)+XU1*R23(J,K)) 
                  UU2 = UU2 +Z*( XI1*R31(J,K)+XQ1*R32(J,K)+XU1*R33(J,K))
		  
  961          CONTINUE
		
C* Luminance montante au niveau du sol (ordre IG)
               I1(SOS_OS_NT,K)=II2*RRMU +XR(K)
               Q1(SOS_OS_NT,K)=QQ2*RRMU
               U1(SOS_OS_NT,K)=UU2*RRMU
  960       CONTINUE

         ENDIF


C* Reflexion de Fresnel sur mer plate
C------------------------------------
      IF (IFRESNEL.EQ.1) THEN
 
          DO 1510 K=1,NBMU
	  
	     I1(SOS_OS_NT,K) = I1(SOS_OS_NT,K) 
     &            + F11(K)*I1(SOS_OS_NT,-K) + F12(K)*Q1(SOS_OS_NT,-K)
     
             Q1(SOS_OS_NT,K) = Q1(SOS_OS_NT,K) 
     &            + F12(K)*I1(SOS_OS_NT,-K) + F11(K)*Q1(SOS_OS_NT,-K)
     
             U1(SOS_OS_NT,K) = U1(SOS_OS_NT,K) + F33(K)*U1(SOS_OS_NT,-K) 
	     
1510 	  CONTINUE 

      ENDIF

C* Integration du champ de luminance montant (mu > 0) du sol au TOA.
C* et du champ de luminance descendant (mu < 0) du TOA au sol.
C-------------------------------------------------------------------
         CALL SOS_INTEGR_EPOPT(NBMU,RMU,H,I2,Q2,U2,I1,Q1,U1)

C* Stockage du champ d'ordre IG.
C------------------------------
         DO K=-NBMU,-1   
C* Au sol pour le test de convergence :
            G1(K)=I1(SOS_OS_NT,K)
            H1(K)=Q1(SOS_OS_NT,K)
            P1(K)=U1(SOS_OS_NT,K)  	
	
C* Au niveau de sortie des resultats : 	
	    G1OUT(K)=I1(NDOWN_OUT,K)
            H1OUT(K)=Q1(NDOWN_OUT,K)
            P1OUT(K)=U1(NDOWN_OUT,K)
         ENDDO

         DO K=1,NBMU     
C* Au TOA pour le test de convergence :
	    G1(K)=I1(0,K)
            H1(K)=Q1(0,K)
            P1(K)=U1(0,K)
	
C* Au niveau de sortie des resultats :
	    G1OUT(K)=I1(NUP_OUT,K)
	    H1OUT(K)=Q1(NUP_OUT,K)
            P1OUT(K)=U1(NUP_OUT,K)
         ENDDO
 

C*  Test de convergence en serie geometrique et stockage de l'ordre ig-2
C-------------------------------------------------------------------------
         IF (IG.EQ.2) GO TO 506
      
         CALL SOS_PARAM_CONV(NBMU,A1,B1,C1,D1,E1,F1,G1,H1,P1,
     &                       I3,Q3,U3,Z1)

C* Si on excede le seuil SOS_PH_SEUIL_CV_SG pour au moins un angle mu, 
C* on n'est pas en condition de convergence. On passe a l'ordre suivant 
C* de diffusion.
         IF (Z1.GT.SOS_PH_SEUIL_CV_SG) GO TO 506


C* Si on est en regime de convergence, on ajoute la queue de la
C* serie geometrique a la somme des diffusions precedentes.
C-------------------------------------------------------
         CALL SOS_AJOUT_QUEUE(NBMU,D1,E1,F1,G1,H1,P1,
     &                     D1OUT,E1OUT,F1OUT,
     &                     G1OUT,H1OUT,P1OUT,
     &                     I3,Q3,U3,I3OUT,Q3OUT,U3OUT)

C* Fichier trace 
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(99,*,err=921) 
     &       '    Convergence en serie geometrique a partir de IG = ',IG
            WRITE(99,*,err=921) ' '	   
         ENDIF

C* On quitte la boucle sur les diffusions.
C* Passage au test sur l'ordre de la decomposition en series de Fourier.
C---------------------------------------------------------------------
         GO TO 505


C* On reste dans la boucle sur les diffusions.
C* Stockage de l'ordre IG-1 dans A1, B1 et C1.
C* Stockage de l'ordre IG dans D1, E1 et F1.
C* pour le test au passage de l'ordre suivant IG+1.
C-------------------------------------------------------
  506    DO 86 K=-NBMU,NBMU
            A1(K)=D1(K)
            B1(K)=E1(K)
            C1(K)=F1(K)
            D1(K)=G1(K)
            E1(K)=H1(K)
            F1(K)=P1(K)
            D1OUT(K)=G1OUT(K)
            E1OUT(K)=H1OUT(K)
            F1OUT(K)=P1OUT(K)
   86    CONTINUE

C* Calcul de la somme des diffusions
C------------------------------------
17       DO 400 J=1,NBMU
            I3(J)=I3(J)+I1(0,J)
            Q3(J)=Q3(J)+Q1(0,J)
            U3(J)=U3(J)+U1(0,J)
            I3(-J)=I3(-J)+I1(SOS_OS_NT,-J)
            Q3(-J)=Q3(-J)+Q1(SOS_OS_NT,-J)
            U3(-J)=U3(-J)+U1(SOS_OS_NT,-J)

            I3OUT(J)=I3OUT(J)+I1(NUP_OUT,J)
            Q3OUT(J)=Q3OUT(J)+Q1(NUP_OUT,J)
            U3OUT(J)=U3OUT(J)+U1(NUP_OUT,J)
            I3OUT(-J)=I3OUT(-J)+I1(NDOWN_OUT,-J)
            Q3OUT(-J)=Q3OUT(-J)+Q1(NDOWN_OUT,-J)
            U3OUT(-J)=U3OUT(-J)+U1(NDOWN_OUT,-J)
  400    CONTINUE

C* Test d'arret sur les diffusions quand les termes suivant 
C* deviennent tres petits
C-------------------------------------------------------
         CALL SOS_ARRET_DIFFUS(NBMU,I1,Q1,U1,I3,Q3,U3,Z1)
   
         IF (Z1.GT.SOS_PH_SEUIL_SUMDIF) goto 508
      
C* Fichier trace 
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(99,*,err=921) 
     &        '    Sortie de la boucle sur les diffusions par le test',
     &        '    sur l''importance de l''ordre IG devant la somme',
     &        '    des diffusions precedentes'
            WRITE(99,*,err=921) ' '	   
         ENDIF
	 
	 
         GO TO 505

C* Test d'arret sur le nombre maximal de diffusions.
C---------------------------------------------------
  508    IF(IG.LT.IGMAX) GO TO 503
  
C* Fichier trace 
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(99,*,err=921) 
     &       '    Sortie de la boucle sur les diffusions pour ordre',
     &       '    maximal atteint : IGMAX = ',IGMAX
            WRITE(99,*,err=921) ' '	   
         ENDIF
	 
  505    CONTINUE

C* Suppression du terme de reflexion direct sur la surface (BRDF / BPDF).
C-----------------------------------------------------------------------	
        IF (IMAT_SURF.EQ.1) THEN
	 
            DO 5200 J=1,NBMU
               I3(J)=I3(J)-RII(J)
               Q3(J)=Q3(J)-RQQ(J)
               U3(J)=U3(J)-RUU(J)
               I3OUT(J)=I3OUT(J)-RIIOUT(J)
               Q3OUT(J)=Q3OUT(J)-RQQOUT(J)
               U3OUT(J)=U3OUT(J)-RUUOUT(J)
 5200       CONTINUE

         ENDIF

C* Calcul de l'eclairement diffus en haut et en bas de l'atmosphere,
C* normalise par le flux solaire incident au TOA.
C*  EMOINS ==> transmission diffuse td(thetas)
C*        si la surface est non reflechissante.
C*  EPLUS ==> inclus les differentes reflexions surface - atmosphere.
C-------------------------------------------------------------------
         IF(IS.EQ.0) THEN
           EMOINS=0.
           EPLUS=0.
           DO J=1,NBMU
             EMOINS = EMOINS + RMU(J)*GA(J)*I3(-J)
             EPLUS  = EPLUS  + RMU(J)*GA(J)*I3(J)
           ENDDO
	   EMOINS = -EMOINS*2/TAB
	   EPLUS  = -EPLUS *2/TAB
         ENDIF	

C* Calcul des fonctions de mu I4, Q4, U4, I5, Q5 et U5:	
C-------------------------------------------------------
         COEF=2.
         IF (IS.EQ.0) COEF=1.

         DO 67 J=-NBMU,NBMU
            IF (J.EQ.0) GOTO 67
	
	    I4(J)=I4(J)+COEF*I3(J)
            Q4(J)=Q4(J)+COEF*Q3(J)
            U4(J)=U4(J)+COEF*U3(J)

	    I5(J)=I5(J)+COEF*I3(J)*SIGN
            Q5(J)=Q5(J)+COEF*Q3(J)*SIGN
            U5(J)=U5(J)+COEF*U3(J)*SIGN
67       CONTINUE

C* Fichier trace 
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(99,*,err=921) '    Cumul de diffusion : '
	    WRITE(99,*,err=921)	    
     &        '       Niveau de sortie du champ montant : ',NUP_OUT
            WRITE(99,*,err=921) 
     &        '       Niveau de sortie du champ descendant : ',NDOWN_OUT        
            WRITE(99,*,err=921) '  '	   	    
            DO K=1,NBMU
              WRITE(99,20,err=921) RMU(K),
     &                             '  ID ',I3OUT(-K),
     &                             '  QD ',Q3OUT(-K),
     &                             '  UD ',U3OUT(-K),	      
     &                             '  IM ',I3OUT(K),
     &                             '  QM ',Q3OUT(K),
     &                             '  UM ',U3OUT(K)     
 	    ENDDO
            WRITE(99,*,err=921) '  '
	    WRITE(99,*,err=921) ' '	   	    	   
         ENDIF


C* Enregistrement des tableaux de luminances a l'ordre IS de 
C* la decomposition en series de Fourier, pour le niveau du profil 
C* requis par l'utilisateur.
C* Bornes d'enregistrement des tableaux limitees au dimensionnement 
C* utile -NBMU:NBMU
C--------------------------------------------------------------
         WRITE(10,err=923)(Q3OUT(K),K=-NBMU,NBMU),
     & 			  (U3OUT(K),K=-NBMU,NBMU),
     &			  (I3OUT(K),K=-NBMU,NBMU)
      


C* Test d'arret sur la decomposition en serie fourier
C-------------------------------------------------------
         CALL SOS_ARRET_FOURIER(NBMU,I3,Q3,U3,I4,Q4,U4,I5,Q5,U5,Z1)
 
C* Si au moins un terme est superieur au seuil, on passe
C* a l'ordre IS suivant.
C-------------------------------------------------------
         IF (Z1.GT.SOS_PH_SEUIL_SF) GO TO 24
	 
C* Sortie de la boucle sur IS
C---------------------------
         GOTO 243

C* FIN DE BOUCLE SUR IS
C----------------------
   24 CONTINUE
   
  243 CONTINUE   
   

C* Fichier trace 
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(99,*,err=921) ' Arret de la decomposition en serie',
     &                          ' de Fourier pour l''ordre : ',IS
            WRITE(99,*,err=921) ' '
	    WRITE(99,*,err=921) ' '
	    WRITE(99,*,err=921) '    Eclairement normalise de la couche'
            WRITE(99,*,err=921) 
     &            '      Flux diffus eclairant le sol ', EMOINS
            WRITE(99,*,err=921) 
     &      '        --> inclus les multiples iterations surface / atm.'
            WRITE(99,*,err=921)
     &      '        --> transmission diffuse td(thetas) uniquement si '
            WRITE(99,*,err=921)
     &      '            la surface est non reflechissante.'
            WRITE(99,*,err=921) 
     &            '      Flux diffus emergeant au TOA ',EPLUS   
            WRITE(99,*,err=921) 
     &      '        --> inclus les multiples iterations surface / atm.'
            WRITE(99,*,err=921) ' '	   
         ENDIF

C* Fin nominale 
C-------------------
      goto 9999
      
C* Cas d'erreur  retour du status -1 au programme appelant
C----------------------------------------------------------
  910 WRITE(6,*) '  ERROR on SURFACE file opening for SOS_OS'
      IER=-1  
      goto 9999 
  913 WRITE(6,*) '  ERROR on SOS binary result file opening for SOS_OS'
      IER=-1
      goto 9999    
  920 WRITE(6,*) '  ERROR on SURFACE file reading for SOS_OS'
      IER=-1
      goto 9999
  921 WRITE(6,*) '  ERROR on logfile writing for SOS_OS'
      IER=-1
      goto 9999
  922 WRITE(6,*) '  ERROR on angles file reading for SOS_OS'
      IER=-1
      goto 9999
  923 WRITE(6,*) '  ERROR on SOS binary result file writing for SOS_OS'
      IER=-1
      goto 9999
  996 WRITE(6,*) '  Inconsistency output level in the ',
     &           '  atmospheric profile'
      IER=-1
      goto 9999

C* Cas particulier
C--------------------     
  995 WRITE(6,*) '  SOS computation are stopped ',
     &           '  because of a limb incidence'
      goto 9999

      
 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(21)
      IF (IMAT_SURF.EQ.1) CLOSE(12)
      CLOSE(10)

C* Format 
C----------  
   15 FORMAT(2(E21.14))
   10 FORMAT(F9.6,A5,E12.5,A5,E12.5,A5,E12.5) 
   20 FORMAT(F9.6,6(A5,E12.5)) 
   25 FORMAT(1X,I3,4(2X,E13.5))
   35 FORMAT(3X,A1,3X,A6,11X,A7,10X,A8,9X,A7)
   40 FORMAT(A17,F6.2,3(A5,E13.7))          
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_OS




C*===========================================================================
C PROCEDURE: SOS_MAT_FRESNEL_PLAN_REFL
C ==========
C      Cette procedure calcule les elements de la matrice de reflexion
C      de Fresnel pour des angles de reflexion equivalents aux angles
C      de Gauss.
C
C      Pour une incidence -mu, la matrice de reflexion de Fresnel s'ecrit :
C           F11(mu)    F12(mu)     0
C           F12(mu)    F11(mu)     0
C              0         0     F33(mu)
C
C
C Description des parametres
C -------------------------- 
C
C     NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C
C     RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)  Tableau des mu de Gauss.
C
C     IND_SURF (double)  (E) Indice de refraction surface / atmosphere.
C
C     IPOLAR  (I4)   (E) : Option de simulation sans polarisation.
C      			   _ 0 : calculs sans polarisation.
C			   _ 1 : calculs normaux avec polarisation.
C	
C     F11(0:SOS_OS_NBMU_MAX)  (double)  (S) 
C     F12(0:SOS_OS_NBMU_MAX)  (double)  (S)
C     F33(0:SOS_OS_NBMU_MAX)  (double)  (S)
C           Elements de la matrice de reflexion de Fresnel.
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
      SUBROUTINE SOS_MAT_FRESNEL_PLAN_REFL(NBMU,RMU,IND_SURF,IPOLAR,
     &                                     F11,F12,F33)
     
      IMPLICIT NONE


C* Definition des variables    
C*-----------------------------------------------------------------------

      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
                                               ! Tableau des mu de Gauss.

      DOUBLE PRECISION IND_SURF  ! Indice de refraction surface / atmosphere.
      
      DOUBLE PRECISION F11(0:SOS_OS_NBMU_MAX)  ! Element de la matrice de Fresnel.
      DOUBLE PRECISION F12(0:SOS_OS_NBMU_MAX)  ! Element de la matrice de Fresnel.
      DOUBLE PRECISION F33(0:SOS_OS_NBMU_MAX)  ! Element de la matrice de Fresnel.


      DOUBLE PRECISION MU   ! Valeur d'un mu de Gauss (direction de reflexion).
      DOUBLE PRECISION MU2  ! Carre de mu.
      DOUBLE PRECISION IND2 ! Carre de IND_SURF.
      DOUBLE PRECISION X    
      DOUBLE PRECISION RL   ! Coefficient de reflexion parallele.
      DOUBLE PRECISION RR   ! Coefficient de reflexion perpendiculaire.

      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations
      INTEGER*4 IPOLAR	! Option de simulation avec / sans polarisation.
			      
      INTEGER J  ! Indice 

C* Traitements    
C*------------

      DO 10 J=0,NBMU

C* Cas J=0 : reflexion du faisceau solaire direct d'incidence rmu(0) < 0.
C* Cas J>0 : reflexion d'un rayon d'incidence -rmu(j) < 0.         
	 IF(J.EQ.0) THEN
	    MU = -RMU(0)
	 ELSE
	    MU = RMU(J)
	 ENDIF   
      
         IND2 = IND_SURF*IND_SURF
	 MU2  = MU*MU	 
	 X = DSQRT(IND2-1.D+00+MU2)
	 RL = (IND2*MU-X)/(IND2*MU+X)
	 RR = (MU-X)/(MU+X)
	 
	 F11(J) = (RL*RL+RR*RR)/2.
	 
	 IF(IPOLAR.EQ.1) THEN
	    F12(J) = (RL*RL-RR*RR)/2.
	    F33(J) =  RL*RR
	 ELSE
            F12(J) = 0.D+00
	    F33(J) = 0.D+00
	 ENDIF
	  
	    	     
10    CONTINUE

      END      !FIN DE LA PROCEDURE SOS_MAT_FRESNEL_PLAN_REFL
       
       
       
       
C*===========================================================================
C PROCEDURE: SOS_NOYAUX
C ==========
C      Cette procedure calcule les elements de la matrice de phase
C      des aerosols a l'ordre IS du developpement en series de Fourier.
C      Les calculs utilisent les polynomes et fonctions de Legendre.
C      Ces fonctions sont sorties pour un ordre L=2. Elles sont utiles
C      pour calculer les elements de la matrice de phase des molecules.
C
C      La matrice suivante P(mu,mup) s'ecrit :
C        BP(mu,mup)    GR(mu,mup)    -GT(mu,mup)
C        GR(mup,mu)    ARR(mu,mup)   -ART(mup,mu)
C       -GT(mup,mu)   -ART(mu,mup)    ATT(mu,mup)
C
C      On calcule BP(mu,mup), GR(mu,mup), GT(mu,mup), ARR(mu,mup),
C      ART(mu,mup) et ATT(mu,mup) pour -1 < mu < 1 et pour -1 < mup < 1.
C
C
C Description des parametres
C -------------------------- 
C
C     IS  (I4)  (E)   Ordre de la decomposition en series de Fourier.
C     NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C     RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)  Tableau des mu de Gauss.
C     OS_NB (I4) : (E) Ordre max des developpements des fonctions de phase 
C		       en Polynomes de Legendre a utiliser.
C
C     ALPHA(0:SOS_OS_NB_MAX) (double)  (E)    Parametres alpha(k).
C     BETA(0:SOS_OS_NB_MAX)  (double)  (E)    Parametres beta(k).
C     GAMMA(0:SOS_OS_NB_MAX) (double)  (E)    Parametres gamma(k).
C     ZETA(0:SOS_OS_NB_MAX)  (double)  (E)    Parametres zeta(k).
C
C     XPL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (S) 
C      		        Fonction de Legendre PSL à l'ordre L=2
C   		        et S=IS en fonction de l'angle mu.
C     XRL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (S)
C      		        Fonction RSL a l'ordre L=2 et S=IS.
C     XTL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (S)
C      		        Fonction TSL a l'ordre L=2 et S=IS
C	             Fonctions de Legendre et fonctions derivees pour le calcul
C	             de la matrice de phase.
C
C	
C     BP(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (S) 
C                                                        Element P11.
C     GR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (S) 
C                                                        Element P12 ou P21.
C     GT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (S) 
C                                                        Element -P13 ou -P31.
C     ARR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (S) 
C                                                        Element P22.
C     ART(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (S)
C                                                        Element -P23 ou -P32.
C     ATT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (S) 
C                                                        Element P33.
C             Element Pij de la matrice de phase des aerosols pour un ordre
C             IS de la decomposition en series de Fourier.
C             Soit mup la direction d'incidence ( -1 < mup < 1)
C             et mu celle de diffusion ( -1 < mu < 1),
C             la diffusion est exprimee par Pij(mu,mup).
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
      SUBROUTINE SOS_NOYAUX(IS,NBMU,RMU,OS_NB,ALPHA,BETA,GAMMA,ZETA,
     &                  XPL,XRL,XTL,BP,GR,GT,ARR,ART,ATT)
     
      IMPLICIT NONE


C* Definition des variables    
C*-----------------------------------------------------------------------

      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				!  Tableau des mu de Gauss.

      DOUBLE PRECISION ALPHA(0:SOS_OS_NB_MAX)	!  Parametres alpha(k).
      DOUBLE PRECISION BETA(0:SOS_OS_NB_MAX)	!  Parametres beta(k).
      DOUBLE PRECISION GAMMA(0:SOS_OS_NB_MAX)	!  Parametres gamma(k).
      DOUBLE PRECISION ZETA(0:SOS_OS_NB_MAX)	!  Parametres zeta(k).

	! Fonctions de Legendre et fonctions derivees pour le calcul
	! de la matrice de phase.
      DOUBLE PRECISION XPL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				!   Fonction de Legendre PSL à l'ordre L=2
      				!   et S=IS en fonction de l'angle mu.
      DOUBLE PRECISION XRL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				!   Fonction RSL a l'ordre L=2 et S=IS.
      DOUBLE PRECISION XTL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				!   Fonction TSL a l'ordre L=2 et S=IS
							
	! Element Pij de la matrice de phase des aerosols pour un ordre
	! IS de la decomposition en series de Fourier.
	! Soit mup la direction d'incidence (-1 < mup < 1)
	! et mu celle de diffusion (-1 < mu < 1),
	! la diffusion est exprimee par Pij(mu,mup).	
      DOUBLE PRECISION 
     &    BP(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! P12 ou P21.
      DOUBLE PRECISION 
     &    GT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! -P13 ou -P31.
      DOUBLE PRECISION 
     &    ARR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! P22.
      DOUBLE PRECISION 
     &    ART(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! -P23 ou -P32.
      DOUBLE PRECISION 
     &    ATT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! P33.



      DOUBLE PRECISION 
     &    PSL(-1:SOS_OS_NB_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Fonction de Legendre PSL pour l'ordre S=IS
      				! en fonction de l'ordre L et de l'angle mu.
      DOUBLE PRECISION 
     &    RSL(-1:SOS_OS_NB_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Fonction RSL.
      DOUBLE PRECISION 
     &    TSL(-1:SOS_OS_NB_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				  ! Fonction TSL.

      DOUBLE PRECISION C	! Valeur ponctuelle du mu de Gauss.

      DOUBLE PRECISION A,B,D,E,F     ! Coefficients pour la recurrence 
      				     ! sur PSL, RSL et TSL.

      DOUBLE PRECISION X	! Variable intermediaire de calcul.
      DOUBLE PRECISION XX	! Variable intermediaire de calcul.
      DOUBLE PRECISION YY	! Variable intermediaire de calcul.

      DOUBLE PRECISION SBP	! Variable de sommation de L = IS a 
                                ! OS_NB pour BP.
      DOUBLE PRECISION SATT	! Variable de sommation de L = IS a 
                                ! OS_NB pour SATT.
      DOUBLE PRECISION SARR	! Variable de sommation de L = IS a 
                                ! OS_NB pour SARR.
      DOUBLE PRECISION SGR	! Variable de sommation de L = IS a 
                                ! OS_NB pour SGR.
      DOUBLE PRECISION SGT	! Variable de sommation de L = IS a 
                                ! OS_NB pour SGT.
      DOUBLE PRECISION SART	! Variable de sommation de L = IS a 
                                ! OS_NB pour SART.
      DOUBLE PRECISION R1		
      DOUBLE PRECISION R2	

      DOUBLE PRECISION RAC3	! Racine de 3.
      DOUBLE PRECISION X26	! 2*racine(6)

      INTEGER*4 OS_NB 	! Ordre max des developpements des fonctions de phase 
			! en Polynomes de Legendre.
      INTEGER IS	! Ordre de la decomposition en series de Fourier.

      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations
      INTEGER*4 I
      INTEGER*4 K
      INTEGER*4 IG		! Variable a 1 ou -1.
      INTEGER*4 L
      INTEGER*4 LP		! Valeur L+1
      INTEGER*4 LM		! Valeur L-1
      INTEGER*4 J	! Indice du tableau de RMU.


C* Initialisation 
C----------------
      RAC3=DSQRT(3.D+00)
      X26=2.*DSQRT(6.D+00)

C-----------------------------------------------------
C* Initialisation des fonctions d'ordre IS pour les 
C* ordres L = 0, 1 et 2.
C-----------------------------------------------------

C* Initialisation pour l'ordre IS = 0
C-------------------------------------
      IF (IS.NE.0) GO TO 700
      
      DO 25 J=0,NBMU
         C=RMU(J)
         PSL(0,-J)=1.
         PSL(0,J)=1.
         PSL(1,J)=C
         PSL(1,-J)=-C
         X=(3.*C*C-1.)*0.5
         PSL(2,-J)=X
         PSL(2,J)=X
         RSL(1,J)=0.
         RSL(1,-J)=0.
         X=3.*(1.-C*C)/X26
         RSL(2,-J)=X
         RSL(2,J)=X
         TSL(1,J)=0.
         TSL(1,-J)=0.
         TSL(2,J)=0.
         TSL(2,-J)=0.
   25 CONTINUE
   
      PSL(1,0)=RMU(0)
      RSL(1,0)=0.
      
      GOTO 501

C* Initialisation pour l'ordre IS = 1
C-------------------------------------
  700 IF (IS.NE.1) GO TO 701
  
      DO 26 J=0,NBMU
         C=RMU(J)
         X=1.-C*C
         PSL(0,J)=0.
         PSL(0,-J)=0.
         PSL(1,-J)=DSQRT(X*0.5)
         PSL(1,J)=DSQRT(X*0.5)
         PSL(2,J)=C*PSL(1,J)*RAC3
         PSL(2,-J)=-PSL(2,J)
         RSL(1,-J)=0.
         RSL(1,J)=0.
         RSL(2,J)=-C*DSQRT(X)*0.5
         RSL(2,-J)=-RSL(2,J)
         TSL(1,-J)=0.
         TSL(1,J)=0.
         TSL(2,J)=-DSQRT(X)*0.5
         TSL(2,-J)=-DSQRT(X)*0.5
   26 CONTINUE
   
      PSL(2,0)=-PSL(2,0)
      RSL(2,0)=-RSL(2,0)
      RSL(1,0)=0.
      TSL(1,0)=0.
      
      GOTO 501

C* Initialisation pour l'ordre IS >= 2
C----------------------------------------
  701 A=1.
  
      DO 27 I=1,IS
        X=I
        A=A*DSQRT((I+IS)/X)*0.5
 27   CONTINUE
 
      B=A*DSQRT(IS/(IS+1.D+00))*DSQRT((IS-1.D+00)/(IS+2.)) 
      
      DO 28 J=0,NBMU
         C=RMU(J)
         XX=1.-C*C
         YY=IS*0.5-1.
         PSL(IS-1,J)=0.
         RSL(IS-1,J)=0.
         TSL(IS-1,J)=0.
         X=A*XX**(IS*0.5)
         PSL(IS,-J)=X
         PSL(IS,J)=X
         X=B*(1.+C*C)*XX**YY
         RSL(IS,-J)=X
         RSL(IS,J)=X
         X=2.*B*C*XX**YY
         TSL(IS,-J)=-X
         TSL(IS,J)=X
   28 CONTINUE


C* Calcul des fonctions PSL(L,mu), RSL(L,mu) et TSL(L,mu) 
C* par recurrence sur l'ordre L pour IS fixe.
C------------------------------------------------------
  501 K=2

      IF (IS.GT.2) K=IS
      
      IF (K.EQ.OS_NB) GOTO 502
      
      IG=-1
      IF (IS.EQ.1) IG=1

      DO 30 L=K,(OS_NB-1)
         LP=L+1
         LM=L-1

C* Coefficients pour la recurrence sur PSL
         A=(2*L+1.)/DSQRT((L+IS+1.D+00)*(L-IS+1.))
         B=DSQRT(DFLOAT((L+IS)*(L-IS)))/(2.*L+1.)

C* Coefficients pour la recurrence sur RSL et TSL
         D=(L+1.)*(2*L+1.)/DSQRT((L+3.D+00)*(L-1.)*(L+IS+1.)  
     &     *(L-IS+1.))
         E=DSQRT((L+2.D+00)*(L-2.)*(L+IS)*(L-IS))/(L*(2.*L+1.))
         F=2.*IS/(L*(L+1.))

C* Application de la recurrence pour chaque mu de Gauss
C------------------------------------------------------
         DO 31 J=0,NBMU
            C=RMU(J)
            X=A*(C*PSL(L,J)-B*PSL(LM,J))
            PSL(LP,J)=X
            X=D*(C*RSL(L,J)-F*TSL(L,J)-E*RSL(LM,J))
            RSL(LP,J)=X
            X=D*(C*TSL(L,J)-F*RSL(L,J)-E*TSL(LM,J))
            TSL(LP,J)=X
	    
            IF(J.EQ.0) GO TO 31
	    
            PSL(LP,-J)=IG*PSL(LP,J)
            RSL(LP,-J)=IG*RSL(LP,J)
            TSL(LP,-J)=-IG*TSL(LP,J)
   31    CONTINUE

         IG=-IG
   30 CONTINUE

  502 CONTINUE


C* Stockage de PSL(L=2,mu), RSL(L=2,mu) et TSL(L=2,mu)
C------------------------------------------------------
      DO 1005 J=-NBMU,NBMU
         XPL(J)=PSL(2,J)
         XRL(J)=RSL(2,J)
         XTL(J)=TSL(2,J)
 1005 CONTINUE


C------------------------------------------------------
C* Calcul des elements de matrice de phase des aerosols
C------------------------------------------------------

C* Boucles sur les directions d'incidence (K = mup, 
C* avec -1 < mup < 1) et d'emergence (J = mu > 0)
C------------------------------------------------------
      DO 32 J=-NBMU,NBMU
         DO 32 K=-NBMU,NBMU
            SBP=0.
            SATT=0.
            SARR=0.
            SGR=0.
            SGT=0.
            SART=0.
	    
            IF(IS.GT.OS_NB) GOTO 1

C* Boucle de sommation sur L = IS a OS_NB.
C-----------------------------------------------
           DO 33 L=IS,OS_NB
              R1=TSL(L,J)*TSL(L,K)
              R2=RSL(L,J)*RSL(L,K)
              SBP=SBP+BETA(L)*PSL(L,J)*PSL(L,K)
              SATT=SATT+ALPHA(L)*R1+ZETA(L)*R2
              SARR=SARR+ZETA(L)*R1+ALPHA(L)*R2
              SGR=SGR+GAMMA(L)*PSL(L,J)*RSL(L,K)
              SGT=SGT+GAMMA(L)*PSL(L,J)*TSL(L,K)
              SART=SART+ALPHA(L)*RSL(L,K)*TSL(L,J)+
     &             ZETA(L)*RSL(L,J)*TSL(L,K)
  33       CONTINUE
  
 1         CONTINUE

            BP(J,K)=SBP
            ATT(J,K)=SATT
            ARR(J,K)=SARR
            GR(J,K)=SGR
            GT(J,K)=SGT
            ART(J,K)=SART
   32 CONTINUE
   
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_NOYAUX









C*===========================================================================
C PROCEDURE: SOS_INTEGR_EPOPT
C ==========
C      Cette procedure realise l'integration sur l'epaisseur optique
C      du champ de rayonnement :
C        - champ montant    : Integration du sol au TOA.
C        - champ descendant : Integration du TOA au sol.
C
C      L'integration est faite de proche en proche du niveau de depart
C      au niveau final.
C      Ceci implique dans la resolution de l'Equation de Transfert Radiatif
C      de faire des integrations sur l'epaisseur optique entre deux niveaux
C      consecutifs. Ces integrations font intervenir les fonctions sources
C      de diffusion (pour I, Q et U - en fonction du niveau d'epaisseur
C      optique et de la direction de propagation) pour les ordres en cours
C      de diffusion et de decomposition en series de Fourier.
C      L'integration entre deux niveaux est faite par linearisation de la
C      fonction source sur la couche.
C
C       Les couches atmospheriques sont supposees homogenes.
C
C
C
C Description des parametres
C -------------------------- 
C
C     NBMU (I4)    (E)  Nombre d'angles (positifs) effectivement utiles
C     RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double) (E) Tableau des mu de Gauss.
C     H(0:SOS_OS_NT) (double) (E)  Epaisseur optique a chaque niveau.
C
C     I2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double) (E) 
C                                                  Fonction source pour I
C     Q2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double) (E) 
C                                                  Fonction source pour Q
C     U2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double) (E) 
C                                                  Fonction source pour U
C                          Vecteur fonction source en fonction du niveau du 
C                          profil et de la direction.
C    	
C     I1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double) (S) Champ I.
C     Q1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double) (S) Champ Q.
C     U1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double) (S) Champ U.
C                          Champ de luminance  en fonction du niveau
C                          du profil et de la direction. 
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
      SUBROUTINE SOS_INTEGR_EPOPT(NBMU,RMU,H,I2,Q2,U2,I1,Q1,U1)
      
      IMPLICIT NONE
      

C* Definition des variables     
C*------------------------------------------------------------------

C* Directions de propagation :
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				!  Tableau des mu de Gauss.

C* Profil d'epaisseur optique :
      DOUBLE PRECISION H(0:SOS_OS_NT)	!  Epaisseur optique pour chaque niveau.
     	
C* Vecteur fonction source en fonction du niveau du profil et de la direction.
      DOUBLE PRECISION I2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      					!  Fonction source pour I.
      DOUBLE PRECISION Q2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      					!  Fonction source pour Q.
      DOUBLE PRECISION U2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      					!  Fonction source pour U.

C* Champ de luminance en fonction du niveau du profil et de la direction.       	
      DOUBLE PRECISION I1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	!  Champ I.
      DOUBLE PRECISION Q1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	!  Champ Q.
      DOUBLE PRECISION U1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	!  Champ U.

      DOUBLE PRECISION RMUK 	! Valeur K du tableau des mu de Gauss 

      DOUBLE PRECISION ZI1	! Valeur intermediaire pour le calcul de I1
      				! a l'ordre IG.
      DOUBLE PRECISION ZQ1	! Valeur intermediaire pour le calcul de Q1.
      DOUBLE PRECISION ZU1	! Valeur intermediaire pour le calcul de U1.

      DOUBLE PRECISION DTAU   	! Epaisseur optique d'une couche.

      DOUBLE PRECISION ATT	! Attenuation entre deux niveaux consecutifs.

      DOUBLE PRECISION MATT	! 1 - att.

      DOUBLE PRECISION ATTDTAU	! att * dtau

      DOUBLE PRECISION A,B	! Coefficients d'interpolation lineaire

      INTEGER*4 I,JJ		! Indices de niveau.
      INTEGER*4 K		! Indice de direction.
      INTEGER*4 NBMU            ! Nombre d'angles effectifs des simulations
		
C------------------------------------------------------	
C* Integration du champ de luminance montant (mu > 0) du sol au TOA.
C*
C* Integration de proche en proche du sol au TOA par linearisation
C* de la fonction source entre deux niveaux successifs.
C------------------------------------------------------
      DO 108 K=1,NBMU

         RMUK=RMU(K)

C* Initialisation par la valeur au sol.
         ZI1=I1(SOS_OS_NT,K)
         ZQ1=Q1(SOS_OS_NT,K)
         ZU1=U1(SOS_OS_NT,K)
		
         DO 108 I=SOS_OS_NT-1,0,-1
            JJ=I+1
            DTAU = H(JJ)-H(I)
            ATT = DEXP(-DTAU/RMUK)
            MATT = 1.D+00 - ATT
            ATTDTAU =  ATT * DTAU

            B = I2(I,K)
            A = (I2(JJ,K) - B) / DTAU
            ZI1 = ZI1*ATT + MATT*(A*RMUK+B) - A*ATTDTAU
            I1(I,K)=ZI1

            B = Q2(I,K)
            A = (Q2(JJ,K) - B) / DTAU
            ZQ1 = ZQ1*ATT + MATT*(A*RMUK+B) - A*ATTDTAU
            Q1(I,K)=ZQ1

            B = U2(I,K)
            A = (U2(JJ,K) - B) / DTAU
            ZU1 = ZU1*ATT + MATT*(A*RMUK+B) - A*ATTDTAU
            U1(I,K)=ZU1
108   CONTINUE        !Boucles sur I et K


C------------------------------------------------------
C* Integration du champ de luminance descendant (mu < 0) du TOA au sol
C*
C* Integration de proche en proche du TOA au sol par linearisation
C* de la fonction source entre deux niveaux successifs.
C* Initialisation par un champ diffus descendant nul au TOA.
C------------------------------------------------------
      DO 109 K=-NBMU,-1

         RMUK=RMU(K)

C* Initialisation par la valeur au TOA.
         I1(0,K)=0.
         Q1(0,K)=0.
         U1(0,K)=0.
         ZI1=0.D+00
         ZQ1=0.D+00
         ZU1=0.D+00

         DO 109 I=1,SOS_OS_NT
            JJ=I-1
            DTAU = H(I)-H(JJ)
            ATT = DEXP(DTAU/RMUK)
            MATT = 1.D+00 - ATT
            ATTDTAU =  ATT * DTAU

            B = I2(I,K)
            A = (B - I2(JJ,K)) / DTAU
            ZI1 = ZI1*ATT + MATT*(A*RMUK+B) + A*ATTDTAU
            I1(I,K)=ZI1

            B = Q2(I,K)
            A = (B - Q2(JJ,K)) / DTAU
            ZQ1 = ZQ1*ATT + MATT*(A*RMUK+B) + A*ATTDTAU
            Q1(I,K)=ZQ1

            B = U2(I,K)
            A = (B- U2(JJ,K)) / DTAU
            ZU1 = ZU1*ATT + MATT*(A*RMUK+B) + A*ATTDTAU
            U1(I,K)=ZU1
  109 CONTINUE  !Boucles sur I et K

      RETURN
      END  	!FIN DE LA PROCEDURE SOS_INTEGR_EPOPT






C*============================================================================
C PROCEDURE: SOS_FSOURCE_ORDRE1
C ==========
C      Cette procedure calcule le vecteur fonction source de la diffusion 
C      primaire, a partir du champ de rayonnement solaire incident et des
C      proprietes de diffusion de l'atmosphere (fonction de phase).
C      
C      Calcul pour chaque niveau du profil et pour chaque direction de 
C      propagation.
C
C
C Description des parametres
C -------------------------- 
C
C       IS   (I4)    (E)  Ordre en cours du developpement en series de Fourier.
C       NBMU (I4)    (E)  Nombre d'angles (positifs) effectivement utiles
C       JK   (I4)    (E)  Numero de stockage de mus dans le tableau RMU
C                        (0 par defaut).
C       XDEL(0:SOS_OS_NT)  (double)  (E)  Pourcentage d'aerosols par niveau.
C       YDEL(0:SOS_OS_NT)  (double)  (E)  Pourcentage de molecules par niveau.
C
C       BETA0	  (double)  (E) 
C       BETA2     (double)  (E) 	.
C       GAMMA2	  (double)  (E) 	
C                   Coefficients du developpement de la matrice de phase
C                   moleculaire en fonctions de Legendre
C
C       XPL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 	
C      		     Fonction de Legendre PSL à l'ordre L=2 et S=IS
C      		     en fonction de l'angle mu.     		
C       XRL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 	
C     		     Fonction RSL a l'ordre L=2 et S=IS en fonction de mu.
C       XTL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 	
C     		     Fonction TSL a l'ordre L=2 et S=IS en fonction de mu.
C	
C       BP(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                                                        Element P11.
C       GR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                                                        Element P12 ou P21.
C       GT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                                                        Element -P13 ou -P31.
C               Element Pij de la matrice de phase des aerosols pour un ordre
C               IS de la decomposition en series de Fourier.
C               Soit mup la direction d'incidence (-1 < mup < 1)
C               et mu celle de diffusion (-1 < mu < 1),
C               la diffusion est exprimee par Pij(mu,mup).
C
C       CH(0:SOS_OS_NT)  (double)  (E) 	Attenuation du faisceau solaire 
C                                       direct/4.
C
C       I2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (S)	
C       Q2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (S)	
C       U2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (S)	
C	         Vecteur fonction source I,Q,U en fonction du niveau du profil
C	         et de la direction.
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
C==============================================================================
      SUBROUTINE SOS_FSOURCE_ORDRE1(IS,NBMU,JK,XDEL,YDEL,
     &                              BETA0,BETA2,GAMMA2,
     &                              XPL,XRL,XTL,BP,GR,GT,CH,I2,Q2,U2)
      
      IMPLICIT NONE
      

C* Definition des variables     
C*------------------------------------------------------------------

      DOUBLE PRECISION XDEL(0:SOS_OS_NT)   !Pourcentage d'aerosols par niveau.
      DOUBLE PRECISION YDEL(0:SOS_OS_NT)   !Pourcentage de molecules par niveau.

      DOUBLE PRECISION BETA0	! Coefficients du developpement de la matrice 
      DOUBLE PRECISION BETA2	! de phase moleculaire en fonctions de Legendre.
      DOUBLE PRECISION GAMMA2

      DOUBLE PRECISION XPL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			      ! Fonction de Legendre PSL a l'ordre L=2 et S=IS
      			      ! en fonction de l'angle mu.
      		
      DOUBLE PRECISION XRL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			! Fonction RSL a l'ordre L=2 et S=IS en fonction de mu
      DOUBLE PRECISION XTL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			! Fonction TSL a l'ordre L=2 et S=IS en fonction de mu

	! Element Pij de la matrice de phase des aerosols pour un ordre
	! IS de la decomposition en series de Fourier.
	! Soit mup la direction d'incidence (-1 < mup < 1)
	! et mu celle de diffusion (-1 < mu < 1),
	! la diffusion est exprimee par Pij(mu,mup).	
      DOUBLE PRECISION 
     &    BP(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! P12 ou P21.
      DOUBLE PRECISION 
     &    GT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! -P13 ou -P31.
    

      DOUBLE PRECISION CH(0:SOS_OS_NT)	! Attenuation du faisceau solaire 
                                        ! direct / 4.

	! Vecteur fonction source en fonction du niveau du profil
	! et de la direction.
      DOUBLE PRECISION I2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour I.
      DOUBLE PRECISION Q2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour Q.
      DOUBLE PRECISION U2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour U.

      DOUBLE PRECISION SPL   ! Valeur de la fonction de XPL pour mu = mus.
      DOUBLE PRECISION SA1   ! Element P11 de la matrice de phase moleculaire.
      DOUBLE PRECISION SB1   ! Element P21 de la matrice de phase moleculaire.
      DOUBLE PRECISION SC1   ! Element P31 de la matrice de phase moleculaire.
      DOUBLE PRECISION SA2   ! Element P11 de la matrice de phase des aerosols.
      DOUBLE PRECISION SB2   ! Element P21 de la matrice de phase des aerosols.
      DOUBLE PRECISION SC2   ! Element P31 de la matrice de phase des aerosols.
 
      DOUBLE PRECISION PCRAY	! Valeur du pourcentage de molecules.
      DOUBLE PRECISION PCAER	! Valeur du pourcentage d'aerosols.

      DOUBLE PRECISION ATTDIR	! Valeur ponctuelle du tableau CH.

          
      INTEGER*4 IS	! Ordre en cours du developpement en series de Fourier.
      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations
      
      INTEGER JK	! Numero de stockage de mus dans le tableau RMU
      			! (0 par defaut).

      INTEGER*4 J,K


C-------------------------------------------------------------      
C* Calcul des elements P11 = SAi , P21 = SBi  et P31 = SCi
C* de la matrice de phase, selon l'ordre IS du developpement 
C* en series de Fourier (i = 1 : molecules, i = 2 : aerosols).
C* La direction d'incidence est la direction d'eclairement 
C* solaire (mus <=> JK = 0).
C-------------------------------------------------------------      


C* Boucle sur les directions de diffusion.
C------------------------------------------
      DO 100 J=-NBMU,NBMU

        IF (IS-2) 200,200,201
      
C* Calculs pour IS <= 2.
C-----------------------
 200    SPL=XPL(JK)
        SA1=BETA0+BETA2*XPL(J)*SPL
        SA2=BP(JK,J)
        SB1=GAMMA2*XRL(J)*SPL
        SB2=GR(JK,J)
        SC1=GAMMA2*XTL(J)*SPL
        SC2=GT(JK,J)

        GOTO 202

C* Calculs pour IS > 2. La matrice de phase des 
C* molecules est nulle pour IS > 2
C------------------------------------------------
 201    SA2=BP(JK,J)
        SA1=0.
        SB2=GR(JK,J)
        SB1=0.
        SC2=GT(JK,J)
        SC1=0.

202     CONTINUE


C------------------------------------------------------------------------
C* Calculs de la fonction source de diffusion primaire pour chaque niveau
C* du profil et chaque direction.
C------------------------------------------------------------------------
        DO 101 K=0,SOS_OS_NT
           ATTDIR=CH(K)
           PCRAY=YDEL(K)
           PCAER=XDEL(K)
           I2(K,J)=ATTDIR*(SA2*PCAER+SA1*PCRAY)
           Q2(K,J)=ATTDIR*(SB2*PCAER+SB1*PCRAY)
           U2(K,J)=-ATTDIR*(SC2*PCAER+SC1*PCRAY)
  101    CONTINUE
  
  100 CONTINUE
  

      END  	!FIN DE LA PROCEDURE SOS_FSOURCE_ORDRE1







C*============================================================================
C PROCEDURE: SOS_FSOURCE_ORDREIG
C ==========
C      Cette procedure calcule le vecteur fonction source de la diffusion
C      d'ordre IG>1, a partir du champ de rayonnement d'ordre IG-1 et des
C      proprietes de diffusion de l'atmosphere (fonction de phase).
C
C      Elle realise le calcul des elements de la matrice de phase, 
C      selon l'ordre IS du developpement en series de Fourier.
C      On construit indirectement la matrice suivante P(mu,mup) :
C         BP(mu,mup)    GR(mu,mup)    -GT(mu,mup)
C         GR(mup,mu)    ARR(mu,mup)   -ART(mup,mu)
C        -GT(mup,mu)   -ART(mu,mup)    ATT(mu,mup)
C 
C      Les termes reellement calcules sont equivalents en suivant les regles
C      de permutation des indices. 
C	
C      Une integration spatiale est appliquee sur les directions d'incidence.
C
C      Calcul pour chaque niveau du profil et pour chaque direction de 
C      propagation.
C
C
C
C Description des parametres
C -------------------------- 
C   
C       IS   (I4)    (E)  Ordre en cours du developpement en series de Fourier.
C       NBMU (I4)    (E)  Nombre d'angles (positifs) effectivement utiles
C       XDEL(0:SOS_OS_NT)  (double)  (E)  Pourcentage d'aerosols par niveau.
C       YDEL(0:SOS_OS_NT)  (double)  (E)  Pourcentage de molecules par niveau.
C
C       BETA0	  (double)  (E) 
C       BETA2     (double)  (E) 	.
C       GAMMA2	  (double)  (E) 
C       ALPHA2    (double)  (E) 	
C                   Coefficients du developpement de la matrice de phase
C                   moleculaire en fonctions de Legendre
C
C       XPL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 	
C      		     Fonction de Legendre PSL à l'ordre L=2 et S=IS
C      		     en fonction de l'angle mu.     		
C       XRL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 	
C     		     Fonction RSL a l'ordre L=2 et S=IS en fonction de mu.
C       XTL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 	
C     		     Fonction TSL a l'ordre L=2 et S=IS en fonction de mu.
C
C      	
C       I1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)    Champ I.
C       Q1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)    Champ Q.
C       U1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)    Champ U.
C         	     Champ de luminance d'ordre  IG en fonction du 
C                    niveau du profil et de la direction. 
C	
C       BP(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                                                        Element P11.
C       GR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                                                        Element P12 ou P21.
C       GT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C                                                        Element -P13 ou -P31.
C       ARR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (E) 
C                                                        Element P22.
C       ART(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (E) 
C                                                        Element -P23 ou -P32.
C       ATT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (E)
C                                                        Element P33.
C               Element Pij de la matrice de phase des aerosols pour un ordre
C               IS de la decomposition en series de Fourier.
C               Soit mup la direction d'incidence ( -1 < mup < 1)
C               et mu celle de diffusion ( -1 < mu < 1),
C               la diffusion est exprimee par Pij(mu,mup).
C
C       GA(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) Tableau des poids de Gauss.
C
C       I2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (S)	
C       Q2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (S)	
C       U2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (S)	
C	         Vecteur fonction source I,Q,U en fonction du niveau du profil
C	         et de la direction.
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
      SUBROUTINE SOS_FSOURCE_ORDREIG(IS,NBMU,XDEL,YDEL,
     &                               BETA0,BETA2,GAMMA2,ALPHA2,
     &                               XPL,XRL,XTL,I1,Q1,U1,
     &                               BP,GR,GT,ARR,ART,ATT,GA,
     &                               I2,Q2,U2)
      
      IMPLICIT NONE
      

C* Definition des variables     
C*------------------------------------------------------------------

      DOUBLE PRECISION XDEL(0:SOS_OS_NT)  ! Pourcentage d'aerosols par niveau.
      DOUBLE PRECISION YDEL(0:SOS_OS_NT)  ! Pourcentage de molecules par niveau.

      DOUBLE PRECISION BETA0	! Coefficients du developpement de la matrice 
      DOUBLE PRECISION BETA2	! de phase moleculaire en fonctions de Legendre.
      DOUBLE PRECISION GAMMA2	
      DOUBLE PRECISION ALPHA2	

      DOUBLE PRECISION XPL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			      ! Fonction de Legendre PSL à l'ordre L=2 et S=IS
      			      ! en fonction de l'angle mu.     		
      DOUBLE PRECISION XRL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Fonction RSL a l'ordre L=2 et S=IS en fonction de mu.
      DOUBLE PRECISION XTL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Fonction TSL a l'ordre L=2 et S=IS en fonction de mu.


      	! Champ de luminance d'ordre IG-1 en fonction du niveau du profil
	! et de la direction.       	
      DOUBLE PRECISION I1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Champ I.
      DOUBLE PRECISION Q1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Champ Q.
      DOUBLE PRECISION U1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Champ U.

	! Vecteur fonction source en fonction du niveau du profil
	! et de la direction.
      DOUBLE PRECISION I2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour I.
      DOUBLE PRECISION Q2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour Q.
      DOUBLE PRECISION U2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour U.


	! Element Pij de la matrice de phase des aerosols pour un ordre
	! IS de la decomposition en series de Fourier.
	! Soit mup la direction d'incidence (-1 < mup < 1)
	! et mu celle de diffusion (-1 < mu < 1),
	! la diffusion est exprimee par Pij(mu,mup).	
      DOUBLE PRECISION 
     &    BP(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! P12 ou P21.
      DOUBLE PRECISION 
     &    GT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! -P13 ou -P31.
      DOUBLE PRECISION 
     &    ARR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! P22.
      DOUBLE PRECISION 
     &    ART(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! -P23 ou -P32.
      DOUBLE PRECISION 
     &    ATT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! P33.


      DOUBLE PRECISION GA(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! Tableau des poids de Gauss

      DOUBLE PRECISION XPK	! Valeur de XPL pour les directions de
      				! diffusion mu > 0 (indice K).	
      DOUBLE PRECISION XRK	! Valeur de XRL pour mu > 0.
      DOUBLE PRECISION XTK	! Valeur de XTL pour mu > 0.
      DOUBLE PRECISION YPK	! Valeur de XPL pour -mu.	
      DOUBLE PRECISION YRK	! Valeur de XRL pour -mu.
      DOUBLE PRECISION YTK	! Valeur de XTL pour -mu.	
      	
      DOUBLE PRECISION XPJ	! Valeur de XPL pour les directions incidentes
      				! mu = mup > 0 (indice J).
      DOUBLE PRECISION XRJ	! Valeur de XRL pour mu = mup > 0.
      DOUBLE PRECISION XTJ	! Valeur de XTL pour mu = mup > 0.
      DOUBLE PRECISION YPJ	! Valeur de XPL pour mu = -mup.	
      DOUBLE PRECISION YRJ	! Valeur de XRL pour mu = -mup.
      DOUBLE PRECISION YTJ	! Valeur de XTL pour mu = -mup.		

      DOUBLE PRECISION XI1	! Valeur de I du champ montant d'ordre IG-1.
      DOUBLE PRECISION XQ1	! Valeur de Q du champ montant d'ordre IG-1.
      DOUBLE PRECISION XU1	! Valeur de U du champ montant d'ordre IG-1.
      DOUBLE PRECISION XI2	! Valeur de I du champ descendant d'ordre IG-1.
      DOUBLE PRECISION XQ2	! Valeur de Q du champ descendant d'ordre IG-1.
      DOUBLE PRECISION XU2	! Valeur de U du champ descendant d'ordre IG-1.

      DOUBLE PRECISION II1	! Valeur intermediaire pour le calcul de I2 
      				! (mu < 0).
      DOUBLE PRECISION QQ1	! Valeur intermediaire pour le calcul de Q2 
      				! (mu < 0).
      DOUBLE PRECISION UU1	! Valeur intermediaire pour le calcul de U2 
      				! (mu < 0).
      DOUBLE PRECISION II2 	! Valeur intermediaire pour le calcul de I2 
      				! (mu > 0).
      DOUBLE PRECISION QQ2 	! Valeur intermediaire pour le calcul de Q2 
      				! (mu > 0).
      DOUBLE PRECISION UU2 	! Valeur intermediaire pour le calcul de U2 
      				! (mu > 0).

      DOUBLE PRECISION Z   ! Valeur ponctuelle du tableau des poids de Gauss.

      DOUBLE PRECISION PCRAY	! Valeur du pourcentage de molecules.
      DOUBLE PRECISION PCAER	! Valeur du pourcentage d'aerosols.
 
      DOUBLE PRECISION BPJK	! Valeur des elements de la matrice de 
      DOUBLE PRECISION BPJMK    ! phase du melange molecules + aerosols.
      DOUBLE PRECISION GTJMK
      DOUBLE PRECISION GTJK
      DOUBLE PRECISION GTKJ
      DOUBLE PRECISION GTKMJ
      DOUBLE PRECISION GRJK
      DOUBLE PRECISION GRKJ
      DOUBLE PRECISION GRKMJ
      DOUBLE PRECISION GRJMK
      DOUBLE PRECISION ARRJK
      DOUBLE PRECISION ARRJMK
      DOUBLE PRECISION ARTJK
      DOUBLE PRECISION ARTJMK
      DOUBLE PRECISION ARTKJ
      DOUBLE PRECISION ARTKMJ
      DOUBLE PRECISION ATTJK
      DOUBLE PRECISION ATTJMK

      INTEGER*4 IS	! Ordre en cours du developpement en series de Fourier.
      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations

      INTEGER*4 I
      INTEGER*4 K
      INTEGER*4 J
      


C------------------------------------------------------
C* Calcul des elements de la matrice de phase, selon 
C* l'ordre IS du developpement en series de Fourier.
C------------------------------------------------------

      IF (IS-2) 210,210,211
      
C--------------------------------
C* Calculs pour IS <= 2.
C-------------------------------

C* Boucle sur les directions de diffusion mu (mu > 0).
C------------------------------------------------------
 210  DO 455 K=1,NBMU

C* Valeur des fonctions PSL, RSL et TSL à l'ordre L=2 et S=IS
C* en fonction de l'angle mup et -mup
C* (pour le calcul de la matrice de phase moleculaire 
C*  par le programme appelant).
C------------------------------------------------------
         XPK=XPL(K)
         XRK=XRL(K)
         XTK=XTL(K)
         YPK=XPL(-K)
         YRK=XRL(-K)
         YTK=XTL(-K)

C* Boucle sur le profil.
C-----------------------
         DO 455 I=0,SOS_OS_NT
            II1=0.
            II2=0.
            QQ1=0.
            QQ2=0.
            UU1=0.
            UU2=0.
            PCAER=XDEL(I)
            PCRAY=YDEL(I)

C* Boucle sur les directions d'incidence mup.
C------------------------------------------------------
            DO 477 J=1,NBMU

C* Calcul des elements de la matrice de phase P(mu,mup) 
C* a l'ordre IS pour le melange molecules + aerosols.
C------------------------------------------------------
               XPJ=XPL(J)
               XRJ=XRL(J)
               XTJ=XTL(J)
               YPJ=XPL(-J)
               YRJ=XRL(-J)
               YTJ=XTL(-J)

               BPJK=BP(J,K)*PCAER+PCRAY*(BETA0+BETA2*XPJ*XPK)
               BPJMK=BP(J,-K)*PCAER+PCRAY*(BETA0+BETA2*XPJ*YPK)
               GTJMK=GT(J,-K)*PCAER+PCRAY*(GAMMA2*XPJ*YTK)
               GTJK=GT(J,K)*PCAER+PCRAY*(GAMMA2*XPJ*XTK)
               GTKMJ=GT(K,-J)*PCAER+PCRAY*(GAMMA2*XPK*YTJ)
               GTKJ=GT(K,J)*PCAER+PCRAY*(GAMMA2*XPK*XTJ)
               GRJK=GR(J,K)*PCAER+PCRAY*(GAMMA2*XPJ*XRK)
               GRJMK=GR(J,-K)*PCAER+PCRAY*(GAMMA2*XPJ*YRK)
               GRKJ=GR(K,J)*PCAER+PCRAY*(GAMMA2*XPK*XRJ)
               GRKMJ=GR(K,-J)*PCAER+PCRAY*(GAMMA2*XPK*YRJ)
               ARRJK=ARR(J,K)*PCAER+PCRAY*(ALPHA2*XRJ*XRK)
               ARRJMK=ARR(J,-K)*PCAER+PCRAY*(ALPHA2*XRJ*YRK)
               ARTJK=ART(J,K)*PCAER+PCRAY*(ALPHA2*XTJ*XRK)
               ARTJMK=ART(J,-K)*PCAER+PCRAY*(ALPHA2*XTJ*YRK)
               ARTKJ=ART(K,J)*PCAER+PCRAY*(ALPHA2*XTK*XRJ)
               ARTKMJ=ART(K,-J)*PCAER+PCRAY*(ALPHA2*XTK*YRJ)
               ATTJMK=ATT(J,-K)*PCAER+PCRAY*(ALPHA2*XTJ*YTK)
               ATTJK=ATT(J,K)*PCAER+PCRAY*(ALPHA2*XTJ*XTK)

C* Calculs de la fonction source d'ordre IG, pour chaque 
C* niveau du profil et chaque direction.
C------------------------------------------------------
               Z=GA(J)

C* Recuperation des valeurs de luminance d'ordre IG-1.
C------------------------------------------------------
               XI1=I1(I,J)
               XI2=I1(I,-J)
               XQ1=Q1(I,J)
               XQ2=Q1(I,-J)
               XU1=U1(I,J)
               XU2=U1(I,-J)

C* Integration spatiale sur mup.
C---------------------------------
               II2=II2+Z*(XI1*BPJK+XI2*BPJMK+XQ1*GRKJ
     &             +XQ2*GRKMJ-XU1*GTKJ-XU2*GTKMJ)
               II1=II1+Z*(XI1*BPJMK+XI2*BPJK+XQ1*GRKMJ
     &             +XQ2*GRKJ+XU1*GTKMJ+XU2*GTKJ)
               QQ2=QQ2+Z*(XI1*GRJK+XI2*GRJMK+XQ1*ARRJK
     &             +XQ2*ARRJMK+XU2*ARTJMK-XU1*ARTJK)
               QQ1=QQ1+Z*(XI1*GRJMK+XI2*GRJK+XQ1*ARRJMK
     &             +XQ2*ARRJK-XU1*ARTJMK+XU2*ARTJK)
               UU2=UU2-Z*(XI1*GTJK-XI2*GTJMK+XQ1*ARTKJ
     &             +XQ2*ARTKMJ-XU1*ATTJK-XU2*ATTJMK)
               UU1=UU1-Z*(XI1*GTJMK-XI2*GTJK-XQ1*ARTKMJ
     &             -XQ2*ARTKJ-XU1*ATTJMK-XU2*ATTJK)
 477        CONTINUE

C* Stockage des resultats
C--------------------------
            I2(I,K)=II2 * 0.5
            I2(I,-K)=II1 * 0.5
            Q2(I,K)=QQ2 * 0.5
            Q2(I,-K)=QQ1 * 0.5
            U2(I,K)=UU2 * 0.5
            U2(I,-K)=UU1 * 0.5
 455  CONTINUE   ! Boucle sur I et K
 
      GOTO 212
      
      
C-----------------------------
C* Calculs pour IS > 2.
C-----------------------------

C* Boucle sur les directions de diffusion mu (mu > 0).
C------------------------------------------------------
 211  DO 45 K=1,NBMU

C* Boucle sur le profil.
C------------------------
         DO 45 I=0,SOS_OS_NT
            II1=0.
            II2=0.
            QQ1=0.
            QQ2=0.
            UU1=0.
            UU2=0.
            PCAER=XDEL(I)

C* Boucle sur les directions d'incidence mup
C-------------------------------------------
            DO 47 J=1,NBMU

C* Calcul des elements de la matrice de phase aerosols P(mu,mup)
C* a l'ordre IS, ponderee par la proportion en aerosols.
C------------------------------------------------------
               BPJK=BP(J,K)*PCAER
               BPJMK=BP(J,-K)*PCAER
               GTJMK=GT(J,-K)*PCAER
               GTJK=GT(J,K)*PCAER
               GTKMJ=GT(K,-J)*PCAER
               GTKJ=GT(K,J)*PCAER
               GRJK=GR(J,K)*PCAER
               GRJMK=GR(J,-K)*PCAER
               GRKJ=GR(K,J)*PCAER
               GRKMJ=GR(K,-J)*PCAER
               ARRJK=ARR(J,K)*PCAER
               ARRJMK=ARR(J,-K)*PCAER
               ARTJK=ART(J,K)*PCAER
               ARTJMK=ART(J,-K)*PCAER
               ARTKJ=ART(K,J)*PCAER
               ARTKMJ=ART(K,-J)*PCAER
               ATTJMK=ATT(J,-K)*PCAER
               ATTJK=ATT(J,K)*PCAER
	       
C* Calculs de la fonction source d'ordre IG, pour chaque 
C* niveau du profil et chaque direction.
C------------------------------------------------------
               Z=GA(J)

C* Recuperation des valeurs de luminance d'ordre IG-1.
C------------------------------------------------------
               XI1=I1(I,J)
               XI2=I1(I,-J)
               XQ1=Q1(I,J)
               XQ2=Q1(I,-J)
               XU1=U1(I,J)
               XU2=U1(I,-J)

C* Integration spatiale sur mup.
C-------------------------------
               II2=II2+Z*(XI1*BPJK+XI2*BPJMK+XQ1*GRKJ
     &             +XQ2*GRKMJ-XU1*GTKJ-XU2*GTKMJ)
               II1=II1+Z*(XI1*BPJMK+XI2*BPJK+XQ1*GRKMJ
     &             +XQ2*GRKJ+XU1*GTKMJ+XU2*GTKJ)
               QQ2=QQ2+Z*(XI1*GRJK+XI2*GRJMK+XQ1*ARRJK
     &             +XQ2*ARRJMK+XU2*ARTJMK-XU1*ARTJK)
               QQ1=QQ1+Z*(XI1*GRJMK+XI2*GRJK+XQ1*ARRJMK
     &             +XQ2*ARRJK-XU1*ARTJMK+XU2*ARTJK)
               UU2=UU2-Z*(XI1*GTJK-XI2*GTJMK+XQ1*ARTKJ
     &             +XQ2*ARTKMJ-XU1*ATTJK-XU2*ATTJMK)
               UU1=UU1-Z*(XI1*GTJMK-XI2*GTJK-XQ1*ARTKMJ
     &             -XQ2*ARTKJ-XU1*ATTJMK-XU2*ATTJK)
     
   47       CONTINUE


C* Stockage des resultats.
C--------------------------
            I2(I,K)=II2 * 0.5
            I2(I,-K)=II1 * 0.5
            Q2(I,K)=QQ2  * 0.5
            Q2(I,-K)=QQ1 * 0.5
            U2(I,K)=UU2  * 0.5
            U2(I,-K)=UU1 * 0.5
	    
   45 CONTINUE    ! Boucle sur I et K

  212 CONTINUE

      RETURN
      END  	!FIN DE LA PROCEDURE SOS_FSOURCE_ORDREIG







C*===========================================================================
C PROCEDURE: SOS_FSOURCE_DIFF_FRESNEL1
C ==========
C      Cette procedure calcule la fonction source de la diffusion 
C      du rayonnement solaire remontant apres une reflexion directe
C      sur mer plate, selon les lois de Fresnel (diffusion du champ
C      de reflexion d'ordre 1).
C
C
C Description des parametres
C -------------------------- 
C
C      F11sun  (double)  (E)     Element F11 de la matrice de reflexion de 
C				 Fresnel pour l'incidence solaire.    
C      F12sun  (double)  (E)     Idem pour l'element F12.   
C
C      XDEL(0:SOS_OS_NT)  (double)  (E) Pourcentage d'aerosols par niveau.
C      YDEL(0:SOS_OS_NT)  (double)  (E) Pourcentage de molecules par niveau. 
C
C      BETA0   (double)  (E)  
C      BETA2   (double)  (E) 
C      GAMMA2  (double)  (E)	
C      ALPHA2  (double)  (E)
C             Coefficients du developpement de la matrice 
C             de phase moleculaire en fonctions de Legendre.
C
C      BP(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                                                        Element P11.
C      GR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                                                        Element P12 ou P21.
C      GT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                                                        Element -P13 ou -P31.
C      ARR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (E) 
C                                                        Element P22.
C      ART(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (E)
C                                                        Element -P23 ou -P32.
C             Element Pij de la matrice de phase des aerosols pour un ordre
C             IS de la decomposition en series de Fourier.
C             Soit mup la direction d'incidence ( -1 < mup < 1)
C             et mu celle de diffusion ( -1 < mu < 1),
C             la diffusion est exprimee par Pij(mu,mup).	
C
C      XPL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (E) 
C      		        Fonction de Legendre PSL à l'ordre L=2
C   		        et S=IS en fonction de l'angle mu.
C      XRL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (E)
C      		        Fonction RSL a l'ordre L=2 et S=IS.
C      XTL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (double)  (E)
C      		        Fonction TSL a l'ordre L=2 et S=IS
C	      Fonctions de Legendre et fonctions derivees pour le calcul
C	      de la matrice de phase.
C
C      IS (I4)   (E)  Ordre de la decomposition en series de Fourier.
C
C      NBMU (I4)    (E)  Nombre d'angles (positifs) effectivement utiles
C
C      MUS (double)  (E)   Cosinus de l'angle zenithal solaire.
C
C      H(0:SOS_OS_NT) (double) (E)  Epaisseur optique a chaque niveau.
C
C      I2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double) (S) 
C                                                  Fonction source pour I
C      Q2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double) (S) 
C                                                  Fonction source pour Q
C      U2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double) (S) 
C                                                  Fonction source pour U
C             Vecteur fonction source, en fonction du niveau du 
C             profil et de la direction, pour la diffusion de
C             la lumiere obtenue par reflexion de Fresnel du faisceau 
C             solaire direct.
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
      SUBROUTINE SOS_FSOURCE_DIFF_FRESNEL1(F11sun,F12sun,XDEL,YDEL,
     &                                     BETA0,BETA2,GAMMA2,ALPHA2,
     &                                     BP,GR,GT,ARR,ART,
     &                                     XPL,XRL,XTL,IS,NBMU,
     &					   MUS,H,I2,Q2,U2)
     
      IMPLICIT NONE


C* Definition des variables    
C*-----------------------------------------------------------------------     
	! Vecteur fonction source en fonction du niveau du profil
	! et de la direction.
      DOUBLE PRECISION I2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour I.
      DOUBLE PRECISION Q2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour Q.
      DOUBLE PRECISION U2(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      						! Fonction source pour U.
						
      DOUBLE PRECISION F11sun   ! Element de la matrice de Fresnel
      			        ! avec pour incidence la direction
			        ! d'eclairement solaire.
      DOUBLE PRECISION F12sun   ! Element de la matrice de Fresnel.

      DOUBLE PRECISION XDEL(0:SOS_OS_NT)  !Pourcentage d'aerosols par niveau.
      DOUBLE PRECISION YDEL(0:SOS_OS_NT)  !Pourcentage de molecules par niveau.

      DOUBLE PRECISION MUS               ! Cosinu de l'angle zenithal solaire.
      DOUBLE PRECISION H(0:SOS_OS_NT)	 !Epaisseur optique pour chaque niveau.
           
      DOUBLE PRECISION YR   ! Pourcentage de molecules au niveau K du profil.
      DOUBLE PRECISION YYR  ! Pourcentage de molecules au niveau K+1 du profil.
      DOUBLE PRECISION XP   ! Pourcentage d'aerosols au niveau K du profil.
      DOUBLE PRECISION XXP  ! Pourcentage d'aerosols au niveau K+1 du profil.      


      DOUBLE PRECISION XPL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Fonction de Legendre PSL à l'ordre L=2 et S=IS
      			! en fonction de l'angle mu.     		
      DOUBLE PRECISION XRL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Fonction RSL a l'ordre L=2 et S=IS en fonction de mu.
      DOUBLE PRECISION XTL(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      			! Fonction TSL a l'ordre L=2 et S=IS en fonction de mu.

      DOUBLE PRECISION SPL   ! Valeur de la fonction de XPL pour mu = mus
           	
      DOUBLE PRECISION BETA0	! Coefficients du developpement de la matrice 
      DOUBLE PRECISION BETA2	! de phase moleculaire en fonctions de Legendre
      DOUBLE PRECISION GAMMA2	
      DOUBLE PRECISION ALPHA2	

	! Element Pij de la matrice de phase des aerosols pour un ordre
	! IS de la decomposition en series de Fourier.
	! Soit mup la direction d'incidence (-1 < mup < 1)
	! et mu celle de diffusion (-1 < mu < 1),
	! la diffusion est exprimee par Pij(mu,mup).	
      DOUBLE PRECISION 
     &    BP(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! P12 ou P21.
      DOUBLE PRECISION 
     &    GT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &       -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! -P13 ou -P31.
      DOUBLE PRECISION 
     &    ARR(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! P22.
      DOUBLE PRECISION 
     &    ART(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX,
     &        -SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) ! -P23 ou -P32.
     
           
      DOUBLE PRECISION BP0MJ
      DOUBLE PRECISION BP0J
      DOUBLE PRECISION GRJ0
      DOUBLE PRECISION GR0J
      DOUBLE PRECISION GR0MJ
      DOUBLE PRECISION GRMJ0
      DOUBLE PRECISION GT0MJ
      DOUBLE PRECISION GT0J
      DOUBLE PRECISION ARR0MJ
      DOUBLE PRECISION ARR0J
      DOUBLE PRECISION ARTJ0
      DOUBLE PRECISION ARTMJ0
       
      
      DOUBLE PRECISION COEFNT
      DOUBLE PRECISION COEFK
      DOUBLE PRECISION COEFKP1 
      
      INTEGER IS	! Ordre de la decomposition en series de Fourier.     
      INTEGER K		! Indice de niveau du profil.
      INTEGER J		! Indice d'angle d'incidence.
      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations
      
       
C* Traitements    
C*------------

C* Initialisation

      DO 10 J=-NBMU,NBMU
         DO 10 K=0,SOS_OS_NT
	    I2(K,J)=0.
	    Q2(K,J)=0.
	    U2(K,J)=0.
10    CONTINUE      


      COEFNT  = DEXP(2.*H(SOS_OS_NT)/MUS) / 4.  

      SPL=XPL(0)

      
      DO 1500 K= 0 , SOS_OS_NT-1
      
         YR=YDEL(K)
         XP=XDEL(K)
         
         YYR=YDEL(K+1)
         XXP=XDEL(K+1)
	 
         DO 1502 J=1,NBMU

C* Elements de matrice de phase pour le melange molecules + aerosols
            IF (IS.LE.2) THEN	 
	    
                BP0MJ=BP(0,-J)*XP+(BETA0+BETA2*XPL(-J)*SPL)*YR
                BP0J=BP(0,J)*XXP+(BETA0+BETA2*XPL(J)*SPL)*YYR
            
                GRJ0=GR(J,0)*XXP+YYR*XRL(0)*XPL(J)*GAMMA2
	        GR0J=GR(0,J)*XXP+YYR*XRL(J)*XPL(0)*GAMMA2
                GR0MJ=GR(0,-J)*XP+YR*XRL(-J)*SPL*GAMMA2
	        GRMJ0 =GR(-J,0)*XP+YR*GAMMA2*XRL(0)*XPL(-J)	
	    
                GT0MJ=GT(0,-J)*XP+YR*GAMMA2*SPL*XTL(-J)
                GT0J=GT(0,J)*XXP+YYR*GAMMA2*SPL*XTL(J)
	                
                ARR0MJ=ARR(0,-J)*XP+ALPHA2*YR*XRL(0)*XRL(-J)
                ARR0J=ARR(0,J)*XXP+ALPHA2*YYR*XRL(0)*XRL(J)            
            
                ARTJ0=ART(J,0)*XXP+YYR*ALPHA2*XTL(J)*XRL(0)
                ARTMJ0=ART(-J,0)*XP+YR*ALPHA2*XTL(-J)*XRL(0)	
		
            ELSE
	    
                BP0MJ=BP(0,-J)*XP
                BP0J=BP(0,J)*XXP
            
                GRJ0=GR(J,0)*XXP
	        GR0J=GR(0,J)*XXP
                GR0MJ=GR(0,-J)*XP
	        GRMJ0 =GR(-J,0)*XP
	    
                GT0MJ=GT(0,-J)*XP
                GT0J=GT(0,J)*XXP
	                
                ARR0MJ=ARR(0,-J)*XP
                ARR0J=ARR(0,J)*XXP          
            
                ARTJ0=ART(J,0)*XXP
                ARTMJ0=ART(-J,0)*XP	  
		  
	    ENDIF
   
	     
	     
C* Fonction source pour le champ montant
            COEFK   = COEFNT * DEXP(-H(K)/MUS)
	    
            I2(K,J)= COEFK * ( F11sun*BP0MJ + F12sun*GRMJ0 )		
            Q2(K,J)= COEFK * (  F11sun*GR0MJ + F12sun*ARR0MJ )
            U2(K,J)= COEFK * (  F11sun*GT0MJ + F12sun*ARTMJ0 )		
       
C* Fonction source pour le champ descendant
            COEFKP1 = COEFNT * DEXP(-H(K+1)/MUS) 
	    
            I2(K+1,-J)=  COEFKP1 * ( F11sun*BP0J + F12sun*GRJ0 )
            Q2(K+1,-J)=  COEFKP1 * ( F11sun*GR0J + F12sun*ARR0J )
            U2(K+1,-J)=  COEFKP1 * ( F11sun*GT0J + F12sun*ARTJ0 )
	     
 1502    CONTINUE
 1500 CONTINUE
      

      END      !FIN DE LA PROCEDURE SOS_FSOURCE_DIFF_FRESNEL1
       
       
       
       
       
       

C*===========================================================================
C PROCEDURE: SOS_PARAM_CONV
C ==========
C      Cette procedure calcule le parametre de convergence pour le test
C      de convergence en serie geometrique. 
C
C      Principe:
C      --------
C      On analyse la convergence en serie geometrique des diffusions 
C      successives a partir des ordres IG-2, IG-1 et IG.
C      Le test est realise pour le champ de luminance emergeant au TOA et pour
C      le champ eclairant le sol. S'il y a convergence a ces limites, on admet
C      qu'il y a convergence dans toute l'atmosphere.
C
C      Supposons avoir un comportement de serie geometrique convergente a 
C      partir de l'ordre IG-2 de diffusion.
C      La raison de la serie est estimee entre les termes IG-2 et IG-1, ainsi
C      qu'entre les termes IG-1 et IG :
C 		_ q  = D1 / A1  (rapport des ordres IG-1 / IG-2 pour I)
C		_ q' = G1 / D1  (rapport des ordres IG / IG-1 pour I)
C     On teste que q et q' sont voisins (condition pour une approximation de
C     serie geometrique).
C     On teste en meme temps le rapport de la somme des termes de la serie
C     geometrique convergente a partir de l'ordre IG (ordre IG / (1-q) )
C     devant la somme des termes precedents
C     (de 1 a IG-1, stockee dans I3, Q3 ou U3).
C
C     La routine sort le parametre de convergence Z1 qui est teste par le
C     programme appelant.
C
C Description des parametres
C -------------------------- 
C
C    NBMU (I4)    (E)  Nombre d'angles (positifs) effectivement utiles
C
C    I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)		
C                         Parametre de Stokes I pour le cumul des
C      		          ordres de diffusion precedent IG en cours.
C    Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)  Parametre de Stokes Q.
C    U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)  Parametre de Stokes U.
C            			
C    A1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C       	          Parametre de Stokes I pour l'ordre IG-2. 
C    B1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			  Parametre de Stokes Q pour l'ordre IG-2.
C    C1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			  Parametre de Stokes U pour l'ordre IG-2.
C
C    D1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			  Parametre de Stokes I pour l'ordre IG-1.	
C    E1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			  Parametre de Stokes Q pour l'ordre IG-1.
C    F1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			  Parametre de Stokes U pour l'ordre IG-1.
C
C    G1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      		          Parametre de Stokes I pour l'ordre IG.	
C    H1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			  Parametre de Stokes Q pour l'ordre IG.
C    P1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			  Parametre de Stokes U pour l'ordre IG.
C
C    Z1	   (double)   (S)    Parametre de convergence.
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
      SUBROUTINE SOS_PARAM_CONV(NBMU,A1,B1,C1,D1,E1,F1,G1,H1,P1,
     &                          I3,Q3,U3,Z1)
      
      IMPLICIT NONE
      

C* Definition des variables     
C*------------------------------------------------------------------

      DOUBLE PRECISION Z1	! Parametre de convergence.

	! Fonctions pour le test de convergence en series geometrique,
	! champ emergeant au TOA (mu > 0) et descendant au sol (mu < 0).
	! Les fonctions sont etablies pour un ordre S du developpement
	! en series de Fourier, en fonction de l'angle mu.		
      DOUBLE PRECISION I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)		
                                ! Parametre de Stokes I pour le cumul des
      				! ordres de diffusion precedent IG en cours.
      DOUBLE PRECISION Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes Q.
      DOUBLE PRECISION U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes U.
            			
      DOUBLE PRECISION A1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
       				! Parametre de Stokes I pour l'ordre IG-2. 
      DOUBLE PRECISION B1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG-2.
      DOUBLE PRECISION C1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes U pour l'ordre IG-2.

      DOUBLE PRECISION D1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG-1.
      DOUBLE PRECISION E1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG-1.
      DOUBLE PRECISION F1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG-1.

      DOUBLE PRECISION G1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG.
      DOUBLE PRECISION H1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG.
      DOUBLE PRECISION P1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG.

      DOUBLE PRECISION Y	! Variable pour les tests d'arret.

      INTEGER*4 K
      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations

C* Determination de la valeur maximale de Z1 pour tous les angles mu.
C----------------------------------------------------------------------------
      Z1=0.D+00

      DO 99 K=-NBMU,NBMU
         IF (K) 56,99,56
   56    IF (A1(K).EQ.0.D+00) GO TO 51
         IF (D1(K).EQ.0.D+00) GO TO 51
         IF (I3(K).EQ.0.D+00) GO TO 51
         Y=((G1(K)/D1(K)-D1(K)/A1(K))
     &       /((1-G1(K)/D1(K))**2)*(G1(K)/I3(K)))
         Z1=DMAX1(Z1,DABS(Y))
	 
   51    CONTINUE
   
         IF (B1(K).EQ.0.D+00) GO TO 53
         IF (E1(K).EQ.0.D+00) GO TO 53
         IF (Q3(K).EQ.0.D+00) GO TO 53
         Y=((H1(K)/E1(K)-E1(K)/B1(K))
     &       /((1-H1(K)/E1(K))**2)*(H1(K)/Q3(K)))
         Z1=DMAX1(Z1,DABS(Y))
	 
   53    CONTINUE
   
         IF (C1(K).EQ.0.D+00) GO TO 55
         IF (F1(K).EQ.0.D+00) GO TO 55
         IF (U3(K).EQ.0.D+00) GO TO 55
         Y=((P1(K)/F1(K)-F1(K)/C1(K))
     &       /((1-P1(K)/F1(K))**2)*(P1(K)/U3(K)))
         Z1=DMAX1(Z1,DABS(Y))
	 
   55    CONTINUE
   
   99 CONTINUE

      RETURN
      END  	!FIN DE LA PROCEDURE SOS_PARAM_CONV





C*============================================================================
C PROCEDURE: SOS_ARRET_DIFFUS
C ==========
C      Cette procedure prepare le test d'arret sur les diffusions quand les 
c      termes suivant deviennent tres petits
C
C
C Description des parametres
C -------------------------- 
C    
C    NBMU (I4) (E)  Nombre d'angles (positifs) effectivement utiles
C
C    I1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E)
C    Q1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E)
C    U1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E)
C      	 Champ de luminance d'ordre IG en fonction du niveau du 
C	 profil et de la direction.       	
C
C    I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)	
C            Parametre de Stokes I pour le cumul des ordres de diffusion 
C            precedent IG en cours.
C    Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) Parametre de Stokes Q.
C    U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) Parametre de Stokes U
C
C    Z1	   (double)   (S)    Variable seuil pour les tests d'arret.
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
      SUBROUTINE SOS_ARRET_DIFFUS(NBMU,I1,Q1,U1,I3,Q3,U3,Z1)
      
      IMPLICIT NONE
      

C* Definition des variables     
C*------------------------------------------------------------------

      DOUBLE PRECISION Z1	! Variable seuil pour les tests d'arret.

      	! Champ de luminance d'ordre IG en fonction du niveau du profil
	! et de la direction.       	
      DOUBLE PRECISION I1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Champ I.
      DOUBLE PRECISION Q1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Champ Q.
      DOUBLE PRECISION U1(0:SOS_OS_NT,-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Champ U.

	! Fonctions pour le test de convergence en series geometrique,
	! champ emergeant au TOA (mu > 0) et descendant au sol (mu < 0).
	! Les fonctions sont etablies pour un ordre S du developpement
	! en series de Fourier, en fonction de l'angle mu.	
      DOUBLE PRECISION I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)		
                                ! Parametre de Stokes I pour le cumul des
      				! ordres de diffusion precedent IG en cours.
      DOUBLE PRECISION Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes Q.
      DOUBLE PRECISION U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes U.
 
      DOUBLE PRECISION Y	! Variable pour les tests d'arret.

      INTEGER*4 IND		! Numero de couche pour les mu > 0 ou mu < 0
      INTEGER*4 K 
      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations

C------------------------------------------------------ 
C* On compare la valeur du terme d'ordre IG a la somme des 
C* termes precedents.
C* Determination de la valeur maximale de leur rapport
C* pour le champ montant au TOA et descendant au sol.
C------------------------------------------------------
      Z1=0.D+00
      DO 88 K=-NBMU,NBMU
      
         IF(K) 91,88,93
	 
  91     IND=SOS_OS_NT		! Au sol
         GOTO 95
	 
  93     IND=0			! Au TOA
  
  95     IF (I3(K).EQ.0.D+00) GO TO 92
  
         Y=I1(IND,K)/I3(K)
         Z1=DMAX1(Z1,DABS(Y))
	 
   92    IF (Q3(K).EQ.0.D+00) GO TO 94
   
         Y=Q1(IND,K)/Q3(K)
         Z1=DMAX1(Z1,DABS(Y))
	 
   94    IF (U3(K).EQ.0.D+00) GO TO 88
   
         Y=U1(IND,K)/U3(K)
         Z1=DMAX1(Z1,DABS(Y))
   88 CONTINUE
   
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_ARRET_DIFFUS




C*============================================================================
C PROCEDURE: SOS_ARRET_FOURIER
C ==========
C      Cette procedure prepare le test d'arret sur la decomposition en 
c      serie fourier.
c
C
C Description des parametres
C -------------------------- 
C
C    NBMU (I4) (E)  Nombre d'angles (positifs) effectivement utiles
C
C    I4(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) 	
C    Q4(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) 
C    U4(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) 
C        Tableaux de sommation des ordres S de la decomposition
C        en series de Fourier, en fonction de l'angle mu.
C
C    I5(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) 	
C    Q5(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) 
C    U5(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) 	
C        Tableaux de sommation des ordres S de la decomposition
C        en series de Fourier, avec permutation du signe des ordres S,
C	 en fonction de l'angle mu.
C
C    I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)	
C            Parametre de Stokes I pour le cumul des ordres de diffusion 
C            precedent IG en cours.
C    Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) Parametre de Stokes Q.
C    U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E) Parametre de Stokes U
C
C    Z1	   (double)   (S)    Variable seuil pour les tests d'arret.
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
C==============================================================================
      SUBROUTINE SOS_ARRET_FOURIER(NBMU,I3,Q3,U3,I4,Q4,U4,I5,Q5,U5,Z1)
      
      IMPLICIT NONE
      

C* Definition des variables     
C*------------------------------------------------------------------

      DOUBLE PRECISION Z1	! Variable seuil pour les tests d'arret.

	! Tableaux de sommation des ordres S de la decomposition
        ! en series de Fourier, en fonction de l'angle mu.
      DOUBLE PRECISION I4(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      DOUBLE PRECISION Q4(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      DOUBLE PRECISION U4(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)

        ! Tableaux de sommation des ordres S de la decomposition
        ! en series de Fourier, avec permutation du signe des ordres S,
	! en fonction de l'angle mu.
      DOUBLE PRECISION I5(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      DOUBLE PRECISION Q5(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      DOUBLE PRECISION U5(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	

	! Champ emergeant au TOA (mu > 0) et descendant au sol (mu < 0).
	! Les fonctions sont etablies pour un ordre S du developpement
	! en series de Fourier, en fonction de l'angle mu.
      DOUBLE PRECISION I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)		
                                ! Parametre de Stokes I pour le cumul des
      				! ordres de diffusion precedent IG en cours.
      DOUBLE PRECISION Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes Q.
      DOUBLE PRECISION U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes U.
 

      DOUBLE PRECISION Y	! Variable pour les tests d'arret.
      
      INTEGER*4 J
      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations

C------------------------------------------------------------------------- 
C* On teste l'arret de la decomposition en series de Fourier par
C* comparaison de l'ordre IS (pour chaque direction mu) par rapport
C* a la somme des contributions precedentes.
C* On sort la valeur maximale Z1 de ce rapport pour le champ emergant 
C* au TOA et descendant au niveau du sol.
C* Si la contribution relative (en valeur absolue) est inferieure
C* au seuil SOS_PH_SEUIL_SF on arrete la decomposition en series de Fourier
C* dans le programme appelant.
C----------------------------------------------------------------------------
      Z1=0.D+00
      DO 169 J=-NBMU,NBMU
         IF (J.EQ.0) GOTO 169

         IF (Q4(J).EQ.0.D+00) GO TO 71
	 
         Y=DABS(Q3(J)/Q4(J))
         Z1=DMAX1(Z1,Y)
	 
   71    IF (I4(J).EQ.0.D+00) GO TO 73
   
         Y=DABS(I3(J)/I4(J))
         Z1=DMAX1(Z1,Y)
	 
   73    IF (U4(J).EQ.0.D+00) GO TO 74
   
         Y=DABS(U3(J)/U4(J))
         Z1=DMAX1(Z1,Y)

   74    IF (Q5(J).EQ.0.D+00) GO TO 771
   
         Y=DABS(Q3(J)/Q5(J))
         Z1=DMAX1(Z1,Y)
	 
  771    IF (U5(J).EQ.0.D+00) GO TO 772
  
         Y=DABS(U3(J)/U5(J))
         Z1=DMAX1(Z1,Y)
	 
  772    IF (I5(J).EQ.0.D+00) GO TO 169
  
         Y=DABS(I3(J)/I5(J))
         Z1=DMAX1(Z1,Y)

  169 CONTINUE

      RETURN
      END  	!FIN DE LA PROCEDURE SOS_ARRET_FOURIER





C*========================================================================
C PROCEDURE: SOS_AJOUT_QUEUE
C ==========
C      Cette procedure realise, si on est en regime de convergence, l'ajout 
C      de la queue de la serie geometrique a la somme des diffusions
C      precedentes.
C
C
C Description des parametres
C -------------------------- 
C
C    NBMU (I4) (E)  Nombre d'angles (positifs) effectivement utiles
C
C    D1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                            Parametre de Stokes I pour l'ordre IG-1.	
C    E1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      		             Parametre de Stokes Q pour l'ordre IG-1.
C    F1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			     Parametre de Stokes U pour l'ordre IG-1.
C
C    G1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			     Parametre de Stokes I pour l'ordre IG.	
C    H1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			     Parametre de Stokes Q pour l'ordre IG.
C    P1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			     Parametre de Stokes U pour l'ordre IG.
C
C    D1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E) 
C                            Parametre de Stokes I pour l'ordre IG-1
C                            (fonctions pour le niveau de sortie 
C                             des resultats)	
C    E1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      		             Parametre de Stokes Q pour l'ordre IG-1.
C    F1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			     Parametre de Stokes U pour l'ordre IG-1.
C
C    G1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			     Parametre de Stokes I pour l'ordre IG.
C                            (fonctions pour le niveau de sortie 
C                             des resultats)		
C    H1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			     Parametre de Stokes Q pour l'ordre IG.
C    P1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E)
C      			     Parametre de Stokes U pour l'ordre IG.
C
C
C    I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E/S)	
C            Parametre de Stokes I pour le cumul des ordres de diffusion 
C            precedent IG en cours.
C    Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E/S) Parametre de Stokes Q.
C    U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)   (double)  (E/S) Parametre de Stokes U
C
C    I3OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E/S)	
C            Parametre de Stokes I pour un ordre S du developpement en 
C            series de Fourier, en fonction de l'angle mu (niveau de sortie 
C            des resultats).
C    Q3OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E/S) Parametre de Stokes Q.
C    U3OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)  (E/S) Parametre de Stokes U.
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
      SUBROUTINE SOS_AJOUT_QUEUE(NBMU,D1,E1,F1,G1,H1,P1,
     &                           D1OUT,E1OUT,F1OUT,
     &                           G1OUT,H1OUT,P1OUT,
     &                           I3,Q3,U3,I3OUT,Q3OUT,U3OUT)
      
      IMPLICIT NONE
      

C* Definition des variables     
C*------------------------------------------------------------------

	! Champ emergeant au TOA (mu > 0) et descendant au sol (mu < 0).
	! Les fonctions sont etablies pour un ordre S du developpement
	! en series de Fourier, en fonction de l'angle mu.	
      DOUBLE PRECISION I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)		
                                ! Parametre de Stokes I pour le cumul des
      				! ordres de diffusion precedent IG en cours.
      DOUBLE PRECISION Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes Q.
      DOUBLE PRECISION U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
      				! Parametre de Stokes U
				
      DOUBLE PRECISION D1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG-1.
      DOUBLE PRECISION E1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG-1.
      DOUBLE PRECISION F1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG-1.

      DOUBLE PRECISION G1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG.
      DOUBLE PRECISION H1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG.
      DOUBLE PRECISION P1(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG.

	! Memes fonctions pour le niveau de sortie des resultats:
      DOUBLE PRECISION D1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG-1.
      DOUBLE PRECISION E1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG-1.
      DOUBLE PRECISION F1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG-1.

      DOUBLE PRECISION G1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes I pour l'ordre IG.
      DOUBLE PRECISION H1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      				! Parametre de Stokes Q pour l'ordre IG.
      DOUBLE PRECISION P1OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)
      			  	! Parametre de Stokes U pour l'ordre IG.


      DOUBLE PRECISION I3OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				!(S)   Parametre de Stokes I pour un ordre S
      				!      du developpement en series de Fourier,
				!      en fonction de l'angle mu.
      DOUBLE PRECISION Q3OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				!(S)   Parametre de Stokes Q.
      DOUBLE PRECISION U3OUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	
      				!(S)   Parametre de Stokes U.


      DOUBLE PRECISION I1QUEUE	
      		! Valeur estimee de la queue de serie geometrique pour I.
      DOUBLE PRECISION Q1QUEUE
      		! Valeur estimee de la queue de serie geometrique pour Q.
      DOUBLE PRECISION U1QUEUE
      		! Valeur estimee de la queue de serie geometrique pour U.

      INTEGER*4 J
      INTEGER*4 NBMU    ! Nombre d'angles effectifs des simulations

C* Pour le champ emergeant au TOA et descendant au sol
C------------------------------------------------------
      DO 610 J=-NBMU,NBMU
      
         IF (J.EQ.0) GOTO 610
	
         IF (D1(J).EQ.0.) THEN
            I1QUEUE = 0.
         ELSE	
            I1QUEUE = G1(J) / ( 1-G1(J)/D1(J) )
         ENDIF
	 	
         IF (E1(J).EQ.0.) THEN
            Q1QUEUE = 0.
         ELSE	
            Q1QUEUE = H1(J) / ( 1-H1(J)/E1(J) )
         ENDIF	
	
         IF (F1(J).EQ.0.) THEN
            U1QUEUE = 0.
         ELSE	
            U1QUEUE = P1(J) / ( 1-P1(J)/F1(J) )
         ENDIF	

         I3(J)=I3(J)+I1QUEUE
         Q3(J)=Q3(J)+Q1QUEUE
         U3(J)=U3(J)+U1QUEUE
	
  610 CONTINUE

C* Pour le champ montant et descendant au niveau de sortie des resultats
C------------------------------------------------------------------------
      DO 612 J=-NBMU,NBMU
      
         IF (J.EQ.0) GOTO 612
	
         IF (D1OUT(J).EQ.0.) THEN
            I1QUEUE = 0.
         ELSE	
            I1QUEUE = G1OUT(J) / ( 1-G1OUT(J)/D1OUT(J) )
         ENDIF
	 	
         IF (E1OUT(J).EQ.0.) THEN
            Q1QUEUE = 0.
         ELSE	
            Q1QUEUE = H1OUT(J) / ( 1-H1OUT(J)/E1OUT(J) )
         ENDIF	
	
         IF (F1OUT(J).EQ.0.) THEN
            U1QUEUE = 0.
         ELSE	
            U1QUEUE = P1OUT(J) / ( 1-P1OUT(J)/F1OUT(J) )
         ENDIF	

         I3OUT(J)=I3OUT(J)+I1QUEUE
         Q3OUT(J)=Q3OUT(J)+Q1QUEUE
         U3OUT(J)=U3OUT(J)+U1QUEUE
	   	
  612 CONTINUE

      RETURN
      END  	!FIN DE LA PROCEDURE SOS_AJOUT_QUEUE
