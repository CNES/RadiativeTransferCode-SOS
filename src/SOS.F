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
C* FICHIER: SOS.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Module principal du code des OS
C*           - lit le profil atmospherique et l'adapte a une eventuelle 
C*             troncature de la fonction de phase des aérosols,
C*           - resoud l'equation de transfert radiatif par ordre 
C*             successif d'interaction (diffusion/reflexion),
C*           - recombine les paramètres de Stokes à partir de leur 
C*             decomposition en series de Fourier.
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, B.Lafrance, X.Lenot
C* DATE: 30/04/2002
C*
C* MOD:VERSION:1.0: Mise a niveau du code des OS du LOA
C*             2.0: Lecture des paramètres aérosols et profils, 
C*                  Mise à jour des commentaires.
C*             2.1: 07/07/2008 
C*                  Modification du format de lecture des paramètres WIND et TETAS
C*                  pour correction d'une erreur de décodage si valeur entière.
C*			WIND passe de F4.1 en F4
C*                   TETAS passe de F9.5 en F9
C*
C* MOD:VERSION:3.0: 24/11/2009 
C*       - Modification du passage des arguments : 
C*       - passage par couples (mot clef, valeur). 
C*       - adaptation de la gestion des cas d'erreurs. 
C*       - Conversion des messages d'erreurs en anglais.
C*
C*
C* MOD:VERSION:3.1: 22/12/2009 
C*                  Modification du format de lecture des paramètres Alpha et Beta du modele Nadal
C*                  pour correction d'une erreur de décodage si valeur entière.
C*			ALPHA passe de F6.4 en F6
C*			BETA  passe de F4.1 en F4
C*
C*                  
C* MOD:VERSION:4.0: 10/01/2010 
C*
C*   - Suppression de la constante CTE_OS_FICANGLE (de SOS.h) 
C*     au profit d'un passage par argument du nom du fichier des angles (produit par SOS_ANGLES).
C*
C*   - Modification de la taille des tableaux de luminance :
C*     ==> Passage de CTE_OS_NBMU à CTE_OS_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'à LUM_NBMU égal au nombre effectif d'angles utilises.
C*
C*   - Suppression de la constante  CTE_OS_NB (de SOS.h)  
C*     au profit : d'un passage par argument de la valeur OS_NB (lue dans le fichier des angles)
C*
C*   - Modification de la taille des tableaux de decomposition en polynomes de Legendre :
C*     ==> Passage de CTE_OS_NB à CTE_OS_NB_MAX.
C*     ==> Boucle sur les ordre jusqu'à OS_NB. 
C* 
C*   - Gestion des angles : 
C*         * Suppression de la lecture d'un fichier d'angles de Gauss
C*         * Ajout d'un appel a la routine SOS_ANGLES pour calculer les angles de Gauss et leurs poids
C*           et pour ajouter des angles utilisateurs (de poids nuls) : angles utiles = LUM_NBMU.
C*         * Suppression de la lecture de SOS.Imus et SOS.thetas : la récupération de ces paramètres
C*           est faite via la lecture du contenu du fichier des angles.
C*         * Ajout des paramètres d'entrée optionnels SOS.ResFileUp.UserAng et SOS.ResFileDown.UserAng
C*           pour l'enregistrement des luminances de sortie uniquement pour les angles "utilisateur"
C*           --> Adaptation de la routine SOS_TRPHI_OPTION.                                  
C*
C*   - Modification de la taille des tableaux de luminance :
C*     ==> Passage de CTE_OS_NBMU à CTE_OS_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'à LUM_NBMU égal au nombre effectif d'angles utilises.
C*
C*   - Routine SOS_TRPHI_OPTION : 
C*         * Modification de la ligne de lecture des tableaux de
C*           la  matrice de reflexion a l'ordre IS de la decomposition en series de Fourier :
C*           ==> limitation aux bornes utiles -LUM_NBMU:LUM_NBMU
C*         * Passage du tableau IND_ANGOUT definissant les angles pour lesquels une sortie
C*           est requise dans les fichiers optionnels SOS.ResFileUp.UserAng 
C*           et SOS.ResFileDown.UserAng
C*
C*   - Adaptation des routines SOS_OS & SOS_TRPHI_OPTION au passage du paramètre
C*     LUM_NBMU (nb d'angles utiles) et de OS_NB (ordre max des développements en SF).
C*
C*   - Adaptation au passage du paramètre FICTRANS (-SOS.Trans) pour sortir les valeurs
C*     de transmission.
C*
C* MOD:VERSION:4.1: 24/06/2010 
C*   -  Les etiquettes de FORMAT en ,X, passent en ,1X, pour eviter des erreurs de compilation.
C*   -  Correction pour effacer le fichier FICOS_TRANS_TMP à la fin des calculs de transmission.
C*
C* MOD:VERSION:6.0: 06/06/2013
C*   -  La fonction SOS.F de la version precedente (V 5.1) a ete scindee en 3 dans
C*	 cette nouvelle version : SOS_PREPAOS.F, SOS.F et SOS_FINOS.F. LA fonction SOS.F
C*	 est executee pour chaque exponentielle CKD, alors que SOS_PREPAOS.F et SOS_FIN_OS.F
C*	 sont lancees respectivement avant et apres la boucle sur les exponentielles (voir SOS_MAIN.F)
C*
C*   -  Les paramètres d'entrée sont fournis sous la forme de liste de paramètres
C*      et non plus sous la forme de couples avec mot cle associe => tous les
C*	 paramètres deviennent donc obligatoires
C*
C*   -  Modification de la valeur de Pi, declaree en constante dans le fichier
C*      (INCTE_PI remplace CTE_PI)
C*
C* MOD:VERSION:6.1: 29/05/2020
C*   -  Ajustement du format du fichier de profil atmosphérique à l'amélioration 
C*      de la précision du profil (étiquette de format 70)
C*
C* MOD:VERSION:6.2: 07/07/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    CTE_LENDIR : Longueur des noms de repertoires.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_OS_NBMU_MAX : Nombre maximal d'angles positifs pour la resolution de l'ETR
C    CTE_OS_NB_MAX :  Valeur maximale pour l'ordre limite du développement en 
C                     polynomes de Legendre.
C    CTE_OS_NT  : Nombre de couches du profil atmospherique.
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROGRAMME: SOS
C ==========      
C      Programme principal du code des OS 
C
C Donnees en entrée du programme 
C ------------------------------
C    -->Les paramètres/arguments du programme
C
C             FICOS (CHARACTER) : (E) Nom complet du fichier resultat des OS.
C                                 (repertoire + nom fichier + extension)
C
C             FICTRANS (CHARACTER): (E) nom complet du fichier des valeurs de transmissions
C                                   (repertoire + nom fichier + extension)
C
C             FICPROFIL(CHARACTER) : (E) nom complet du fichier resultat de PROFIL
C                                    (repertoire + nom fichier + extension)
C
C             NT (I4) : (E) Nombre de couches dans le profil atmospherique
C
C      	ZOUT (F9) : (E) Altitude du niveau de sortie
C			         -1 =>  Sortie standard:
C				         Luminance émergeante au TOA
C				         + descendante au niveau du sol.
C			          z =>  Luminance montante et descendante
C				         pour l'altitude z (0 <= z <= CTE_TOA_ALT).
C			
C             IGMAX (I4) : (E) Ordre maximal d'interaction.
C		 
C             IPOLAR (I4) : (E) Option de coupure de la polarisation.
C                               0 : simulation sans polarisation.
C			           1 : simulation normale avec polarisation.
C
C             RON  (DOUBLE) : (E) Facteur de depolarisation (moleculaire) 
C                             (format de decodage en réel F9.5)
C
C             IND_SURF (DOUBLE) : (E) Indice de réfraction de la surface / air.		 	
C
C             RHO (DOUBLE) : (E) Albédo de surface     	 	 	
C
C             IMAT_SURF (I4) :(E) Indice precisant si les simulations
C      			      utilisent un modele de reflexion BRDF ou BPDF :
C			             1 : pour calcul de BRDF ou BPDF (il faut utiliser N0 > 0).
C			             sinon : pas de simulation de BRDF ou BPDF. 	
C
C             IFRESNEL (I4) : (E) Indice precisant si les simulations
C      			      introduisent la reflexion de Fresnel sur mer plate :
C			             1 : pour un calcul de la reflexion de Fresnel.
C			             Sinon : pas de prise en compte.    
C        
C             FICSURF (CHARACTER) : ==> Associé à : ISURF = 1 3 4 5 ou 6 (voir SOS_PREPA_OS)
C					nom complet du fichier contenant la 
C                   			matrice de reflexion de surface  (fichier Glitter ou autres BRDF/BPDF).
C                   			(repertoire + nom fichier + extension)
C                    		> DEFAULT si on utilise un fichier de la base des BRDF / BPDF.	 
C
C      	N0 (I4) : (E) Numéro de mu de Gauss pour l'angle zenithal solaire.      			
C			      (angle complementaire : mus = -rmu(N0))
C			      (N0 > 0 => pour en tenir compte).
C
C      	PIZ (DOUBLE) : (E) Albédo de simple diffusion.
C
C      	PIZTR (DOUBLE) : (E) Albédo de simple diffusion.apres troncature
C
C      	A (DOUBLE) : (E) Coefficient de troncature
C
C      	RMU (DOUBLE) : (E) Tableau 1D de 2*CTE_OS_NBMU_MAX + 1 elements
C			            contenant le cosinus des angles de Gauss et des
C			            angles fournis par l'utilisateur
C
C      	GA (DOUBLE) : (E) Tableau 1D de 2*CTE_OS_NBMU_MAX + 1 elements
C			           contenant le poids associe a chaque angle de Gauss et chaque
C			           angle fourni par l'utilisateur
C 
C             TETAS (DOUBLE) : (E) Angle zenithal solaire en degrés (0 < tetas < 90).
C      		               (angle complementaire : mus= -cos(tetas))
C
C             OS_NB (I4) : (E) Ordre des développements des fonctions de phase 
C			          en Polynomes de Legendre.      			               
C
C             LUM_NBMU (I4) : (E) Nombre d'angles effectifs pour les calculs de luminances      			               
C
C             ALPHA (DOUBLE) : (E) Tableau 1D de CTE_OS_NB_MAX + 1 elements    
C				       contenant les coefficients Alpha de la decomposition en fonctions 
C				       de Legendre des elements Pij de la fonction de phase  			               
C
C             BETA (DOUBLE) : (E) Tableau 1D de CTE_OS_NB_MAX + 1 elements    
C				      contenant les coefficients Beta de la decomposition en fonctions 
C				      de Legendre des elements Pij de la fonction de phase  			               
C
C             GAMMA (DOUBLE) : (E) Tableau 1D de CTE_OS_NB_MAX + 1 elements    
C				       contenant les coefficients Gamma de la decomposition en fonctions 
C				       de Legendre des elements Pij de la fonction de phase  			               
C
C             ZETA (DOUBLE) : (E) Tableau 1D de CTE_OS_NB_MAX + 1 elements    
C				      contenant les coefficients Zeta de la decomposition en fonctions 
C				      de Legendre des elements Pij de la fonction de phase  			               
C
C             TTOT_TRONC (DOUBLE) : (S) Epaisseur optique totale au sol (tronquée)    
C
C             TTOT_VRAI (DOUBLE) : (S) Epaisseur optique totale au sol (non tronquée)    
C
C             TAUOUT (DOUBLE) : (S) Epaisseur optique au niveau de sortie pour le champ montant
C			               (tronquée)    
C
C             EMOINS  (DOUBLE): (S) Flux diffus eclairant le sol, normalisé par
C                   			 le flux solaire incident au TOA.
C             
C             EPLUS  (DOUBLE) : (S) Flux diffus emergeant au TOA, normalisé par
C                   			 le flux solaire incident au TOA.
C
C             TDIFMUS (DOUBLE) : (S) Transmission diffuse pour la direction d'incidence solaire
C      				 de l'atmosphere equivalente   
C
C             TDIFMUG (DOUBLE) : (S) Tableau 1D de 2*CTE_OS_NBMU_MAX + 1 elements
C			     	         contenant la transmission diffuse de l'atmosphere equivalente
C				         pour les angles de Gauss et les angles fournis par l'utilisateur   
C
C	       TRACE (LOGICAL) : (E) Si vrai => ecriture des fichiers trace
C
C	       IDLOG (I4) : (E) Identifiant du fichier log
C
C             IER (I4) : (S) code d'erreur =0 si pas d'erreur, =-1 sinon
C
C    -->Fichiers
C       Le fichier Matrice de reflexion de surface
C       Le fichier Profil
C       Le fichier des transmissions
C
C
C Resultats fournis par le programme
C ----------------------------------
C   Fichier resultat des OS (decomposition en series de Fourier)
C   Transmissions diffuses
C   Epaisseurs optiques (sol et au niveau de sortie)
C   Fichier de trace (optionnel) (No d'unite logique du fichier trace =IDLOG)
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
C   -->Contenu du fichier resultat des OS (FICOS)
C       Fichier binaire non formatte, cree par SOS_OS
C       Il contient les paramètres de Stokes Q , U et I pour un ordre S
C       du développement en series de Fourier, en fonction de l'angle mu
C       (valeurs en double précision).
C
C   -->Contenu du fichier profil
C	Fichier ASCII contenant pour chaque niveau du profil atmospherique :
C	L'indice du niveau, l'épaisseur optique molecules + aérosols, 
C	l'épaisseur optique d'absorption, l'épaisseur optique totale,
C	le pourcentage d'aérosols et le pourcentage de molecules
C
C Common utilise:
C --------------
C  Aucun
C
C
C Cas d'erreur :
C ------------        
C     -Decodage des paramètres incorrect
C     -Erreur dans un sous-programme
C     -Erreur a l'ouverture du fichier Trace
C     -Paramètres non valables
C
C     Affichage d'un message à l'écran, arret du programme et 
C     retour du status 1 au shell
C   
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel aux routines:
C     - SOS_OS
C
C============================================================================
       SUBROUTINE SOS(FICOS,FICTRANS,FICPROFIL,
     &		      NT,ZOUT,IGMAX,IPOLAR,RON,           
     &	             IND_SURF,RHO,IMAT_SURF,IFRESNEL,FICSURF,N0,PIZ,
     &		      PIZTR,A,RMU,GA,TETAS,OS_NB,LUM_NBMU,ALPHA,BETA,
     &		      GAMMA,ZETA,TTOT_TRONC,TTOT_VRAI,TAUOUT,TDIFMUS,
     &		      TDIFMUG,EMOINS,EPLUS,TRACE,IDLOG,IER)      
      
      IMPLICIT NONE

C* Definition des variables                         
C*-----------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                         !Cosinus des angles de Gauss (et utilisateur).
      DOUBLE PRECISION GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                         !Poids des angles de Gauss (et utilisateur).				
 
 
      DOUBLE PRECISION TETAS ! Angle zenithal solaire en degrés (0 < tetas < 90).
      				 ! (angle complementaire : mus= -cos(tetas))

      
      DOUBLE PRECISION IND_SURF  ! Indice de réfraction de la surface / air.
           
      DOUBLE PRECISION RON	! Facteur de depolarisation (moleculaire).

      DOUBLE PRECISION RHO	! Albédo de surface.
      
      DOUBLE PRECISION A	! Coefficient de troncature.
      DOUBLE PRECISION PIZ	! Albédo de simple diffusion.
      DOUBLE PRECISION PIZTR	! Albédo de simple diffusion apres troncature.

      DOUBLE PRECISION ALPHA(0:CTE_OS_NB_MAX)	! Paramètres alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NB_MAX)	! Paramètres beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NB_MAX)	! Paramètres gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Paramètres zeta(k).

      DOUBLE PRECISION VA,VATR
      DOUBLE PRECISION VR
      DOUBLE PRECISION VG
      
      DOUBLE PRECISION H(0:CTE_OS_NT)	! Epaisseur optique pour chaque niveau
      DOUBLE PRECISION Hma(0:CTE_OS_NT)	! Epaisseur optique molecules + aérosols pour chaque niveau
      DOUBLE PRECISION Ha(0:CTE_OS_NT)	! Epaisseur optique absorption pour chaque niveau
      DOUBLE PRECISION HTR(0:CTE_OS_NT)	! Epaisseur optique tronquée pour chaque niveau
      DOUBLE PRECISION XDEL(0:CTE_OS_NT)  ! Pourcentage d'aérosols par niveau
      DOUBLE PRECISION YDEL(0:CTE_OS_NT)  ! Pourcentage de molecules par niveau
      DOUBLE PRECISION ZPROF(0:CTE_OS_NT) ! Altitudes des niveaux du profil (en km)

      DOUBLE PRECISION TTOT_TRONC ! Epaisseur optique totale au sol (tronquée).
      DOUBLE PRECISION TTOT_VRAI  ! Epaisseur optique totale au sol (non tronquée).
      DOUBLE PRECISION TAUOUT     ! Epaisseur optique au niveau de
      				      ! sortie des resultats pour le champ montant (tronquée).

      DOUBLE PRECISION EMOINS   ! Flux diffus eclairant le sol, normalisé par
      				    ! le flux solaire incident au TOA.

      DOUBLE PRECISION EPLUS    ! Flux diffus emergeant au TOA, normalisé par
      				    ! le flux solaire incident au TOA.

      DOUBLE PRECISION EPLUS_SOLNOIR    ! Flux diffus emergeant au TOA, normalisé par
      				            ! le flux solaire incident au TOA, pour le cas d'un sol noir.
								
      DOUBLE PRECISION TDIFMUS  !Transmission diffuse pour la direction d'incidence solaire
      				    !de l'atmosphere equivalente
				               
      DOUBLE PRECISION TDIFMUG(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				!Transmission diffuse pour les angles de Gauss de l'atmosphere equivalente

      DOUBLE PRECISION TDIFMUS_TMP  !Transmission diffuse pour chaque angle de Gauss
      				        !de l'atmosphere equivalente 

      DOUBLE PRECISION  ZZ      ! Coefficient d'interpolation lineaire entre 2 altitudes	    

      INTEGER*4 NT      ! Nombre de couches dans le profil
     
      INTEGER*4 N0	! Numéro de mu de Gauss pour l'angle
      			! zenithal solaire.
			! (angle complementaire : mus = -rmu(N0))
			! (N0 > 0 => pour en tenir compte).

				
      DOUBLE PRECISION ZOUT	! Type de sortie :
				!      -1 =>  Sortie standard:
				!             Luminance émergeante au TOA
				!           + descendante au niveau du sol.
				!       z =>  Luminance montante et descendante
				!             pour l'altitude z (0 <= z <= CTE_TOA_ALT).
			

      INTEGER*4 IGMAX	! Ordre maximal de diffusion.
      
      INTEGER*4 IPOLAR	! Option de simulation sans polarisation.
      			       ! 0 : calculs sans polarisation.
			       ! 1 : calculs normaux avec polarisation.

      INTEGER*4 ITRPHI	! Indice de type de sortie.
			       !   1 =>  Plan de visee pour une valeur d'azimut.
			       !   2 =>  Diagramme polaire en fonction de la direction
			       !         d'azimut et de visee.
			
      INTEGER*4 PAS_PHI	! Pas sur l'azimut en degrés (pour ITRPHI=2)

      INTEGER*4 IGLI	! Indice precisant si les simulations
      			! introduisent le Glitter :
			!    1 => pour calcul du Glitter (il faut alors utiliser N0 > 0).
			!    Sinon => pas de simulation du Glitter.            
	
      INTEGER*4 IFRESNEL ! Indice precisant si les simulations
      			    ! introduisent la reflexion de Fresnel sur mer plate :
			    !    1 => pour un calcul de la reflexion de Fresnel.
			    !    Sinon => pas de prise en compte.            
	
		 
	
      INTEGER*4 IMAT_SURF ! Indice precisant si les simulations
      			     ! utilisent un modele de reflexion BRDF ou BPDF :
			     !    1 => pour calcul de BRDF ou BPDF (il faut utiliser N0 > 0).
			     !    Sinon => pas de simulation de BRDF ou BPDF. 	
			  			  			 											
      INTEGER*4 IBORM ! Ordre maximal du développement en series de Fourier
      			 ! en cours de traitement :
			 !   _ 2 pour une atmosphere purement moleculaire,
			 !   _ NB en presence d'aérosols.

      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon

      INTEGER*4 K       ! Indice de tableau
      INTEGER*4 I,J     ! Indice de tableau
      
      INTEGER*4 LDIRSURF ! longueur de la chaine de caractere REPSURF.
      INTEGER*4 LFICSURF ! longueur de la chaine de caractere FICSURF.
      
      INTEGER*4 NBMU_GAUSS      ! Nombre d'angles de Gauss utilises       
      INTEGER*4 LUM_NBMU        ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB 	 ! Ordre des développements des fonctions de phase 
			        ! en Polynomes de Legendre.      			               
      
      LOGICAL TRACE     ! = vrai,si ecriture dans le fichier trace
      LOGICAL LTA       ! = vrai,si épaisseur optique des aérosols nulle
      LOGICAL EX	! Verifie l'existence d'un fichier

      CHARACTER*CTE_LENFIC2 FICOS       ! Nom du fichier resultat des OS.
      CHARACTER*CTE_LENFIC2 FICOS_TRANS_TMP   ! Nom du fichier temporaire des OS (pour les transmissions)
      
      CHARACTER*CTE_LENFIC2 FICSURF ! Fichier contenant la matrice de reflexion de
      				        ! surface (fichier Glitter ou autres BRDF/BPDF).
      
      CHARACTER*CTE_LENFIC2 FICPROFIL   ! nom du fichier d'entrée PROFIL

      CHARACTER*CTE_LENFIC2 FICTRANS    ! Nom du fichier optionnel de sortie des transmissions.

             
      INTEGER*4 IDLOG	        ! Numéro identifiant du fichier Trace	    
     	  
      	        		 

C* Initialisation
C-----------------
      IER=0
      LTA=.TRUE.			!épaisseur optique des aérosols nulle
      VA=0.
      VATR=0.
      VR=0.


	 
      
C* Lecture du fichier PROFIL 
C  et stockage de l'épaisseur optique non tronquée
C  -----------------------------------------------
      
      OPEN(2,FILE=FICPROFIL,STATUS='OLD',ERR=912)
      
      DO I=0,NT
c         READ(2,70,err=922)J,ZPROF(I),Hma(I),Ha(I),H(I),XDEL(I),YDEL(I)
         READ(2,70,err=922)J,ZPROF(I),H(I),XDEL(I),YDEL(I)
      ENDDO

      TTOT_VRAI=H(NT)
      
            
C* Adaptation du profil a la troncature si A non nul
C  -------------------------------------------------
      HTR(0)=H(0)

      IF(A.NE.0.) THEN
         DO I=1,NT	
	    VA = XDEL(I)*(H(I)-H(I-1))
	    VATR = VA*(1-PIZ*0.5*A)
	    VR = YDEL(I)*(H(I)-H(I-1))
	    VG = (1-XDEL(I)-YDEL(I))*(H(I)-H(I-1))
	    
	    HTR(I) = (VATR+VR+VG) + HTR(I-1)
	    
	    XDEL(I) = VATR/(VATR+VR+VG)
	    YDEL(I) = VR/(VATR+VR+VG)
         ENDDO
      ENDIF

      DO I=0,NT
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
      CALL SOS_OS(LUM_NBMU,RMU,GA,OS_NB,NT,
     &            FICSURF,FICOS,
     &            N0,TETAS,RHO,IMAT_SURF,IFRESNEL,IND_SURF,
     &            H,XDEL,YDEL,ZPROF,RON,
     &            ALPHA,BETA,GAMMA,ZETA,
     &            ZOUT,IGMAX,IBORM,IPOLAR,TRACE,IDLOG,
     &            EMOINS,EPLUS,IER)
      IF (IER.NE.0) GOTO 995   ! cas d'erreur


C* Definition de l'épaisseur optique au niveau du profil 
C* pour la sortie du champ montant 
C-------------------------------------------------------
      IF (ZOUT.EQ.-1) THEN
	       TAUOUT=H(0)		!Toa.
      ELSE
      
c Recherche du niveau d'altitude inferieur a ZOUT
	    J=1
	    DO WHILE(ZOUT.LT.ZPROF(J))
		   J=J+1
	    ENDDO

C Position de l'altitude par rapport aux niveaux précendents et suivants
	    ZZ = (ZOUT-ZPROF(J-1))/(ZPROF(J)-ZPROF(J-1))

C* Calcul des champs descendants et montants
            TAUOUT=(1-ZZ)*H(J-1)+ZZ*H(J)
	    
      ENDIF

C* Epaisseur optique au niveau du sol (éventuellement tronquée)
      TTOT_TRONC=H(NT)
      
C* Calcul des transmissions directes et diffuses
C-----------------------------------------------

C     Pas d'écriture des resultats (I,Q,U) dans un fichier : FICOS_TRANS_TMP="NO_OUTPUT"
C     Respect des conditions pour simuler une transmission diffuse :
C        RHO = 0 : sol noir (pas de reflexion lambertienne)
C        IMAT_SURF = 0 : sol noir (pas de reflexion speculaire)
C        IFRESNEL = 0 : sol noir (pas de reflexion speculaire)
C        ZOUT = -1 : sortie pour champ down = sol
C        IBORM = 0 : développement limité à l'ordre IS=0
C        TRACE = .FALSE. : on ne complete pas le fichier TRACE 
C
C     Les transmissions diffuses en sortie de SOS_OS sont fournies pour 
C     l'atmosphere equivalente (suite a la troncature de la fonction de phase
C     des aérosols). Elles sont ramenees a la valeur correspondant a 
C     l'atmosphere vraie.

      IF (FICTRANS.NE.'NO_OUTPUT') THEN

	  FICOS_TRANS_TMP = "NO_OUTPUT"	!Fichier non créé

C     --> Transmission diffuse TOA --> Sol pour la direction d'incidence solaire
	  	 	  
	  CALL SOS_OS(LUM_NBMU,RMU,GA,OS_NB,NT,FICSURF,
     &                FICOS_TRANS_TMP,
     &                N0,TETAS,0.D+00,0,0,IND_SURF,
     &                H,XDEL,YDEL,ZPROF,RON,ALPHA,BETA,GAMMA,ZETA,
     &                -1.D+00,IGMAX,0,IPOLAR,.FALSE.,
     &		        IDLOG,TDIFMUS,EPLUS_SOLNOIR,IER)
          IF (IER.NE.0) GOTO 995   ! cas d'erreur	  
	
	  
C     --> Transmission diffuse Sol --> TOA : 
C         Boucle sur les valeurs IMUS : N0 = J	pour calcul TOA --> Sol (réciprocité)  
	  DO J=1,LUM_NBMU
	  
	  
             CALL SOS_OS(LUM_NBMU,RMU,GA,OS_NB,NT,FICSURF,
     &                   FICOS_TRANS_TMP,
     &                   J,TETAS,0.D+00,0,0,IND_SURF,
     &                   H,XDEL,YDEL,ZPROF,RON,ALPHA,BETA,GAMMA,ZETA,
     &                   -1.D+00,IGMAX,0,IPOLAR,.FALSE.,
     &			    IDLOG,TDIFMUS_TMP,EPLUS_SOLNOIR,IER)
             IF (IER.NE.0) GOTO 995   ! cas d'erreur
	      
	     TDIFMUG(J)= TDIFMUS_TMP
	     
          ENDDO	    
	  
      ENDIF !Fin test si fichier transmission à créer      
      
C* Fermeture fichier
C-------------------
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

            
  912 WRITE(6,*) '  ERROR on PROFILE file opening for SOS'
       IER=-1
      GOTO 9999 
      
  913 WRITE(6,*) '  ERROR on AEROSOLS file opening for SOS'
       IER=-1
      GOTO 9999    
      
  922 WRITE(6,*) '  ERROR on PROFILE file reading for SOS'
       IER=-1
      GOTO 9999
      
  923 WRITE(6,*) '  ERROR on AEROSOLS file reading for SOS'
       IER=-1
      GOTO 9999
      
  995 WRITE(6,*) '  ERROR on subroutine SOS_OS'
       IER=-1
      GOTO 9999
                                                                                         
                               
C* Format
C---------   
   50 FORMAT(2X,F9.5)
   60 FORMAT(4(2X,E16.8))
c   70 FORMAT(2X,I5,F10.5,5(F13.7))
   70 FORMAT(2X,I5,F10.5,3(E15.8))
  223 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X,I2)
    
 9999 RETURN
 
      END	!FIN DU PROGRAMME SOS




