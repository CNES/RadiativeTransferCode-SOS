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
C* FICHIER: SOS_PROFIL.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Definition du profil atmosphérique.
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique atmosphérique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, B.Lafrance, X.Lenot
C* DATE: 30/04/02
C* 
C* MOD:VERSION:1.0: Mise a niveau du code des OS du LOA
C*             2.0: Correction d'une application de couche moléculaire sous-jacente 
C*                  pour un profil a 2 couches (iprofil=2 et Zmin=0) : NBSC_C3.
C*             2.1: 07/07/2008 
C*                  Modification du format de lecture des paramètres HR et HA
C*                  pour correction d'une erreur de dï¿œcodage si valeur entière.
C*		      HR et HA passent de F9.5 en F9
C*
C* MOD:VERSION:3.0: 19/11/2009 
C*       - Modification du passage des arguments : 
C*           - passage par couples (mot clef, valeur). 
C*           - adaptation de la gestion des cas d'erreurs. 
C*       - Conversion des messages d'erreurs en anglais.
C*
C* MOD:VERSION:6.0: 06/06/2013
C*   - Prise en compte de l'absorption (cas 3: diffusion moléculaire et aérosols)
C*
C*   - Modification de CTE_DISC pour prendre en compte l'épaisseur optique d'absorption
C*
C*   - Modification du calcul des épaisseurs optiques des couches (calcul vrai N - (N-1))
C*     au lieu du calcul par derivees
C*
C*   - Les paramètres d'entrée sont fournis sous la forme de liste de paramètres
C*     et non plus sous la forme de couples avec mot cle associe => tous les
C*	paramètres deviennent donc obligatoires
C*
C*
C* MOD:VERSION:6.1: 19/02/2016
C*   -  Calcul de l'altitude à chaque niveau du profil:
C*      Initialisation de l'altitude du profil à ZPROF(0)=300.D+00
C*      pour éviter un résultat en Inf qui pose souci à la lecture du fichier. 
C*
C*   -  Ajout d'un test sur la valeur de TA (AOT) par rapport aux nombres
C*      de couches aérosols pour les CAS 2 (profil aérosols seuls) et CAS 4
C*      (profil aérosols homogène entre deux niveaux), pour s'assurer qu'il
C*      est possible de définir un profil AOT variable d'un niveau à l'autre,
C*      par rapport à la précision en 0.00001 du profil d'épaisseurs optiques.
C*
C* MOD:VERSION:6.2: 29/05/2020
C*   -  Modification du format du fichier de profil atmosphérique pour améliorer 
C*      la précision du profil (étiquette de format 20)
C*
C* MOD:VERSION:7.0: 28/09/2020
C*   -  Adaptation pour optimiser la définition du profil en fonction de 
C*      l'épaisseur optique d'absorption. 
C*      Utilisation d'un seuil limite (CTE_SEUIL_TAUABS) jusqu'auquel le profil
C*      est défini par couches fines (d'épaisseur CTE_TCOUCHE).
C*
C* MOD:VERSION:7.1: 31/08/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Conversion des traces en Anglais.
C*
C* MOD:VERSION:7.2: 05/04/2023
C*    - Usage de la constante CTE_TOA_ALT pour définir l'altitude max (TOA)
C*
C*
C* MOD:VERSION:8.0: 27/07/2023
C*    - Refonte complète de la routine pour permettre de gérer les cas de faible
C*      épaisseur optique de diffusion devant celle d'absorption gazeuse.
C*    - Introduction des constantes CTE_OS_NT_MIN et CTE_TOA_FIRST_LAYER_OPT_THICKNESS
C*      (définies dans SOS.h)
C*
C* MOD:VERSION:8.1: 26/09/2023
C*    - Modification de la syntaxe du test sur (TI-TZMOY) dans la routine SOS_DISC pour
C*      compatibilité avec Fortran 2018 pour compilateur gcc/10.2.0
C*
C* MOD:VERSION:8.2: 12/10/2023
C*    - Correction d'une erreur d'affectation des valeurs d'altitude du profil
C*      pour un profil simulant une couche de mélange d’aérosols et molécules, 
C*      entre deux couches moléculaires pures (cas IPROFIL =  2), 
C*    - Correction du test de proximité de l'altitude du dernier niveau calculé
C*      à l'altitude limite (ZLIM) du profil, visant à supprimer ce niveau si très voisin de ZLIM.
C*     
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisées :
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_OS_NT : Nombre max de couches du profil atmosphérique. 
C    CTE_OS_NT_MIN : Nombre min de couches du profil atmosphérique. 
C    CTE_TOA_ALT : Altitude du TOA (en km).
C    CTE_TOA_FIRST_LAYER_OPT_THICKNESS : Epaisseur optique de la première couche atmosphérique, sous TOA.
C    CTE_DELTA_Z : Pas en altitude (en km) pour la recherche du premier niveau 
C                  d'épaisseur optique CTE_TOA_FIRST_LAYER_OPT_THICKNESS.
C    CTE_THRESHOLD_DZ : Seuil (en km) pour la comparaison des altitudes lors de la définition du profil.
C    CTE_ABS_NBLEV : Nombre de niveaux du profil initial en absorption.    
C    CTE_TCOUCHE : Epaisseur optique maximale d'une couche du profil atmosphérique.
C    CTE_THRESHOLD_TAUABS : Seuil sur l'épaisseur optique d'absorption pour optimisation 
C                           du profil atmosphérique.  
C    CTE_DZTRANSI : Epaisseur d'une couche de transition (en km).
C    CTE_PROFIL_MIN_NBC : Nombre minimal de sous-couches des couches 1 et 3
C                         pour un profil d'aérosols entre deux altitudes.       
C----------------------------------------------------------------------------
#include "SOS.h"
 
      


C==============================================================================
C PROGRAMME: SOS_PROFIL
C ==========
C      Ce programme definit le profil atmosphérique en épaisseur optique
C      et en proportion de diffusants (molécules , aérosols) en fonction
C      du niveau, avec prise en compte de l'absorption atmosphérique. 
C      Il est également possible de définir le profil en trois
C      couches: molécules, molécules + aérosols (mélange homogène) et 
C      molécules. La couche intermediaire est repérée par ses altitudes 
C      minimale et maximale.  
C      
C      NB : Le profil n'est pas ajusté pour une troncature de la fonction
C           de phase. Il le sera par le code des OS.
C
C
C Donnees en entrée du programme
C ------------------------------
C     -->Les paramètres/arguments du programme
C
C
C      IPROFIL (I2) : (E) Type de profil des aérosols simulé.
C      		       1 : Profil d'aérosols selon une decroissance exponentielle 
C				    par l'echelle de hauteur
C				2 : Couche homogene d'aérosols entre deux altitudes
C
C      TR (DOUBLE) : (E) Epaisseur optique rayleigh.
C	       
C      HR (DOUBLE) : (E) Echelle de hauteur moléculaire (km).
C                          Associée a TR >= 0.
C  
C      TA (DOUBLE) : (E) Epaisseur optique aérosols (non tronquée).
C	       	       
C      HA (DOUBLE) : (E) Echelle de hauteur du profil d'aérosols (km).
C                          Associée a TA >= 0. et  IPROFIL = 1
C 
C      ZMIN (DOUBLE) : (E) Altitude minimale de la couche d'aérosols (km).
C                          Associée a IPROFIL = 2
C 	   	 
C      ZMAX (DOUBLE) : (E) Altitude maximale de la couche d'aérosols (km).
C                          Associée a IPROFIL = 2
C
C      ABSPROFIL (I2) : (E) Type de profil atmosphérique pour l'absorption
C			         0:User, 1:TROPICAL, 2:MLS, 3:MLW, 4:SAS, 5:SAW,
C                             6:USST62, 7:Pas d'absorption 	 	 
C
C      ALTABS (double) : (E)  Tableau de (CTE_ABS_NBLEV) contenant les niveaux d'altitude du profil
C
C      TABS (double) : (E) Tableau de (CTE_ABS_NBLEV-1) contenant les épaisseurs optiques
C			      d'absorption (continumm + CKD) pour chaque couche atmosphérique
C
C	TRACE (LOGICAL) : (E) Si vrai => écriture des fichiers trace
C
C	IDLOG (I4) : (E) Identifiant du fichier log
C
C      FICPROFIL (CHARACTER) : nom complet du fichier resultat PROFIL
C                              (répertoire + nom fichier + extension)
C
C      NT (I4) : (S) Nombre de couches dans le profil atmosphérique
C
C      IER (I4) : (S) code d'erreur =0 si pas d'erreur, =-1 sinon
C
C
C Resultats fournis par le programme
C ----------------------------------
C   Fichier resultat PROFIL contenant :
C	- l'indice du niveau
C	- l'épaisseur optique molécules + aérosols pour chaque niveau
C      - l'épaisseur optique d'absorption pour chaque niveau
C      - l'épaisseur optique totale pour chaque niveau
C	- le pourcentage d'aérosols par niveau
C	- le pourcentage de molécules par niveau.
C
C   Le nombre de niveaux dans le profil 
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Lecture ou écriture dans un fichier
C     -Paramètres non valables
C
C     Affichage d'un message à l'écran, arret du programme et
C     retour du status 1 au shell
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel aux routines:
C     - SOS_DISC
C
C==============================================================================

      SUBROUTINE SOS_PROFILE(IPROFIL, TR, HR, TA, HA, 
     &		               ZMIN, ZMAX, ABSPROFIL, ALTABS, TABS,
     &		               TRACE, IDLOG, FICPROFIL, NT, IER) 

      
      IMPLICIT NONE


C* Definition des variables                                                 
C--------------------------------------------


C*   Profils
      INTEGER*2 IPROFIL	 ! (E) Type du profil d'aérosol
      DOUBLE PRECISION TTOT  ! épaisseur optique totale .
      INTEGER*4 NT           ! Nombre de couches dans le profil
      INTEGER*4 NT_NG        ! Nombre de couches dans le profil sans absorption gazeuse
      
C*   molécules :
      DOUBLE PRECISION TR    ! (E) Epaisseur optique rayleigh.
      DOUBLE PRECISION HR    ! (E) Echelle de hauteur du profil moléculaire (km).

C*   Aérosols :
      DOUBLE PRECISION TA    ! (E) Epaisseur optique des aérosols (non tronquée).
      DOUBLE PRECISION HA    ! (E) Echelle de hauteur du profil d'aérosols (km).
      DOUBLE PRECISION ZMIN  ! (E) Altitude minimale de la couche molécules + aérosols.
      DOUBLE PRECISION ZMAX  ! (E) Altitude maximale de la couche molécules + aérosols.
      
      DOUBLE PRECISION H(0:CTE_OS_NT)     ! Epaisseur optique à chaque niveau.
      DOUBLE PRECISION PCAER(0:CTE_OS_NT) ! Pourcentage d'aérosols par niveau.
      DOUBLE PRECISION PCMOL(0:CTE_OS_NT) ! Pourcentage de molécules par niveau.

      DOUBLE PRECISION H_NG(0:CTE_OS_NT)     ! Epaisseur optique à chaque niveau du profil sans absorption gazeuse.
      DOUBLE PRECISION PCAER_NG(0:CTE_OS_NT) ! Pourcentage d'aérosols par niveau du profil sans absorption gazeuse.
      DOUBLE PRECISION PCMOL_NG(0:CTE_OS_NT) ! Pourcentage de molécules par niveau du profil sans absorption gazeuse.
      	
C*   Absorption gazeuse
      DOUBLE PRECISION TABS(1:CTE_ABS_NBLEV)     ! Epaisseur optique absorption
      DOUBLE PRECISION ALTABS(1:CTE_ABS_NBLEV)   ! Altitudes du profil d'absorption
      DOUBLE PRECISION ZLIM    ! Valeur minimale d'altitude de définition du profil (non nulle si épaisseur optique d'absorption très forte).
      DOUBLE PRECISION TG_ZLIM ! Epaisseur optique d'absorption à l'altitude ZLIM
      DOUBLE PRECISION ALIN, BLIN ! Coefficients d'interpolation linéaire de l'ep opt gazeuse sur l'altitude
      INTEGER*2 ABSPROFIL    ! (E) Type du profil en absorption


      DOUBLE PRECISION TTOT_ZLIM ! Epaisseur optique totale à l'altitude ZLIM
      
C*   Variables specifiques au profil pour une couche d'aérosols entre deux altitudes :     
      DOUBLE PRECISION VR_C1	! Epaisseur optique moléculaire de la couche 1.	
      DOUBLE PRECISION VR_C2	! Epaisseur optique moléculaire de la couche 2.
      DOUBLE PRECISION VR_C3	! Epaisseur optique moléculaire de la couche 3.
      DOUBLE PRECISION VR_SC	! Epaisseur optique moléculaire entre deux niveaux.

      INTEGER*4 NBSC_C1		! Nombre de sous couches dans la couche 1
      INTEGER*4 NBSC_C2		! Nombre de sous couches dans la couche 2
      INTEGER*4 NBSC_C3 	       ! Nombre de sous couches dans la couche 3
      INTEGER*4 NB_TR	 	! Nombre de couches de transition


      DOUBLE PRECISION DELTA_Z	! Ecart d'altitude entre deux niveaux (en km).
      
C*   Autres variables 
      DOUBLE PRECISION Z 		! Altitude d'un niveau (en km).
             
      DOUBLE PRECISION Hmol(0:CTE_OS_NT)   ! Profil d'épaisseur optique moléculaire.
      DOUBLE PRECISION Haer(0:CTE_OS_NT)   ! Profil d'épaisseur optique des aérosols.
      DOUBLE PRECISION Habs(0:CTE_OS_NT)   ! Profil d'épaisseur optique des gaz absorbants

      DOUBLE PRECISION Hmol_NG(0:CTE_OS_NT)   ! Profil d'épaisseur optique moléculaire sans absorption gazeuse.
      DOUBLE PRECISION Haer_NG(0:CTE_OS_NT)   ! Profil d'épaisseur optique des aérosols sans absorption gazeuse.

      
      DOUBLE PRECISION ZPROF(0:CTE_OS_NT)  ! Altitudes des niveaux du profil (en km)
      DOUBLE PRECISION ZPROF_NG(0:CTE_OS_NT)   ! Altitudes des niveaux du profil sans absorption gazeuse (en km)
      DOUBLE PRECISION Zing			 ! Valeur d'altitude du profil sans gaz
      INTEGER*4 ing                        ! Indice de ZPROF_NG donnant Zing

      DOUBLE PRECISION VA,VR,VG
            
      DOUBLE PRECISION  ZZ  ! Coefficient d'interpolation linéaire entre 2 altitudes	    
      INTEGER*4 I,J  	! Niveau.
      INTEGER*4 II          ! Indice de niveau
      INTEGER*4 K  		! Eléments du profil.en gaz absorbants
      
      INTEGER*4 IER     	! code d'erreur =0 si pas d'erreur, =-1 sinon     
      
      LOGICAL TRACE        	! = vrai,si écriture dans le fichier trace

      INTEGER*4 IDLOG	! Numéro identifiant du fichier Trace	    

      CHARACTER*CTE_LENFIC2 FICPROFIL       ! chemin et nom du fichier resultat PROFIL 

      DOUBLE PRECISION DTAU			!Epaisseur optique du TOA à une altitude donnée

      DOUBLE PRECISION T_FIRST_LAYER 	!Epaisseur optique de la première couche au TOA
      DOUBLE PRECISION T_LAYER 		!Epaisseur optique des couches hors première

      DOUBLE PRECISION T_FIRST_LAYER_NG 	!Epaisseur optique de la couche au TOA (sans absorption gazeuse)
      DOUBLE PRECISION T_LAYER_NG 		!Epaisseur optique des couches hors première (sans absorption gazeuse)

      CHARACTER*100 CAS_TXT_NG		!Texte précisant le cas sans gaz retenu
      CHARACTER*100 CAS_TXT 		!Texte précisant le cas avec gaz retenu
      INTEGER*2 LTXT				!Longueur chaine de caractères

C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0
  
      IF((IPROFIL.NE.1).AND.(IPROFIL.NE.2)) GOTO 940
               
      IF(IPROFIL.EQ.2) GOTO 1000	!On passe les traitements ci-dessous du CAS IPROFIL=1


C===============================================================================
C Etape 1 : Calcul du profil SANS absorption gazeuse 
C ==> Il est calculé également si on définit un profil avec absorption
C ==> On estime les profils: ZPROF_NG, Hmol_NG, Haer_NG
C===============================================================================
     
      !-------------------------------------------------------------------------
      ! Détermination du nombre de couches et des épaisseurs optiques:     
      !  - Cas sans absorption gazeuse : NT_NG  
      !  - Estimation de l'ep. opt. de la première couche au TOA : T_FIRST_LAYER_NG 
      !  - Estimation de l'ep. opt. des couches sous jacentes : T_LAYER_NG
      !-------------------------------------------------------------------------
      IF (CTE_TOA_FIRST_LAYER_OPT_THICKNESS.GT.CTE_TCOUCHE) GOTO 950

      ! Sans absorption gazeuse
      ! =======================
      TTOT = TR+TA
      IF ((TTOT/CTE_OS_NT_MIN).LE. 
     &     CTE_TOA_FIRST_LAYER_OPT_THICKNESS) THEN
         !Cas d'un milieu très peu diffusant: toutes les couches sont prises à même ep. opt même la première  
         CAS_TXT_NG = 
     &   "NO GAS case 1: very low scattering optical thickness"      
         NT_NG = CTE_OS_NT_MIN
         T_LAYER_NG = TTOT / NT_NG
         T_FIRST_LAYER_NG = T_LAYER_NG
         
      ELSE
         IF ((TTOT/CTE_OS_NT_MIN).LT.CTE_TCOUCHE) THEN
            !Cas d'un milieu peu diffusant: 
            !    on garde la première couche distincte avec une ep. opt valant CTE_TOA_FIRST_LAYER_OPT_THICKNESS 
            !    on prend CTE_OS_NT_MIN couches à même ep. opt en dessous (< CTE_TCOUCHE)
            CAS_TXT_NG="NO GAS case 2: low scattering optical thickness"  
            NT_NG = CTE_OS_NT_MIN + 1
            T_FIRST_LAYER_NG = CTE_TOA_FIRST_LAYER_OPT_THICKNESS
            T_LAYER_NG = (TTOT-T_FIRST_LAYER_NG) / CTE_OS_NT_MIN
            
         ELSE
            !Cas TTOT/CTE_TCOUCHE > CTE_OS_NT_MIN
            !    on définit une première couche distincte avec une ep. opt valant CTE_TOA_FIRST_LAYER_OPT_THICKNESS 
            !    on prend (NT_NG-1) couches à même ep. opt en dessous, T_LAYER_NG (proche de CTE_TCOUCHE)
            CAS_TXT_NG="NO GAS case 3: more layers than CTE_OS_NT_MIN-"
            LTXT=INDEX(CAS_TXT_NG,'-')-1 
            CAS_TXT_NG = 
     &      CAS_TXT_NG(1:LTXT)//" for layers of opt. thick.-"
            LTXT=INDEX(CAS_TXT_NG,'-')-1 
            CAS_TXT_NG =  CAS_TXT_NG(1:LTXT)//" CTE_TCOUCHE"    

            T_FIRST_LAYER_NG = CTE_TOA_FIRST_LAYER_OPT_THICKNESS	!Ep. Opt de la 1ère couche

            NT_NG = INT((TTOT-T_FIRST_LAYER_NG)/CTE_TCOUCHE)
            T_LAYER_NG =  (TTOT-T_FIRST_LAYER_NG) / NT_NG	!Ep. Opt des couches sous la couche ultra fine au TOA
	
            NT_NG = NT_NG + 1	!Ajout de la couche ultra fine au TOA

         ENDIF
      ENDIF

    
      !-------------------------------------------------------------------------
      ! Calcul du profil sans absorption gazeuse 
      ! ==> Il est calculé également si on définit un profil avec absorption
      !     pour stocker le profil sans gaz ZPROF_NG qui sera utilisé 
      !     pour la définition du profil avec gaz
      ! ==> Sort le nombre de couches NT_NG 
      !     et les profils Hmol_NG, Haer_NG, H_NG, ZPROF_NG, PCMOL_NG, PCAER_NG
      !-------------------------------------------------------------------------

      TTOT = TR+TA

      IF (TA.EQ.0.0D+00) THEN
          ! ---------------------------------------------------
          ! CAS 1: Cas d'une couche sans diffusion aérosols : 
          !        Pas d'absorption ==> le profil est homogène. 
          ! ---------------------------------------------------
  
          Hmol_NG(0) = 0.D+00
          Hmol_NG(1) = T_FIRST_LAYER_NG

          DO I=2,NT_NG
             Hmol_NG(I) = (I-1)*T_LAYER_NG + T_FIRST_LAYER_NG
          ENDDO

          DO I=0,NT_NG
             PCMOL_NG(I) = 1.000
             PCAER_NG(I) = 0.000
             Haer_NG(I) = 0.D+00
             H_NG(I) = Hmol_NG(I)
          ENDDO

          ZPROF_NG(0) = CTE_TOA_ALT
          DO I=1,NT_NG         
             ZPROF_NG(I)=HR*DLOG(TR/Hmol_NG(I))
          ENDDO
 
      ELSE
          ! ---------------------------------------------------
          ! CAS 2: Cas d'une couche avec diffusion aérosols  
          !        et molécules: Pas d'absorption 
          ! ---------------------------------------------------

          ZPROF_NG(0) = CTE_TOA_ALT
          Hmol_NG(0) = 0.D+00
          Haer_NG(0) = 0.D+00
          H_NG(0) = 0.D+00

          !On cherche l'altitude pour le niveau 1 dont l'épaisseur optique (mol+aer) sera T_FIRST_LAYER_NG
          DTAU = 0.
          Z = CTE_TOA_ALT 
          DO WHILE(DTAU.LT.T_FIRST_LAYER_NG)
               Z = Z - CTE_DELTA_Z !km
               DTAU = TR*DEXP(-Z/HR) + TA*DEXP(-Z/HA) 
	   ENDDO
          ZPROF_NG(1) = Z 
          VR = TR*DEXP(-Z/HR) 
          VA = TA*DEXP(-Z/HA)
          Hmol_NG(1) = VR
          Haer_NG(1) = VA
          H_NG(1) = DTAU
          PCMOL_NG(1) = VR/DTAU
          PCAER_NG(1) = VA/DTAU

          !Duplication des proportions sur le niveau 0 au TOA (couche sous-jacente ultra fine)
          PCMOL_NG(0) = PCMOL_NG(1)
          PCAER_NG(0) = PCAER_NG(1)

          !On cherche l'altitude des niveaux 2 à NT_NG-1
          !et calcule les paramètres de profil pour chaque niveau
          DO I=2,NT_NG-1
               CALL SOS_DISC(T_LAYER_NG,TA,HA,TR,HR,TABS,ALTABS,
     &                       H_NG(I-1),ZPROF_NG(1),0.D+00,0.D+00,Z)	
               ZPROF_NG(I) = Z  
               VR = TR*DEXP(-Z/HR) 
               VA = TA*DEXP(-Z/HA)
               Hmol_NG(I) = VR
               Haer_NG(I) = VA
               H_NG(I) =  VR + VA

               VR = VR - Hmol_NG(I-1)
               VA = VA - Haer_NG(I-1)
               PCMOL_NG(I) = VR / (VR+VA)
               PCAER_NG(I) = VA / (VR+VA)
	   ENDDO 
 
          !Gestion du niveau en surface (à NT_NG)
          ZPROF_NG(NT_NG) = 0.D+00 
          Hmol_NG(NT_NG) = TR
          Haer_NG(NT_NG) = TA
          H_NG(NT_NG) =  TR +TA

          VR = TR - Hmol_NG(NT_NG-1)
          VA = TA - Haer_NG(NT_NG-1)
          PCMOL_NG(NT_NG) = VR / (VR+VA)
          PCAER_NG(NT_NG) = VA / (VR+VA)

      ENDIF !Fin test TA = 0.


      IF ((ABSPROFIL.EQ.7).OR.(TABS(CTE_ABS_NBLEV).EQ.0.D+00)) THEN
	  !En absence d'absorption on stocke les resultats du profil sans gaz 
         !en tant que résultats finaux

        NT = NT_NG
              
         DO I=0,NT
            ZPROF(I) = ZPROF_NG(I)   
	     Hmol(I) = Hmol_NG(I)
	     Haer(I) = Haer_NG(I)
	     Habs(I) = 0.D+00
	     H(I) = Hmol(I) + Haer(I)    
            PCAER(I) = PCAER_NG(I)
            PCMOL(I) = PCMOL_NG(I)
         ENDDO

      ELSE
C===============================================================================
C Etape 2 : Calcul du profil AVEC absorption gazeuse 
C ==> Si présence d'absorption gazeuse
C ==> On estime les profils: ZPROF, Hmol, Haer, PCMOL, PCAER 
C===============================================================================
       


         IF (CTE_TCOUCHE.GT.CTE_THRESHOLD_TAUABS) GOTO 960
    
         !-------------------------------------------------------------------------
         ! Détermination des épaisseurs optiques avec absorption gazeuse:  
         !  - Estimation de l'ep. opt. de la première couche au TOA : T_FIRST_LAYER 
         !  - Estimation de l'ep. opt. des couches sous jacentes : T_LAYER
         !-------------------------------------------------------------------------
  
         ! Epaisseur optique totale sur le trajet
         TTOT = TR + TA + TABS(CTE_ABS_NBLEV)

         IF (TABS(CTE_ABS_NBLEV).GT.CTE_THRESHOLD_TAUABS) THEN   
            !Cas ep. opt. d'absorption forte > seuil limite 
            !Le profil d'absorption est pris en compte jusqu'à l'altitude ZLIM à estimer
            CAS_TXT =  
     &      "WITH GAS case 0 : Strong absorption ==> Limit layer used"
 
            !Détermination de l'altitude pour laquelle l'épaisseur optique d'absorption gazeuse
            !dépasse la valeur seuil
            I=1
            DO WHILE (TABS(I).LT.CTE_THRESHOLD_TAUABS) 
               I=I+1
            ENDDO

            ALIN = (TABS(I) - TABS(I-1)) / (ALTABS(I) - ALTABS(I-1))
            BLIN = TABS(I) - ALIN*ALTABS(I)
          	
            TG_ZLIM = CTE_THRESHOLD_TAUABS
            ZLIM = ( TG_ZLIM - BLIN ) / ALIN


            T_FIRST_LAYER = CTE_TOA_FIRST_LAYER_OPT_THICKNESS

            TTOT_ZLIM = TA*DEXP(-ZLIM/HA)+TR*DEXP(-ZLIM/HR)+TG_ZLIM 
            T_LAYER  = (TTOT_ZLIM-T_FIRST_LAYER) / (CTE_OS_NT-NT_NG-2) !Nb max de couches - celles du profil sans gaz - celle TOA - celle sous limite ep. opt. 
            T_LAYER  = MAX(T_LAYER, CTE_TCOUCHE)
 

         ELSE
            !Cas ep. opt. d'absorption inférieure au seuil limite 
            !Le profil d'absorption est pris en compte jusqu'au sol

            ZLIM = 0.		
            TG_ZLIM = TABS(CTE_ABS_NBLEV)

    
            IF ((TTOT/CTE_OS_NT_MIN).LE. 
     &           CTE_TOA_FIRST_LAYER_OPT_THICKNESS) THEN
               !Cas d'un milieu à très faible ep. opt : toutes les couches sont prises à même ep. opt même la première  
               CAS_TXT =  "WITH GAS case 1: No limit layer -"
               LTXT=INDEX(CAS_TXT,'-')-1 
               CAS_TXT = 
     &         CAS_TXT(1:LTXT)//" + very low optical thickness"

               NT  = CTE_OS_NT_MIN
               T_LAYER  = TTOT / NT 
               T_FIRST_LAYER = T_LAYER
            ELSE
               IF ((TTOT/CTE_OS_NT_MIN).LT.CTE_TCOUCHE) THEN
                  CAS_TXT =  "WITH GAS case 2: No limit layer -"
                  LTXT=INDEX(CAS_TXT,'-')-1 
                  CAS_TXT = CAS_TXT(1:LTXT)//" + low optical thickness"

                  !Cas d'un milieu à faible ep. opt : 
                  !      on garde la première couche distincte avec une ep. opt valant CTE_TOA_FIRST_LAYER_OPT_THICKNESS  
                  !      on prend CTE_OS_NT_MIN couches à même ep. opt en dessous (< CTE_TCOUCHE)
                  NT = CTE_OS_NT_MIN + 1
                  T_FIRST_LAYER = CTE_TOA_FIRST_LAYER_OPT_THICKNESS
                  T_LAYER = (TTOT-T_FIRST_LAYER) / CTE_OS_NT_MIN

               ELSE
                  !Cas TTOT/CTE_TCOUCHE > CTE_OS_NT_MIN
                  CAS_TXT= 
     &            "WITH GAS case 3: more layers than CTE_OS_NT_MIN-"
                  LTXT=INDEX(CAS_TXT,'-')-1 
                  CAS_TXT = 
     &            CAS_TXT(1:LTXT)//" for layers of opt. thick.-"
                  LTXT=INDEX(CAS_TXT,'-')-1 
                  CAS_TXT =  CAS_TXT(1:LTXT)//" CTE_TCOUCHE"

                  T_FIRST_LAYER = CTE_TOA_FIRST_LAYER_OPT_THICKNESS	!Ep. Opt de la 1ère couche

                  NT = INT((TTOT-T_FIRST_LAYER)/CTE_TCOUCHE)	
                  T_LAYER =  (TTOT-T_FIRST_LAYER) / NT		!Ep. Opt des couches sous la couche ultra fine au TOA

                  NT = NT + 1	!+1 pour ajout de la couche ultra fine au TOA

               ENDIF
            ENDIF ! Fin test sur valeur TTOT

         ENDIF !Fin test ep. opt. d'absorption  > seuil limite 
		 

    
         !-------------------------------------------------------------------------
         ! Calcul du profil AVEC absorption gazeuse 
         !
         !    Calcul des altitudes du profil ZPROF(I)
         !    Calcul des épaisseurs optiques Hmol(I), Haer(I) et Habs(I) 
         !    et des proportions de molécules et aérosols (PCMOL(I) et PCAER(I))
         !
         !    Calcul du nombre de couches atmosphériques NT:
         !       Nombre de couches max = CTE_OS_NT
         !       Nombre de couches min = CTE_OS_NT_MIN
         !       Nombre de couches dépendant d'un possible dépassement du seuil d'épaisseur optique d'absorption
         !  ------------------------------------------

         !-------------------------------------------------------------------------
 

         ! Calcul de l'épaisseur optique et de la contribution à la 
         ! diffusion des molécules et des aérosols du TOA (300 km) au sol.
         ! ----------------------------------------------------------------
         NT = 1
         Z = CTE_TOA_ALT
         Zing = ZPROF_NG(1)
         Hmol(0) = 0.D+00
         Haer(0) = 0.D+00
         Habs(0) = 0.D+00
         H(0) = 0.D+00
	  I = 0

         TTOT_ZLIM = TA*DEXP(-ZLIM/HA)+TR*DEXP(-ZLIM/HR)+TG_ZLIM  
		
         DO WHILE ((TTOT_ZLIM-H(NT-1)).GT.T_LAYER)
 		
            I=NT
            IF (I.EQ.1) THEN
               !On cherche l'altitude pour le niveau 1 dont l'épaisseur optique (mol+aer+gaz) sera T_FIRST_LAYER

               DTAU = 0.
               DO WHILE(DTAU.LT.T_FIRST_LAYER)

                  Z = Z - CTE_DELTA_Z

                  !Recherche de l'indice J tel que ALTABS(J-1) > z >= (ALTABS(J)
	           J=2
	           DO WHILE(Z.LT.ALTABS(J))				
		       J=J+1
	           ENDDO
    
                  ! Interpolation linéaire de l'épaisseur optique d'absorption gazeuse VG entre les deux niveaux du profil gaz
                  IF (Z.LE.ALTABS(1)) THEN
                     ZZ=(Z-ALTABS(J-1))/(ALTABS(J)-ALTABS(J-1))
	              VG=(1-ZZ)*TABS(J-1)+ZZ*TABS(J)
                  ELSE
                     VG=0.
	           ENDIF
 

                  ! Epaisseur optique moléculaire et aérosols 
                  VR=TR*DEXP(-Z/HR)
                  VA=TA*DEXP(-Z/HA)       

                  ! Epaisseur optique totale   
	           DTAU = VR + VA + VG
 
	        ENDDO
               ZPROF(1) = Z
               H(1) = DTAU
               ing = 1
         
            ELSE ! Détermination de l'altitude pour l'épaisseur optique au niveau I (pour I > 2)

               CALL SOS_DISC(T_LAYER,TA,HA,TR,HR,TABS,ALTABS,H(I-1),
     &                       ZPROF(1),TG_ZLIM,ZLIM,Z)	      


            ENDIF !Fin distinction cas I=1 ou I>1


	     IF (Z.LE.Zing) THEN
                 Z=Zing	        !Si Z est sous un niveau Zing du profil sans gaz : on force son introduction dans le profil avec gaz
                 ing = ing + 1      !On considère ensuite le niveau Zing suivant (plus bas)
                 Zing = ZPROF_NG(ing) 
            ELSE
                IF ((Z-Zing).LE.CTE_THRESHOLD_DZ) THEN
                   ing = ing + 1           !Z très proche de Zing, on ne va pas forcer Z suivant à valeur Zing: 
                   Zing = ZPROF_NG(ing)    !on modifie Zing en prenant la valeur suivante (plus basse)
                ENDIF               
            ENDIF  



            ZPROF(I)=Z

            ! Calcul de l'épaisseur optique d'absorption par interpolation linéaire
            ! Recherche des niveaux d'altitude du profil en absorption entourant l'altitude Z
	     J=2
	     DO WHILE(Z.LT.ALTABS(J))
		 J=J+1
	     ENDDO
 
            ! Interpolation linéaire entre les deux couches
            IF (Z.GT.ALTABS(1)) THEN
               ! Gestion du cas ou Z est en dehors du profil (Z > 120 km)
		 !print*,"!!!Cas Z=",Z," > 	ALTABS(1) =",ALTABS(1)     
	        VG=TABS(J-1)
	     ELSE
	        ZZ=(Z-ALTABS(J-1))/(ALTABS(J)-ALTABS(J-1))
	        VG=(1-ZZ)*TABS(J-1)+ZZ*TABS(J)
	     ENDIF

	    
            ! Calcul des proportions moléculaires et aérosols.
            VR=TR*DEXP(-Z/HR)
            VA=TA*DEXP(-Z/HA)
	     Hmol(I)=VR
	     Haer(I)=VA
	     Habs(I)=VG
            H(I)=VA+VR+VG

	     VA=VA-Haer(I-1)
            VR=VR-Hmol(I-1)
	     VG=VG-Habs(I-1)
	     PCAER(I)=VA/(VA+VR+VG)
            PCMOL(I)=VR/(VA+VR+VG)

            NT = NT + 1

         ENDDO !Fin boucle estimation du niveau Z tant que Z > T_LAYER 

         !Gestion du niveau limite de profondeur (ZLIM!=0. ou surface)
         IF ((ZPROF(NT-1)-ZLIM).LE.CTE_THRESHOLD_DZ) THEN
            NT=NT-1	!Si l'écart est faible entre Z et ZLIM, on ne tient pas compte du dernier niveau calculé
         ENDIF
         ZPROF(NT) = ZLIM  
         VR = TR*DEXP(-ZLIM /HR)
         VA = TA*DEXP(-ZLIM /HA)
         VG = TG_ZLIM
         Hmol(NT) = VR
         Haer(NT) = VA
         Habs(NT) = VG
         H(NT) =  VR + VA + TG_ZLIM

         VA=VA-Haer(NT-1)
         VR=VR-Hmol(NT-1)
	  VG=VG-Habs(NT-1)
	  PCAER(NT)=VA/(VA+VR+VG)
         PCMOL(NT)=VR/(VA+VR+VG)

	  ! Définition des valeurs du profil au TOA 
         ZPROF(0)= CTE_TOA_ALT
         PCAER(0) = PCAER(1)
         PCMOL(0) = PCMOL(1)
	  Hmol(0) = 0.
	  Haer(0) = 0.
	  Habs(0) = 0.
         H(0) = 0.

	  

         IF (TABS(CTE_ABS_NBLEV).GT.CTE_THRESHOLD_TAUABS) THEN   
            !On complète le profil du niveau dont l'épaisseur optique gazeuse atteint le seuil
            !au sol en ajoutant une couche supplémentaire.
            !On affecte à cette dernière couche (du niveau limite au sol) la totalité des contenus sous le niveau limite.
            !==> suppose une atténuation telle que cette couche n'est pas visible.
            NT = NT + 1
            Hmol(NT)=TR
            Haer(NT)=TA
            Habs(NT)=TABS(CTE_ABS_NBLEV)
            H(NT)=Hmol(NT)+Haer(NT)+Habs(NT)

            VR=Hmol(NT)-Hmol(NT-1)
            VA=Haer(NT)-Haer(NT-1)
            VG=Habs(NT)-Habs(NT-1)

            PCAER(NT)=VA/(VA+VR+VG)
            PCMOL(NT)=VR/(VA+VR+VG)

            ZPROF(NT)=0.D+00

         ENDIF

      ENDIF !Fin test cas sans/avec absorption gazeuse


      GOTO 7777 !On passe le cas IPROFIL=2



     
1000  CONTINUE
C======================================================
C* CAS IPROFIL=2 : 
C*   Diffusion des aérosols et des molécules.
C*   Couche uniforme d'aérosols entre deux altitudes
C*   Cas sans absorption gazeuse.
C======================================================
C*      
C*		T O A
C*	___________________________ 
C*
C*		couche 1
C*	--------------------------- Zmax+DZtransi
C*	___________________________ Zmax
C*
C*
C*		couche 2
C*
C*	___________________________ Zmin	
C*	--------------------------- Zmin-DZtransi
C*
C*		couche 3
C*	___________________________

	         
	 
C* Controle de la validite des arguments 
C---------------------------------------
      IF ((ZMIN.LT.0.D-00).OR.(ZMAX.LE.ZMIN)) GOTO 1010


C* Nombre de couches du profil 
C---------------------------------------
      TTOT = TR + TA 
      NT = int(TTOT / CTE_TCOUCHE) 
      IF (NT.GT.CTE_OS_NT) NT = CTE_OS_NT
          
 
C* Epaisseur optique moléculaire des couches
C  -----------------------------------------
      VR_C1=TR*DEXP(-(ZMAX+CTE_DZTRANSI)/HR)
      VR_C2=TR*( DEXP(-ZMIN/HR)-DEXP(-(ZMAX+CTE_DZTRANSI)/HR) )
      IF(ZMIN.EQ.0.D-00) THEN
	  VR_C3 = 0.D+00
	  NB_TR = 1
      ELSE   
	  VR_C3=TR*( 1.D+00-DEXP(-(ZMIN-CTE_DZTRANSI)/HR) )
	  NB_TR = 2     
      ENDIF 
	   
C* Nombre de sous-couches dans chaque couche
C  -----------------------------------------

      NBSC_C1 = (NT-NB_TR)*VR_C1/(TR+TA) 
      NBSC_C1 = MAX( CTE_PROFIL_MIN_NBC, NBSC_C1 )
	 
      IF(ZMIN.EQ.0.D+00) THEN
	  NBSC_C3 = 0	    
      ELSE
	  NBSC_C3 = (NT-NB_TR)*VR_C3/(TR+TA)	 
	  NBSC_C3 = MAX( CTE_PROFIL_MIN_NBC, NBSC_C3 )
      ENDIF
	 	 
      NBSC_C2 = (NT-NB_TR)-NBSC_C1-NBSC_C3

      IF((TA/NBSC_C2).LT.0.00001) GOTO 1020

C* Profil moléculaire : couche 1
C  -----------------------------
      VR_SC=VR_C1/NBSC_C1

      DO I=1,NBSC_C1
	  Hmol(I)=Hmol(I-1)+VR_SC
	  Haer(I)=0.
      ENDDO
	       
C* Transition : ZMAX+CTE_DZTRANSI -> ZMAX
C  -------------------------------------
      I=NBSC_C1+1
      Hmol(I)=TR*DEXP(-ZMAX/HR)
      VR_SC=Hmol(I)-Hmol(I-1)
      Haer(I)=Haer(I-1)+(TA*VR_SC/VR_C2)
	      	       
C* Melange homogène: aérosols + molécules: couche 2
C  ------------------------------------------------
      DELTA_Z=(ZMAX-ZMIN)/NBSC_C2
	       
      !Initialisation
      Z=ZMAX
	       
      DO I = NBSC_C1+2 , NBSC_C1+NBSC_C2+1
	         
	  Z=Z-DELTA_Z
	  Hmol(I)=TR*DEXP(-Z/HR)
	  VR_SC=Hmol(I)-Hmol(I-1)
	  Haer(I)=Haer(I-1)+(TA*VR_SC/VR_C2)
	        	  
      ENDDO

	
      IF (ZMIN.NE.0.D+00) THEN       	       
	    
C* Transition : ZMIN -> ZMIN-CTE_DZTRANSI
C  -------------------------------------
	  I=NBSC_C1+NBSC_C2+2
	  Hmol(I)=TR*DEXP(-(ZMIN-CTE_DZTRANSI)/HR)
	  Haer(I)=Haer(I-1)  

C* Profil moléculaire : couche 3
C  -----------------------------
	       
	  VR_SC=VR_C3/NBSC_C3
	       
	  DO I = NBSC_C1+NBSC_C2+3 , NT
            Hmol(I)=VR_SC+Hmol(I-1)
            Haer(I)=Haer(I-1)
	  ENDDO
	       
      ENDIF !(Zmin.ne.0)


	       	       	       
C* Calcul de H, PCAER, PCMOL
C  -----------------------	
      ZPROF(0)=CTE_TOA_ALT
      Hmol(0)=TR*DEXP(-CTE_TOA_ALT/HR)
      Haer(0)=0.
      H(0)=Hmol(0)+Haer(0)     
      PCAER(0)=0.
      PCMOL(0)=1.
         
	
      DO I=1,NT
	       	  
         H(I)=Hmol(I)+Haer(I)
         ZPROF(I)=HR*DLOG(TR/Hmol(I))

         IF(Haer(I).EQ.Haer(I-1)) THEN
            PCAER(I)=0.
            PCMOL(I)=1
         ELSE
            PCMOL(I)=1/(1+(TA/VR_C2))
            PCAER(I)=1-PCMOL(I)
         ENDIF
		    
      ENDDO
	    
C FIN DU CAS IPROFIL = 2   
C------------------------------------------------------ 
 7777 CONTINUE



C* Ecriture du fichier trace
C-------------------------------------------
      
      IF (TRACE) THEN      

         WRITE(IDLOG,*,err=921) "****************"
         WRITE(IDLOG,*,err=921) "SOS_PROFIL log :"
         WRITE(IDLOG,*,err=921) ' '
  
         IF(IPROFIL.EQ.1) THEN
            
	     WRITE(IDLOG,*,err=921) "  Profile without gaz"
	     WRITE(IDLOG,*,err=921) "  -------------------"                   
	     WRITE(IDLOG,*,err=921) "    Case : ", CAS_TXT_NG   
	     WRITE(IDLOG,*,err=921) 
     &         "    Number of layers of the profile: ",NT_NG         
	     WRITE(IDLOG,*,err=921)
     &      "    First layer opt. thickness : ",T_FIRST_LAYER_NG        
	     WRITE(IDLOG,*,err=921)
     &      "    Min Opt. thickness of layers below: ",T_LAYER_NG

  	     WRITE(IDLOG,*,err=921) ' '
  	     WRITE(IDLOG,*,err=921)   
     &      "    Total optical thickness of molecules: ",Hmol_NG(NT_NG)
  	     WRITE(IDLOG,*,err=921)   
     &      "    Total optical thickness of aerosols: ",Haer_NG(NT_NG)
  	     WRITE(IDLOG,*,err=921)   
     &      "    Total optical thickness: ",H_NG(NT_NG)
  	     WRITE(IDLOG,*,err=921) ' '
	     WRITE(IDLOG,10,err=921)  
     &             'I','ALT','Hmol(I)','Haer(I)', 
     &             'Habs(I)','H(I)','PCAER(I)','PCMOL(I)'
	     DO I=0,NT_NG          
	        WRITE(IDLOG,16,err=921)  
     &               I,ZPROF_NG(I),Hmol_NG(I),Haer_NG(I),0.D+00, 
     &               H_NG(I),PCAER_NG(I),PCMOL_NG(I)
	     ENDDO

            IF ((ABSPROFIL.NE.7).AND. 
     &          (TABS(CTE_ABS_NBLEV).NE.0.D+00)) THEN

  	        WRITE(IDLOG,*,err=921) ' '
  	        WRITE(IDLOG,*,err=921) ' '
	        WRITE(IDLOG,*,err=921) 
     &         "  Profile with gaz (used for RT calculation)"
	        WRITE(IDLOG,*,err=921) 
     &         "  ------------------------------------------"                 
	        WRITE(IDLOG,*,err=921) "    Case : ", CAS_TXT         
	        WRITE(IDLOG,*,err=921) 
     &         "    Number of layers  of the profile: ",NT       
	        WRITE(IDLOG,*,err=921)
     &         "    First layer opt. thickness : ",T_FIRST_LAYER      
	        WRITE(IDLOG,*,err=921)
     &         "    Min Opt. thickness of layers below: ",T_LAYER

   	        WRITE(IDLOG,*,err=921) ' '
  	        WRITE(IDLOG,*,err=921)   
     &         "    Total optical thickness of molecules: ",Hmol(NT)
  	        WRITE(IDLOG,*,err=921)   
     &         "    Total optical thickness of aerosols: ",Haer(NT)
  	        WRITE(IDLOG,*,err=921)   
     &         "    Total optical thickness of gas absorption: ",  
     &         Habs(NT)
  	        WRITE(IDLOG,*,err=921)   
     &         "    Total optical thickness: ",H(NT)
 
   	        WRITE(IDLOG,*,err=921) ' '
  	        WRITE(IDLOG,*,err=921)   
     &         "    Min altitude (km) for profile definition: ",ZLIM
  	        WRITE(IDLOG,*,err=921)   
     &         "    (If not surface : an additional layer is added"
  	        WRITE(IDLOG,*,err=921)   
     &         "     to complete the profile)"
 
  	        WRITE(IDLOG,*,err=921) ' '
	        WRITE(IDLOG,10,err=921)  
     &             'I','ALT','Hmol(I)','Haer(I)', 
     &             'Habs(I)','H(I)','PCAER(I)','PCMOL(I)'
	        DO I=0,NT          
	           WRITE(IDLOG,16,err=921)  
     &               I,ZPROF(I),Hmol(I),Haer(I),Habs(I), 
     &               H(I),PCAER(I),PCMOL(I)
	        ENDDO
            ENDIF !Fin cas avec absorption gazeuse



         ELSE !Cas IPROFIL=2

	     WRITE(IDLOG,*,err=921)     
     &      '  Profile for a variable altitude of the aerosol layer'
	     WRITE(IDLOG,*,err=921)     
     &      '  Molecular layer + aerosols : Minimal altitude  : ',ZMIN
	     WRITE(IDLOG,*,err=921)     
     &      '  Molecular layer + aerosols : Maximal altitude : ',ZMAX
            WRITE(IDLOG,*,err=921) ' '
	     WRITE(IDLOG,*,err=921) 
     &	     '  Number of transition layers      : ',NB_TR
	     WRITE(IDLOG,*,err=921) 
     &	     '  Number of sub-layers of the layer 1 : ',NBSC_C1
	     WRITE(IDLOG,*,err=921) 
     &	     '  Number of sub-layers of the layer 2 : ',NBSC_C2
	     WRITE(IDLOG,*,err=921)
     &	     '  Number of sub-layers of the layer 3 : ',NBSC_C3

   	     WRITE(IDLOG,*,err=921) ' '
	     WRITE(IDLOG,10,err=921)'I','ALT','Hmol(I)','Haer(I)',
     &      'Habs(I)','H(I)','PCAER(I)','PCMOL(I)'
	     DO I=0,NT	        
	        WRITE(IDLOG,15,err=921) I,ZPROF(I),Hmol(I),Haer(I),
     &                                  0.D+00,H(I),PCAER(I),PCMOL(I)
	     ENDDO

         ENDIF       
      
         WRITE(IDLOG,*,err=921) ' '
         WRITE(IDLOG,*,err=921) "****************"
	 	   	   	   	   	   
      ENDIF ! Trace





      
C* Enregistrement des paramètres du profil atmosphérique dans un fichier
C  ---------------------------------------------------------------------
      OPEN(2,FILE=FICPROFIL,ERR=991)
      
      DO I=0,NT
         WRITE(2,20,err=992)I,ZPROF(I),H(I),PCAER(I),PCMOL(I)
      ENDDO

C* Fermeture fichier
C-------------------
      CLOSE(2)		! Fermeture du fichier resultat PROFIL

C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------

  940 WRITE(6,*) '  ERROR 940 in SOS_PROFIL :'
      WRITE(6,*) '  IPROFIL must be 1 or 2' 
      WRITE(6,*) '  Current value :', IPROFIL  
      GOTO 9998

  950 WRITE(6,*) '  ERROR 950 in SOS_PROFIL :'
      WRITE(6,*) '  CTE_TOA_FIRST_LAYER_OPT_THICKNESS must be ' 
      WRITE(6,*) '  lower than CTE_TCOUCHE (see inc/SOS.h file)'
      WRITE(6,*) '     CTE_TOA_FIRST_LAYER_OPT_THICKNESS:', 
     &            CTE_TOA_FIRST_LAYER_OPT_THICKNESS    
      WRITE(6,*) '     CTE_TCOUCHE:', CTE_TCOUCHE  
      GOTO 9998

  960 WRITE(6,*) '  ERROR 960 in SOS_PROFIL :'
      WRITE(6,*) 
     &'  CTE_TCOUCHE must be lower than CTE_THRESHOLD_TAUABS'
      WRITE(6,*) '  (see inc/SOS.h file)'
      WRITE(6,*) '     CTE_TCOUCHE:', CTE_TCOUCHE  
      WRITE(6,*) 
     &'     CTE_THRESHOLD_TAUABS:', CTE_THRESHOLD_TAUABS 
      GOTO 9998

  991 WRITE(6,*) '  ERROR on SOS_PROFIL result file opening'
      GOTO 9998

  992 WRITE(6,*) '  ERROR on SOS_PROFIL result file writing'
      GOTO 9998

  912 WRITE(6,*) '  ERROR on logfile opening for SOS_PROFIL'
      GOTO 9998

  921 WRITE(6,*) '  ERROR on logfile writing for SOS_PROFIL'
      GOTO 9998

 1010 WRITE(6,*) '  ERROR on boundary altitudes of the aerosol layer :'
      WRITE(6,*) '  --> check -AP.AerLayer.Zmin and AP.AerLayer.Zmax'
      GOTO 9998

 1020 WRITE(6,*) '  AOT = ',TA, ' is lower than the minimal value',
     &           ' for the profile definition (0.00001)'
      WRITE(6,*) '  ==> Use AER.AOTref=0 if a very low value is wanted'

 9998 IER=-1
      GOTO 9999

                   
C* Format
C---------   
c   10 FORMAT(5X,A1,4X,A3,7X,A7,10X,A7,10X,A7,10X,A4,12X,2(A7,11X))
   10 FORMAT(5X,A1,4X,A3,7X,A7,3X,A7,3X,A7,3X,A4,6X,2(A8,2X))
   15 FORMAT(2X,I4,F11.5,6(F10.5))
   16 FORMAT(2X,I4,F10.5,1X,6(E15.8,2X))	
   20 FORMAT(2X,I5,F10.5,3(E15.8))
   21 FORMAT(2X,I4,6(E15.8,2X))	
   30 FORMAT(5X,A1,2X,A7,2X,A7,2X,A7,2X,A4,5X,2(A7,2X))
   35 FORMAT(2X,I5,6(F10.5))
   
 9999 RETURN
      END	!FIN DU PROGRAMME SOS_PROFIL




C==============================================================================
C PROCEDURE: SOS_DISC
C ==========
C      Cette procédure estime l'altitude à attribuer au niveau I du profil
C      atmosphérique.
C
C      A) On estime l'épaisseur optique du niveau I en fonction de celle du niveau I-1 
C         et du pas en épaisseur optique.
C
C      B) On évalue l'altitude correspondante à ce niveau I par une méthode de
C         dichotomie, avec le calcul de l'épaisseur optique d'absorption gazeuse et
C         de celles moléculaires et aérosols pour les altitudes Z testées.
C         On sort de la boucle quand l'épaisseur optique de l'altitude Z 
C         est très proche de celle souhaitée au niveau I.
C
C Description des paramètres
C --------------------------
C   molécules :
C      TR  (double)  (E)   Epaisseur optique rayleigh.
C      HR  (double)  (E)   Echelle de hauteur du profil moléculaire (km).
C
C   Aérosols :
C      TA  (double)  (E)   Epaisseur optique des aérosols.
C      HA  (double)  (E)   Echelle de hauteur du profil d'aérosols (km).
C
C   Gaz :
C      TABS(1:CTE_ABS_NBLEV)   (double)  (E)   Profil initial en altitude de l'épaisseur optique absorption
C      ALTABS(1:CTE_ABS_NBLEV) (double)  (E)   Altitudes du profil d'absorption initial
C      TG_ZLIM (double)  (E) Epaisseur optique d'absorption à l'altitude ZLIM
C      ZLIM (double)  (E) Valeur minimale d'altitude à considérer (en km).
C      ZMAX_INIT  (double)  (E)  Altitude max à considérer (en km).
C
C   Informations sur le niveau I-1 :
C      I     (I4)      (E)  Numéro du niveau traite.
C      TIM1  (double)  (E)  Epaisseur optique totale du niveau I-1.
C      DT    (double)  (E)  Pas en épaisseur optique
C
C      ZX    (double)  (S)  Altitude estimée pour le niveau I.
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
      SUBROUTINE SOS_DISC(DT,TA,HA,TR,HR,TABS,ALTABS,TIM1,
     &                    ZMAX_INIT, TG_ZLIM, ZLIM, ZX)

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------


C*   molécules :
      DOUBLE PRECISION TR	! Epaisseur optique rayleigh.
      DOUBLE PRECISION HR	! Echelle de hauteur du profil moléculaire.

C*   Aérosols :
      DOUBLE PRECISION TA	! Epaisseur optique des aérosols (non tronquée).
      DOUBLE PRECISION HA	! Echelle de hauteur du profil d'aérosols.

C*   Gaz absorbants
      DOUBLE PRECISION TABS(1:CTE_ABS_NBLEV)   ! Profil en altitude de l'épaisseur optique absorption
      DOUBLE PRECISION ALTABS(1:CTE_ABS_NBLEV) ! Altitudes du profil d'absorption
      DOUBLE PRECISION TG      ! Epaisseur optique d'absorption
      DOUBLE PRECISION ZLIM    ! Valeur minimale d'altitude pour la recherche dichotomique.
      DOUBLE PRECISION TG_ZLIM ! Epaisseur optique d'absorption à l'altitude ZLIM

       
C*   Informations sur le niveau I-1 :
      DOUBLE PRECISION TIM1	! Epaisseur optique totale du niveau I-1.

      DOUBLE PRECISION ZX	! Altitude estimée pour le niveau I.
      	
      DOUBLE PRECISION ZZ	! Coefficient d'interpolation linéaire entre 2 altitudes
      	
      DOUBLE PRECISION DT	! Pas en épaisseur optique.
      DOUBLE PRECISION TI	! Epaisseur optique estimee au niveau I.
      DOUBLE PRECISION ZMAX	! Altitude (km).
      DOUBLE PRECISION ZMAX_INIT	! Altitude (km).
      DOUBLE PRECISION ZMIN


      DOUBLE PRECISION ZMOY	! Altitude moyenne.
      DOUBLE PRECISION TZMOY ! Epaisseur optique totale pour l'altitude Zmoy.
      DOUBLE PRECISION XD	! Ecart absolu.
      
      INTEGER*4 J		 ! Numéro du niveau traité.


C------------------------------------------------------
C* A) Estimation de l'épaisseur optique du niveau I en fonction
C*    de celle du niveau I-1 et du pas DT.
C------------------------------------------------------     
      TI=TIM1+DT

C------------------------------------------------------
C* B) Estimation de l'altitude correspondant à l'épaisseur 
C*    optique TI du niveau I-1.
C------------------------------------------------------

C*   Initialisation des altitudes .
C--------------------------------------
      ZMAX=ZMAX_INIT
      ZMIN=ZLIM

C*   Initialisation
C---------------------------
  706 ZMOY=(ZMAX+ZMIN)/2.
     

C*   Calcul de l'épaisseur optique d'absorption à l'altitude ZMOY
C-----------------------------------------------------------------
      IF (TG_ZLIM.GT.0.0D+00) THEN
         J=2
         DO WHILE(ZMOY.LT.ALTABS(J))
	     J=J+1
         ENDDO
	    
         ! Interpolation linéaire entre les deux couches
         IF (ZMOY.GT.ALTABS(1)) THEN !Gestion du cas ou Z est en dehors du profil (Z > 120 km)	    
	     ZZ=0
         ELSE
	     ZZ=(ZMOY-ALTABS(J-1))/(ALTABS(J)-ALTABS(J-1))
         ENDIF
         TG=(1-ZZ)*TABS(J-1)+ZZ*TABS(J)
      
      ELSE
         TG=0.0D+00
      ENDIF	  

	
  
C*   Epaisseur optique totale pour l'altitude Zmoy.
C----------------------------------------------------
      TZMOY=TA*DEXP(-ZMOY/HA)+TR*DEXP(-ZMOY/HR)+TG
 
	        
C*   Ecart absolu avec l'épaisseur optique estimee du niveau I.
C-------------------------------------------------------------
      XD=DABS(TI-TZMOY)
  
C*   Test de sortie si Zmoy correspond au niveau I.
C------------------------------------------------------
      IF(XD.LT..000001) GOTO 705
      IF(ZMOY.EQ..00D+00) GOTO 705
      
C*   Affinage de l'encadrement de l'altitude du niveau I.
C------------------------------------------------------
      IF( (TI-TZMOY).LT.0.D+00 ) THEN
         ZMIN=ZMOY
      ELSE
         ZMAX=ZMOY
      ENDIF
      
      GOTO 706    
       
C*   Affectation du résultat
C------------------------------------------------------
 705  ZX=ZMOY
 
 
      RETURN
      END	!FIN DE LA PROCEDURE SOS_DISC



