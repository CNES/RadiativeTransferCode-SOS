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
C* FICHIER: SOS_PREPA_ABSPROFILE.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Préparation pour le calcul du profil atmosphérique en absorption
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, X.Lenot, B.Lafrance
C* 
C*
C* MOD:VERSION:1.0: DATE: 06/06/2013
C*    Creation du fichier pour adapter les OS a la méthode CKD (absorption atmosphérique)
C*
C* MOD:VERSION:1.1: DATE: 24/02/2016
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*    - Tableau NEXP passe de Double Precision à Integer
C*
C* MOD:VERSION:1.2: DATE: 06/04/2019
C*    - Correction de l'information sur le CH4 dans la table ATMOCM
C       ==> ATMOCM(6) passe de la valeur erronée 13.837E+22 à 3.837E+22
C*
C* MOD:VERSION:2.0: DATE: 02/09/2021
C*    - Ajustements pour la suppression de la distinction des raies et continuums
C*      d'absorption avec l'utilisation de nouvelles données CKD mixant les deux.
C*      ==> Suppression du tableur TAUCONT qui donnait les épaisseurs optiques
C*          de continuum d'absorption.
C*      ==> Extension de la table des espèces du tableur IABS pour inclure les gaz
C*          dont l'absorption est par continuums uniquement (NO2 seulement)
C*          et suppression de la table ICONT.
C*      ==> Suppression des tables AA et BB de correction de la dépendance à la 
C*          température des coefficients CKD.
C*      ==> Le tableau XKI est remplacé par KDIS_KI avec la variable Température en plus.
C*      ==> Le tableau AL est remplacé par KDIS_AI
C*      ==> Ajustement des paramètres d'E/S.
C*      ==> Calcul de l'indice de longueur d'onde des tableaux de coefficients CKD
C*          en fonction du pas spectral retenu et des bornes spectrales des données CKD
C*
C*    - Correction de la pondération des profils par les quantités utilisateur de H2O et O3.
C*
C*    - Ajout de la possibilité de pondérer les quantités de CH4 et CO2 
C*      ==> Ajout des paramètres CH4 et CO2.
C*      ==> Pondération des profils CH4 et CO2 par le rapport entre la quantité en surface 
C*          définie par l'utilisateur et celle du profil par défaut.
C*
C* MOD:VERSION:2.1: 07/07/2022
C*    - Suppression du passage par la variable d'environnement SOS_ABS_FIC.
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Conversion des traces en Anglais.
C*
C* MOD:VERSION:2.2: 16/06/2023
C*    - Ajout des paramètres KDIS_KI_H2O,TAB_CONC_H2O et NB_CONC_H2O à l'appel de la routine READ_CKD_COEFF
C******************************************************************************

C----------------------------------------------------------------------------
C Définition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_ABS_NBLEV : nombre de niveaux du profil en absorption
c    CTE_NBABS : nombre d'absorbants     
C    CTE_ABS_NBCOL : nombre de colonnes dans un profil d'absorption fourni
C		       en entrée  
C    CTE_ABS_H2O  : Indicateur de prise en compte du H2O (0 non / 1 oui)
C    CTE_ABS_CO2  : Indicateur de prise en compte du CO2 (0 non / 1 oui)
C    CTE_ABS_O3   : Indicateur de prise en compte du O3  (0 non / 1 oui)
C    CTE_ABS_N2O  : Indicateur de prise en compte du N2O (0 non / 1 oui)
C    CTE_ABS_CO   : Indicateur de prise en compte du CO  (0 non / 1 oui)
C    CTE_ABS_CH4  : Indicateur de prise en compte du CH4 (0 non / 1 oui)
C    CTE_ABS_O2   : Indicateur de prise en compte du O2  (0 non / 1 oui)
C    CTE_CONT_NO2 : Indicateur de prise en compte du NO2 (0 non / 1 oui)
C
C

C----------------------------------------------------------------------------
#include "SOS.h"



C==============================================================================
C PROGRAMME: SOS_PREPA_ABSPROFILE
C ==========
C      Ce programme prépare les calculs pour le profil en absorption:
C
C        - Il lit le fichier de profil : fichier utilisateur ou fichiers prédéfinis
C          NB: les fichiers prédéfinis distinguent le profil du continuum NO2 des autres gaz
C              --> fichier SO2-NO2 spécifique, indépendant du type de profil.
C
C        - Il convertit les données du profil en nombre de particules/cm2 par couche,
C          compatible avec la dimension des coefficients ki en cm2/particule.
C
C        - Il surcharge les profils H2O et O3 des contenus intégrés fournis par l'utilisateur.
C
C        - Il surcharge les profils CH4 et CO2 par le rapport entre la quantité en surface 
C          définie par l'utilisateur et celle du profil par défaut.
C
C        - Il lit, pour toutes les longueurs d'onde du domaine couvert par les fichiers CKD,
C          les coefficients des fichiers CKD :
C             * Nombre d'exponentielles de la décomposition CKD de chaque absorbant
C               pour la gamme de longueurs d'onde couverte.
C             * coefficients Ki et Ai pour chaque exponentielle.
C
C        - Il fournit l'indice de tableau correspondant à la longueur d'onde traitée
C
C
C Données en entrée du programme
C ------------------------------
C     -->Les paramètres/arguments du programme
C
C
C     WA (DOUBLE) : (E) Longueur d'onde (microns).
C
C     PSURF (DOUBLE)  : (E) Pression atmosphérique de surface (mbar).
C
C     H2O   (DOUBLE)  : (E) Contenu intégré en vapeur d'eau (g/cm2).
C
C     O3   (DOUBLE)  : (E) Contenu intégré en ozone (dobson).
C
C     CO2  (DOUBLE)  : (E) Concentration du dioxyde de carbone en surface (ppmv).
C
C     CH4  (DOUBLE)  : (E) Concentration du méthane en surface (ppmv).
C
C     ABSPROFIL (I2) : (E) Type de profil atmosphérique pour l'absorption
C		         0:User, 1:TROPICAL, 2:MLS, 3:MLW, 4:SAS, 5:SAW,
C                      6:USST62, 7:Pas d'absorption 	 	 
C
C     FICABSPROFIL (CHARACTER) : (E) associe a ABSPROFIL = 0
C				  nom complet du fichier du profil d'absorption atmosphérique
C				  (répertoire + nom fichier + extension)
C
C     TRACE (LOGICAL) : (E) Si vrai => écriture des fichiers trace
C
C     IDLOG (I4) : (E) Identifiant du fichier log
C
C     NU  (double) : (S)   Nombre d'onde en cm-1
C
C     LAMB1  (I4) : (S)   Indice de l'intervalle CKD
C
C     IABS  (I2) : (S) Tableau de CTE_NBABS éléments. Chaque élément correspond à un 
C			  gaz absorbant (Absorbants: H2O CO2 O3 N2O CO CH4 O2 NO2).
C	 		  Pour chaque absorbant, vaut 1 si traité avec méthode CKD, 0 sinon 
C			   
C
C     USERPROFIL (double) (S)  Tableau de (CTE_ABS_NBLEV * CTE_ABS_NBCOL) éléments
C			    contenant pour chaque gaz absorbant le profil exprimé en ppmv
C			    Il contient dans l'ordre : 
C                  		USERPROFIL(J, 1) : altitude (km)
C                  		USERPROFIL(J, 2) : pression (mbar)
C                  		USERPROFIL(J, 3) : température (K)
C                  		USERPROFIL(J, 4) : H2O (ppmv)
C                  		USERPROFIL(J, 5) : CO2 (ppmv)
C                  		USERPROFIL(J, 6) : O3  (ppmv)
C                  		USERPROFIL(J, 7) : N2O (ppmv)
C                  		USERPROFIL(J, 8) : CO  (ppmv)
C                  		USERPROFIL(J, 9) : CH4 (ppmv)
C                  		USERPROFIL(J, 10): O2  (ppmv)
C                 		USERPROFIL(J, 11): densité 		(not used)
C                  		USERPROFIL(J, 12): NO2 (ppmv)
C                  		USERPROFIL(J, 13): SO2 (ppmv)	(not used)
C
C
C     ALTABS (double) : (S)  Tableau de (CTE_ABS_NBLEV) contenant les niveaux d'altitude du profil
C
C     RO (double) (S) : Tableau de (CTE_NBABS * CTE_ABS_NBLEV) éléments
C		          contenant, pour chaque gaz absorbant de la méthode CKD, le profil 
C			   exprimé en particules/cm2 pour chaque couche.
C			   Il contient dans l'ordre : H2O, CO2, 03, N2O, CO, CH4, O2, NO2
C                       --> Gaz qui ont une absorption par raies fines (pas uniquement par continnums)
C
C     NEXP (I4) :   (S)  Tableau de  (CTE_NBABS * CTE_CKD_NWVL_MAX) éléments contenant pour chaque absorbant
C		           et chaque intervalle spectral le nombre d'exponentielles de la méthode CKD
C
C
C     KDIS_AI (double) : (S) Tableau de (CTE_CKD_NAI_MAX * CTE_NBABS * CTE_CKD_NWVL_MAX) éléments contenant 
C			        pour chaque longueur d'onde, chaque espèce et chaque exponentielle de la décomposition, 
C			        la valeur du coefficient ai de pondération de l'exponentielle.
C
C     KDIS_KI(double) : (S) Tableau des coefficients ki des exponentilles de la décomposition, contenant
C                           (CTE_CKD_NT_MAX * CTE_CKD_NP_MAX * CTE_CKD_NAI_MAX * CTE_NBABS * CTE_CKD_NWVL_MAX)
C                           éléments couvrant les longueurs d'onde, espèces, ordres i de la décomposition en
C                           exponentielles, pressions et températures
C
C     KDIS_KI_H2O(double) : (S) Tableau des coefficients ki des exponentilles de la décomposition de la vapeur d'eau, 
C                               contenant (CTE_CKD_NT_MAX * CTE_CKD_NP_MAX * CTE_CKD_NC_MAX * CTE_CKD_NAI_MAX *  CTE_CKD_NWVL_MAX)
C                               éléments couvrant les longueurs d'onde, ordres i de la décomposition en exponentielles,
C                               concentrations en H2O, pressions et températures
C
C     TAB_PRES (double): (S) Tableau des valeurs de pression associées aux valeurs des coefficients ki(P,T) 
C                            de CTE_CKD_NP_MAX éléments.
C     NB_PRES (integer): (S) Nombre de valeurs de pression effectives 
C
C     TAB_TEMP (double): (S) Tableau des valeurs de température associées aux valeurs des coefficients ki(P,T) 
C                            de CTE_CKD_NT_MAX éléments.
C     NB_TEMP (integer): (S) Nombre de valeurs de température effectives 
C
C     TAB_CONC_H2O (double): (S) Tableau des valeurs de concentration en H2O associées aux valeurs des coefficients ki(P,T,C) 
C                                de la vapeur d'eau, à CTE_CKD_NC_MAX éléments. (en ppmv)
C     NB_CONC_H2O (integer): (S) Nombre de valeurs de concentrations effectives 
C
C     IER (I4) : (S) code d'erreur =0 si pas d'erreur, =-1 sinon
C
C
C
C Résultats fournis par le programme
C ----------------------------------
C	- Profil de gaz absorbants
C	- Coefficients CKD
C
C Variable d'environnement
C ------------------------
C   SOS_ABS_ROOT : qui définit l'arborescence d'accès au répertoire principal du code SOS_ABS
C   Cette variable d'environnement permet de localiser les fichiers sous $SOS_ABS_ROOT/fic 
C
C Common utilisé:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Ecriture dans un fichier log
C     -Lecture du fichier profil utilisateur
C
C  Sous programmes utilisés:
C --------------------------
C    Fonctions incluses dans SOS_SUB_TRS.F : DATATM, READ_CKD_COEFF
C
C==============================================================================

      SUBROUTINE SOS_PREPA_ABSPROFILE(
     &			      WA, NUSTEP, PSURF, H2O, O3, CO2, CH4,
     &			      ABSPROFIL, FICABSPROFIL, TRACE, IDLOG,
     &		 	      NU, LAMB1, IABS, USERPROFIL, ALTABS,
     &			      RO, NEXP, KDIS_AI, KDIS_KI, KDIS_KI_H2O,
     &                     TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                     TAB_CONC_H2O,NB_CONC_H2O, IER)
      			 
      IMPLICIT NONE

C* Définition des paramètres      
      INTEGER*4  NLEVEL  ! nb de niveaux en altitude du profil      
      PARAMETER(NLEVEL=CTE_ABS_NBLEV)

      INTEGER    MAXKD   ! Nombre max d'exponientielles par gaz absorbant (méthode CKD)
      PARAMETER (MAXKD=CTE_CKD_NAI_MAX)
      
      INTEGER    MAXP    ! Nombre max de niveaux de pression atmosphérique  (méthode CKD)
      PARAMETER (MAXP=CTE_CKD_NP_MAX)

      INTEGER    MAXT    ! Nombre max de niveaux de température atmosphérique  (méthode CKD)
      PARAMETER (MAXT=CTE_CKD_NT_MAX)




      
C* Definition des variables                                                 
C--------------------------------------------

C*   Simu
      DOUBLE PRECISION WA       ! (E) Longueur d'onde (microns).	
      DOUBLE PRECISION NU       ! Nombre d'onde (cm-1)      
      DOUBLE PRECISION NUSTEP   ! (E) Pas du nombre d'onde (cm-1)
      DOUBLE PRECISION NUMIN, NUMAX ! Nombre d'onde (cm-1) min et max sur les données spectrales

C*   Profils
      INTEGER*2 ABSPROFIL    ! (E) Type du profil en absorption
      CHARACTER*CTE_LENFIC2 FICABSPROFIL ! (E) Fichier profil en absorption fourni par utilisateur
      DOUBLE PRECISION ALTC(1:CTE_ABS_NBLEV)	! Altitudes du profil d'absorption (en ordre croissant)      
      DOUBLE PRECISION ALTABS(1:CTE_ABS_NBLEV)   ! Altitudes du profil d'absorption (ordre decroissant)
      DOUBLE PRECISION USERPROFIL(CTE_ABS_NBLEV,CTE_ABS_NBCOL)  ! profil en absorption fourni par utilisateur
      DOUBLE PRECISION RO(CTE_NBABS,NLEVEL),DENS(NLEVEL)	 ! Nb de particules/cm2 par espèce et couche
      DOUBLE PRECISION DP(NLEVEL),P(NLEVEL),T(NLEVEL)

C*   Coefficients CKD
      INTEGER*4  NEXP(CTE_NBABS, CTE_CKD_NWVL_MAX) 
                        ! (S) Nombre d'exponentielles par longueur d'onde et espèce 

      DOUBLE PRECISION  KDIS_AI(CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX) 
                        ! (S) Valeur des coefficients ai
                        !     par longueur d'onde, espèce 
                        !     et ordre i de la décomposition en exponentielles
                                          
      DOUBLE PRECISION  KDIS_KI(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                          CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX )
                        ! (S) Valeur des coefficients ki
                        !     par longueur d'onde, espèce, 
                        !     ordre i de la décomposition en exponentielles,
                        !     pression et température

      DOUBLE PRECISION  KDIS_KI_H2O(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                              CTE_CKD_NC_MAX, CTE_CKD_NAI_MAX, 
     &                              CTE_CKD_NWVL_MAX )
                        ! (S) Valeur des coefficients ki
                        !     par longueur d'onde,  
                        !     ordre i de la décomposition en exponentielles,
                        !     concentration, pression et température

      INTEGER NB_TEMP, NB_PRES  ! (S) Nombre de valeurs de température et pression dans le fichier CKD
      DOUBLE PRECISION TAB_TEMP(CTE_CKD_NT_MAX) ! (S) Tableau des valeurs de température lues
      DOUBLE PRECISION TAB_PRES(CTE_CKD_NP_MAX) ! (S) Tableau des valeurs de pression lues            
 
      INTEGER NB_CONC_H2O  ! (S) Nombre de valeurs de concentrations H2O dans le fichier CKD H2O
      DOUBLE PRECISION TAB_CONC_H2O(CTE_CKD_NC_MAX) ! (S) Tableau des valeurs de concentrations H2O lues (ppmv)

          
C*   Continuum
      DOUBLE PRECISION WVCONT   ! Longueur d'onde 
      
C*   Absorption gazeuse
      DOUBLE PRECISION PSURF  ! (E) Pression atmosphérique de surface (mbar).
      DOUBLE PRECISION H2O    ! (E) Contenu intégré en vapeur d'eau (g/cm2).
      DOUBLE PRECISION O3     ! (E) Contenu intégré en ozone (Dobson = 1000 cm.atm).
      DOUBLE PRECISION CO2    ! (E) Concentration du dioxyde de carbone en surface (ppmv).
      DOUBLE PRECISION CH4    ! (E) Concentration du méthane en surface (ppmv).

      INTEGER*2 IABS(1:CTE_NBABS)   ! Tableau identifiant les bandes d'absorption simulée 
                                    ! (H2O, CO2, O3, N2O, CO, CH4, O2, NO2)

      DOUBLE PRECISION QO3      ! contenu en ozone ajusté
      DOUBLE PRECISION QH2O     ! contenu en vapeur d'eau ajusté

      DOUBLE PRECISION CO2_DEFAULT_PPMV_SURF  ! Concentration de CO2 en surface du profil par défaut (en ppmv)
      DOUBLE PRECISION CH4_DEFAULT_PPMV_SURF  ! Concentration de CH4 en surface du profil par défaut (en ppmv)

      INTEGER*2 NABS            ! indice du corps absorbant
      INTEGER*2 JABS		    ! Statut du corps absorbant (1 pour traité, 0 sinon)
      
C*    Autres variables
      INTEGER*4 I,J,II 	! Niveau.
      INTEGER*4 K

      INTEGER*4 LAMB1

      INTEGER*4 IER     	! code d'erreur =0 si pas d'erreur, =-1 sinon     

      LOGICAL TRACE        	! = vrai,si écriture dans le fichier trace

      INTEGER*4 IDLOG	        ! Numéro identifiant du fichier Trace	    

      CHARACTER*CTE_LENDIR DIRFIC	! chemin d'acces au fichier de profil SO2-NO2
      INTEGER*4 IFIN       		! nb caractères significatifs d'une chaine                
      
C*    Données atmo      
      DOUBLE PRECISION ATMOCM(CTE_NBABS) !Coef par espèce gazeuse : H2O, CO2, O3, N2O, CO, CH4, O2, NO2
                                         ! = (10*Nb Avogadro) /(Masse molaire * accélération de la pesanteur)
      DATA (ATMOCM(I),I=1,CTE_NBABS)/3.410E+22,1.395E+22,1.279E+22,
     &1.395E+22,2.192E+22,3.837E+22,1.918E+22,1.3340E+22/
     
      DOUBLE PRECISION PDSMOL(CTE_NBABS)
      DATA (PDSMOL(I),I=1,CTE_NBABS)/18.,44.,48.,44.,28.,16.,32.,46./



C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0
      

 	
C* Initialisation des tableaux pour la prise en compte des absorbants
C--------------------------------------------------------------------

C
C-------- ABSORBANTS      --------------------------------------
C LINES AND/OR CONTINUUM
C 1-H2O  2-CO2  3-O3  4-N2O  5-CO  6-CH4  7-O2  8-NO2 
      
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C Absorbants pris en compte (définis dans SOS.h)
      IABS(1)=CTE_ABS_H2O
      IABS(2)=CTE_ABS_CO2
      IABS(3)=CTE_ABS_O3
      IABS(4)=CTE_ABS_N2O
      IABS(5)=CTE_ABS_CO
      IABS(6)=CTE_ABS_CH4
      IABS(7)=CTE_ABS_O2
      IABS(8)=CTE_CONT_NO2
      
C Initialisation du profil d'absorption en altitude
      DO I=1,CTE_ABS_NBLEV
         ALTABS(I)=0.
      ENDDO

C Si Psurf, H2O, O3, CO2 ou CH4 ne sont pas définis, on prend les
C valeurs par défaut des profils
      
      IF (TRACE) THEN      

         IF (H2O.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
            WRITE(IDLOG,*,err=921) 'H2O (g/cm2) :', H2O
         ENDIF            
         
         IF (O3.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
            WRITE(IDLOG,*,err=921) 'O3 (cm.atm) :', O3
         ENDIF
       
         IF (CO2.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
            WRITE(IDLOG,*,err=921) 'CO2 at surface (ppmv) :', CO2
         ENDIF
         
         IF (CH4.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
            WRITE(IDLOG,*,err=921) 'CH4 at surface (ppmv) :', CH4
         ENDIF
         
      ENDIF         
   
   
C Lecture du profil utilisateur si il est fourni en entrée
c USERPROFIL(J,K), for each level J, enter (K=1,CTE_ABS_NBCOL):
C Z(km),P(hPa),T(K), H2O, CO2, O3, N2O, CO, CH4, O2, density (not used), NO2, SO2
c with gaz concentrations in ppmv.
C
C Ou lecture des fichiers prédéfinis.

      IF(ABSPROFIL.EQ.0) THEN !Cas fichier utilisateur
         write(6,*) ' Reading user file for gas profile:',FICABSPROFIL
         OPEN(2,FILE=FICABSPROFIL,err=1012)
         DO I=1,CTE_ABS_NBLEV
            READ(2,*,err=1021)II,(USERPROFIL(I,K),K=1,CTE_ABS_NBCOL)
         ENDDO
         CLOSE(2)

      ELSE !Cas utilisation de fichiers de profils prédéfinis

         !Nom du répertoire avec les données
         CALL GETENV('SOS_ABS_ROOT',DIRFIC)
         IFIN=INDEX(DIRFIC,' ')
         IFIN=IFIN-1
         IF (IFIN.LE.0) GOTO 925
      
         IF ((IFIN+4).GT.CTE_LENDIR) GOTO 926
      
         DIRFIC(IFIN+1:IFIN+4)='/fic'
         IFIN=IFIN+4
         IF (IFIN.LE.0) IFIN=CTE_LENDIR
      
         !Lecture du profil de SO2,NO2 et calcule du contenu NO2 en ppm (mg/kg) sur le profil
         !--> Renseigne l'information sur le contenu en NO2 dans RO(8,*)
         WRITE(6,*) '  SOS_PREPA_ABSPROFILE: read profile SO2,NO2'
         OPEN(1,FILE=DIRFIC(1:IFIN)//'/SO2-NO2')
         DO J=1,NLEVEL
            READ(1,*,ERR=927)USERPROFIL(J,13),USERPROFIL(J,12)	!Fichier donne SO2 en 1ère colonne
								       !NO2 stocké en indice 12, SO2 en indice 13
         ENDDO
         CLOSE(1)

      ENDIF

C     Récupération des quantités de gaz par niveau en ppm (mg/kg) sur le profil (tableau RO)
C     pour les gaz 1 à 7 : 1-H2O  2-CO2  3-O3  4-N2O  5-CO  6-CH4  7-O2
C     --> Usage du profil en ppmv*10^6 du tableau USERPROFIL, 
C         renseigné en amont par l'option d'utilisation d'un fichier utilisateur (ABSPROFIL.EQ.0)
C         ou renseigné directement par la routine DATATM si ABSPROFIL > 0.
C     --> A ce stade, le contenu du tableau RO donne le profil par niveau (pas par couche atmosphérique).
      CALL DATATM(RO,P,T,ALTC,ABSPROFIL,DENS,USERPROFIL,NLEVEL,PSURF)

C     Conversion de la quantité de NO2, non gérée par la routine DATATM,
C     en ppm (mg/kg) sur le profil ==> Complète le tableau RO
      DO J=1,NLEVEL           
         RO(8,J)=USERPROFIL(J,12)*1.0E-06*46/28.9
      ENDDO

      IF(CO2.GE.0.) CO2_DEFAULT_PPMV_SURF = RO(2,1)*28.97/(44.0E-06)  ! Reconversion en ppmv
      IF(CH4.GE.0.) CH4_DEFAULT_PPMV_SURF = RO(6,1)*28.97/(16.0E-06)


      
c     Inversion des altitudes pour SOS => de la plus haute à la plus basse
      DO J=1,NLEVEL
	     ALTABS(J)=ALTC(NLEVEL-J+1)
      ENDDO
C
C     Conversion du profil de contenus de ppm par niveau en particules/cm2 par couche
C                                             ----------                   ----------
      DO J=1,NLEVEL-1
         DP(J)=P(J)-P(J+1)
         DO K=1,CTE_NBABS
            RO(K,J)=DP(J)*(RO(K,J)+RO(K,J+1))/2.0*ATMOCM(K)	!Expression en particules/cm2 par couche
c            RO(K,J)=RO(K,J)/6.022E+23*PDSMOL(K)			!Expression en g/cm2
         ENDDO
      ENDDO
      
C     Correction H2O
      IF(H2O.GE.0.)THEN
         WRITE(6,*) 
     &   "  Correction of the profile with the user's H2O amount"
         QH2O=0.
         DO J=1,NLEVEL
            QH2O=QH2O+RO(1,J)	! Contenu intégré du profil en particules/cm2
         ENDDO

         QH2O=QH2O/6.022E+23*PDSMOL(1)	!Passage en g/cm2
         DO J=1,NLEVEL
            RO(1,J)=RO(1,J)*H2O/QH2O	! Reste en particules/cm2
            USERPROFIL(J,4)=USERPROFIL(J,4)*H2O/QH2O
         ENDDO
      ENDIF

C
C     Correction O3
      IF(O3.GE.0.)THEN
         WRITE(6,*) 
     &   "  Correction of the profile with the user's O3 amount"
         O3=O3/1000.		! Passage Dobson ==> cm.atm
         QO3=0.         
         DO J=1,NLEVEL
            QO3=QO3+RO(3,J)	! Contenu intégré du profil en particules/cm2
         ENDDO

         QO3=QO3/6.022E+23*PDSMOL(3)	!Passage particules/cm2 ==> g/cm2
         QO3=QO3*466.23			!Passage g/cm2 ==> cm.atm
         DO J=1,NLEVEL
            RO(3,J)=RO(3,J)*O3/QO3				! Reste en particules/cm2
            USERPROFIL(J,6)=USERPROFIL(J,6)*O3/QO3
         ENDDO
      ENDIF
   
C     Correction CO2
      IF(CO2.GE.0.)THEN
         DO J=1,NLEVEL
            RO(2,J)=RO(2,J)*CO2/CO2_DEFAULT_PPMV_SURF		! Reste en particules/cm2
            USERPROFIL(J,5)=USERPROFIL(J,5)*CO2/CO2_DEFAULT_PPMV_SURF
         ENDDO         
      ENDIF

C     Correction CH4
      IF(CH4.GE.0.)THEN
         DO J=1,NLEVEL
            RO(6,J)=RO(6,J)*CH4/CH4_DEFAULT_PPMV_SURF		! Reste en particules/cm2
            USERPROFIL(J,9)=USERPROFIL(J,9)*CH4/CH4_DEFAULT_PPMV_SURF
         ENDDO         
      ENDIF

C*******************************************************
C*******************************************************
C*
C* Lecture des coefficients de la CKD
C* Estimation de l'indice LAMB1 des tables
C*
C*******************************************************
C*******************************************************
C     Détermination du nombre d'onde pour le calcul en absorption
C     (nu = 1 / Wa en cm-1)
      NU=1.0E+4/WA   
      IF ((NU.GT.CTE_CKD_NUMAX).OR.(NU.LT.CTE_CKD_NUMIN)) GOTO 905

c     Lecture des coefficients d'absorption hors H2O
      !WRITE(6,*) '  Read Absorption Coefficients'
      !WRITE(6,*) '  Wave number (cm-1) : ', NU
      DO NABS=1,CTE_NBABS               !Couvre les 8 espèces gazeuses
         JABS=IABS(NABS)
         !print*,"    ===> NABS : ", NABS
         CALL READ_CKD_COEFF(NABS, JABS, NU, NUSTEP,
     &                       NEXP, KDIS_AI, KDIS_KI, KDIS_KI_H2O,
     &                       NUMAX, NUMIN,
     &                       TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &				 TAB_CONC_H2O, NB_CONC_H2O, IER)

         IF (IER.NE.0) GOTO 900   ! cas d'erreur
      ENDDO
       
C     Indice pour les tableaux CKD correspondant au nombre d'onde NU traité
      LAMB1 = 1+INT((NUMAX-NU)/NUSTEP)  


C* Ecriture du fichier trace
C-------------------------------------------
      IF (TRACE) THEN      

         WRITE(IDLOG,*,err=921) "Wavelength (nm): ",WA
	  WRITE(IDLOG,*,err=921) "Wave number (cm-1): ",NU
	  WRITE(IDLOG,*,err=921) 
     &      'Spectral resolution of CKD coefficients (cm-1): ',NUSTEP
	  WRITE(IDLOG,*,err=921) 'Index CKD : ',LAMB1
	  WRITE(IDLOG,*,err=921) 'Type of absorption profile: ',
     &   ABSPROFIL
         WRITE(IDLOG,*,err=921) ' '
	 
	  IF (ABSPROFIL.NE.7) THEN
	 
	      WRITE(IDLOG,*,err=921) 
     &       'Calculation of the absorption profile'
	      WRITE(IDLOG,*,err=921) 
     &       '-------------------------------------'
	      WRITE(IDLOG,*,err=921) ' '
	      WRITE(IDLOG,*,err=921) ' Species considered'
             WRITE(IDLOG,16,err=921)'H2O:  ',IABS(1)
	      WRITE(IDLOG,16,err=921)'CO2:  ',IABS(2)
	      WRITE(IDLOG,16,err=921)'O3:   ',IABS(3)
	      WRITE(IDLOG,16,err=921)'N2O:  ',IABS(4)
	      WRITE(IDLOG,16,err=921)'CO:   ',IABS(5)
	      WRITE(IDLOG,16,err=921)'CH4:  ',IABS(6)
	      WRITE(IDLOG,16,err=921)'O2:   ',IABS(7)
	      WRITE(IDLOG,16,err=921)'NO2:  ',IABS(8) 
	      WRITE(IDLOG,*,err=921) ' '

	      IF (PSURF.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921) 'User Psurf (mbar)', PSURF
	        WRITE(IDLOG,*,err=921) ' '
             ENDIF

             IF (H2O.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921)  
     $         "Correction of the profile with the user's H2O amount:"
               WRITE(IDLOG,*,err=921) 'User H2O content (g/cm2)', H2O
             ENDIF

             !Calcul de la quantité totale de vapeur d'eau sur la colonne atmosphérique
             QH2O=0.		!Recalcul du nouveau contenu intégré
             DO J=1,NLEVEL
                QH2O=QH2O+RO(1,J)	! Contenu intégré du profil en particules/cm2
             ENDDO
             QH2O=QH2O/6.022E+23*PDSMOL(1)	!Passage en g/cm2
             WRITE(IDLOG,17,err=921)   
     $       '==> Used column amount of Water vapor (g/cm2):',QH2O 
	      WRITE(IDLOG,*,err=921) ' '

             IF (O3.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921)  
     $         "Correction of the profile with the user's O3 amount:"
               WRITE(IDLOG,*,err=921)   
     $         'User O3 content (Dobsons):', O3*1000.
               WRITE(IDLOG,*,err=921)  
     $         'Default profile Ozone content (Dobsons):',QO3*1000.
             ENDIF

             !Calcul de la quantité totale d'ozone sur la colonne atmosphérique
             QO3=0.
             DO J=1,NLEVEL
                QO3=QO3+RO(3,J)	
             ENDDO			! 
             QO3=QO3/6.022E+23*PDSMOL(3)	!Passage particules/cm2 ==> g/cm2
             QO3=QO3*466.23			!Passage g/cm2==> cm.atm
             WRITE(IDLOG,17,err=921)   
     $       '==> Used column amount of Ozone (Dobsons):    ',QO3*1000. 
	      WRITE(IDLOG,*,err=921) ' '

             IF (CO2.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921)  
     $         "Correction of the profile with the user's CO2 ratio:"
               WRITE(IDLOG,*,err=921)   
     $         'User CO2 surface concentration (ppmv):', CO2
               WRITE(IDLOG,*,err=921)  
     $         'Default profile CO2 surface concentration (ppmv):',
     $         CO2_DEFAULT_PPMV_SURF
	        WRITE(IDLOG,*,err=921) ' '
             ENDIF

             IF (CH4.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921)  
     $         "Correction of the profile with the user's CH4 ratio:"
               WRITE(IDLOG,*,err=921)   
     $         'User CH4 surface concentration (ppmv):', CH4
               WRITE(IDLOG,*,err=921)  
     $         'Default profile CH4 surface concentration (ppmv):',
     $         CH4_DEFAULT_PPMV_SURF
	        WRITE(IDLOG,*,err=921) ' '
             ENDIF

	      WRITE(IDLOG,*,err=921) ' '
	      WRITE(IDLOG,*,err=921) ' Gas profile'
	      WRITE(IDLOG,*,err=921)'LEVEL ALT(km) P(hPa) T(K), + 
     $       GASES (ppmv), except column 12: air density'
             DO j=1,CTE_ABS_NBLEV
                WRITE(IDLOG,19,err=921)j,(USERPROFIL(j,k),k=1,13)
             ENDDO

	      WRITE(IDLOG,*,err=921) ' '

	  ENDIF ! fin test absprofil

      ENDIF ! fin trace
      

C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
 925  WRITE(6,*) '  SOS_PREPA_ABSPROFILE : ERROR_925: '
      WRITE(6,*) '  => Error while getting SOS_ABS_ROOT variable'
      GOTO 9998

 926  WRITE(6,*) '  SOS_PREPA_ABSPROFILE : ERROR_926: '
      WRITE(6,*) '  => Length of $SOS_ABS_ROOT/fic is too long' 
      WRITE(6,*) '  => CTE_LENDIR has to be adjusted.' 
      GOTO 9998

 927  WRITE(6,*) '  SOS_PREPA_ABSPROFILE : ERROR_927: '
      WRITE(6,*) '  => Error while reading the file fic/SO2-NO2'
      GOTO 9998

  900 WRITE(6,*) '  ERROR in subroutine READ_CKD_COEFF'
      GOTO 9998

  905 WRITE(6,*) '  The simulation wavelength is not included ' //
     &           'in the spectral range of CKD data'
      WRITE(6,*) '  --> NU (cm-1) = ', NU
      WRITE(6,*) '  --> Covered spectral range (cm-1) : ', 
     &           CTE_CKD_NUMIN, CTE_CKD_NUMAX
      GOTO 9998
            
  912 WRITE(6,*) '  ERROR while opening logfile by SOS_PREPAPROFIL'
      GOTO 9998
      
  921 WRITE(6,*) '  ERROR while writing logfile by SOS_PREPAPROFIL'
      GOTO 9998
      
 1012 WRITE(6,*)'  ERROR while opening user profile by SOS_PREPAPROFIL'
      GOTO 9998
      
 1021 WRITE(6,*)'  ERROR while reading user profile by SOS_PREPAPROFIL'
      GOTO 9998

9998  IER=-1
      GOTO 9999

C* Format
C---------   
   16 FORMAT(a6,1x,i2)
   17 FORMAT(a46,1x,f9.3)   
   19 FORMAT(i2.2,1x,f5.1,1x,f10.5,1x,f5.1,1x,11(e11.5,1x))   

 9999 RETURN
	 
      END
      
      
      
  
