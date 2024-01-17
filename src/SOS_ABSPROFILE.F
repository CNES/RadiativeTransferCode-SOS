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
C* FICHIER: SOS_ABSPROFILE.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Calcul du profil atmosphérique en absorption pour une exponentielle CKD
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmosphérique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, B.Lafrance, X.Lenot
C* 
C*
C* MOD:VERSION:1.0: DATE: 06/06/2013 
C*    Création du fichier pour adapter les OS à la méthode CKD (absorption atmosphérique)
C*
C* MOD:VERSION:1.1: DATE: 24/02/2016
C*    - Tableau NEXP passe de Double Précision à Integer
C*    - Modification du type de la variable K (numéro d'indice d'exponentiel) : 
C*      passe de INTEGER*4 à INTEGER*2
C*      pour être compatible avec le type dans la routine coefabs.
C*
C* MOD:VERSION:2.0: DATE: 08/06/2020
C*    - Le tableau XKI est remplacé par KDIS_KI avec la variable Température en plus,
C*      ainsi que les coefficients CKD pour le continuum NO2
C*    - Ajustement des arguments d'entrée (IK8 pour le NO2).
C*    - Suppression des calculs spécifiques aux continuums (car intégrés dans calculs CKD).
C*    - Suppression de la sortie de la transmission totale sur la colonne atmosphérique.
C*
C* MOD:VERSION:2.1: 31/08/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Conversion des traces en Anglais.
C*
C* MOD:VERSION:2.2: 16/06/2023
C*    - Ajout des paramètres KDIS_KI_H2O,TAB_CONC_H2O et NB_CONC_H2O en arguments de la routine.
C*    - Introduction de la prise en compte de la concentration en vapeur d'eau.
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Définition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisées :
C    CTE_LENFIC2 :   Longueur des noms de fichiers avec arborescence.
C    CTE_ABS_NBLEV : Nombre de niveaux du profil en absorption.
c    CTE_NBABS :     Nombre d'absorbants.     
C    CTE_ABS_NBCOL : Nombre de colonnes dans un profil d'absorption fourni
C		       en entrée.  
C----------------------------------------------------------------------------
#include "SOS.h"



C==============================================================================
C PROCEDURE: SOS_ABSPROFILE.F
C ==========
C      Cette procédure calcule le profil en épaisseur optique d'absorption
C      en fonction de l'altitude (sur 50 niveaux). 
C
C      Le profil est calculé pour un nombre d'onde donné à partir de la méthode CKD
C
C      Cette fonction est une adaptation de la fonction TRSKD fournie par le LOA.
C      Cette dernière permettait de calculer la transmission sur un intervalle
C      spectral donné.
C
C Description des paramètres
C --------------------------
C
C     ABSPROFIL (I2) : (E) Type de profil atmosphérique pour l'absorption
C			0:User, 1:TROPICAL, 2:MLS, 3:MLW, 4:SAS, 5:SAW,
C                    6:USST62, 7:Pas d'absorption 	 	 
C
C     NU  (double) : (E)  Nombre d'onde en cm-1
C
C     LAMB1  (I4) : (E)   Indice de l'intervalle CKD
C
C     IABS  (I2) : (E) Tableau de CTE_NBABS éléments. Chaque élément correspond à un 
C			  gaz absorbant (Absorbants: H2O CO2 O3 N2O CO CH4 O2 NO2).
C	 		  Pour chaque absorbant, vaut 1 si traité avec méthode CKD, 0 sinon 
C
C     USERPROFIL (double) (E)  Tableau de (CTE_ABS_NBLEV * CTE_ABS_NBCOL) éléments
C			    contenant pour chaque gaz absorbant le profil exprimé en ppmv*1E+6
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
C     ALTABS (double) : (E)  Tableau de (CTE_ABS_NBLEV) contenant les niveaux d'altitude du profil
C
C     RO (double) (E) : Tableau de (CTE_NBABS * CTE_ABS_NBLEV) éléments
C		          contenant, pour chaque gaz absorbant de la méthode CKD, le profil 
C			   exprimé en particules/cm2 pour chaque couche.
C			   Il contient dans l'ordre : H2O, CO2, 03, N2O, CO, CH4, O2, NO2
C                       --> Gaz qui ont une absorption par raies fines (pas uniquement par continnums)
C
C     NEXP (I4) :   (E)  Tableau de  (CTE_NBABS * CTE_CKD_NWVL_MAX) éléments contenant pour chaque absorbant
C		           et chaque intervalle spectral le nombre d'exponentielles de la méthode CKD
C
C     KDIS_KI(double) : (E) Tableau des coefficients ki des exponentilles de la décomposition, contenant
C                           (CTE_CKD_NT_MAX * CTE_CKD_NP_MAX * CTE_CKD_NAI_MAX * CTE_NBABS * CTE_CKD_NWVL_MAX)
C                           éléments couvrant les longueurs d'onde, espèces, ordres i de la décomposition en
C                           exponentielles, pressions et températures
C
C     KDIS_KI_H2O(double) : (E) Tableau des coefficients ki des exponentilles de la décomposition de la vapeur d'eau, 
C                               contenant (CTE_CKD_NT_MAX * CTE_CKD_NP_MAX * CTE_CKD_NC_MAX * CTE_CKD_NAI_MAX *  CTE_CKD_NWVL_MAX)
C                               éléments couvrant les longueurs d'onde, ordres i de la décomposition en exponentielles,
C                               concentrations en H2O, pressions et températures
C
C     IK1  (I4) : (E) Indice de l'exponentielle pour l'absorption de H2O 
C
C     IK2  (I4) : (E) Indice de l'exponentielle pour l'absorption de CO2 
C
C     IK3  (I4) : (E) Indice de l'exponentielle pour l'absorption de O3 
C
C     IK4  (I4) : (E) Indice de l'exponentielle pour l'absorption de N2O 
C
C     IK5  (I4) : (E) Indice de l'exponentielle pour l'absorption de CO 
C
C     IK6  (I4) : (E) Indice de l'exponentielle pour l'absorption de CH4 
C
C     IK7  (I4) : (E) Indice de l'exponentielle pour l'absorption de O2 
C
C     IK8  (I4) : (E) Indice de l'exponentielle pour l'absorption de NO2 
C
C     TAB_PRES (double): (E) Tableau des valeurs de pression associés aux valeurs des coefficients ki(P,T) 
C                            de CTE_CKD_NP_MAX éléments.
C     NB_PRES (integer): (E) Nombre de valeurs de pression effectives 
C
C     TAB_TEMP (double): (E) Tableau des valeurs de température associés aux valeurs des coefficients ki(P,T) 
C                            de CTE_CKD_NT_MAX éléments.
C     NB_TEMP (integer): (E) Nombre de valeurs de température effectives 
C
C     TAB_CONC_H2O (double): (E) Tableau des valeurs de concentration en H2O associées aux valeurs des coefficients ki(P,T,C) 
C                                de la vapeur d'eau, à CTE_CKD_NC_MAX éléments (ppmv)
C     NB_CONC_H2O (integer): (E) Nombre de valeurs de concentrations effectives.
C
C     TAUABSTOT (double) : (S) Tableau de (CTE_ABS_NBLEV-1) contenant les épaisseurs optiques
C			          d'absorption (continumm + CKD) pour chaque couche atmosphérique
C
C     TRACE (LOGICAL) : (E) Si vrai => écriture des fichiers trace
C
C     IDLOG (I4) : (E) Identifiant du fichier log
C
C     IER (I4) : (S) code d'erreur =0 si pas d'erreur, =-1 sinon
C  
C
C Common utilisé:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Ecriture du Log
C
C==============================================================================

      SUBROUTINE SOS_ABSPROFILE(ABSPROFIL, NU, LAMB1, IABS, 
     &		         USERPROFIL, ALTABS, RO, NEXP,
     &		         KDIS_KI, KDIS_KI_H2O,
     &		         IK1, IK2, IK3, IK4, IK5, IK6, IK7, IK8,
     &                 TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                 TAB_CONC_H2O, NB_CONC_H2O,
     &			  TAUABSTOT, TRACE, IDLOG, IER) 
    
      IMPLICIT NONE

C* Définition des paramètres
     
      INTEGER*4  NLEVEL    ! nb de niveaux en altitude du profil      
      PARAMETER(NLEVEL=CTE_ABS_NBLEV)
      
      INTEGER    MAXABS   ! Nombre de gaz absorbant (méthode CKD)
      PARAMETER (MAXABS=CTE_NBABS)

      INTEGER    MAXKD   ! Nombre max d'exponientielles par gaz absorbant (méthode CKD)
      PARAMETER (MAXKD=CTE_CKD_NAI_MAX)

      INTEGER    MAXT    ! Nombre max de niveaux de température atmosphérique  (méthode CKD)
      PARAMETER (MAXT=CTE_CKD_NT_MAX)
            
      INTEGER    MAXP    ! Nombre max de niveaux de pression atmosphérique  (méthode CKD)
      PARAMETER (MAXP=CTE_CKD_NP_MAX)
      
      INTEGER    MAXNU    ! nombre max de nombres d'ondes pour le calcul de l'absorption
      PARAMETER (MAXNU=CTE_CKD_NWVL_MAX)
  
      integer ij !pour test     
C* Definition des variables                                                 
C--------------------------------------------

C*   Simu
      DOUBLE PRECISION NU    ! Nombre d'onde (cm-1)

C*   Profils
      INTEGER*2 ABSPROFIL    ! (E) Type du profil en absorption
      CHARACTER*CTE_LENFIC2 FICABSPROFIL ! (E) Fichier profil en absorption fourni par utilisateur
      
      DOUBLE PRECISION ALTABS(1:NLEVEL)         ! Altitudes du profil d'absorption (ordre decroissant)
      DOUBLE PRECISION TAUABSTOT(1:NLEVEL)      ! Profil d'épaisseur optique
      DOUBLE PRECISION TRANS(1:NLEVEL)          ! Profil de transmission
      
      DOUBLE PRECISION USERPROFIL(NLEVEL,CTE_ABS_NBCOL) ! profil en absorption fourni en entrée (ppmv)
      DOUBLE PRECISION RO(CTE_NBABS,NLEVEL)	! Nombre de particules/cm2 par espèce gazeuse et par couche du profil atmosphérique
      DOUBLE PRECISION TAUABS1C(NLEVEL-1)        ! Epaisseur optique de chaque couche


C*   Coefficients CKD
      INTEGER*4  NEXP(CTE_NBABS,CTE_CKD_NWVL_MAX) 
                              ! (E) Nombre d'exponentielles par longueur d'onde et espèce 
                              
      
      DOUBLE PRECISION  KDIS_KI(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                          CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX )
                        ! (E) Valeur des coefficients ki
                        !     par longueur d'onde, espèce, 
                        !     ordre i de la décomposition en exponentielles,
                        !     pression et température
                        

      DOUBLE PRECISION  KDIS_KI_H2O(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                              CTE_CKD_NC_MAX, CTE_CKD_NAI_MAX, 
     &                              CTE_CKD_NWVL_MAX )
                        ! (E) Valeur des coefficients ki
                        !     par longueur d'onde,  
                        !     ordre i de la décomposition en exponentielles,
                        !     concentration, pression et température

      INTEGER NB_TEMP, NB_PRES  ! (E) Nombre de valeurs de température et pression dans le fichier CKD
      DOUBLE PRECISION TAB_TEMP(CTE_CKD_NT_MAX) ! (E) Tableau des valeurs de température lues
      DOUBLE PRECISION TAB_PRES(CTE_CKD_NP_MAX) ! (E) Tableau des valeurs de pression lues            

      INTEGER NB_CONC_H2O  ! (E) Nombre de valeurs de concentrations H2O dans le fichier CKD H2O
      DOUBLE PRECISION TAB_CONC_H2O(CTE_CKD_NC_MAX) ! (E) Tableau des valeurs de concentrations H2O lues
           
                                         
      DOUBLE PRECISION  PRS     ! Pression moyenne d'une couche du profil (mbar)
      DOUBLE PRECISION  TMP     ! Température moyenne d'une couche du profil (K)
      DOUBLE PRECISION  CONC_H2O! Concentration H2O moyenne d'une couche du profil (en ppmv)
      DOUBLE PRECISION  XK

c*   Transmissions
      DOUBLE PRECISION TRS
      
C*   Absorption gazeuse
      INTEGER*2 IABS(1:CTE_NBABS)   ! Tableau identifiant les bandes d'absorption simulée 
                                    ! (H2O, CO2, O3, N2O, CO, CH4, O2, NO2)
                                    ! Inclus raies et/ou continuums d'absorption.

      INTEGER*2 NABS        ! Indice du corps absorbant
      INTEGER*2 JABS		! Statut du corps absorbant (1 pour traite, 0 sinon)
       
C*    Autres variables
      INTEGER*4 I,J  					! Niveau.
      INTEGER*4 IK,IK1,IK2,IK3,IK4,IK5,IK6,IK7,IK8  	! Indice d'exponentielle      
      INTEGER*2 K	

      INTEGER*4 LAMB1

      INTEGER*4 IER     	! code d'erreur =0 si pas d'erreur, =-1 sinon     

      LOGICAL TRACE        	! = vrai,si écriture dans le fichier trace

      INTEGER*4 IDLOG	! Numéro identifiant du fichier Trace	    
      



C* Initialisation 
C------------------------------------------------------

C* Code d'erreur
      IER=0

C* Initialisation des tableaux
      DO J=1,NLEVEL
         TAUABSTOT(J)=0.D+00
      ENDDO
        
C* Pas de traitement si ABSPROFIL=7
C----------------------------------

      IF (ABSPROFIL.EQ.7) RETURN
	
C* Calcul des épaisseur optiques d'absorptions pour un profil en altitude
C* Méthode CKD 
C----------------------------------------------------------------------

     
c La transmission au niveau supérieur vaut 1     
      TRANS(1)=1.D+00	
      TAUABSTOT(1)=0.D+00
            	   
c Calcul de la transmission + épaisseur optique pour chaque couche du profil
c Il y a N niveaux donc N-1 couches
       !Init de la transmission
       TRS=1.0
       DO J=1,NLEVEL-1	!Boucle sur les couches du TOA au sol
          TAUABS1C(J)=0.		! Ep opt d'absorption de la couche J (profil couche haute=1, couche basse=NLEVEL-1)
          PRS=(USERPROFIL(NLEVEL-J,2)+USERPROFIL(NLEVEL-J+1,2))/2.  !USERPROFIL : couche haute=NLEVEL-1, couche basse=NLEVEL-1
	   TMP=(USERPROFIL(NLEVEL-J,3)+USERPROFIL(NLEVEL-J+1,3))/2.
          CONC_H2O=(USERPROFIL(NLEVEL-J,4)+USERPROFIL(NLEVEL-J+1,4))/2.
          CONC_H2O=CONC_H2O*1.D-06 !Conversion en ppmv

          DO K=1,CTE_NBABS
             IF(K.EQ.1)IK=IK1
             IF(K.EQ.2)IK=IK2
             IF(K.EQ.3)IK=IK3
             IF(K.EQ.4)IK=IK4
             IF(K.EQ.5)IK=IK5
             IF(K.EQ.6)IK=IK6
             IF(K.EQ.7)IK=IK7
             IF(K.EQ.8)IK=IK8
             XK=0.
             
             IF (NEXP(K,LAMB1).GE.1) THEN       ! Si nb exponentiel = 1 : Cas continuum seul (pas de raie fine)
                
     		  CALL COEFF_ABS_CKD(K,LAMB1,IK,PRS,TMP,CONC_H2O,
     &                             TAB_PRES,NB_PRES,TAB_TEMP,NB_TEMP,
     &                             TAB_CONC_H2O, NB_CONC_H2O,
     &                             KDIS_KI,KDIS_KI_H2O,XK,IER)
		
                IF(IER.NE.0) GOTO 910
             ENDIF

             !Epaisseur optique d'absorption totale (raies + continuums)
             !pour la couche atmosphérique J et pour le composant K
             TAUABS1C(J)=TAUABS1C(J)+XK*RO(K,NLEVEL-J)

          ENDDO !Fin de boucle sur l'espèce gazeuse pour le niveau J
  
c calcul de la transmission entre le sommet et le niveau nLevel - J	 
          TRS=TRS*exp(-TAUABS1C(J))
	   TRANS(J+1)= TRS		!Transmission bas de couche J = niveau J+1


c calcul du profil d'épaisseur optique (bas de couche)
          IF (TRS.GT.0.0D+00) THEN
	      TAUABSTOT(J+1)=-DLOG(TRS)
          ELSE
             !Epaisseur optique extrême non définie: on la pose valant CTE_TAUABS_MAX
	      TAUABSTOT(J+1)=CTE_TAUABS_MAX
          ENDIF
       ENDDO !Fin boucle sur couches

       

      
C* Ecriture du fichier trace
C-------------------------------------------
      IF (TRACE) THEN      

	     WRITE(IDLOG,*,ERR=921) ' '	 
	     WRITE(IDLOG,*,ERR=921) ' ***** Absorption profile ****'
	     WRITE(IDLOG,*,ERR=921) ' IK1 (H2O): ',IK1 
	     WRITE(IDLOG,*,ERR=921) ' IK2 (CO2): ',IK2
	     WRITE(IDLOG,*,ERR=921) ' IK3  (O3): ',IK3
	     WRITE(IDLOG,*,ERR=921) ' IK4 (N2O): ',IK4
	     WRITE(IDLOG,*,ERR=921) ' IK5  (CO): ',IK5
	     WRITE(IDLOG,*,ERR=921) ' IK6 (CH4): ',IK6
	     WRITE(IDLOG,*,ERR=921) ' IK7  (O2): ',IK7
	     WRITE(IDLOG,*,ERR=921) ' IK8 (NO2): ',IK8
	     WRITE(IDLOG,*)'LEVEL ALT(km)  OT   TRANS'
	     DO j=1,NLEVEL
               WRITE(IDLOG,18)j,ALTABS(j),TAUABSTOT(j),TRANS(j)
            ENDDO
	     WRITE(IDLOG,*,ERR=921) ' '

      ENDIF
 


C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------

  910 WRITE(6,*) '  SOS_ABSPROFILE : ERROR_910 : '
      WRITE(6,*) '  --> Error in subroutine COEFF_ABS_CKD'
      GOTO 9998
      
  921 WRITE(6,*) '  ERROR writing logfile for SOS_PREPAPROFIL'
      GOTO 9998

 9998 IER=-1
      GOTO 9999
      
C* Format
C---------   
  
!!   18 FORMAT(i2.2,1x,f7.2,2f9.5)   
   18 FORMAT(i2.2,1x,f7.2,2e15.5)   
   
 9999 RETURN
	 
      END
