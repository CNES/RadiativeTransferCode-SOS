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
C* FICHIER: SOS_PROFIL.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Definition du profil atmospherique.
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et developpement : CS SI, B.Lafrance, C. Casteras,
C*					  C.Ruffel.
C* DATE: 30/04/02
C* 
C* MOD:VERSION:1.0: Mise a niveau du code des OS du LOA
C*             2.0: Correction d'une application de couche moleculaire sous-jacente 
C*                  pour un profil a 2 couches (iprofil=2 et Zmin=0) : NBSC_C3.
C*             2.1: 07/07/2008 
C*                  Modification du format de lecture des parametres HR et HA
C*                  pour correction d'une erreur de décodage si valeur entière.
C*			HR et HA passent de F9.5 en F9
C*
C* MOD:VERSION:3.0: 19/11/2009 
C*       - Modification du passage des arguments : 
C*           - passage par couples (mot clef, valeur). 
C*           - adaptation de la gestion des cas d'erreurs. 
C*       - Conversion des messages d'erreurs en anglais. 
C*
C* MOD:VERSION:3.1: 08/12/2015
C*    - Suppression des formats à la lecture des paramètres d'entrée en faveur 
C*      de READ(PROG_ARG,*,ERR=NNN)
C*      ==> Permet une compilation sous gfortran et évite les mauvaises conversions 
C*          si la valeur d'entrée ne respecte pas strictement le format attendu 
C*          (par exemple, saisie d'un entier N au lieu d'un réel N. 
C*           aboutissant à une valeur interprétée fausse). 
C*
C*    - Ajout d'un chiffre significatif de plus aux valeurs du profil pour éviter
C*      d'avoir des épaisseurs optiques identiques sur deux niveaux consécutifs
C*      dans le cas d'une atmosphère très ténue et d'un grand nombre de couches
C*      simulées.
C*      Etiquette 20 : Passage du format FORMAT(2X,I4,3(F9.5)) 
C*                     au format FORMAT(2X,I4,3(F10.6))
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_LENKEYWORD : Taille maximale des Keywords pour le passage des arguments.
C    SOS_OS_NT : Nombre de couches du profil atmospherique. 
C    SOS_DZTRANSI : Epaisseur d'une couche de transition (en km).
C    SOS_PROFIL_MIN_NBC : Nombre minimal de sous-couches des couches 1 et 3
C                         pour un profil d'aerosols entre deux altitudes.               
C----------------------------------------------------------------------------
#include "SOS.h"



C==============================================================================
C PROGRAMME: SOS_PROFIL
C ==========
C      Ce programme definit le profil atmospherique en epaisseur optique
C      et en proportion de diffusants (molecules , aerosols) en fonction
C      du niveau. Il est egalement possible de definir le profil en trois
C      couches: molecules, molecules + aerosols (melange homogene) et 
C      molecules. La couche intermediaire est reperee par ses altitudes 
C      minimale et maximale.  
C      
C      NB : Le profil n'est pas ajuste pour une troncature de la fonction
C           de phase. Il le sera par le code des OS.
C
C
C Donnees en entree du programme
C ------------------------------
C     -->Les parametres/arguments du programme
C
C        Les arguments sont passes sous la forme d'une liste de couples (-Keyword Value)
C
C        Definition des mots clef, statut obligatoire / option, format de la valeur associee :
C         -AP.ResFile ==> Obligatoire 
C      	     FICPROFIL (CHARACTER*SOS_LENFIC2) : nom complet du fichier resultat PROFIL
C                                       (repertoire + nom fichier + extension)
C 
C         -AP.Log ==> Obligatoire  
C            FICLOGPROFIL (CHARACTER*SOS_LENFIC2) : Fichier trace du programme PROFIL 
C						  (0 ou chemin complet au fichier Trace)
C                 (repertoire + nom fichier + extension)
C
C         -AP.MOT ==> Obligatoire 
C  	      PROG_ARG (CHARACTER*10)  : Epaisseur optique rayleigh.
C                                       (format de décodage : libre --> Float)
C	       
C         -AP.HR ==> Optionnel : Associe a AP.MOT >= 0.0001   
C 	      PROG_ARG (CHARACTER*10)  : Echelle de hauteur moleculaire (km).
C                                       (format de décodage : libre --> Float)
C  
C         -AP.AOT ==> Obligatoire 
C 	      PROG_ARG (CHARACTER*10)  : Epaisseur optique aerosols (non tronquee).
C                                       (format de décodage : libre --> Float)
C
C         -AP.Type ==> Obligatoire 
C 	      PROG_ARG (CHARACTER*10)  : Type de profil des aerosols simule.
C                                       (format de décodage : libre --> Integer) 
C      		        1 : Profil d'aerosols selon une decroissance exponentielle par échelle de hauteur
C			2 : Couche homogene d'aerosols entre deux altitudes
C	       	       
C         -AP.AerHS.HA ==> Optionnel : Associe a -AP.Type 1  et AP.AOT >= 0.0001
C 	      PROG_ARG (CHARACTER*10)  : Echelle de hauteur du profil d'aerosols (km).
C                                       (format de décodage : libre --> Float)
C  	 	 	       	       
C         -AP.AerLayer.Zmin ==> Optionnel : Associe a -AP.Type 2  
C 	      PROG_ARG (CHARACTER*10)  : Altitude minimale de la couche d'aerosols (km).
C                                       (format de décodage : libre --> Float)
C 	   	 
C         -AP.AerLayer.Zmax ==> Optionnel : Associe a -AP.Type 2  
C 	      PROG_ARG (CHARACTER*10)  : Altitude maximale de la couche d'aerosols (km).
C                                       (format de décodage : libre --> Float)
C 	 	 
C
C
C Resultats fournis par le programme
C ----------------------------------
C   Fichier resultat PROFIL contenant :
C	- l'indice du niveau
C	- l'epaisseur optique pour chaque niveau
C	- le pourcentage d'aerosols par niveau
C	- le pourcentage de molecules par niveau.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Lecture ou ecriture dans un fichier
C     -Parametres non valables
C
C     Affichage d'un message a l'ecran, arret du programme et
C     retour du status 1 au shell
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel aux routines:
C     - SOS_DISC
C
C==============================================================================
      PROGRAM SOS_PROFIL
      
      IMPLICIT NONE

C Constantes PROFIL
C -----------------

      INTEGER SOS_NB_PROFIL	!Valeur max du parametre profil 
      PARAMETER (SOS_NB_PROFIL=2)

C* Definition des variables                                                 
C--------------------------------------------

      INTEGER*2 IPROFIL	     ! (E) Type du profil.
      
C*   Molecules :
      DOUBLE PRECISION TR    ! (E) Epaisseur optique rayleigh.
      DOUBLE PRECISION HR    ! (E) Echelle de hauteur du profil moleculaire (km).

C*   Aerosols :
      DOUBLE PRECISION TA    ! (E) Epaisseur optique des aerosols (non tronquee).
      DOUBLE PRECISION HA    ! (E) Echelle de hauteur du profil d'aerosols (km).
      DOUBLE PRECISION ZMIN  ! (E) Altitude minimale de la couche molecules + aerosols.
      DOUBLE PRECISION ZMAX  ! (E) Altitude maximale de la couche molecules + aerosols.
      
      DOUBLE PRECISION H(0:SOS_OS_NT)	 ! Epaisseur optique a chaque niveau.
      DOUBLE PRECISION XDEL(0:SOS_OS_NT) ! Pourcentage d'aerosols par niveau.
      DOUBLE PRECISION YDEL(0:SOS_OS_NT) ! Pourcentage de molecules par niveau.
 
C*   Variables specifiques au profil pour une couche d'aerosols entre deux altitudes :     
      DOUBLE PRECISION VR_C1		! Epaisseur optique moleculaire de la couche 1.	
      DOUBLE PRECISION VR_C2		! Epaisseur optique moleculaire de la couche 2.
      DOUBLE PRECISION VR_C3		! Epaisseur optique moleculaire de la couche 3.
      DOUBLE PRECISION VR_SC		! Epaisseur optique moleculaire entre deux niveaux.

      INTEGER*4 NBSC_C1		! Nombre de sous couches dans la couche 1
      INTEGER*4 NBSC_C2		! Nombre de sous couches dans la couche 2
      INTEGER*4 NBSC_C3 	! Nombre de sous couches dans la couche 3
      INTEGER*4 NB_TR	 	! Nombre de couches de transition

      INTEGER*2 NBARGS		! Nombre d'arguments en entrée
      INTEGER*2 DEF_FICPROFIL	! Indicateur d'affectation d'une valeur a FICPROFIL 
      INTEGER*2 DEF_FICLOGPROFIL ! Indicateur d'affectation d'une valeur a FICLOGPROFIL 
      INTEGER*2 DEF_TR 		! Indicateur d'affectation d'une valeur a TR
      INTEGER*2 DEF_HR 		! Indicateur d'affectation d'une valeur a HR
      INTEGER*2 DEF_TA 		! Indicateur d'affectation d'une valeur a TA
      INTEGER*2 DEF_IPROFIL	! Indicateur d'affectation d'une valeur a IPROFIL
      INTEGER*2 DEF_HA 		! Indicateur d'affectation d'une valeur a HA
      INTEGER*2 DEF_ZMIN 	! Indicateur d'affectation d'une valeur a ZMIN
      INTEGER*2 DEF_ZMAX	! Indicateur d'affectation d'une valeur a ZMAX

      DOUBLE PRECISION DELTA_Z	! Ecart d'altitude entre deux niveaux (en km).
      
C*   Autres variables 
      DOUBLE PRECISION Z 		! Altitude d'un niveau (en km).
            
      DOUBLE PRECISION Hmol(0:SOS_OS_NT)   ! Profil d'epaisseur optique moleculaire.
      DOUBLE PRECISION Haer(0:SOS_OS_NT)   ! Profil d'epaisseur optique des aerosols.
    
      DOUBLE PRECISION VA,VR
            
      INTEGER*4 I  		! Niveau.
      INTEGER*4 IARGC           ! Type de la fonction IARGC
      INTEGER*4 NPARAM 		! Numero du parametre en entrée (type INTEGER*4 obligatoire
                                ! pour le f77)
      
      
      LOGICAL TRACE        	! = vrai,si ecriture dans le fichier trace

      CHARACTER*10 PROG_ARG		    ! Arguments du programme

      CHARACTER*SOS_LENFIC2 FICPROFIL       ! chemin et nom du fichier resultat PROFIL 
      CHARACTER*SOS_LENFIC2 FICLOGPROFIL    ! chemin et nom du fichier Trace
      CHARACTER*SOS_LENKEYWORD KEYWORD	    ! Mot clef du passage d'arguments.


C* Lecture des parametres/arguments d'entree et decodage 
C-------------------------------------------------------

C --- Initialisation des indicateurs de definition des parametres (0 pour indefinis)
      DEF_FICPROFIL=0 
      DEF_FICLOGPROFIL=0 
      DEF_TR=0
      DEF_HR=0
      DEF_TA=0 
      DEF_IPROFIL=0
      DEF_HA=0
      DEF_ZMIN=0
      DEF_ZMAX=0
      

C --- Nombre d'arguments en entree
      NBARGS = IARGC()

C --- Boucle de lecture des couples d'arguments       
      DO NPARAM=1,NBARGS, 2
         CALL GETARG(NPARAM,KEYWORD)
	 
	 IF (KEYWORD.EQ."-AP.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICPROFIL)	!Fichier resultat
	     DEF_FICPROFIL=1
	 ENDIF	 
	 	 
	 IF (KEYWORD.EQ."-AP.Log") THEN
	     CALL GETARG (NPARAM+1,FICLOGPROFIL) !Fichier Trace du programme
	     DEF_FICLOGPROFIL=1
	 ENDIF	 	  	 
	 	 
	 IF (KEYWORD.EQ."-AP.MOT") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Epaisseur optique moleculaire
	     READ(PROG_ARG,*,ERR=990)TR
	     DEF_TR=1
	 ENDIF	
	 	 
	 IF (KEYWORD.EQ."-AP.HR") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Echelle de hauteur moleculaire (km)
	     READ(PROG_ARG,*,ERR=990)HR
	     DEF_HR=1
	 ENDIF	

	 IF (KEYWORD.EQ."-AP.AOT") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Epaisseur optique aerosols
	     READ(PROG_ARG,*,ERR=990)TA
	     DEF_TA=1
	 ENDIF	

	 IF (KEYWORD.EQ."-AP.Type") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Type de profil des aerosols simule
	     READ(PROG_ARG,*,ERR=990)IPROFIL 
	     DEF_IPROFIL=1
	 ENDIF	
	 	 
	 IF (KEYWORD.EQ."-AP.AerHS.HA") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Echelle de hauteur des aerosols (km)
	     READ(PROG_ARG,*,ERR=990)HA
	     DEF_HA=1
	 ENDIF	
	 	 
	 IF (KEYWORD.EQ."-AP.AerLayer.Zmin") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Altitude minimale de la couche d'aerosols (km)
	     READ(PROG_ARG,*,ERR=990)ZMIN
	     DEF_ZMIN=1
	 ENDIF		 
	 	 
	 IF (KEYWORD.EQ."-AP.AerLayer.Zmax") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Altitude maximale de la couche d'aerosols (km)
	     READ(PROG_ARG,*,ERR=990)ZMAX
	     DEF_ZMAX=1
	 ENDIF	
	 
      ENDDO	!Fin de boucle de lecture des arguments
	  	


	 
C --- Gestion de la coherence des arguments 

C     Argument FICPROFIL obligatoire
      IF (DEF_FICPROFIL.EQ.0) GOTO 801
      IF (DEF_FICLOGPROFIL.EQ.0) GOTO 802    
      IF (DEF_TR.EQ.0) GOTO 803
      IF (DEF_TA.EQ.0) GOTO 805
      IF ((DEF_HR.EQ.0).AND.(TR.GT.0.0001)) GOTO 804
      IF (DEF_IPROFIL.EQ.0) GOTO 806

      IF (IPROFIL.EQ.1) THEN
          IF ((DEF_HA.EQ.0).AND.(TA.GT.0.0001)) GOTO 807
      ENDIF
      
      IF (IPROFIL.EQ.2) THEN
	  IF ((DEF_ZMIN.EQ.0).OR.(DEF_ZMAX.EQ.0)) GOTO 808
      ENDIF
      
      
          
C* Controle de la validite des arguments 
C---------------------------------------
      IF ((TR.LT.0.D-00).OR.(TA.LT.0.D-00)) GOTO 1005
      
      !HR ne doit etre definie que si TR >= 0.0001
      IF ((HR.LE.0.D-00).AND.(TR.GE.0.0001)) GOTO 1007
      !HA ne doit etre definie que si TA >= 0.0001 (pour cas IPROFIL = 1)
      IF ((IPROFIL.EQ.1).AND.(HA.LE.0.D-00).AND.(TA.GE.0.0001))GOTO 1008
      
      IF ((IPROFIL.LT.1).OR.(IPROFIL.GT.SOS_NB_PROFIL)) GOTO 1009

C* Ouverture du fichier Trace 
C-----------------------------
 ! Creation d'un fichier TRACE ou non
      IF ((FICLOGPROFIL(1:1).EQ.'0').AND.
     &    (FICLOGPROFIL(2:2).EQ.' ')) THEN
	TRACE=.FALSE.
      ELSE
        TRACE=.TRUE.
	OPEN (88,FILE=FICLOGPROFIL,ERR=912)
      END IF
      
	
C* Initialisation
C  --------------
       DO I=0,SOS_OS_NT
	  Hmol(I)=0.
	  Haer(I)=0.
	  H(I)=0.
	  XDEL(I)=0.
	  YDEL(I)=0.
       ENDDO
     
C======================================================
C* CAS 1: Cas d'une couche sans diffusion moleculaire : 
C*        le profil est homogene.
C======================================================
       IF (TR.LT.0.0001) THEN
          DO 433 I=0,SOS_OS_NT
             YDEL(I) = 0.000
             XDEL(I) = 1.000
             H(I) = I*TA/SOS_OS_NT
	     Haer(I)= H(I)
 433      CONTINUE
          GOTO  7777   
       ENDIF

C======================================================       
C* CAS 2: Cas d'une couche sans diffusion aerosol : 
C*        le profil est homogene
C======================================================
       IF (TA.LT.0.0001) THEN
          DO 633 I=0,SOS_OS_NT
             YDEL(I) = 1.000
             XDEL(I) = 0.000
             H(I) = I*TR/SOS_OS_NT
	     Hmol(I)= H(I)
 633      CONTINUE
          GOTO  7777
       ENDIF
       
C======================================================
C* CAS 3: diffusion des aerosols et des molecules.
C======================================================

      IF(IPROFIL.EQ.1) THEN


C* Initialisation du niveau d'altitude a 300 km.
C------------------------------------------------------
         Z=300.D+00

C* Calcul de l'epaisseur optique et de la contribution a la 
C* diffusion des molecules et des aerosols du TOA (300 km) au sol.
C----------------------------------------------------------------
         DO 14 I=0,SOS_OS_NT

C* Determination de l'altitude pour l'epaisseur optique au niveau I.
C------------------------------------------------------------------
            IF (I.NE.0) THEN
               CALL SOS_DISC(TA,HA,TR,HR,I,H(I-1),YDEL(I-1),Z)
	    ENDIF

C* Calcul des proportions moleculaires et aerosols.
C------------------------------------------------------
            VR=TR*DEXP(-Z/HR)
            VA=TA*DEXP(-Z/HA)
	    Hmol(I)=VR
	    Haer(I)=VA
            H(I)=VA+VR
            VR=VR/HR
            VA=VA/HA
            VR=VR/(VA+VR)
            XDEL(I)=(1.00-VR)
            YDEL(I)=VR
  14     CONTINUE


C======================================================
C* CAS 4: diffusion des aerosols et des molecules.
C======================================================

      ELSE IF(IPROFIL.EQ.2) THEN
      
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
 
 
C* Epaisseur optique moleculaire des couches
C  -----------------------------------------
	 VR_C1=TR*DEXP(-(ZMAX+SOS_DZTRANSI)/HR)
	 VR_C2=TR*( DEXP(-ZMIN/HR)-DEXP(-(ZMAX+SOS_DZTRANSI)/HR) )
	 IF(ZMIN.EQ.0.D-00) THEN
	    VR_C3 = 0.D+00
	    NB_TR = 1
	 ELSE   
	    VR_C3=TR*( 1.D+00-DEXP(-(ZMIN-SOS_DZTRANSI)/HR) )
	    NB_TR = 2     
	 ENDIF 
	   
C* Nombre de sous-couches dans chaque couche
C  -----------------------------------------

         NBSC_C1 = (SOS_OS_NT-NB_TR)*VR_C1/(TR+TA) 
         NBSC_C1 = MAX( SOS_PROFIL_MIN_NBC, NBSC_C1 )
	 
	 IF(ZMIN.EQ.0.D+00) THEN
	    NBSC_C3 = 0	    
	 ELSE
	    NBSC_C3 = (SOS_OS_NT-NB_TR)*VR_C3/(TR+TA)	 
	    NBSC_C3 = MAX( SOS_PROFIL_MIN_NBC, NBSC_C3 )
	 ENDIF
	 	 
	 NBSC_C2 = (SOS_OS_NT-NB_TR)-NBSC_C1-NBSC_C3



C* Profil moleculaire : couche 1
C  -----------------------------
	 VR_SC=VR_C1/NBSC_C1

	 DO I=1,NBSC_C1
	      Hmol(I)=Hmol(I-1)+VR_SC
	      Haer(I)=0.
	 ENDDO
	       
C* Transition : ZMAX+SOS_DZTRANSI -> ZMAX
C  -------------------------------------
	I=NBSC_C1+1
	Hmol(I)=TR*DEXP(-ZMAX/HR)
	VR_SC=Hmol(I)-Hmol(I-1)
	Haer(I)=Haer(I-1)+(TA*VR_SC/VR_C2)
	      	       
C* Melange homogene: aerosols + molecules: couche 2
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
	    
C* Transition : ZMIN -> ZMIN-SOS_DZTRANSI
C  -------------------------------------
	    I=NBSC_C1+NBSC_C2+2
	    Hmol(I)=TR*DEXP(-(ZMIN-SOS_DZTRANSI)/HR)
	    Haer(I)=Haer(I-1)  

C* Profil moleculaire : couche 3
C  -----------------------------
	       
	    VR_SC=VR_C3/NBSC_C3
	       
	    DO I = NBSC_C1+NBSC_C2+3 , SOS_OS_NT
	          Hmol(I)=VR_SC+Hmol(I-1)
		  Haer(I)=Haer(I-1)
	    ENDDO
	       
	ENDIF !(Zmin.ne.0)


	       	       	       
C* Calcul de H, XDEL, YDEL
C  -----------------------	       
        XDEL(0)=0.
	YDEL(0)=1
	
	DO I=1,SOS_OS_NT
	       	  
	   H(I)=Hmol(I)+Haer(I)
	   
	   IF(Haer(I).EQ.Haer(I-1)) THEN
		XDEL(I)=0.
		YDEL(I)=1
	   ELSE
		YDEL(I)=1/(1+(TA/VR_C2))
		XDEL(I)=1-YDEL(I)
	   ENDIF
		    
	ENDDO
	    
	    
      ENDIF !Iprofil
      
      
C* Etiquette pour terminer les CAS 1 et 2 
C* (equivalent d'un SINON)
C------------------------------------------------------ 
 7777 CONTINUE

C* Ecriture du fichier trace
C-------------------------------------------
      IF (TRACE) THEN      
         WRITE(88,*,err=921) 'ROUTINE PROFIL'
	 WRITE(88,*,err=921) '=============='
	 WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) '  Profil sans ajustement de  troncature '
         WRITE(88,*,err=921) '      et sans application sur XDEL '// 
     &      'de l''albedo de simple diffusion aerosols.'     
         WRITE(88,*,err=921) ' '
         WRITE(88,*,err=921) 'PARAMETRES DE SIMULATIONS'
	 WRITE(88,*,err=921) '------------------------- '
	 WRITE(88,*,err=921) 
     &   '  Epaisseur optique moleculaire :  ',TR
	 IF (TR.GT.0.0001)  WRITE(88,*,err=921) 
     &   '  Echelle de hauteur du profil moleculaire (km) : ',HR
	 WRITE(88,*,err=921) 
     &   '  Epaisseur optique des aerosols : ',TA
     
      IF(TA.LT.0.0001) WRITE(88,*,err=921) '  Profil moleculaire pur'
      IF(TR.LT.0.0001) WRITE(88,*,err=921) '  Profil d''aerosols pur'
      
      IF((IPROFIL.NE.2).AND.(TA.GE.0.0001).AND.(TR.GE.0.0001)) THEN
	 WRITE(88,*,err=921) 
     &   '  Echelle de hauteur du profil d''aerosols (km) : ',HA
         WRITE(88,*,err=921)     
     &      '  Profil par echelles de hauteur '    
      ENDIF
      
      IF((TA.GE.0.0001).AND.(TR.GE.0.0001)) THEN
      
         IF(IPROFIL.EQ.2) THEN
	    
	    WRITE(88,*,err=921)     
     &      '  Profil pour une couche d''aerosols d''altitude variable '    
	    WRITE(88,*,err=921)     
     &      '  Couche molecules + aerosols : Altitude minimale : ',ZMIN
	    WRITE(88,*,err=921)     
     &      '  Couche molecules + aerosols : Altitude maximale : ',ZMAX
            WRITE(88,*,err=921) ' '
            WRITE(88,*,err=921) 'RESULTATS'
            WRITE(88,*,err=921) '--------- '
	    WRITE(88,*,err=921) 
     &	    '  Nombre de couches de transition       : ',NB_TR
	    WRITE(88,*,err=921) 
     &	    '  Nombre de sous-couches de la couche 1 : ',NBSC_C1
	    WRITE(88,*,err=921) 
     &	    '  Nombre de sous-couches de la couche 2 : ',NBSC_C2
	    WRITE(88,*,err=921)
     &	    '  Nombre de sous-couches de la couche 3 : ',NBSC_C3
	    
	 ENDIF ! Iprofil=2 
	 
   	    WRITE(88,*,err=921) ' '
	    WRITE(88,10,err=921) 'I','ALT','Hmol(I)','Haer(I)',
     &                           'H(I)','XDEL(I)','YDEL(I)'
	    DO I=0,SOS_OS_NT
	       Z=HR*DLOG(TR/Hmol(I))
	       WRITE(88,15,err=921) I,Z,Hmol(I),Haer(I),
     &                              H(I),XDEL(I),YDEL(I)
	    ENDDO
	    
      ELSE ! Cas sans molecules ou sans aerosols
     
	    WRITE(88,*,err=921) ' '	    
	    WRITE(88,30,err=921) 'I','Hmol(I)','Haer(I)',
     &                           'H(I)','XDEL(I)','YDEL(I)'
	    DO I=0,SOS_OS_NT
	       WRITE(88,35,err=921) I,Hmol(I),Haer(I),
     &                              H(I),XDEL(I),YDEL(I)
	    ENDDO      
      
      ENDIF  ! TA >= 0.0001 et TR >= 0.0001
         
   	 WRITE(88,*,err=921) ' '
	 	   	   	   	   	   
      END IF ! Trace


C* Enregistrement des parametres du profil atmospherique dans un fichier
C  ---------------------------------------------------------------------
      OPEN(2,FILE=FICPROFIL,ERR=991)
      
      DO I=0,SOS_OS_NT
         WRITE(2,20,err=992)I,H(I),XDEL(I),YDEL(I)
      ENDDO

C* Fermeture fichier
C-------------------
      CLOSE(88)		! Fermeture du fichier trace PROFIL
      CLOSE(2)		! Fermeture du fichier resultat PROFIL

C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  801 WRITE(6,*) '  ERROR on parameters for SOS_PROFIL'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AP.ResFile Value'
      CALL EXIT(1)
      GOTO 9999 
    
  802 WRITE(6,*) '  ERROR on parameters for SOS_PROFIL'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AP.Log Value'
      CALL EXIT(1)
      GOTO 9999  
    
  803 WRITE(6,*) '  ERROR on parameters for SOS_PROFIL'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AP.MOT Value   (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999 
    
  804 WRITE(6,*) '  ERROR on parameters for SOS_PROFIL'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AP.HR Value   (format : F9)'
      CALL EXIT(1)
      GOTO 9999    
   
  805 WRITE(6,*) '  ERROR on parameters for SOS_PROFIL'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AP.AOT Value   (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999  
              
  806 WRITE(6,*) '  ERROR on parameters for SOS_PROFIL'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AP.Type Value   (format : I2)'
      CALL EXIT(1)
      GOTO 9999  

  807 WRITE(6,*) '  ERROR on parameters for SOS_PROFIL'
      WRITE(6,*) '  --> Height scale aerosol profil'
      WRITE(6,*) '      requires (-Keyword Value): ', 
     &           '     -AP.AerHS.HA Value   (format : F9)'
      CALL EXIT(1)
      GOTO 9999        

  808 WRITE(6,*) '  ERROR on parameters for SOS_PROFIL'
      WRITE(6,*) '  --> Aerosol layer between Zmin and Zmax'
      WRITE(6,*) '      requires (-Keyword Value): ' 
      WRITE(6,*) '     -AP.AerLayer.Zmin Value   (format : F9.5)'
      WRITE(6,*) '     -AP.AerLayer.Zmax Value   (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999      
             
  990 WRITE(6,*) '  ERROR on argument loading for SOS_PROFIL'
      WRITE(6,*) '  --> Check required formats.'
      CALL EXIT(1)
      goto 9999
  991 WRITE(6,*) '  ERROR on SOS_PROFIL result file opening'
      CALL EXIT(1)
      goto 9999
  992 WRITE(6,*) '  ERROR on SOS_PROFIL result file writing'
      CALL EXIT(1)
      goto 9999
  912 WRITE(6,*) '  ERROR on logfile opening for SOS_PROFIL'
      CALL EXIT(1)
      goto 9999
  921 WRITE(6,*) '  ERROR on logfile writing for SOS_PROFIL'
      CALL EXIT(1)
      goto 9999
 1005 WRITE(6,*) '  ERROR on optical thickness value :'
      WRITE(6,*) '  --> -AP.MOT and -AP.AOT have to be positive'
      CALL EXIT(1)
      goto 9999
 1007 WRITE(6,*) '  ERROR on height scale value :'
      WRITE(6,*) '  --> -AP.HR has to be positive'
      CALL EXIT(1)
      goto 9999
 1008 WRITE(6,*) '  ERROR on height scale value :'
      WRITE(6,*) '  --> -AP.AerHS.HA has to be positive'
      CALL EXIT(1)
      goto 9999
 1009 WRITE(6,*) '  ERROR -AP.Type : wrong profile option value'
      CALL EXIT(1)
      goto 9999 
 1010 WRITE(6,*) '  ERROR on boundary altitudes of the aerosol layer :'
      WRITE(6,*) '  --> check -AP.AerLayer.Zmin and AP.AerLayer.Zmax'
      CALL EXIT(1)
      goto 9999
                 
C* Format
C---------   
   10 FORMAT(5X,A1,4X,A3,6X,A7,2X,A7,2X,A4,5X,2(A7,2X))
   15 FORMAT(2X,I4,F11.5,6(F9.5))
   20 FORMAT(2X,I4,3(F10.6))
   30 FORMAT(5X,A1,2X,A7,2X,A7,2X,A4,5X,2(A7,2X))
   35 FORMAT(2X,I4,6(F9.5))
   
 9999 STOP
      END	!FIN DU PROGRAMME SOS_PROFIL




C==============================================================================
C PROCEDURE: SOS_DISC
C ==========
C      Cette procedure estime l'altitude a attribuer au niveau I du profil
C      atmospherique.
C      L'estimation depend de l'epaisseur optique du niveau I-1 et de sa
C      contribution moleculaire.
C
C      A) La couche atmospherique du sol au niveau I-1 est decoupee en
C      sous couches de memes epaisseurs optiques. On estime ainsi une premiere
C      valeur de l'epaisseur optique du niveau I.
C
C      B) On evalue l'altitude correspondante a ce niveau I par une methode de
C      dichotomie.
C
C      C) La proportion moleculaire est calculee pour le niveau I. Si elle
C      s'ecarte de la proportion du niveau I-1 de plus de 75%, la variation
C      est trop importante. On reduit d'un facteur 2 le pas en epaisseur
C      optique entre les niveaux I-1 et I. On reprend l'estimation de
C      l'altitude correspondante au point B.
C
C Description des parametres
C --------------------------
C   Molecules :
C      TR  (double)  (E)   Epaisseur optique rayleigh.
C      HR  (double)  (E)   Echelle de hauteur du profil moleculaire (km).
C
C   Aerosols :
C      TA  (double)  (E)   Epaisseur optique des aerosols.
C      HA  (double)  (E)   Echelle de hauteur du profil d'aerosols (km).
C
C   Informations sur le niveau I-1 :
C      I     (I4)      (E)  Numero du niveau traite.
C      TIM1  (double)  (E)  Epaisseur optique totale du niveau I-1.
C      DD    (double)  (E)  Pourcentage moleculaire au niveau I-1.
C
C      ZX    (double)  (S)  Altitude estimee pour le niveau I.
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
      SUBROUTINE SOS_DISC(TA,HA,TR,HR,I,TIM1,DD,ZX)

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------

C*   Molecules :
      DOUBLE PRECISION TR	! Epaisseur optique rayleigh.
      DOUBLE PRECISION HR	! Echelle de hauteur du profil moleculaire.

C*   Aerosols :
      DOUBLE PRECISION TA	! Epaisseur optique des aerosols (non tronquee).
      DOUBLE PRECISION HA	! Echelle de hauteur du profil d'aerosols.

C*   Informations sur le niveau I-1 :
      DOUBLE PRECISION TIM1	! Epaisseur optique totale du niveau I-1.
      DOUBLE PRECISION DD	! Pourcentage moleculaire au niveau I-1.

      DOUBLE PRECISION ZX	! Altitude estimee pour le niveau I.
      	
      DOUBLE PRECISION DT	! Pas en epaisseur optique.
      DOUBLE PRECISION TI	! Epaisseur optique estimee au niveau I.
      DOUBLE PRECISION ZMAX	! Altitude (km).
      DOUBLE PRECISION ZMIN
      DOUBLE PRECISION ZMOY	! Altitude moyenne.
      DOUBLE PRECISION TZMOY	! Epaisseur optique totale pour l'altitude Zmoy.
      DOUBLE PRECISION XD	! Ecart absolu.
      DOUBLE PRECISION DELTA
      DOUBLE PRECISION ECART

      INTEGER*4 I		 ! Numero du niveau traite.


C------------------------------------------------------
C* A) Decoupe de l'atmosphere entre le sol et le niveau I-1
C*    en NT-(I-1) couches de memes epaisseurs optiques DT.
C------------------------------------------------------
      DT=2.*(TA+TR-TIM1)/(SOS_OS_NT-I+1.)
      
  99  DT=DT/2.

C* Estimation de l'epaisseur optique du niveau I en fonction
C* de celle du niveau I-1 et de l'estimation du pas DT.
C------------------------------------------------------
      TI=TIM1+DT
      
C------------------------------------------------------
C* B) Estimation de l'altitude correspondant a l'epaisseur 
C*    optique TI du niveau I-1.
C------------------------------------------------------

C*   Initialisation des altitudes .
C--------------------------------------
      ZMAX=300.
      ZMIN=0.

C*   Initialisation
C---------------------------
  706 ZMOY=(ZMAX+ZMIN)/2.
  
C*   Epaisseur optique totale pour l'altitude Zmoy.
C----------------------------------------------------
      TZMOY=TA*DEXP(-ZMOY/HA)+TR*DEXP(-ZMOY/HR)
      
C*   Ecart absolu avec l'epaisseur optique estimee du niveau I.
C-------------------------------------------------------------
      XD=DABS(TI-TZMOY)
      
C*   Test de sortie si Zmoy correspond au niveau I.
C------------------------------------------------------
      IF(XD.LT..000001) GOTO 705
      
C*   Affinage de l'encadrement de l'altitude du niveau I.
C------------------------------------------------------
      IF(TI-TZMOY) 701,703,703
      
 701  ZMIN=ZMOY
      GOTO 706
      
 703  ZMAX=ZMOY
      GOTO 706
      
C*   Affectation du resultat
C------------------------------------------------------
 705  ZX=ZMOY
 
C------------------------------------------------------
C* C) Comparaison des proportions moleculaires entre 
C*    le niveau I-1 et celui en cours.
C------------------------------------------------------
      IF(TR.EQ.0.) GOTO 112
      
C*   Evaluation de la contribution moleculaire au niveau I 
C*   pour l'altitude ZX estimee.
C------------------------------------------------------
      DELTA=1./(1.+TA*HR/TR/HA*DEXP(ZX*(1./HR-1./HA)))
      
C*   Comparaison avec celle du niveau I-1. On reduit l'ecart 
C*   entre le niveau I-1 et le niveau I si les contributions 
C*   moleculaires different de plus de 75%.
C-----------------------------------------------------------
      ECART=DABS((DD-DELTA)/DD)
      
      IF(ECART.GT.0.75) GOTO 99
      
 112  CONTINUE
 
      RETURN
      END	!FIN DE LA PROCEDURE SOS_DISC



