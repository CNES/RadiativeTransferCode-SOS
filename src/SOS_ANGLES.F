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
C* FICHIER: SOS_ANGLES.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Definition des angles utilises pour la discretisation des proprietes 
C*       radiatives simulees : fonctions de phase, BRDF/BPDF, champs de luminance. 
C*       - Pour les calculs de fonctions de phase.
C*       - Pour les calculs de matrice de BRDF/BPDF et champs de luminance.
C*       ==> Estime les angles de Gauss et les poids associes pour les integrations spatiales.
C*       ==> Ajoute des angles utilisateurs de poids nuls (sans impact sur les valeurs integrees).
C*       ==> Ajoute l'angle zenithal solaire (de poids nuls).
C*
C* AUTEUR: 
C*   Code de Gauss initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et developpements complementaires : CS SI, B.Lafrance.
C*
C*
C* MOD:VERSION:1.0: 11/01/2010 Creation
C*
C*                                   
C*    - Ajout de la constante SOS_DEFAULT_NBMU_LUM (de SOS.h) pour affectation par
C*      defaut au nombre d'angles de Gauss a utiliser NBMU_GAUSS s'il n'est pas precise
C*      par l'utilisateur 
C*
C* MOD:VERSION:1.1: 24/06/2010 
C*     Les etiquettes de FORMAT en ,X, passent en ,1X, pour eviter des erreurs de compilation.
C*   
C* MOD:VERSION:1.2: 08/12/2015
C*    - Suppression des formats à la lecture des paramètres d'entrée en faveur 
C*      de READ(PROG_ARG,*,ERR=NNN)
C*      ==> Permet une compilation sous gfortran et évite les mauvaises conversions 
C*          si la valeur d'entrée ne respecte pas strictement le format attendu 
C*          (par exemple, saisie d'un entier N au lieu d'un réel N. 
C*           aboutissant à une valeur interprétée fausse). 
C*
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_LENKEYWORD : Taille maximale des Keywords pour le passage des arguments.
C    SOS_NBMAX_USER_ANGLES : Nombre maximal d'angles utilisateurs
C    SOS_PI : Valeur de PI
C    SOS_MIE_NBMU_MAX : Nombre maximal d'angles de Gauss positifs 
C                       pour les calculs de MIE
C    SOS_OS_NBMU_MAX  : Nombre maximal d'angles de Gauss positifs 
C                       pour les calculs de luminance
C    SOS_OS_NB_MAX :    Valeur maximale pour l'ordre limite du developpement en 
C                       polynomes de Legendre.
C    SOS_OS_NS_MAX :    Valeur maximale pour l'ordre limite du developpement en polynomes 
C                       de Legendre des elements de la matrice de Fresnel 
C                       et en series de Fourier des luminances.
C    SOS_OS_NM_MAX :    Valeur maximale pour l'ordre limite  du developpement 
C                       en series de Fourier de la fonction G.
C    SOS_NBANGLES_MAX : Nombre maximal d'angles positifs pour la 
C                       declaration des tableaux XMU et XHR de la routine    
C    SOS_DEFAULT_NBMU_MIE : Nombre d'angles de Gauss positifs par defaut pour les 
C                           calculs de fonctions de phase (si non renseigne par l'utilisateur).
C    SOS_DEFAULT_OS_NB : Ordre par defaut des developpements des fonctions de phase 
C                        en Polynomes de Legendre.
C    SOS_DEFAULT_NBMU_LUM : Nombre d'angles de Gauss positifs par defaut pour les 
C                           calculs de luminance (si non renseigne par l'utilisateur).
C    SOS_DEFAULT_OS_NS    : Ordre maximal par defaut du developpement en polynomes de Legendre   
C                           des elements de la matrice de Fresnel.
C    SOS_DEFAULT_OS_NM    : Ordre maximal par defaut du developpement en series de Fourier  
C                           de la fonction G des calculs de matrice de reflexion.
C    SOS_SEUIL_ECART_MUS  : Seuil pour le test d'ecart entre cos(teta(j)) et cos(tetas).
C----------------------------------------------------------------------------
#include "SOS.h"



C==============================================================================
C PROGRAMME: SOS_ANGLES
C ==========
C      Programme de generation des fichiers d'angles utilises par les codes SOS_AEROSOLS,
C      SOS_SURFACE et SOS.
C
C Donnees en entree du programme 
C ------------------------------
C    -->Les parametres/arguments du programme
C
C        Les arguments sont passes sous la forme d'une liste de couples (-Keyword Value)
C
C        Definition des mots clef, statut obligatoire / option, format de la valeur associee :
C         -ANG.Rad.NbGauss  ==> Optionnel : si non defini, on utilise par defaut la valeur SOS_DEFAULT_NBMU_LUM 
C            PROG_ARG (CHARACTER*10)  : Nombre d'angles de Gauss pour les calculs de luminance.
C                                       (format de décodage : libre --> Integer) 
C
C         -ANG.Rad.UserAngFile  ==> Optionnel : 
C             FICANGLES_USER_LUM (CHARACTER*SOS_LENFIC2) : nom complet du fichier des angles "utilisateurs"
C                                                          pour les calculs de luminance et BRDF
C                 (repertoire + nom fichier + extension)
C
C         -ANG.Rad.Thetas ==> Obligatoire 
C 	      PROG_ARG (CHARACTER*10)  : Angle zenithal solaire en degres (0 < tetas < 90).
C      		                        (angle complementaire : mus= -cos(tetas))
C                                       (format de décodage : libre --> Float)
C
C         -ANG.Rad.ResFile  ==> Obligatoire (chemin complet au fichier ANGLES_LUM) 
C             FICRES_LUM (CHARACTER*SOS_LENFIC2) : nom complet du fichier des angles utilises 
C                                                  pour les calculs de luminance et BRDF
C                 (repertoire + nom fichier + extension)
C
C         -ANG.Aer.NbGauss  ==> Optionnel : si non defini, on utilise par defaut la valeur SOS_DEFAULT_NBMU_MIE 
C            PROG_ARG (CHARACTER*10)  : Nombre d'angles de Gauss pour les calculs de fonctions de phase.
C                                       (format de décodage : libre --> Integer) 
C
C         -ANG.Aer.UserAngFile  ==> Optionnel : 
C             FICANGLES_USER_MIE (CHARACTER*SOS_LENFIC2) : nom complet du fichier des angles "utilisateurs"
C                                                          pour les calculs de fonctions de phase.
C                 (repertoire + nom fichier + extension)
C
C         -ANG.Aer.ResFile  ==> Obligatoire (chemin complet au fichier ANGLES_LUM) 
C             FICRES_MIE (CHARACTER*SOS_LENFIC2) : nom complet du fichier des angles utilises 
C                                                  pour les calculs de fonctions de phase.
C                 (repertoire + nom fichier + extension)
C
C         -ANG.Log ==> Obligatoire (0 ou chemin complet au fichier Trace)
C             FICLOG (CHARACTER*SOS_LENFIC2): nom complet du fichier trace
C                                            (repertoire + nom fichier + extension)
C                                            ou =0, si pas de  trace	 
C
C
C Resultat fourni par la procedure
C --------------------------------
C   Les fichiers decrivant les angles a utiliser pour les calculs de fonction de phase
C   et pour les calculs de luminance + matrices de reflexion.
C   Fichiers contenant egalement les parametres precisant les ordres de developpement
C   des diverses fonctions : OS_NB (fonctions de phase), 
C                            OS_NS (ordre max des developpements des
C                                   fonctions de Fresnel et des luminances)
C   			     OS_NM (ordre maximal du developpement en series de Fourier  
C                                   de la fonction G des calculs de matrice de reflexion
C 
C Description des fichiers utilises
C ---------------------------------
C
C    Fichier utilisateur FICANGLES_USER: 
C           contient une liste d'angles (en degres) a ajouter aux angles de Gauss
C           pour lesquels le poids attribues sera nuls (=> sans impact sur les 
C           integrations spatiales)
C
C        Format : Une valeur d'angle par ligne d'enregistrement.
C                 Les angles doivent etre compris entre 0 et 90 degres.
C                 Nombre d'enregistrements maximal : SOS_NBMAX_USER_ANGLES
C
C    Fichier trace : identifiant = 99 
C
C    Fichier resultat FICRES_MIE:
C           contient les informations suivantes
C		NB_TOTAL_ANGLES :  nb d'angles total a utiliser
C               NB_GAUSS_ANGLES :  nb d'angles de Gauss
C               ANGLES_USERFILE :  Nom du fichier des angles utilisateur
C               INTERNAL_OS_NB :   Ordre des developpements des fonctions de phase 
C				   en Polynomes de Legendre a utiliser.
C               INDEX COS_ANGLE WEIGHT : lignes donnant le numero d'angle, son cosinus et son poids
C                                        au format I4,1X,2D21.14
C
C    Fichier resultat FICRES_LUM:
C           contient les informations suivantes
C		NB_TOTAL_ANGLES :    nb d'angles total a utiliser
C               NB_GAUSS_ANGLES :    nb d'angles de Gauss
C               ANGLES_USERFILE :  Nom du fichier des angles utilisateur
C	        SOLAR ZENITH ANGLE : Angle zenithal solaire (degres)
C		INTERNAL_IMUS :      Numero d'indice de table correspondant
C               INTERNAL_OS_NB :   Ordre des developpements des fonctions de phase 
C				   en Polynomes de Legendre a utiliser.
C               INTERNAL_OS_NS :     Ordre maximal a utiliser pour le developpement en
C                                    polynomes de Legendre des elements de la matrice de Fresnel 
C                                    et en series de Fourier des calculs de luminance.
C               INTERNAL_OS_NM :     Ordre maximal du developpement en series de Fourier  
C                                    de la fonction G des calculs de matrice de reflexion.
C               INDEX COS_ANGLE WEIGHT OUTPUT: lignes donnant le numero d'angle, son cosinus 
C                                              et son poids, s'il est donne en sortie "utilisateur"
C                                              au format I4,1X,2D21.14,1X,I4
C         
C
C Common utilise:
C --------------
C  Aucun
C
C Cas d'erreur :
C ------------
C     - Decodage des parametres incorrect
C     - Erreur dans un sous-programme
C     - Erreur a l'ouverture, lecture et ecriture du fichier Trace 
C       et des fichiers E/S :  fichier des angles utilisateurs et fichiers resultats
C     - Parametres non valables
C     - Incoherence de la dimension SOS_NBANGLES_MAX des tableaux d'angles de la routine 
C       par rapport a celle des tableaux utilisees dans le code SOS :
C          calculs de MIE : SOS_MIE_NBMU_MAX et de luminance : SOS_OS_NBMU_MAX.
C     - Incoherence d'un angle utilisateur par rapport au domaine attendu (0 a 90 degres)
C     - Nombre d'angles utilisateurs trop important devant SOS_NBMAX_USER_ANGLES
C     - Nombre d'angles utiles (de Gauss + utilisateurs) trop important devant ???????  
C
C     Affichage d'un message a l'ecran, arret du programme et 
C     retour du status 1 au shell
C
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel aux routines:
C      - SOS_ANGLES_GAUSS_USER
C      - SOS_GAUSS  
C
C==============================================================================

      PROGRAM SOS_ANGLES
      
      IMPLICIT NONE


C Constantes ANGLES
C -----------------
      DOUBLE PRECISION CONVDEGRAD       !Coefficient de conversion des degres en radians.
      PARAMETER (CONVDEGRAD=SOS_PI/180.D+00)

C* Definition des variables                                                 
C--------------------------------------------
      
      
      DOUBLE PRECISION TETAS	! Angle zenithal solaire en degres (0 < tetas < 90).
      				! (angle complementaire : mus= -cos(tetas))
      DOUBLE PRECISION XMUS     ! Cosinus de l'angle zenithal solaire
      DOUBLE PRECISION XMU(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX) !Cosinus des angles (de Gauss et utilisateur)
      DOUBLE PRECISION XHR(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX) !Poids des angles (de Gauss et utilisateur)
      INTEGER*2 IND_ANGOUT(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX) !Indices des tableaux d'angles pour lesquels 
      							       !il doit y avoir une sortie dans le fichier
							       !résultat des angles utilisateur 
							       !(1 : pour une sortie requise)
      DOUBLE PRECISION XMU_copy(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX) !copie de XMU avant ajout de TETAS
      DOUBLE PRECISION XHR_copy(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX) !copie de XHR avant ajout de TETAS
      INTEGER*2 IND_ANGOUT_copy(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX) !copie de IND_ANGOUT avant ajout de TETAS
      
      INTEGER*2 NBARGS		        ! Nombre d'arguments en entrée
      INTEGER*2 DEF_NBMU_GAUSS_LUM      ! Indicateur d'affectation d'une valeur a NBMU_GAUSS_LUM
      INTEGER*2 DEF_FICANGLES_USER_LUM  ! Indicateur d'affectation d'une valeur a FICANGLES_USER_LUM
      INTEGER*2 DEF_FICRES_LUM  	! Indicateur d'affectation d'une valeur a FICRES_LUM
      INTEGER*2 DEF_NBMU_GAUSS_MIE      ! Indicateur d'affectation d'une valeur a NBMU_GAUSS_MIE  
      INTEGER*2 DEF_FICANGLES_USER_MIE  ! Indicateur d'affectation d'une valeur a FICANGLES_USER_MIE
      INTEGER*2 DEF_FICRES_MIE  	! Indicateur d'affectation d'une valeur a FICRES_MIE
      INTEGER*2 DEF_TETAS		! Indicateur d'affectation d'une valeur a TETAS
      INTEGER*2 DEF_FICLOG	        ! Indicateur d'affectation d'une valeur a FICLOG  
      
      INTEGER*4 IARGC     	        ! Type de la fonction IARGC
      INTEGER*4 NPARAM 		        ! Numero du parametre en entrée (type INTEGER*4 obligatoire
                                        ! pour le f77)
      INTEGER*4 NBMU_GAUSS_LUM          ! Nombre d'angles de Gauss utilises pour les calculs de luminance      
      INTEGER*4 NBMU_LUM                ! Nombre d'angles effectifs des simulations de luminance
      INTEGER*4 NBMU_GAUSS_MIE          ! Nombre d'angles de Gauss utilises pour les calculs de Mie      
      INTEGER*4 NBMU_MIE                ! Nombre d'angles effectifs des calculs de Mie
      INTEGER*4 IER     		! code d'erreur =0 si pas d'erreur, =-1 sinon
      INTEGER*4 J          		! Indices des tableaux d'angles
      INTEGER*4 IMUS          		! Indice de l'angle zenithal solaire dans la table des angles
      INTEGER*4 OS_NB 			! Ordre des developpements des fonctions de phase 
					! en Polynomes de Legendre.
      INTEGER*4 OS_NS    		! Ordre maximal du developpement en polynomes de Legendre   
					! des elements de la matrice de Fresnel.
      INTEGER*4 OS_NM    		! Ordre maximal du developpement en series de Fourier  
					! de la fonction G des calculs de matrice de reflexion.
               
      CHARACTER*SOS_LENFIC2 FICLOG              ! nom du fichier trace
      			!Calculs de luminance :
      CHARACTER*SOS_LENFIC2 FICANGLES_USER_LUM  ! chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*SOS_LENFIC2 FICRES_LUM          ! chemin et nom du fichier des angles utilises 
      			!Calculs de fonctions de phase :
      CHARACTER*SOS_LENFIC2 FICANGLES_USER_MIE  ! chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*SOS_LENFIC2 FICRES_MIE          ! chemin et nom du fichier des angles utilises 
      
      CHARACTER*SOS_LENKEYWORD KEYWORD	    ! Mot clef du passage d'arguments.

      CHARACTER*10 PROG_ARG		    ! Arguments du programme

      LOGICAL TRACE     ! = vrai,si ecriture dans le fichier trace
      

C* Lecture des parametres/arguments d'entree et decodage 
C-------------------------------------------------------

C --- Initialisation des indicateurs de definition des parametres (0 pour indefinis)
      
      DEF_NBMU_GAUSS_LUM=0
      DEF_FICANGLES_USER_LUM=0   
      DEF_TETAS=0
      DEF_FICRES_LUM=0
      DEF_NBMU_GAUSS_MIE=0
      DEF_FICANGLES_USER_MIE=0  
      DEF_FICRES_MIE=0
      DEF_FICLOG=0                      


C --- Nombre d'arguments en entree
      NBARGS = IARGC()

C --- Boucle de lecture des couples d'arguments       
      DO NPARAM=1,NBARGS, 2
         CALL GETARG(NPARAM,KEYWORD)
	 
	 IF (KEYWORD.EQ."-ANG.Rad.NbGauss") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Nb d'angles de Gauss pour les simulations de luminance
	     READ(PROG_ARG,*,ERR=900) NBMU_GAUSS_LUM
	     DEF_NBMU_GAUSS_LUM=1
	 ENDIF	 
	 	 
	 IF (KEYWORD.EQ."-ANG.Rad.UserAngFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_USER_LUM)	!Nom du fichier des angles "utilisateurs" (chemin complet)
	     DEF_FICANGLES_USER_LUM=1			!pour les calculs de luminance
	 ENDIF	  	 
	 		 
	 IF (KEYWORD.EQ."-ANG.Rad.Thetas") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)			!Angle zenithal solaire (degres)
	     READ(PROG_ARG,*,ERR=910) TETAS 
             DEF_TETAS=1
	 ENDIF	

	 IF (KEYWORD.EQ."-ANG.Rad.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICRES_LUM)	        !Nom du fichier resultat (chemin complet)
	     DEF_FICRES_LUM=1				!pour les calculs de luminance
	 ENDIF	  	
				
	 IF (KEYWORD.EQ."-ANG.Aer.NbGauss") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Nb d'angles de Gauss pour les calculs de Mie
	     READ(PROG_ARG,*,ERR=905) NBMU_GAUSS_MIE
	     DEF_NBMU_GAUSS_MIE=1
	 ENDIF	 
	 	 
	 IF (KEYWORD.EQ."-ANG.Aer.UserAngFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_USER_MIE)	!Nom du fichier des angles "utilisateurs" (chemin complet)
	     DEF_FICANGLES_USER_MIE=1			!pour les calculs de Mie
	 ENDIF	  	 
	 	 	 
	 IF (KEYWORD.EQ."-ANG.Aer.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICRES_MIE)	        !Nom du fichier resultat (chemin complet)
	     DEF_FICRES_MIE=1				!pour les calculs de Mie
	 ENDIF	  	
	 	 
	 IF (KEYWORD.EQ."-ANG.Log") THEN
	     CALL GETARG (NPARAM+1,FICLOG)		!Nom du fichier Trace.
	     DEF_FICLOG=1
	 ENDIF		 
	 
	 
      ENDDO	!Fin de boucle de lecture des arguments
	  	

	 
C --- Gestion de la coherence des arguments 

C     Angle solaire : obligatoire           
      IF (DEF_TETAS.EQ.0) GOTO 920
      IF ((TETAS.LT.0.D-00).OR.(TETAS.GE.90.D-00)) GOTO 925
      
C     Fichier resultat pour les luminances : obligatoire           
      IF (DEF_FICRES_LUM.EQ.0) GOTO 930

C     Fichier resultat pour les fonctions de phase : obligatoire           
      IF (DEF_FICRES_MIE.EQ.0) GOTO 935

C     Argument FICLOG obligatoire          
      IF (DEF_FICLOG.EQ.0) THEN
         GOTO 940
      ELSE !Determination de la requete de suivi des traitements
         IF ((FICLOG(1:1).EQ.'0').AND.(FICLOG(2:2).EQ.' ')) THEN
	    TRACE=.FALSE.
         ELSE
            TRACE=.TRUE.
         ENDIF
      ENDIF

C     Argument optionnels
C     --> Affectation de la valeur par defaut si NBMU_GAUSS_MIE non renseigne
      IF (DEF_NBMU_GAUSS_MIE.EQ.0) THEN
          NBMU_GAUSS_MIE = SOS_DEFAULT_NBMU_MIE
	  OS_NB = SOS_DEFAULT_OS_NB
      ELSE
          !Pour une fonction de phase discretisee en NBMU_GAUSS_MIE angles de Gauss (positifs), 
	  !sa decomposition en polynomes de Legendre ne peut etre correcte que si elle est poussee
	  !jusqu'a un ordre valant au moins 2 fois le nb d'angles de Gauss utilises dans la 
	  !discretisation de la fonction.
          OS_NB = 2*NBMU_GAUSS_MIE	  
      ENDIF 

C     --> Affectation de la valeur par defaut si NBMU_GAUSS_LUM non renseigne
      IF (DEF_NBMU_GAUSS_LUM.EQ.0) THEN
          NBMU_GAUSS_LUM = SOS_DEFAULT_NBMU_LUM
	  OS_NS = SOS_DEFAULT_OS_NS
	  OS_NM = SOS_DEFAULT_OS_NM
      ELSE
          !Pour un champs de luminance discretise en NBMU_GAUSS_LUM angles de Gauss (positifs), 
	  !la decomposition en polynomes de Legendre des fonctions de reflexion de Fresnel
	  !et des luminances atmospherique n'est  correcte que si elle est poussee jusqu'a 
	  !un ordre valant au moins 2 fois le nb d'angles de Gauss utilises pour discretiser
	  !la luminance.
          OS_NS = 2*NBMU_GAUSS_LUM
	  
	  !Le calcul de la fonction G des calculs de BRDF (glitter, ...) impose un ordre maximal 
	  !des developpements en series de Fourier tel que que OS_NM >= OS_NB + OS_NS
	  OS_NM = OS_NB + OS_NS
      ENDIF 
      

C     --> Controle des ordres limites des developpements pour les calculs 
C         des coefficients des matrices de phase et de reflexion
      IF (OS_NB.LT.2*NBMU_GAUSS_MIE) GOTO 945
      IF (OS_NS.LT.2*NBMU_GAUSS_LUM) GOTO 946
      IF (OS_NM.LT.(OS_NB+OS_NS))    GOTO 947
       
      IF (OS_NB.GT.SOS_OS_NB_MAX) GOTO 9451
      IF (OS_NS.GT.SOS_OS_NS_MAX) GOTO 9461
      IF (OS_NM.GT.SOS_OS_NM_MAX) GOTO 9471
            
C     --> Utilisation d'un fichier d'angles "utilisateur"  
      IF (DEF_FICANGLES_USER_LUM.EQ.0) THEN
          FICANGLES_USER_LUM="NO_USER_ANGLES"
      ENDIF
      
C     --> Utilisation d'un fichier d'angles "utilisateur"  
      IF (DEF_FICANGLES_USER_MIE.EQ.0) THEN
          FICANGLES_USER_MIE="NO_USER_ANGLES"
      ENDIF           
      


C* Ouverture du fichier Trace 
C------------------------------------------------------
      IF (TRACE) THEN
         OPEN (99,FILE=FICLOG,ERR=950)
      ENDIF
      

C* Preparation des angles pour les calculs de fonction de phase 
C--------------------------------------------------------------      
      CALL SOS_ANGLES_GAUSS_USER("MIE",NBMU_GAUSS_MIE,
     &				 FICANGLES_USER_MIE,TRACE,
     &                           NBMU_MIE, XMU, XHR, IER) 
      
      IF (IER.EQ.-1) GOTO 955
      

C* Ecriture du fichier d'angles resultat pour les calculs de fonction
C* de phase 
C-------------------------------------------------------------------         
      OPEN(UNIT=10,FILE=FICRES_MIE,FORM='FORMATTED',ERR=960)
      WRITE(10,100,ERR=961) NBMU_MIE      
      WRITE(10,110,ERR=961) NBMU_GAUSS_MIE 
      WRITE(10,120,ERR=961) FICANGLES_USER_MIE
      WRITE(10,200,ERR=961) OS_NB
      WRITE(10,250,ERR=961) 
      DO J=1,NBMU_MIE
         WRITE(10,300,ERR=961) J, XMU(J), XHR(J)
      ENDDO
      
C* Preparation des angles de Gauss et Utilisateur
C* pour les calculs de luminance et BRDF/BPDF 
C-------------------------------------------------------------------      
      CALL SOS_ANGLES_GAUSS_USER("LUM",NBMU_GAUSS_LUM,
     &                           FICANGLES_USER_LUM,TRACE,
     &                           NBMU_LUM, XMU, XHR, IER) 
      
      IF (IER.EQ.-1) GOTO 955      

C* Preparation de la table des indices a retenir pour une sortie des
C* angles utilisateurs 
C-------------------------------------------------------------------         
      DO J=1,NBMU_LUM
      
         IF (XHR(J).EQ.0.D+00) THEN
	    IND_ANGOUT(J)=1
	 ELSE
	    IND_ANGOUT(J)=0
         ENDIF
	 
      ENDDO

C* Ajout de l'angle d'incidence solaire 
C-------------------------------------------------------------------         
      IMUS=-1
      XMUS=DCOS(TETAS*CONVDEGRAD)
      
      !Recherche d'une correspondance avec un angle deja defini
      !-->On releve juste la valeur de IMUS, les tables ne sont pas a modifier
      J=1
      DO J=1,NBMU_LUM 
         IF (DABS(XMUS-XMU(J)).LT.SOS_SEUIL_ECART_MUS) IMUS=J
      ENDDO
      
      !Si TETAS n'est pas un angle deja defini
      !-->Ajout dans la table des angles, sans sortie requise pour les angles utilisateur
      IF (IMUS.EQ.-1) THEN 
         !On estime l'indice IMUS dans la nouvelle table de NBMU_LUM+1 elements.
	 
	 !-->Archivage des tables sans TETAS
         DO J=1,NBMU_LUM
	    XMU_copy(J) = XMU(J)
	    XHR_copy(J) = XHR(J)
	    IND_ANGOUT_copy(J) = IND_ANGOUT(J)
	 ENDDO
	 
	 !La table XMU_copy est ordonnee du plus grand au plus petit element
	 !-->On cherche si XMUS doit etre le premier element de la nouvelle table
	 IF (XMUS.GT.XMU_copy(1)) THEN 
	     IMUS=1
	 ELSE
	     !-->On cherche si XMUS doit etre le dernier element de la nouvelle table
	     IF (XMUS.LT.XMU_copy(NBMU_LUM)) THEN 
	        IMUS=NBMU_LUM+1
	     ELSE
	        !-->On cherche quel element J de la table initiale doit etre remplace par XMUS
	        DO J=1,NBMU_LUM-1
		   IF ( (XMUS.GT.XMU_copy(J+1)).
     &                   AND.(XMUS.LT.XMU_copy(J)) ) IMUS=J+1
                ENDDO 
	     ENDIF
	 ENDIF
	 
	 !-->Mise a jour de la taille des tables
	 NBMU_LUM=NBMU_LUM+1
	 IF (NBMU_LUM.GT.SOS_OS_NBMU_MAX) GOTO 980
	 
	 !-->Positionnement de TETAS dans les tables
	 XMU(IMUS) = XMUS
	 XHR(IMUS) = 0.D+00
	 IND_ANGOUT(IMUS)= 0
	 
	 !-->Completude des tables avec decalage des indices
	 DO J=IMUS+1,NBMU_LUM  
	    XMU(J) = XMU_copy(J-1)
	    XHR(J) = XHR_copy(J-1)
	    IND_ANGOUT(J) = IND_ANGOUT_copy(J-1)
	 ENDDO
	    
        !-->Completudes des angles definis sur J=-NBMU_UTILES,0
         DO J=1,NBMU_LUM
	        XMU(-J) = -XMU(J)
	        XHR(-J) =  XHR(J)
		IND_ANGOUT(-J) = IND_ANGOUT(J)
         ENDDO 
         XMU(0) = 0.D+00
         XHR(0) = 0.D+00
	 IND_ANGOUT(0) = 0
	 	 
      ENDIF !Test si IMUS non defini
 
      
C* Information de Trace 
C-------------------------------------------------------------------         
      IF (TRACE) THEN
         WRITE (99,*,ERR=990) "Ajout de l'angle zenithal solaire"
	 WRITE (99,*,ERR=990) "---------------------------------"
	 WRITE (99,*,ERR=990) " "
	 WRITE(99,*,ERR=990) "TETAS : ", TETAS
	 WRITE(99,*,ERR=990) "IMUS  : ", IMUS
	 WRITE(99,*,ERR=990) "Nombre final d'angles : ", NBMU_LUM 
         WRITE (99,*,ERR=990) " " 
	 WRITE (99,*,err=990) "Liste des angles utilises :"
	 WRITE (99,*,err=990) "   Indice, mu, angle (deg.), poids, user"
	 DO J=1,NBMU_LUM
	     WRITE(99,223,ERR=990) J, XMU(J),DACOS(XMU(J))/CONVDEGRAD,
     &	                           XHR(J),IND_ANGOUT(J)  
	 ENDDO
	 WRITE (99,*,ERR=990) " " 
	 WRITE (99,*,ERR=990) " "  	
      ENDIF !Fichier Trace	       
		

C* Ecriture du fichier d'angles resultat pour les calculs de fonction
C* de phase 
C-------------------------------------------------------------------         
      OPEN(UNIT=20,FILE=FICRES_LUM,FORM='FORMATTED',ERR=960)
      WRITE(20,100,ERR=961) NBMU_LUM
      WRITE(20,110,ERR=961) NBMU_GAUSS_LUM
      WRITE(20,120,ERR=961) FICANGLES_USER_LUM
      WRITE(20,410,ERR=961) TETAS
      WRITE(20,420,ERR=961) IMUS
      WRITE(20,200,ERR=961) OS_NB
      WRITE(20,500,ERR=961) OS_NS
      WRITE(20,510,ERR=961) OS_NM
      WRITE(20,550,ERR=961) 
      DO J=1,NBMU_LUM
         WRITE(20,600,ERR=961) J, XMU(J), XHR(J), IND_ANGOUT(J)
      ENDDO
     
     
C* Fermeture fichier
C-------------------
      CLOSE(10)
      CLOSE(20)
      		
C* Fin nominale 
C-------------------
      GOTO 9999     

C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  900 WRITE(6,*) '  ERROR on argument loading for SOS_ANGLES'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -ANG.Rad.NbGauss value : Number of Gauss'
      WRITE(6,*) '        angles for radiance computations (I4).' 
      CALL EXIT(1)
      GOTO 9999     

  905 WRITE(6,*) '  ERROR on argument loading for SOS_ANGLES'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -ANG.Aer.NbGauss value : Number of Gauss'
      WRITE(6,*) '        angles for radiance computations (I4).' 
      CALL EXIT(1)
      GOTO 9999     
      
  910 WRITE(6,*) '  ERROR on argument loading for SOS_ANGLES'
      WRITE(6,*) '  --> Check required formats for -ANG.Thetas value :'
      WRITE(6,*) '        Solar zenithal angle (F9, degrees).'
      CALL EXIT(1)
      GOTO 9999  

  920 WRITE(6,*) '  ERROR on parameters for SOS_ANGLES'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -ANG.Rad.Thetas Value'
      CALL EXIT(1)
      GOTO 9999 

  925 WRITE(6,*) '  ERROR on -ANG.Rad.Thetas value : ',
     &           'out of available range'
      CALL EXIT(1)
      GOTO 9999 
      
      
  930 WRITE(6,*) '  ERROR on parameters for SOS_ANGLES'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -ANG.Rad.ResFile Value'
      CALL EXIT(1)
      GOTO 9999 

  935 WRITE(6,*) '  ERROR on parameters for SOS_ANGLES'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -ANG.Aer.ResFile Value'
      CALL EXIT(1)
      GOTO 9999        
                 
  940 WRITE(6,*) '  ERROR on parameters for SOS_ANGLES'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -ANG.Log Value'
      CALL EXIT(1)
      GOTO 9999 

  945 WRITE(6,*) '  ERROR : SOS_DEFAULT_OS_NB <= 2*SOS_DEFAULT_NBMU_MIE'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
      CALL EXIT(1)
      GOTO 9999 

  946 WRITE(6,*) '  ERROR : SOS_DEFAULT_OS_NS <= 2*SOS_DEFAULT_NBMU_LUM'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
      CALL EXIT(1)
      GOTO 9999   

  947 WRITE(6,*) '  ERROR : ',
     &    'SOS_DEFAULT_OS_NM < SOS_DEFAULT_OS_NB + SOS_DEFAULT_OS_NS'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
      CALL EXIT(1)
      GOTO 9999                

 9451 WRITE(6,*) '  ERROR : OS_NB > SOS_OS_NB_MAX'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
      CALL EXIT(1)
      GOTO 9999 

 9461 WRITE(6,*) '  ERROR : OS_NS > SOS_OS_NS_MAX'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
      CALL EXIT(1)
      GOTO 9999   

 9471 WRITE(6,*) '  ERROR : OS_NM > SOS_OS_NM_MAX'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
      CALL EXIT(1)
      GOTO 9999                
            
  950 WRITE(6,*) '  ERROR on logfile opening for SOS_ANGLES'
      CALL EXIT(1)  
      GOTO 9999  

  955 WRITE(6,*) '  ERROR on subroutine SOS_ANGLES_GAUSS_USER'
      CALL EXIT(1)
      GOTO 9999

  960 WRITE(6,*) '  ERROR on ANGLES file opening for SOS_ANGLES'
      CALL EXIT(1)  
      GOTO 9999  

  961 WRITE(6,*) '  ERROR on ANGLES file writing for SOS_ANGLES'
      CALL EXIT(1)  
      GOTO 9999  
                   
  980 WRITE(6,*) "Error on the number of used radiance angles "
      WRITE(6,*) "(gauss angles + possibly solar angle)  : "
      WRITE(6,*) "  --> greater than the maximal number allowed "
      WRITE(6,*) "       SOS_OS_NBMU_MAX =",SOS_OS_NBMU_MAX
      CALL EXIT(1)
      GOTO 9999 
      
  990 WRITE(6,*) '  ERROR on logfile writing for SOS_ANGLES'
      CALL EXIT(1)
      GOTO 9999
                                              
 9999 STOP
 
C* Format
C---------   

      
  100 FORMAT(17hNB_TOTAL_ANGLES :,I4)
  110 FORMAT(17hNB_GAUSS_ANGLES :,I4)
  120 FORMAT(17hANGLES_USERFILE :,A)
  200 FORMAT(16hINTERNAL_OS_NB :,I4)
  250 FORMAT(35hINDEX   COS_ANGLE            WEIGHT)
  300 FORMAT(I4,1X,2D21.14)
  
  410 FORMAT(20hSOLAR ZENITH ANGLE :,F7.3)
  420 FORMAT(15hINTERNAL_IMUS :,I4)
  500 FORMAT(16hINTERNAL_OS_NS :,I4)
  510 FORMAT(16hINTERNAL_OS_NM :,I4)
  550 FORMAT(55hINDEX   COS_ANGLE            WEIGHT              OUTPUT)
  600 FORMAT(I4,1X,2D21.14,1X,I4)
  223 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X,I2)
  
      END	!FIN DU PROGRAMME SOS_ANGLES


      

C==============================================================================
C PROCEDURE: SOS_ANGLES_GAUSS_USER
C ==========
C
C  La routine calcule les angles et poids de Gauss et ajoute, si requis, 
C  les angles du fichier utilisateur.
C
C  Elle trie les angles
C
C Description des parametres
C ----------------------------		  
C
C       USAGE (CHAR*3)  (E) : indicateur de l'usage qui sera fait du fichier d'angles : 
C                             calculs de MIE ou simulation du champs de rayonnement LUM	
C	NBMU_GAUSS (I4) (E) : nombre d'angles de Gauss souhaites
C       FICANGLES_USER (CHAR*SOS_LENFIC2) (E) : nom du fichier des angles "utilisateur"
C                                               (chemin complet)
C       TRACE (LOGICAL) (E)  : vrai,si ecriture dans le fichier trace
C       NBMU_UTILES (I4) (S) : nombre d'angles apres concatenation des angles de Gauss 
C                              et des angles "utilisateur"
C	XMU (double)     (S) : tableau des cosinus des angles calcules.
C	XHR (double)     (S) : tableau des poids associes aux angles calcules.	     
C       IER (I4)         (S) : indice d'erreur : 0 si pas d'erreur, -1 sinon.
C     
C Resultat fourni par la procedure
C --------------------------------
C   Les tableaux XMU et XHR des cosinus et poids des angles utiles
C   Le nb d'angles NBMU_UTILES
C   Un retour d'erreur IER 
C   
C Description des fichiers utilises
C ---------------------------------
C
C    Fichier utilisateur FICANGLES_USER: 
C           contient une liste d'angles (en degres) a ajouter aux angles de Gauss
C           pour lesquels le poids attribues sera nuls (=> sans impact sur les 
C           integrations spatiales)
C
C        Format : Une valeur d'angle par ligne d'enregistrement.
C                 Les angles doivent etre compris entre 0 et 90 degres.
C                 Nombre d'enregistrements maximal : SOS_NBMAX_USER_ANGLES
C
C    Fichier trace : identifiant = 99 (conforme a l'ouverture par le code appelant)
C
C Common utilise:
C --------------
C  Aucun
C
C Cas d'erreur :
C ------------
C  Aucun
C
C  Sous programmes utilises:
C --------------------------
C  SOS_GAUSS
C
C==============================================================================
      
      SUBROUTINE  SOS_ANGLES_GAUSS_USER(USAGE,NBMU_GAUSS,FICANGLES_USER,
     &                                  TRACE,NBMU_UTILES,XMU,XHR,IER) 

      IMPLICIT NONE


C Constantes ANGLES
C -----------------
      DOUBLE PRECISION CONVDEGRAD       !Coefficient de conversion des degres en radians.
      PARAMETER (CONVDEGRAD=SOS_PI/180.D+00)

C* Definition des variables                                                 
C--------------------------------------------
      INTEGER*4 NBMU_GAUSS      ! Nombre d'angles de Gauss utilises       
      INTEGER*4 NBMU_USER       ! Nombre d'angles contenus dans le fichier utilisateur 
      INTEGER*4 NBMU_UTILES     ! Nombre d'angles effectifs des simulations 
      INTEGER*4 IER             ! code d'erreur =0 si pas d'erreur, =-1 sinon
      INTEGER*4 J, K            ! Indices des tableaux d'angles
			    
      DOUBLE PRECISION XMU(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX) !Cosinus des angles (de Gauss et utilisateur)
      DOUBLE PRECISION XHR(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX) !Poids des angles (de Gauss et utilisateur)
      DOUBLE PRECISION VAL	!Valeur d'angle lue dans le fichier utilisateur					
      DOUBLE PRECISION MU_COMP(SOS_NBMAX_USER_ANGLES)  
                                        !Cosinus des angles (de Gauss et utilisateur)
 
      CHARACTER*SOS_LENFIC2 FICANGLES_USER  ! chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*3 USAGE     ! Indicateur de l'usage qui sera fait du fichier d'angles : 
      			    ! calculs de MIE ou simulation du champs de rayonnement LUM			


      LOGICAL TRACE     ! = vrai,si ecriture dans le fichier trace

C* Initialisation
C-----------------
      IER=0
      NBMU_USER=0 

C* Controle de validite du parametre USAGE
C------------------------------------------------------      
      IF ((USAGE.NE.'MIE').AND.(USAGE.NE.'LUM')) GOTO 950
      
C* Controle de coherence des dimensionnements
C------------------------------------------------------      
      IF (USAGE.EQ.'MIE') THEN
          IF (SOS_NBANGLES_MAX.LT.SOS_MIE_NBMU_MAX) GOTO 960
      ENDIF
      
      IF (USAGE.EQ.'LUM') THEN
          IF (SOS_NBANGLES_MAX.LT.SOS_OS_NBMU_MAX) GOTO 961
      ENDIF
      
      IF (NBMU_GAUSS.GE.SOS_NBANGLES_MAX) GOTO 962
           
      
C* Calcul des angles et poids de Gauss
C  -----------------------------------
      CALL SOS_GAUSS(NBMU_GAUSS+1, XMU, XHR)
	
C* Lecture des angles utilisateurs
C------------------------------------------------------ 	
      IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
      	 
         OPEN(1,FILE=FICANGLES_USER,STATUS='OLD',ERR=970)      
      
          J=0  
   20	  READ(1,*,ERR=971,END=30) VAL
               IF ((VAL.LT.0.).OR.(VAL.GT.90.)) GOTO 972              
	       J=J+1
	       IF (J.GT.SOS_NBMAX_USER_ANGLES) GOTO 973
	       MU_COMP(J)=DCOS(VAL*CONVDEGRAD)
	  GOTO 20
   30	  CLOSE(1)
          NBMU_USER=J
    
      ENDIF

C* Determination du nombre d'angles utiles 
C------------------------------------------------------      
      NBMU_UTILES = NBMU_GAUSS + NBMU_USER

      IF (USAGE.EQ.'MIE') THEN
          IF (NBMU_UTILES.GT.SOS_MIE_NBMU_MAX) GOTO 980
      ENDIF 
      
      IF (USAGE.EQ.'LUM') THEN
          IF (NBMU_UTILES.GT.SOS_OS_NBMU_MAX) GOTO 981
      ENDIF 

C*  Ajout des angles utilisateurs aux angles de Gauss
C------------------------------------------------------           
      IF (NBMU_USER.NE.0) THEN
      
          DO J=1,NBMU_USER
	     XMU(NBMU_GAUSS+J) = MU_COMP(J)
	     XHR(NBMU_GAUSS+J) = 0.D+00
	  ENDDO 

      ENDIF
      
C*  Tri des angles : cas USAGE="MIE" --> Trie de MUmin à MUmax
C*                   cas USAGE="LUM" --> Trie de MUmax à MUmin
C*    Tri sur les angles definis sur J=1,NBMU_UTILES
C------------------------------------------------------           
      IF (USAGE.EQ.'MIE') THEN
 
          DO J=1,NBMU_UTILES	     
	     DO K=J+1,NBMU_UTILES
	        IF (XMU(J).GT.XMU(K)) THEN		    		    
		    !Permutation des valeurs XMU(J) et XMU(K)
		    VAL=XMU(J)
		    XMU(J)=XMU(K)
		    XMU(K)=VAL
		    !Permutation des valeurs XMU(J) et XHR(K)
		    VAL=XHR(J)
		    XHR(J)=XHR(K)
		    XHR(K)=VAL   	    	    
		ENDIF
	     ENDDO 
	  ENDDO 
	  
      ENDIF !USAGE=MIE	  

      IF (USAGE.EQ.'LUM') THEN
 
          DO J=1,NBMU_UTILES	     
	     DO K=J+1,NBMU_UTILES
	        IF (XMU(J).LT.XMU(K)) THEN		    		    
		    !Permutation des valeurs XMU(J) et XMU(K)
		    VAL=XMU(J)
		    XMU(J)=XMU(K)
		    XMU(K)=VAL
		    !Permutation des valeurs XMU(J) et XHR(K)
		    VAL=XHR(J)
		    XHR(J)=XHR(K)
		    XHR(K)=VAL       
		ENDIF
	     ENDDO 
	  ENDDO 
	  
      ENDIF !USAGE=LUM	  
      	  	
	  	  
C*  Completudes des angles definis sur J=-NBMU_UTILES,0
C------------------------------------------------------      	  
      DO J=1,NBMU_UTILES
	     XMU(-J) = -XMU(J)
	     XHR(-J) = XHR(J)
      ENDDO 
      XMU(0) = 0.D+00
      XHR(0) = 0.D+00
	
C*  Information au fichier TRACE
C------------------------------------------------------   
      IF (TRACE) THEN
         IF (USAGE.EQ.'MIE') THEN
	     WRITE (99,*,ERR=990) "CALCULS DE FONCTIONS DE PHASE"
	     WRITE (99,*,ERR=990) "*****************************"
	 ENDIF
	 IF (USAGE.EQ.'LUM') THEN
	     WRITE (99,*,ERR=990) "CALCULS DE LUMINANCE"
	     WRITE (99,*,ERR=990) "********************"
	 ENDIF
	 WRITE (99,*,ERR=990) " "
	 WRITE(99,*,ERR=990) "Nombre d'angles de Gauss : ", NBMU_GAUSS
	 WRITE (99,*,ERR=990) " "
         IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	     WRITE(99,*,ERR=990)  
     &	       "Usage du fichier d'angles complementaires :",
     &         FICANGLES_USER
     	     WRITE(99,*,ERR=990) "Nb d'angles complementaires :",
     &         (NBMU_UTILES-NBMU_GAUSS)
	 ELSE
	     WRITE(99,*,ERR=990)  
     &	       "Pas de fichier d'angles complementaires"
	 ENDIF
         WRITE (99,*,ERR=990) " " 
	 WRITE (99,*,err=990) "Liste des angles utilises :"
	 WRITE (99,*,err=990) "   Indice, mu, angle (degres), poids"
	 DO J=1,NBMU_UTILES
	     WRITE(99,222,ERR=990) J, XMU(J),
     &	                           DACOS(XMU(J))/CONVDEGRAD, XHR(J)  
	 ENDDO
	 WRITE (99,*,ERR=990) " " 
	 WRITE (99,*,ERR=990) " "  	
      ENDIF !Fichier Trace	        
C* Fin nominale 
C-------------------
      GOTO 9999


C* Cas d'erreur et retour du status -1 au programme appelant
C----------------------------------------------------------            
  950  WRITE(6,*) "Error on USAGE parameters : has to be 'MIE' or 'LUM'"
       IER=-1
       GOTO 9999        
       
  960  WRITE(6,*) "Inconsistence on size of tables in SOS.h."
       WRITE(6,*) " --> required : SOS_NBANGLES_MAX >= SOS_MIE_NBMU_MAX"
       IER=-1
       GOTO 9999   

  961  WRITE(6,*) "Inconsistence on size of tables in SOS.h."
       WRITE(6,*) " --> required : SOS_NBANGLES_MAX >= SOS_OS_NBMU_MAX"
       IER=-1
       GOTO 9999        

  962  WRITE(6,*) "Inconsistence on size of tables in SOS.h."
       WRITE(6,*) " --> Number of required Gauss angles is ", NBMU_GAUSS
       WRITE(6,*) "     greater than the limit SOS_NBANGLES_MAX :", 
     &            SOS_NBANGLES_MAX
       IER=-1
       GOTO 9999        
          
  970  WRITE(6,*) "Error while opening the user's angles file : "
       WRITE(6,*)  FICANGLES_USER
       WRITE(6,*) " --> check its localization"
       IER=-1
       GOTO 9999      

  971  WRITE(6,*) "Error while reading the user's angles file :"
       WRITE(6,*)  FICANGLES_USER
       IER=-1
       GOTO 9999   

  972  WRITE(6,*) "Error on an user's angles value in file: "
       WRITE(6,*)  FICANGLES_USER
       WRITE(6,*) "  --> have to defined between 0 and 90°"
       IER=-1
       GOTO 9999   

  973  WRITE(6,*) "Error on the user's angles file : "
       WRITE(6,*)  FICANGLES_USER
       WRITE(6,*) "  --> The number of angles has to be < or = to"
       WRITE(6,*) "       SOS_NBMAX_USER_ANGLES =",SOS_NBMAX_USER_ANGLES
       IER=-1
       GOTO 9999

  980  WRITE(6,*) "Error on the number of used phaze function angles : "
       WRITE(6,*) "  --> greater than the maximal number allowed "
       WRITE(6,*) "       SOS_MIE_NBMU_MAX =",SOS_MIE_NBMU_MAX
       IER=-1
       GOTO 9999 
       
  981  WRITE(6,*) "Error on the number of used radiance angles : "
       WRITE(6,*) "  --> greater than the maximal number allowed "
       WRITE(6,*) "       SOS_OS_NBMU_MAX =",SOS_OS_NBMU_MAX
       IER=-1
       GOTO 9999 
      
  990 WRITE(6,*) '  ERROR on logfile writing for SOS_ANGLES'
      CALL EXIT(1)
      GOTO 9999
              
9999  RETURN


C* Format
C---------   
  222 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X)
   
      END !FIN DE LA PROCEDURE SOS_ANGLES_GAUSS_USER       
             
	     
	     
	     
C==============================================================================
C PROCEDURE: SOS_GAUSS
C ==========
C
C  Routine originale fournie par le LOA (M. Herman, octobre 2009).
C  Pas de documentation.
C
C  La routine calcule les angles et poids de Gauss pour MM-1 valeurs.
C
C  Les tableaux generes couvrent les valeurs XMU(-MM:MM) et XHR(-MM:MM)
C  avec pour XMU(-MM)=-1 et XMU(MM)=1 de poids nuls.
C
C  Une contrainte est d'avoir SOS_NBANGLES_MAX >= MM
C
C Description des parametres
C ----------------------------
C	MM (I4) (E) : nombre d'angles de Gauss souhaites
C	XMU (double) (S) : tableau des cosinus des angles calcules.
C	XHR (double) (S) : tableau des poids associes aux angles calcules.			  
C
C Resultat fourni par la procedure
C --------------------------------
C   Les tableaux XMU et XHR
C   
C Description des fichiers utilises
C ---------------------------------
C  Aucun
C
C Common utilise:
C --------------
C  Aucun
C
C Cas d'erreur :
C ------------
C  Aucun
C
C  Sous programmes utilises:
C --------------------------
C  Aucun
C
C==============================================================================
      
      SUBROUTINE  SOS_GAUSS(MM,AMU,PMU)
C     N : ORDRE DE LA QUADRATURE SUR XL A XU

C* Definition des variables                         
C*-----------------------------------------------------------------
      PARAMETER (IX=2*SOS_NBANGLES_MAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(IX),PA(IX),W(IX),R(IX)
      DOUBLE PRECISION AMU(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX)  !Angles de Gauss
      DOUBLE PRECISION PMU(-SOS_NBANGLES_MAX:SOS_NBANGLES_MAX)  !Poids de Gauss

      TOL = 1.0D-15
      PI = 3.141592653589793D+00
      N=2*MM-2
      XL=-1
      XU=+1
      AA = 2.0D+00/PI**2
      AB = -62.0D+00/(3.0D+00*PI**4)
      AC = 15116.0D+00/(15.0D+00*PI**6)
      AD = -12554474.D+00/(105.0D+00*PI**8)
      PA(1) = 1.D0
      EN = N
      NP1 = N+1
      U= 1.0D+00-(2.0D+00/PI)**2
      D = 1.0D+00/DSQRT((EN+0.5D+00)**2+U/4.0D+00)

      DO 100 I= 1,N
         SM = I
         AZ = 4.0D+00*SM-1.0D+00
         AE = AA/AZ
         AF = AB/AZ**3
         AG = AC/AZ**5
         AH = AD/AZ**7
         Z(I) = 0.25D+00*PI*(AZ+AE+AF+AG+AH)
	 
100   CONTINUE
	
      DO 200 K = 1,N
         X = DCOS(Z(K)*D)
    1    PA(2) = X
    	
         DO 201  NN = 3,NP1
         ENN = NN-1
         PA(NN) =
     &   ((2.0D+00*ENN-1.0D+00)*X*PA(NN-1)-(ENN-1.0D+00)*PA(NN-2))/ENN
  201    CONTINUE
  
         PNP = EN*(PA(N)-X*PA(NP1))/(1.0D+00-X*X)
         XI = X-PA(NP1)/PNP
         XD = DABS(XI-X)
         XDD = XD-TOL
	 
         IF (XDD) 3,3,2
    2       X = XI
            GO TO 1
    3       R(K) = X
            W(K) = 2.0D+00*(1.0D+00-X*X)/(EN*PA(N))**2
  200 CONTINUE
	  
      AP = (XU-XL)/2.D0
      BP = (XU+XL)/2.D0
      DO I=1,MM-1
         K=MM-I
         AMU(K)=BP+AP*R(I)
         PMU(K)=AP*W(I)
         AMU(-K)=-AMU(K)
         PMU(-K)=PMU(K)
      ENDDO
    	
      AMU(-MM)=-1.
      AMU(MM)=1.
      AMU(0)=0.
      PMU(0)=0.
      PMU(-MM)=0.
      PMU(MM)=0.
C
      RETURN
      END !FIN DE LA PROCEDURE SOS_GAUSS
      
