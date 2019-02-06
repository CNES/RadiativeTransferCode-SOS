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
C* FICHIER: SOS_TRPHI.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module genere les fichiers resultats des parametres de Stokes 
C*       I, Q et U .
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et developpement : CS SI, B.Lafrance, C. Casteras.
C* DATE: 30/04/2002
C*
C* MOD:VERSION:1.0: Mise a niveau du code des OS du LOA
C*
C* MOD:VERSION:3.0: 
C*    - Correction d'une erreur pour le cas de la reflexion de Fresnel
C*      sur mer plate : calcul de la reflexion du faisceau solaire direct
C*      affecte a la direction de visee J=N0.
C*
C*    - Introduction de la reflexion directe des modeles:
C*           - Roujean.
C*           - Rondeaux.
C*           - Breon.
C*           - Nadal.
C*
C* MOD:VERSION:3.1:  Amelioration de commentaires.
C*
C* MOD:VERSION:4.0: 10/01/2010 
C*   - Modification de la taille des tableaux de champs de luminance :
C*     ==> Passage de SOS_OS_NBMU a SOS_OS_NBMU_MAX.
C*
C*   - Ajout du parametre d'entree : NBMU (nb d'angles effectifs pour les simulations)
C*
C*   - Modification des lignes de lecture des tableaux de luminances 
C*     a l'ordre IS de la decomposition en series de Fourier 
C*     ==> limitation aux bornes utiles -NBMU:NBMU
C*
C*   - Adaptation de la routine SOS_TRPHI au passage du parametre NBMU.
C*
C*   - Adaptation de la routine SOS_TRPHI :
C*        * au passage  des parametres SOS.UserAng.ResFileUp et 
C*          SOS.UserAng.ResFileDown precisant les noms des fichiers
C*          pour l'enregistrement des luminances de sortie uniquement 
C*          pour les angles "utilisateur".
C*        * au passage de la table IND_ANGOUT precisant les angles
C*          "utilisateur" pour lesquels il faut une sortie de luminance.
C*
C* MOD:VERSION:4.1: 21/01/2016
C*   - Correction du terme de réflexion directe sur mer plate
C*     ==> Introduction de la luminance solaire directe (normalisée à Es=PI)
C*         MAIS pour la distance Terre-Soleil moyenne (selon la valeur
C*         de SOS_SOLAR_DISC_SOLID_ANGLE définie dans SOS.h.
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    SOS_LENDIR : Longueur des noms de repertoires.
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_PI : valeur de PI                           
C    SOS_OS_NBMU_MAX : Nombre maximal d'angles positifs pour la resolution de l'ETR
C    SOS_SEUIL_Z : Valeur seuil pour le calcul des angles de rotation
C    SOS_SEUIL_X : Valeur seuil pour le calcul des elements M21 et M31 de la 
C                  matrice de reflexion exprimee dans un repere lie au plan 
C                  meridien
C    SOS_SOLAR_DISC_SOLID_ANGLE : Valeur de l'angle solide du disque solaire
C                                 apparent.
C----------------------------------------------------------------------------
#include "SOS.h"


C============================================================================
C PROCEDURE: SOS_TRPHI_OPTION
C ==========
C      Cette procedure genere les fichiers resultats des parametres de Stokes 
C      I, Q et U selon le choix de l'utilisateur.
C      
C      Si ITRPHI = 1, on obtient le plan de visees pour PHIos 
C                     (azimut relatif par rapport a l'azimut solaire)
C         Convention : 
C             L'angle de visee est positif pour le demi plan PHI = PHIos.
C	      L'angle de visee est negatif pour le demi plan PHI = PHIos + 180.
C
C         Deux fichiers sont crees : _ un pour le champ montant, 
C				     _ un pour le champ descendant.
C
C
C      Si ITRPHI = 2, on obtient un diagramme polaire en fonction de l'angle 
C      de visee (teta) et de l'angle azimutal (phi) :
C
C         X =  phi    ,   Y = teta    ,   I  ,  Q  ,  U
C
C	  L'angle de visee croit du centre au bord de l'image.
C
C         L'utilisateur choisit le pas d'echantillonnage en azimut (PAS_PHI).
C
C         Le code genere un fichier pour les parametres de Stokes (I, Q, U) 
C         pour les champs montant et descendant.
C
C         ATTENTION : Les fichiers obtenus peuvent etre volumineux.
C
C
C      Pour le cas mer plate, le code ajoute la réflexion de la luminance solaire
C      directe, descendante du TOA à la surface, atténuée par l'atmosphère.
C      Réflexion symétrique à la direction d'incidence solaire par rapport au zénith.
C      NB : La luminance solaire est calculée pour un éclairement solaire au TOA 
C           valant PI, pour une distance Terre-Soleil moyenne,  
C           la constante SOS_SOLAR_DISC_SOLID_ANGLE donnant l'angle solide du disque
C           solaire apparent pour la distance Terre-Soleil moyenne. 
C
C
C Description des parametres
C --------------------------
C     NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C
C     RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (DOUBLE) : (E) Cosinus des angles.
C     GA(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) (DOUBLE)  : (E) Poids des angles.
C
C     FICHOS   (CHARACTER*SOS_LENFIC2)  (E) :  Nom du fichier resultat des OS 
C                                     (repertoire + nom fichier + extension)
C     TAU (double)    (E) : Epaisseur optique totale de l'atmosphere 
C                          (eventuellement tronquee).
C     TAUOUT (double) (E) : Epaisseur optique au niveau de sortie des resultats pour le
C                           champ montant.
C     IGLI  (I4)      (E) : Option de Glitter :  0 : pas de Glitter.
C                                                1 : calcul du Glitter.
C     N0  (I4)        (E) : Numero de mu de Gauss pour l'angle zenithal solaire. 
C                          (angle complementaire : mus = -rmu(N0)) 
C                          (N0 > 0 => pour en tenir compte).
C     WIND (double)   (E) : Vitesse du vent en m/s.
C     IND_SURF (double) (E) : Indice de refraction relatif surface / air.
C     IFRESNEL  (I4)  (E) : Option de calcul de reflexion de Fresnel 
C                            sur mer plate :
C                        = 0 : pas de reflexion sur mer plate.
C                        = 1 : calcul de la reflexion de Fresnel sur mer plate.
C     IROUJEAN (I4)   (E) : Option de calcul de la reflexion de Roujean
C                        = 0 : pas de reflexion par la BRDF de Roujean.
C                        = 1 : calcul de la reflexion par la BRDF de Roujean.
C     K0_ROUJEAN (double) (E) : Parametre k0 de la BRDF de Roujean.
C     K1_ROUJEAN (double) (E) : Parametre k1 de la BRDF de Roujean.
C     K2_ROUJEAN (double) (E) : Parametre k2 de la BRDF de Roujean.
C     IRONDEAUX (I4)  (E) : Option de calcul de la reflexion de Rondeaux
C                        = 0 : pas de reflexion par la BPDF de Rondeaux.
C                        = 1 : calcul de la reflexion par la BPDF de Rondeaux.
C     IBREON (I4)     (E) : Option de calcul de la reflexion de Breon
C                        = 0 : pas de reflexion par la BPDF de Breon.
C                        = 1 : calcul de la reflexion par la BPDF de Breon.
C     INADAL (I4)     (E) : Option de calcul de la reflexion de Nadal
C                        = 0 : pas de reflexion par la BPDF de Nadal.
C                        = 1 : calcul de la reflexion par la BPDF de Nadal.
C     ALPHA_NADAL (double)  (E) : Parametre alpha du modele de Nadal.         	
C     BETA_NADAL (double)  (E) : Parametre beta du modele de Nadal.          
C     ITRPHI (I4)     (E) : Indice de type de sortie.
C     PHIOS (double)  (E) : Azimut relatif en degres (pour ITRPHI=1).
C     PAS_PHI (I4)    (E) : Pas sur l'azimut en degres (pour ITRPHI=2)
C                           Valeur entiere.
C     IPOLAR  (I4)   (E) : Option de simulation sans polarisation.
C      			   _ 0 : calculs sans polarisation.
C			   _ 1 : calculs normaux avec polarisation.
C     FICHUP   (CHARACTER*SOS_LENFIC2)  (E) : 
C                           Nom du fichier resultat pour le champ montant.
C                           (repertoire + nom fichier + extension)
C     FICHDOWN   (CHARACTER*SOS_LENFIC2)  (E) :
C                           Nom du fichier resultat pour le champ descendant.
C                           (repertoire + nom fichier + extension)
C     IND_ANGOUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (I4) : (E) 
C		 Indices des tableaux d'angles pour lesquels il doit y avoir 
C                une sortie dans le fichier résultat des angles utilisateur 
C               (1 : pour une sortie requise, 0 sinon)				
C
C     FICUP_ANGLES_USER (CHARACTER*SOS_LENFIC2)  (E) :
C                           Nom du fichier resultat pour le champ montant.
C                           (repertoire + nom fichier + extension)
C                           --> sortie pour les angles utilisateur uniquement.
C     FICDOWN_ANGLES_USER (CHARACTER*SOS_LENFIC2)  (E) :
C                           Nom du fichier resultat pour le champ descendant.
C                           (repertoire + nom fichier + extension)
C                           --> sortie pour les angles utilisateur uniquement.    
C
C     IER  (I4)       (S) : Code d'erreur =0 si pas d'erreur, =-1 sinon 	
C
C  
C Autre resultat fourni par la procedure
C -------------------------------------
C   Fichier resultat pour le champ montant
C   Fichier resultat pour le champ descendant
C    --> de base pour tous les angles utiles 
C        + optionnellement pour les angles utilisateurs uniquement
C
C Variable d'environnement
C ------------------------
C   Aucune
C
C
C Description des fichiers utilises
C ---------------------------------
C
C   -->Contenu du fichier resultat des OS
C       Fichier binaire non formatte, cree par SOS_OS
C       Il contient les parametres de Stokes Q , U et I pour un ordre S
C       du developpement en series de Fourier, en fonction de l'angle mu
C       (valeurs en double precision).
C
C   -->Contenu du fichier resultat pour le champ montant
C       Fichier Ascii, format d'un enregistrement (F7.2,2X,3(F13.6,2X)) 
C       ou (2(F7.2,2X),3(F13.6,2X))
C       2 cas : ITRPHI = 1, 
C               on obtient le plan de visees pour PHIos .
C               Nombre d'enregistrements : 
C                   un pour chaque valeur de -NBMU a NBMU (angles utiles)
C               Chaque enregistrement contient 4 valeurs :
C                   angle de visee(en degres) et parametres de Stokes I,Q,U 
C
C               ITRPHI = 2, 
C               on obtient un diagramme polaire en fonction de l'angle 
C               de visee (teta) et de l'angle azimutal (phi) 
C               Nombre d'enregistrements :
C                   un pour chaque valeur de 1 a NBMU (angles utiles)
C               Chaque enregistrement contient 5 valeurs :
C                   coordonnees phi, teta et parametres de Stokes I,Q,U 
C
C   -->Contenu du fichier resultat pour le champ descendant
C      Fichier Ascii, format d'un enregistrement (F7.2,2X,3(F13.6,2X)) 
C      ou (2(F7.2,2X),3(F13.6,2X))
C      Meme description que pour le fichier resultat du champ montant
C
C
C Common utilise:
C --------------
C     Aucun
C     
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C       Erreur si le fichier resultat des OS n'existe pas
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C     
C  Routines utilises par la procedure:
C ------------------------------------
C  Ce programme fait appel aux routines:
C      - SOS_TRPHI
C
C==============================================================================
  
      SUBROUTINE SOS_TRPHI_OPTION(NBMU,RMU,GA,
     &                        FICHOS,TAU,TAUOUT,
     &                        IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &                        IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                        IRONDEAUX,IBREON,
     &                        INADAL,ALPHA_NADAL,BETA_NADAL,
     &                        ITRPHI,PHIOS,PAS_PHI,IPOLAR,
     &			      FICHUP,FICHDOWN,IND_ANGOUT,
     &                        FICUP_ANGLES_USER,FICDOWN_ANGLES_USER,
     &                        IER)
 	

C* Definition des variables                 
C*-----------------------------------------
      DOUBLE PRECISION TAU	! Epaisseur optique totale de l'atmosphere 
				! (eventuellement tronquee).

      DOUBLE PRECISION TAUOUT	! Epaisseur optique au niveau de sortie
      				! des resultats pour le champ montant.
						
      DOUBLE PRECISION WIND	! Vitesse du vent en m/s. 			
      DOUBLE PRECISION IND_SURF	! Indice de refraction relatif surface / air.
      DOUBLE PRECISION K0_ROUJEAN ! Parametre k0 du modele de BRDF de Roujean.
      DOUBLE PRECISION K1_ROUJEAN ! Parametre k1 du modele de BRDF de Roujean.
      DOUBLE PRECISION K2_ROUJEAN ! Parametre k2 du modele de BRDF de Roujean.
      DOUBLE PRECISION ALPHA_NADAL ! Parametre alpha du modele de BPDF de Nadal.
      DOUBLE PRECISION BETA_NADAL  ! Parametre beta du modele de BPDF de Nadal.
      DOUBLE PRECISION PHIOS	! Azimut relatif en degres (pour ITRPHI=1).
      
      DOUBLE PRECISION TETA	! Angle de visee (en degres).
      DOUBLE PRECISION PHI	! Azimut relatif (en radians).
      				     
      DOUBLE PRECISION XIT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Parametre de Stokes I 
                                        ! en fonction de mu pour l'azimut phi.
					!    mu < 0 : champ descendant.
					!    mu > 0 : champ montant.
      DOUBLE PRECISION XQT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Parametre de Stokes Q.
      DOUBLE PRECISION XUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)	! Parametre de Stokes U.
            	  
      
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX) 
                                             !Tableau des angles de Gauss (et utilisateur). 
      DOUBLE PRECISION GA(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
                                             !Poids des angles de Gauss (et utilisateur).

      CHARACTER*SOS_LENFIC2 FICHOS    !Nom du fichier resultat des OS
      CHARACTER*SOS_LENFIC2 FICHUP    !Nom du fichier resultat pour le champ montant
      CHARACTER*SOS_LENFIC2 FICHDOWN  !Nom du fichier resultat pour le champ descendant

      CHARACTER*SOS_LENFIC2 FICUP_ANGLES_USER	! Nom du fichier resultat pour le champ
      				                ! montant (pour les angles utilisateur).					
      CHARACTER*SOS_LENFIC2 FICDOWN_ANGLES_USER	! Nom du fichier resultat pour le champ
      				                ! descendant (pour les angles utilisateur).	

      INTEGER*4 IGLI	! Option de Glitter : 0 : pas de Glitter.
			!                     1 : calcul du Glitter.

      INTEGER*4 IFRESNEL ! Indice precisant si les simulations
      			 ! introduisent la reflexion de Fresnel sur mer plate :
			 !    1 => pour un calcul de la reflexion de Fresnel.
			 !    Sinon => pas de prise en compte.      
	
      INTEGER*4 IROUJEAN ! Indice precisant si les simulations
      			 ! introduisent la reflexion par la BRDF de Roujean :
			 !    1 => pour un calcul avec la BRDF de Roujean.
			 !    Sinon => pas de prise en compte.      
      
      INTEGER*4 IRONDEAUX ! Indice precisant si les simulations
      			  ! introduisent la reflexion par la BPDF de Rondeaux :
			  !    1 => pour un calcul avec la BPDF de Rondeaux.
			  !    Sinon => pas de prise en compte.   
			  
      INTEGER*4 IBREON    ! Indice precisant si les simulations
      			  ! introduisent la reflexion par la BPDF de Breon :
			  !    1 => pour un calcul avec la BPDF de Breon.
			  !    Sinon => pas de prise en compte.   
			  
      INTEGER*4 INADAL	  ! Indice precisant si les simulations
      			  ! introduisent la reflexion par la BPDF de Nadal :
			  !    1 => pour un calcul avec la BPDF de Nadal.
			  !    Sinon => pas de prise en compte.   
			             
      INTEGER*4 N0	! Numero de mu de Gauss pour l'angle zenithal solaire. 
			!  (angle complementaire : mus = -rmu(N0)) 
			!   (N0 > 0 => pour en tenir compte).

      INTEGER*4 ITRPHI	! Indice de type de sortie.
			!   1 =>  Plan de visee pour une valeur d'azimut.		           
			!   2 =>  Diagramme polaire en fonction de la direction
			!         d'azimut et de visee.
			      
      INTEGER*4 PAS_PHI	! Pas sur l'azimut en degres (pour ITRPHI=2, valeur entiere)
      
      INTEGER*4 IPOLAR	! Option de simulation sans polarisation.
      			! 0 : calculs sans polarisation.
			! 1 : calculs normaux avec polarisation.
			      
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon
      INTEGER*4 J	! Indice du tableau des angles de Gauss.
      INTEGER*4 IPHI	! Indice de boucle sur l'azimut relatif.
      INTEGER*4 NBMU    ! nombre d'angles utiles					     

      INTEGER*4 IND_ANGOUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  
                        ! Indices des tableaux d'angles pour lesquels il doit y avoir 
      			! une sortie dans le fichier résultat des angles utilisateur 
			! (1 : pour une sortie requise)				

			   
	
C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0

C* Ouverture des fichiers 
C------------------------------------------------------
     
C     Fichiers resultats  pour le champ montant et le champ descendant
C     Si le fichier existe deja, il est ecrase par le nouveau fichier     
      OPEN(UNIT=14,FILE=FICHUP,FORM='FORMATTED',ERR=991)
      OPEN(UNIT=15,FILE=FICHDOWN,FORM='FORMATTED',ERR=991)

C     Fichiers resultats  pour le champ montant et le champ descendant
C     Si le fichier existe deja, il est ecrase par le nouveau fichier  
C     --> Cas des fichiers optionnels contenant les resultats pour les
C         angles utilisateurs uniquement   
      IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
          OPEN(UNIT=24,FILE=FICUP_ANGLES_USER,FORM='FORMATTED',ERR=991)
      ENDIF
      	  
      IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
          OPEN(UNIT=25,FILE=FICDOWN_ANGLES_USER,FORM='FORMATTED',
     &	       ERR=991)
      ENDIF
      
      
C* Sortie d'un plan azimutal
C------------------------------------------------------
      IF(ITRPHI.EQ.1) THEN

C* Demi plan PHIos + 180 degres : Angles de visee negatifs par convention
C------------------------------------------------------------------------       
         PHI=SOS_PI+PHIOS*SOS_PI/180.D+00
         CALL SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &	                IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &                  IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                  IRONDEAUX,IBREON,
     &                  INADAL,ALPHA_NADAL,BETA_NADAL,
     &                  IPOLAR,XIT,XQT,XUT,IER)
         IF (IER.NE.0) GOTO 9999   ! cas d'erreur

         
C* Cas nominal : tous les angles sont a sortir	 
         DO J=NBMU,1,-1
            TETA=DACOS(RMU(J))*180.D+00/SOS_PI
	    
C*          Enregistrement pour le champ montant
            WRITE(14,55,err=992) -TETA,XIT(J),XQT(J),XUT(J) 
	    
C*          Enregistrement pour le champ descendant
	    WRITE(15,55,err=992) -TETA,XIT(-J),XQT(-J),XUT(-J)	    	     
	    
         ENDDO
	
	 
C* Cas optionnel : champ montant pour les angles utilisateurs uniquement
         IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN

            DO J=NBMU,1,-1   
               IF (IND_ANGOUT(J).EQ.1) THEN
	            TETA=DACOS(RMU(J))*180.D+00/SOS_PI
		    WRITE(24,55,err=992) -TETA,XIT(J),XQT(J),XUT(J)
		ENDIF     	     	    
            ENDDO	 
	 
	 ENDIF
	 
C* Cas optionnel : champ descendant pour les angles utilisateurs uniquement
         IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN

            DO J=NBMU,1,-1   
               IF (IND_ANGOUT(J).EQ.1) THEN
	            TETA=DACOS(RMU(J))*180.D+00/SOS_PI
		    WRITE(25,55,err=992) -TETA,XIT(-J),XQT(-J),XUT(-J)
		ENDIF     	     	    
            ENDDO	 
	 
	 ENDIF	 

	    

	    
C* Demi plan PHIos : Angles de visee positifs par convention.
C-----------------------------------------------------------      
         PHI=PHIOS*SOS_PI/180.D+00
         CALL SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &	                IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &                  IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                  IRONDEAUX,IBREON,
     &                  INADAL,ALPHA_NADAL,BETA_NADAL,
     &                  IPOLAR,XIT,XQT,XUT,IER)
         IF (IER.NE.0) GOTO 9999   ! cas d'erreur

         DO J=1,NBMU
            TETA=DACOS(RMU(J))*180.D+00/SOS_PI
	    
C* Cas nominal : tous les angles sont a sortir
C*          Enregistrement pour le champ montant
            WRITE(14,55,err=992) TETA,XIT(J),XQT(J),XUT(J)
	    
C*          Enregistrement pour le champ descendant
	    WRITE(15,55,err=992) TETA,XIT(-J),XQT(-J),XUT(-J)

         ENDDO	    
      
      
      
C* Cas optionnel : champ montant pour les angles utilisateurs uniquement
	 IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
            DO J=1,NBMU
	        IF (IND_ANGOUT(J).EQ.1) THEN
		    TETA=DACOS(RMU(J))*180.D+00/SOS_PI
		    WRITE(24,55,err=992) TETA,XIT(J),XQT(J),XUT(J)
		ENDIF
	    ENDDO
	 ENDIF    
	 
C* Cas optionnel : champ descendant pour les angles utilisateurs uniquement
	 IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
	    DO J=1,NBMU
	        IF (IND_ANGOUT(J).EQ.1) THEN
		    TETA=DACOS(RMU(J))*180.D+00/SOS_PI
		    WRITE(25,55,err=992) TETA,XIT(-J),XQT(-J),XUT(-J)
		ENDIF
            ENDDO
	 ENDIF
	 
	    
      ENDIF  !ITRPHI=1
      
          
C* Sortie d'un diagramme polaire
C---------------------------------
      IF(ITRPHI.EQ.2) THEN        
	 
	 DO IPHI=0,360,PAS_PHI
	    PHI = SOS_PI * IPHI / 180.D+00
	    CALL SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &	                   IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &                     IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                     IRONDEAUX,IBREON,
     &                     INADAL,ALPHA_NADAL,BETA_NADAL,
     &                     IPOLAR,XIT,XQT,XUT,IER)
            IF (IER.NE.0) GOTO 9999   ! cas d'erreur
     
            PHI = dble(IPHI)
            DO J=1,NBMU
               TETA=DACOS(RMU(J))*180.D+00/SOS_PI

C* Cas nominal : tous les angles sont a sortir	      	       
C*             Enregistrement pour le champ montant      
               WRITE(14,56,err=992) PHI,TETA,XIT(J),XQT(J),XUT(J)

C*             Enregistrement pour le champ descendant	 
	       WRITE(15,56,err=992) PHI,TETA,XIT(-J),XQT(-J),XUT(-J)
       
            ENDDO !MU
	 ENDDO !PHI	        
	 
C* Cas optionnel : champ montant pour les angles utilisateurs uniquement	 
         IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
	    DO IPHI=0,360,PAS_PHI
	       PHI = SOS_PI * IPHI / 180.D+00
	       CALL SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &	                      IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &                        IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                        IRONDEAUX,IBREON,
     &                        INADAL,ALPHA_NADAL,BETA_NADAL,
     &                        IPOLAR,XIT,XQT,XUT,IER)
               IF (IER.NE.0) GOTO 9999   ! cas d'erreur
     
               PHI = dble(IPHI)
               DO J=1,NBMU
	          IF (IND_ANGOUT(J).EQ.1) THEN
		      TETA=DACOS(RMU(J))*180.D+00/SOS_PI
		      WRITE(24,56,err=992) PHI,TETA,XIT(J),XQT(J),XUT(J)
		  ENDIF
               ENDDO !MU
	    ENDDO !PHI	 
	 ENDIF
	 
C* Cas optionnel : champ descendant pour les angles utilisateurs uniquement
	 IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
	    DO IPHI=0,360,PAS_PHI
	       PHI = SOS_PI * IPHI / 180.D+00
	       CALL SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &	                      IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &                        IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                        IRONDEAUX,IBREON,
     &                        INADAL,ALPHA_NADAL,BETA_NADAL,
     &                        IPOLAR,XIT,XQT,XUT,IER)
               IF (IER.NE.0) GOTO 9999   ! cas d'erreur
     
               PHI = dble(IPHI)
               DO J=1,NBMU
	          IF (IND_ANGOUT(J).EQ.1) THEN
		      TETA=DACOS(RMU(J))*180.D+00/SOS_PI
		      WRITE(25,56,err=992)
     &		            PHI,TETA,XIT(-J),XQT(-J),XUT(-J)
		  ENDIF
               ENDDO !MU
	    ENDDO !PHI	 
	 ENDIF
	       
      ENDIF !ITRPHI=2
     

C* Fin nominale 
C-------------------
      goto 9999
      
C* Cas d'erreur  retour du status -1 au programme appelant
C----------------------------------------------------------
  991 WRITE(6,*) "SOS_TRPHI_OPTION : Error while opening a file"
      IER=-1
      goto 9999
  992 WRITE(6,*) 
     & "SOS_TRPHI_OPTION : Error while reading or writing on a file"
      IER=-1
      goto 9999
      
      
 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(14)
      CLOSE(15)
      CLOSE(21)
            
C* Format 
C----------
15    FORMAT(2(E21.14))        
55    FORMAT(F7.2,2X,3(F13.6,2X))
56    FORMAT(2(F7.2,2X),3(F13.6,2X))

      RETURN  	
      END  	!FIN DE LA PROCEDURE SOS_OPTION_TRPHI








C==============================================================================
C PROCEDURE: SOS_TRPHI
C ==========
C      Cette procedure recombine la decomposition en series de Fourier des
C      parametres de Stokes I, Q et U pour un azimut relatif fixe.
C      
C      Elle permet de recalculer analytiquement la reflexion directe sur mer
C      (glitter) supprimee des calculs des OS pour reduire l'ordre maximale
C      utile pour la decomposition.
C
C
C Description des parametres
C --------------------------
C     FICHOS   (CHARACTER*SOS_LENFIC2)  (E) : Nom du fichier resultat des OS 
C                                   (repertoire + nom fichier + extension)     
C     NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C     RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)   (E) : Tableau des cosinus des angles.
C     TAU  (double)      (E) : Epaisseur optique totale de l'atmosphere 
C                              (eventuellement tronquee).
C     TAUOUT  (double)   (E) : Epaisseur optique au niveau de sortie
C                              des resultats pour le champ montant.
C     PHI  (double)      (E) : Azimut relatif en radians.     
C     IGLI  (I4)      (E) : Option de Glitter :  0 : pas de Glitter.
C                                                1 : calcul du Glitter.
C     N0  (I4)        (E) : Numero de mu de Gauss pour l'angle zenithal solaire. 
C                          (angle complementaire : mus = -rmu(N0)) 
C                          (N0 > 0 => pour en tenir compte).
C     WIND (double)   (E) : Vitesse du vent en m/s.
C     IND_SURF (double) (E) : Indice de refraction relatif surface / air.
C     IFRESNEL  (I4)  (E) : Option de calcul de reflexion de Fresnel 
C                            sur mer plate :
C                        = 0 : pas de reflexion sur mer plate.
C                        = 1 : calcul de la reflexion de Fresnel sur mer plate.
C     IROUJEAN (I4)   (E) : Option de calcul de la reflexion de Roujean
C                        = 0 : pas de reflexion par la BRDF de Roujean.
C                        = 1 : calcul de la reflexion par la BRDF de Roujean.
C     K0_ROUJEAN (double) (E) : Parametre k0 de la BRDF de Roujean.
C     K1_ROUJEAN (double) (E) : Parametre k1 de la BRDF de Roujean.
C     K2_ROUJEAN (double) (E) : Parametre k2 de la BRDF de Roujean.
C     IRONDEAUX (I4)  (E) : Option de calcul de la reflexion de Rondeaux
C                        = 0 : pas de reflexion par la BPDF de Rondeaux.
C                        = 1 : calcul de la reflexion par la BPDF de Rondeaux.
C     IBREON (I4)     (E) : Option de calcul de la reflexion de Breon
C                        = 0 : pas de reflexion par la BPDF de Breon.
C                        = 1 : calcul de la reflexion par la BPDF de Breon.
C     INADAL (I4)     (E) : Option de calcul de la reflexion de Nadal
C                        = 0 : pas de reflexion par la BPDF de Nadal.
C                        = 1 : calcul de la reflexion par la BPDF de Nadal.
C     ALPHA_NADAL (double)  (E) : Parametre alpha du modele de Nadal.         	
C     BETA_NADAL (double)  (E) : Parametre beta du modele de Nadal.          
C     IPOLAR  (I4)    (E) : Option de simulation sans polarisation.
C      			 = 0 : calculs sans polarisation.
C			 = 1 : calculs normaux avec polarisation.
C     XIT(-SOS_OS_NBMU:SOS_OS_NBMU)  (double)   (S) : 
C                              Parametre de Stokes I en fonction de l'angle mu, 
C                              pour l'azimut phi.
C                                 mu < 0 : champ descendant.
C                                 mu > 0 : champ montant.
C     XQT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)   (S) : Parametre de Stokes Q.
C     XUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  (double)   (S) : Parametre de Stokes U.
C     IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier resultat des OS
C       Fichier binaire non formatte, cree par SOS_OS
C       Il contient les parametres de Stokes Q , U et I pour un ordre S
C       du developpement en series de Fourier, en fonction de l'angle mu
C       (valeurs en double precision).
C
C
C Common utilise:
C --------------
C     Aucun
C     
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C       Erreur si le fichier des resultat des OS n'existe pas
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C     
C  Routines utilises par la procedure:
C ------------------------------------
C  Ce programme fait appel aux routines:
C      - SOS_GLITTE
C      - SOS_ANGLE
C      - SOS_REFLEX
C      - SOS_MATRIC
C
C==============================================================================

      SUBROUTINE SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &                     IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &                     IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                     IRONDEAUX,IBREON,
     &                     INADAL,ALPHA_NADAL,BETA_NADAL,
     &                     IPOLAR,XIT,XQT,XUT,IER)
           
      IMPLICIT NONE

C* Definition des variables       
C*----------------------------------------------------------------------     
      DOUBLE PRECISION RMU(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Tableau des mu de Gauss 

      DOUBLE PRECISION TAU	! Epaisseur optique totale de l'atmosphere 
                                ! (eventuellement tronquee).
				
      DOUBLE PRECISION TAUOUT	! Epaisseur optique au niveau de sortie
      				! des resultats pour le champ montant.	
				
      DOUBLE PRECISION PHI	! Azimut relatif en radians.			
      DOUBLE PRECISION WIND	! Vitesse du vent en m/s.		
      DOUBLE PRECISION IND_SURF	! Indice de refraction relatif surface / air.
      DOUBLE PRECISION K0_ROUJEAN ! Parametre k0 du modele de BRDF de Roujean.
      DOUBLE PRECISION K1_ROUJEAN ! Parametre k1 du modele de BRDF de Roujean.
      DOUBLE PRECISION K2_ROUJEAN ! Parametre k2 du modele de BRDF de Roujean.
      DOUBLE PRECISION PHIRJ	! Azimut relatif en radians selon la convention 
      				! du modele de Roujean.	
      DOUBLE PRECISION F	! BRDF de Roujean * cos(angle incidence).
      	
      DOUBLE PRECISION ALPHA_NADAL ! Parametre alpha du modele de BPDF de Nadal.
      DOUBLE PRECISION BETA_NADAL  ! Parametre beta du modele de BPDF de Nadal.	
      DOUBLE PRECISION F21FRESNEL	! Valeur de la fonction F21 de Fresnel.
      DOUBLE PRECISION F21NADAL		! Valeur de la fonction F21 de Nadal.
            				
      DOUBLE PRECISION XIT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  !Parametre de Stokes I 
      					!    en fonction de l'angle mu,
					!    pour l'azimut phi.
					!    mu < 0 : champ descendant.
					!    mu > 0 : champ montant.
      DOUBLE PRECISION XQT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  !Parametre de Stokes Q.
      DOUBLE PRECISION XUT(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  !Parametre de Stokes U.
     					     
      DOUBLE PRECISION I3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Parametre de Stokes I 
      				! pour un ordre IS du developpement en series
				! de Fourier, en fonction de l'angle mu.
      DOUBLE PRECISION Q3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Parametre de Stokes Q.
      DOUBLE PRECISION U3(-SOS_OS_NBMU_MAX:SOS_OS_NBMU_MAX)  ! Parametre de Stokes U.
            	
      DOUBLE PRECISION XPHI	! Produit de PHI par IS.
      DOUBLE PRECISION C0	!Cosinus de l'angle d'elevation solaire rmu(N0)
      DOUBLE PRECISION C1	!Cosinus de l'angle de visee rmu(J) > 0.
      DOUBLE PRECISION S0	!Sinus de l'angle d'elevation solaire rmu(N0)
      DOUBLE PRECISION S1	!Sinus de l'angle de visee rmu(J) > 0.

      DOUBLE PRECISION AT0	! Transmission du TOA au sol pour la direction 
      				! d'incidence.
					
      DOUBLE PRECISION AT(1:SOS_OS_NBMU_MAX)  ! Transmission totale : 
      				! direction d'incidence (du TOA au sol) +
      				! direction d'emergence (du sol au niveau de
				! sortie pour le champ montant). 

      DOUBLE PRECISION ATJ	! Valeur du tableau AT pour la direction J. 
      
      DOUBLE PRECISION SIGMA2	! Carre du parametre Sigma de la theorie de Cox
      				! et Munck.
      
      DOUBLE PRECISION P        ! Probabilite de reflexion vers la 
                                ! direction (mu, phi).
      DOUBLE PRECISION COSDIF	! Cosinus de l'angle de diffusion.
      DOUBLE PRECISION COSKIP	! Cosinus de l'angle de rotation ki'.
      
      DOUBLE PRECISION R11	! Element de la matrice de reflexion de Fresnel
      DOUBLE PRECISION R12      ! pour un repere lie au plan de reflexion.
      DOUBLE PRECISION R33
        
      DOUBLE PRECISION M11	! Element de la matrice de reflexion dans le 
      DOUBLE PRECISION M21      ! repere lie au plan meridien.
      DOUBLE PRECISION M31

      DOUBLE PRECISION COEF_SUN	! Valeur du rapport de l'éclairement solaire au TOA 
      				! (PI) sur l'angle solide du disque solaire.   

      
          
      CHARACTER*SOS_LENFIC2 FICHOS      ! Nom du fichier resultat des OS.

 
      INTEGER*4 IGLI	! Option de Glitter :
      			! 0 : pas de Glitter / 1 : calcul du Glitter.

      INTEGER*4 IFRESNEL ! Indice precisant si les simulations
      			 ! introduisent la reflexion de Fresnel sur mer plate :
			 !    1 => pour un calcul de la reflexion de Fresnel.
			 !    Sinon => pas de prise en compte. 
			      
      INTEGER*4 IROUJEAN ! Indice precisant si les simulations
      			 ! introduisent la reflexion par la BRDF de Roujean :
			 !    1 => pour un calcul avec la BRDF de Roujean.
			 !    Sinon => pas de prise en compte.      
      
      INTEGER*4 IRONDEAUX ! Indice precisant si les simulations
      			  ! introduisent la reflexion par la BPDF de Rondeaux :
			  !    1 => pour un calcul avec la BPDF de Rondeaux.
			  !    Sinon => pas de prise en compte.   
			  
      INTEGER*4 IBREON    ! Indice precisant si les simulations
      			  ! introduisent la reflexion par la BPDF de Breon :
			  !    1 => pour un calcul avec la BPDF de Breon.
			  !    Sinon => pas de prise en compte.   
			  
      INTEGER*4 INADAL	  ! Indice precisant si les simulations
      			  ! introduisent la reflexion par la BPDF de Nadal :
			  !    1 => pour un calcul avec la BPDF de Nadal.
			  !    Sinon => pas de prise en compte.   

      INTEGER*4 IPOLAR	! Option de simulation sans polarisation.
      			! 0 : calculs sans polarisation.
			! 1 : calculs normaux avec polarisation.
								            			     
      INTEGER*4 N0	! Numero de mu de Gauss pour l'angle zenithal solaire. 
			! (angle complementaire : mus = -rmu(N0)) 
			! (N0 > 0 => pour en tenir compte).

      INTEGER*4 IS	! Ordre de la decomposition en series de Fourier.      
      INTEGER*4 IER     ! code d'erreur =0 si pas d'erreur, =-1 sinon  
          
      INTEGER*4 J	! Indice du tableau des angles de Gauss.
      
      INTEGER*4 NBMU    ! nombre d'angles utiles		
	
C* Ouverture des fichiers 
C------------------------------------------------     
C     Ouverture du fichier resultat des OS
      OPEN(UNIT=16,FILE=FICHOS,FORM='UNFORMATTED',
     &      STATUS='OLD',ERR=991) 
     
C* Recombinaison pour l'azimut PHI des parametres de Stokes I, Q et U 
C* a partir de la decomposition en series de Fourier.  
C------------------------------------------------------------------
C*    Initialisation pour IS = 0    

C*      Lecture du fichier qui a ete enregistre avec des 
C*      bornes de tableaux limitees au dimensionnement 
C*      utile -NBMU:NBMU
      READ(16,err=992) (Q3(J),J=-NBMU,NBMU),
     & 		       (U3(J),J=-NBMU,NBMU),
     &		       (I3(J),J=-NBMU,NBMU)
             
      
      DO 2 J= -NBMU, NBMU
         IF(J.EQ.0) GOTO 2
         XQT(J) = Q3(J)
         XUT(J) = U3(J)
         XIT(J) = I3(J)   
    2 CONTINUE

C*   Boucle sur l'ordre IS.    
      IS=1
      
C*      Lecture du fichier qui a ete enregistre avec des 
C*      bornes de tableaux limitees au dimensionnement 
C*      utile -NBMU:NBMU
   10 READ(16,END=9,err=992) (Q3(J),J=-NBMU,NBMU),
     & 		             (U3(J),J=-NBMU,NBMU),
     &		             (I3(J),J=-NBMU,NBMU)
   
      XPHI= IS*PHI
      
      DO 4 J= -NBMU,NBMU
         IF(J.EQ.0) GOTO 4     
         XQT(J) = XQT(J) +2.*Q3(J)*DCOS(XPHI)
         XUT(J) = XUT(J) +2.*U3(J)*DSIN(XPHI)
         XIT(J) = XIT(J) +2.*I3(J)*DCOS(XPHI)
   4  CONTINUE
   
      IS=IS+1
      GOTO 10
      
 9    CLOSE(16)

C* Ajout de la reflexion directe du faisceau solaire sur la mer. 
C-----------------------------------------------------------------     
      IF (IGLI.EQ.1) THEN

C* Transmission du TOA au sol pour la direction d'incidence
C------------------------------------------------------------ 
         C0 = RMU(N0)
         AT0=DEXP(-TAU/C0)

C* Transmission pour la direction d'incidence et pour chaque direction
C* d'emergence J du sol au niveau de sortie des resultats. 
C------------------------------------------------------------------      
         DO 11 J=1, NBMU
            AT(J) = AT0*DEXP(-(TAU-TAUOUT)/RMU(J))
  11     CONTINUE


C* Carre du parametre Sigma de la theorie de Cox et Munck.
C---------------------------------------------------------- 
         SIGMA2 = .003 + .00512*WIND
         
         DO 101 J=1,NBMU
            C1=RMU(J)

C* Calcul de la probabilite P de reflexion dans la direction J 
C* pour la statistique de pentes de vague  
C----------------------------------------------------------    
            CALL SOS_GLITTE(SIGMA2,C0,C1,PHI,P)
      
C* Calcul des angles de diffusion et de rotation. 
C----------------------------------------------------------     
            CALL SOS_ANGLE(C0,C1,PHI,COSKIP,COSDIF)
      
C* Calcul des elements de matrice de la reflexion de Fresnel 
C* pour un repere lie au plan de reflexion.
C----------------------------------------------------------   
            CALL SOS_REFLEX(COSDIF,IND_SURF,R11,R12,R33)
      
C* Calcul des elements de matrice de la reflexion de Fresnel 
C* pour un repere lie au plan meridien (rotation). 
C----------------------------------------------------------        
            CALL SOS_MATRIC(COSKIP,R11,R12,M11,M21,M31)

C* Ajout de la reflexion directe pour un eclairement solaire 
C* normalisee (PI) attenue (ATJ) avec la probabilite P liee 
C* aux pentes de vagues.
C----------------------------------------------------------   
            ATJ =AT(J)
            XIT(J)= XIT(J)+M11*SOS_PI*ATJ*P
	    
	    IF(IPOLAR.EQ.1) THEN
               XQT(J)= XQT(J)+M21*SOS_PI*ATJ*P
               XUT(J)= XUT(J)+M31*SOS_PI*ATJ*P
	    ENDIF
	
  101    CONTINUE
  
      ENDIF	!Test igli=1





C* Ajout de la reflexion du faisceau solaire direct sur mer plate  
C----------------------------------------------------------   

      IF (IFRESNEL.EQ.1) THEN
		
C* Test si le plan d'incidence correspond au plan principal solaire
C* et l'angle zenithal solaire fait parti des angles de sortie
C* (c'est a dire, s'il s'agit d'un angle de Gauss).     
          IF ((DCOS(PHI).EQ.1.D+00).AND.(N0.GT.0)) THEN

C* Transmission du TOA au sol pour la direction d'incidence 
              C0  = RMU(N0)
              AT0 = DEXP(-TAU/C0)

C* Transmission pour la direction d'incidence et de reflexion
C* au niveau de sortie des resultats     
              ATJ = AT0*DEXP(-(TAU-TAUOUT)/C0)	     

C* Coefficients de la matrice de reflexion de Fresnel pour
C* l'incidence solaire

              COSDIF = 1. - 2.*C0*C0
              CALL SOS_REFLEX(COSDIF,IND_SURF,R11,R12,R33) 

C* Valeur de la luminance solaire au TOA (normalisée à Es=PI)
C* == > Ls = Es/angle_solide_disque_solaire
              COEF_SUN = SOS_PI / SOS_SOLAR_DISC_SOLID_ANGLE

C* Reflexion du faisceau solaire direct, atténué	      
	      XIT(N0) = XIT(N0) + R11*COEF_SUN*ATJ 
	      IF(IPOLAR.EQ.1) XQT(N0) = XQT(N0) + R12*COEF_SUN*ATJ  
	      
	  ENDIF
	  
      ENDIF     !Test ifresnel=1





C* Ajout de la reflexion du faisceau solaire direct pour une BRDF
C* du modele de Roujean  
C----------------------------------------------------------   

      IF (IROUJEAN.EQ.1) THEN
 
C* Transmission du TOA au sol pour la direction d'incidence 
          C0 = RMU(N0)
          AT0=DEXP(-TAU/C0)
	 
	  S0 = DSQRT(1.D+00-C0*C0)

C* Boucle sur les directions d'emergence  
	  DO 103 J=1,NBMU
	  
C*     Transmission pour la direction d'incidence et pour la direction
C*     d'emergence J du sol au niveau de sortie des resultats. 
	     C1 = RMU(J)
	     ATJ = AT0*DEXP(-(TAU-TAUOUT)/C1) 
	     
	     S1 = DSQRT(1.D+00-C1*C1)

C*     Calcul de la BRDF de Roujean (luminance normalisee*cos(tetav))	     
	     PHIRJ = SOS_PI - PHI
             CALL SOS_CALC_F_ROUJEAN(K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &	                             C0,S0,C1,S1,PHIRJ,F)   
     
C*     Ajout de la reflexion directe pour un eclairement solaire 
C*     normalisee (PI) attenue (ATJ) de la reflexion de Roujean.
            XIT(J)= XIT(J) + ATJ*F/C1
            
103       CONTINUE

      ENDIF     !Test iroujean=1





C* Ajout de la reflexion du faisceau solaire direct pour une BPDF
C* du modele de Rondeaux ou de Breon
C----------------------------------------------------------   

      IF( (IRONDEAUX.EQ.1).OR.(IBREON.EQ.1) )  THEN

C* Transmission du TOA au sol pour la direction d'incidence 
          C0 = RMU(N0)
          AT0=DEXP(-TAU/C0)
	  
C* Boucle sur les directions d'emergence  
	  DO 105 J=1,NBMU
	  
C*     Transmission pour la direction d'incidence et pour la direction
C*     d'emergence J du sol au niveau de sortie des resultats. 
	     C1 = RMU(J)
	     ATJ = AT0*DEXP(-(TAU-TAUOUT)/C1) 	     

C***** Calcul de la BPDF de Rondeaux ou Breon

C*        Calcul des angles de diffusion et de rotation.     
             CALL SOS_ANGLE(C0,C1,PHI,COSKIP,COSDIF)
      
C*       Calcul des elements de matrice de la reflexion de Fresnel 
C*       pour un repere lie au plan de reflexion.
             CALL SOS_REFLEX(COSDIF,IND_SURF,R11,R12,R33)
      
C*       Calcul des elements de matrice de la reflexion de Fresnel 
C*       pour un repere lie au plan meridien (rotation).    
             CALL SOS_MATRIC(COSKIP,R11,R12,M11,M21,M31)

C*       Ajout de la reflexion directe pour un eclairement solaire 
C*       normalisee (PI) attenue (ATJ) de la reflexion de Rondeaux 
C*       ou de Breon.
	     IF (IRONDEAUX.EQ.1) THEN
	         P = 1. / ( 4.*(1+C1/C0) )	! Modele de Rondeaux
	     ELSE	
	         P = 1. / ( 4.*C1 )		! Modele de Breon
	     ENDIF 	
            
	    XIT(J)= XIT(J)+M11*ATJ*P
	    IF(IPOLAR.EQ.1) THEN
	       XQT(J)= XQT(J)+M21*ATJ*P
               XUT(J)= XUT(J)+M31*ATJ*P
	    ENDIF         
	           
105       CONTINUE


      ENDIF     !Test irondeaux ou ibreon = 1
      





C* Ajout de la reflexion du faisceau solaire direct pour une BPDF
C* du modele de Nadal  
C----------------------------------------------------------   

      IF (INADAL.EQ.1) THEN

C* Transmission du TOA au sol pour la direction d'incidence 
          C0 = RMU(N0)
          AT0=DEXP(-TAU/C0)
	  
C* Boucle sur les directions d'emergence  
	  DO 107 J=1,NBMU
	  
C*     Transmission pour la direction d'incidence et pour la direction
C*     d'emergence J du sol au niveau de sortie des resultats. 
	     C1 = RMU(J)
	     ATJ = AT0*DEXP(-(TAU-TAUOUT)/C1) 	     

C***** Calcul de la BPDF de Nadal

C*        Calcul des angles de diffusion et de rotation.     
             CALL SOS_ANGLE(C0,C1,PHI,COSKIP,COSDIF)
      
C*       Calcul des elements de matrice de la reflexion de Fresnel 
C*       pour un repere lie au plan de reflexion.
             CALL SOS_REFLEX(COSDIF,IND_SURF,R11,R12,R33)
      
C*       Calcul des elements de matrice de la reflexion de Fresnel 
C*       pour un repere lie au plan meridien (rotation).    
             CALL SOS_MATRIC(COSKIP,R11,R12,M11,M21,M31)

C*       Ajout de la reflexion directe pour un eclairement solaire 
C*       normalisee (PI) attenue (ATJ) de la reflexion de Nadal

C*         F21 de Fresnel (convention de signe de Nadal)                
             F21FRESNEL = -R12     
      
C*         Calcul de F21 de Nadal
             F21NADAL = -BETA_NADAL * F21FRESNEL / (C0 + C1)     
	     F21NADAL = ALPHA_NADAL * ( 1.D+00 - DEXP(F21NADAL) )

C*         Calcul de P               
             IF(F21FRESNEL.LT.1.0D-10) THEN	!Developpement limite au voisinage 
      					        !de 0 pour F21Fresnel
                P = ALPHA_NADAL*BETA_NADAL / (C0 + C1)
             ELSE	 
                P = F21NADAL / F21FRESNEL
             ENDIF

               XIT(J)= XIT(J)+M11*ATJ*P
C*         Ajout aux termes de polarisation              
	     IF(IPOLAR.EQ.1) THEN
	       XQT(J)= XQT(J)+M21*ATJ*P
               XUT(J)= XUT(J)+M31*ATJ*P
	    ENDIF         
	    
	                    
107       CONTINUE
	  
      ENDIF     !Test inadal=1

C* Fin nominale 
C-------------------
      goto 9999
      
C* Cas d'erreur : retour du status 1 au programme appelant
C------------------------------------------------------
  991 WRITE(6,*) "SOS_TRPHI : Error while opening a file"
      IER=-1
      goto 9999
  992 WRITE(6,*) 
     &  "SOS_TRPHI : Error while reading or writing on a file"
      IER=-1 
      goto 9999

 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(16)
            
      RETURN  	!FIN DE LA PROCEDURE SOS_TRPHI
      END



C==============================================================================
C PROCEDURE: SOS_GLITTE
C ==========
C      Cette procedure calcule pour une incidence (mup, phip) sur des vagues
C      la probabilite de reflexion vers la direction (mu, phi).
C      L'agitation de la mer depend de la vitesse du vent. 
C
C      Calculs bases sur la theorie de Cox et Munck (1954).
C
C
C Description des parametres
C --------------------------
C     SIG	  (double)  (E)  Carre du parametre Sigma de la theorie de 
C      				 Cox et Munck.
C     C0	  (double)  (E)  Valeur absolue du cosinus de l'angle zenithal
C      				 d'incidence (mup).
C     C1	  (double)  (E)  Valeur absolue du cosinus de l'angle zenithal 
C     			         d'emergence (mu).
C     PHI	  (double)  (E)  Azimut relatif entre la direction d'emergence
C                                et celle d'incidence phi - phip (en radians).
C     P	          (double)  (S)  Probabilite de reflexion vers la direction 
C                               (mu, phi).
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
      SUBROUTINE SOS_GLITTE(SIG,C0,C1,PHI,P)
       
      IMPLICIT NONE
      
C* Definition des variables   
C*-----------------------------------------     
      DOUBLE PRECISION SIG 	!Carre du parametre Sigma de la theorie de Cox
      				!et Munck.

      DOUBLE PRECISION C0	!Valeur absolue du cosinus de l'angle zenithal
      				!d'incidence (mup).
      DOUBLE PRECISION C1	!Valeur absolue du cosinus de l'angle zenithal 
      				!d'emergence (mu).
      DOUBLE PRECISION PHI	!Azimut relatif entre la direction d'emergence
      				!et celle d'incidence phi - phip (en radians).
      
      DOUBLE PRECISION P	!Probabilite de reflexion vers la 
                                !direction (mu, phi).
      
      DOUBLE PRECISION X1
      DOUBLE PRECISION X2
      DOUBLE PRECISION X3
      DOUBLE PRECISION C0N	! Cosinus de l'angle zenithal de la normale 
                                ! a la facette des vagues.
      DOUBLE PRECISION XXX
      
      DOUBLE PRECISION PP

       
C* Traitement                
C--------------
       X1=DSQRT(1-C1**2)-DCOS(PHI)*DSQRT(1-C0**2)
       X2= DSQRT(1-C0**2)*DSIN(PHI)
       X3=C0+C1
       C0N=(X3/(DSQRT(X1**2+X2**2+X3**2)))
       XXX=(-(1-C0N**2)/(SIG*C0N**2))
       
       IF (XXX.LT.-100) THEN
           P=0.
       ELSE
           PP=(1/(SOS_PI*SIG))*DEXP(XXX)
           P=PP/(4*C1*C0N**4)
       ENDIF     
 
       RETURN
       END  	!FIN DE LA PROCEDURE SOS_GLITTE



C==============================================================================
C PROCEDURE: SOS_ANGLE
C ==========
C      Cette procedure calcule :
C         _ le cosinus de l'angle de diffusion entre la direction d'incidence 
C	    solaire et la direction de visee.
C	  _ le cosinus de l'angle de rotation KI'. 
C
C Description des parametres
C --------------------------
C       C0        (double)  (E)	 Valeur absolue du cosinus de l'angle zenithal
C                                solaire.
C       C1	  (double)  (E)	 Valeur absolue du cosinus de l'angle zenithal 
C                                de visee.
C       PHI	  (double)  (S)	 Azimut relatif entre la direction de visee
C                                et celle du Soleil (en radians).
C       COSKIP	  (double)  (S)  Cosinus de l'angle de rotation Ki'.
C       COSDIF	  (double)  (S)  Cosinus de l'angle de diffusion.
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
      SUBROUTINE SOS_ANGLE(C0,C1,PHI,COSKIP,COSDIF)
      
      IMPLICIT NONE

C* Definition des variables   
C*----------------------------------------
      DOUBLE PRECISION C0	!Valeur absolue du cosinus de l'angle zenithal
      				!solaire.
      DOUBLE PRECISION C1	!Valeur absolue du cosinus de l'angle zenithal 
      				!de visee.
      DOUBLE PRECISION PHI	!Azimut relatif entre la direction de visee
      				!et celle du Soleil (en radians).
      
      DOUBLE PRECISION COSKIP	!(S) Cosinus de l'angle de rotation Ki'.
      DOUBLE PRECISION COSDIF	!(S) Cosinus de l'angle de diffusion.
            				 
      DOUBLE PRECISION S	! 1 ou -1.
      DOUBLE PRECISION Z	! 

C* Traitement                
C--------------         
      S=1.
      IF (DSIN(PHI).GT.0.D+00)  S=-1.
      
      COSDIF=-C0*C1+DSQRT(1-C0**2)*DSQRT(1-C1**2)*DCOS(PHI)
      Z=S*(DSQRT(1-COSDIF**2))*(DSQRT(1-C1*C1))
      COSKIP=0.
      
      IF (DABS(Z).GT. SOS_SEUIL_Z)  COSKIP=(C1*COSDIF+C0)/Z
      
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_ANGLE





C==============================================================================
C PROCEDURE: SOS_REFLEX
C ==========
C      Cette procedure calcule les elements R11, R12 et R33 de la matrice
C      de reflexion de Fresnel exprimee dans un repere lie au plan de reflexion.
C
C      La matrice de Fresnel s'ecrit :
C                  | R11(teta)  R12(teta)    0      |
C         F(teta)= | R12(teta)  R11(teta)    0      |
C                  |    0          0      R33(teta) |
C
C         avec :
C           _ teta : angle de reflexion par rapport a la normale au dioptre,
C           _ R11(teta) = 0.5 * (RL*RL + RR*RR),
C           _ R12(teta) = 0.5 * (RL*RL - RR*RR),
C           _ R33(teta) = RL*RR,
C
C         L'indice de refraction relatif du dioptre, entre le milieu de 
C         reflexion et le milieu de transmission, est IND.
C
C         On calcule RL et RR en fonction de IND et teta, par :
C           _ RL : le coefficient de reflexion parallele
C
C               RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C           _ RR : le coefficient de reflexion perpendiculaire
C
C               RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C
C
C Description des parametres
C --------------------------
C       COSDIF	  (double)  (E)  Cosinus de l'angle de diffusion.
C       IND	  (double)  (E)  indice de refraction relatif 
C     				 (eau / air pour une reflexion sur mer).
C       R11	  (double)  (S) Element de la matrice de reflexion de Fresnel.
C       R12	  (double)  (S) idem R11
C       R33	  (double)  (S) idem R11
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
      SUBROUTINE SOS_REFLEX(COSDIF,IND,R11,R12,R33)
      
      IMPLICIT NONE      

C* Definition des variables                
C*---------------------------    
      DOUBLE PRECISION COSDIF	! Cosinus de l'angle de diffusion.
      
      DOUBLE PRECISION IND	!indice de refraction relatif 
      				!(eau / air pour une reflexion sur mer).
      
      DOUBLE PRECISION R11	!Element de la matrice de reflexion de Fresnel.
      DOUBLE PRECISION R12
      DOUBLE PRECISION R33
      
      
      DOUBLE PRECISION IND2	! Carre de l'indice de refraction relatif.
				
      DOUBLE PRECISION COSW	! Cosinus de l'angle de reflexion W
      				! (angle d'incidence par rapport a la 
				!  normale du dioptre)
      
      DOUBLE PRECISION V	! Sinus carre de W.
      DOUBLE PRECISION X	! 
      
      DOUBLE PRECISION RL	! Coefficient de reflexion parallele.
      DOUBLE PRECISION RR	! Coefficient de reflexion perpendiculaire.
      

C* Traitement                
C--------------     
      IND2=IND*IND
      
      COSW=DSQRT(.5*(1-COSDIF))
      V=.5*(1+COSDIF)
      X=DSQRT(IND2-V)     
      RL=(IND2*COSW-X)/(IND2*COSW+X)
      RR=(COSW-X)/(COSW+X)
      R11=(RL**2+RR**2)/2.
      R12=(RL**2-RR**2)/2.
      R33=RR*RL
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_REFLEX




C==============================================================================
C PROCEDURE: SOS_MATRIC
C ==========
C      Cette procedure calcule les elements M11, M21 et M31 de la matrice
C      de reflexion exprimee dans un repere lie au plan meridien (1ere colonne).
C      
C
C Description des parametres
C --------------------------
C       COSKIP	  (double)  (E) Cosinus de l'angle de rotation KI'      
C       R11	  (double)  (E)	Element de la matrice de reflexion de Fresnel
C     				dans le repere lie au plan de reflexion.
C       R12	  (double)  (E)	idem R11
C       M11	  (double)  (S) Element de la matrice de reflexion dans le 
C                               repere lie au plan meridien.
C       M21	  (double)  (S) idem M11
C       M31	  (double)  (S) idem M11
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
      SUBROUTINE SOS_MATRIC(COSKIP,R11,R12,M11,M21,M31)
      
      IMPLICIT NONE

C* Definition des variables   
C*------------------------------------------------   
      DOUBLE PRECISION COSKIP	!Cosinus de l'angle de rotation KI'.
      
      DOUBLE PRECISION R11	!Element de la matrice de reflexion de Fresnel
      				!dans le repere lie au plan de reflexion.
      DOUBLE PRECISION R12
      
      DOUBLE PRECISION M11	!Element de la matrice de reflexion dans le 
      				!repere lie au plan meridien.
      DOUBLE PRECISION M21
      DOUBLE PRECISION M31
   
      DOUBLE PRECISION C2	! Cosinus de 2*KI'.
      DOUBLE PRECISION S2	! Sinus de 2*KI'.
      DOUBLE PRECISION X

C* Traitement                
C--------------               
      X=1.-DABS(COSKIP)
      C2=1.
      S2=0.
      IF (X.GE.SOS_SEUIL_X) THEN
         C2=2.*COSKIP*COSKIP-1.
         S2=2.*COSKIP*DSQRT(1.-COSKIP*COSKIP)
      ENDIF
      
      IF (COSKIP.EQ.0.D+00) R12=0.
      
      M11=R11
      M21=C2*R12
      M31=S2*R12
      
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_MATRIC




