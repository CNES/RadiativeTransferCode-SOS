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
C* FICHIER: SOS_TRPHI.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module génère les fichiers résultats des paramètres de Stokes 
C*       I, Q et U .
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmosphérique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, B.Lafrance
C* DATE: 30/04/2002
C*
C* MOD:VERSION:1.0: Mise à niveau du code des OS du LOA
C*
C* MOD:VERSION:3.0: 
C*    - Correction d'une erreur pour le cas de la réflexion de Fresnel
C*      sur mer plate : calcul de la réflexion du faisceau solaire direct
C*      affecté à la direction de visée J=N0.
C*
C*    - Introduction de la réflexion directe des modèles:
C*           - Roujean.
C*           - Rondeaux.
C*           - Bréon.
C*           - Nadal.
C*
C* MOD:VERSION:3.1:  Amélioration de commentaires.
C*
C* MOD:VERSION:4.0: 10/01/2010 
C*   - Modification de la taille des tableaux de champs de luminance :
C*     ==> Passage de CTE_OS_NBMU à CTE_OS_NBMU_MAX.
C*
C*   - Ajout du paramètre d'entrée : NBMU (nb d'angles effectifs pour les simulations)
C*
C*   - Modification des lignes de lecture des tableaux de luminances 
C*     à l'ordre IS de la décomposition en séries de Fourier 
C*     ==> limitation aux bornes utiles -NBMU:NBMU
C*
C*   - Adaptation de la routine SOS_TRPHI au passage du paramètre NBMU.
C*
C*   - Adaptation de la routine SOS_TRPHI :
C*        * au passage  des paramètres SOS.UserAng.ResFileUp et 
C*          SOS.UserAng.ResFileDown précisant les noms des fichiers
C*          pour l'enregistrement des luminances de sortie uniquement 
C*          pour les angles "utilisateur".
C*        * au passage de la table IND_ANGOUT précisant les angles
C*          "utilisateur" pour lesquels il faut une sortie de luminance.
C*
C* MOD:VERSION:6.0: 06/06/2013
C*
C*   -  Modification de la valeur de Pi, déclarée en constante dans le fichier
C*      (INCTE_PI remplace CTE_PI)
C*
C* MOD:VERSION:6.1: 19/02/2016
C*   - Correction du terme de réflexion directe sur mer plate
C*     ==> Introduction de la luminance solaire directe (normalisée à Es=PI)
C*         MAIS pour la distance Terre-Soleil moyenne (selon la valeur
C*         de CTE_SOLAR_DISC_SOLID_ANGLE définie dans SOS.h.
C*
C* MOD:VERSION:7.0: 07/08/2020
C*    - Introduction de la réflexion directe du modèle de BPDF de Maignan.
C*    
C* MOD:VERSION:7.1: 23/08/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Ajout d'une entête aux fichiers de sortie et ajout des données 
C*      d'angle de diffusion, angle de polarisation, luminance polarisée 
C*      et taux de polarisation.
C*   
C*  MOD:VERSION:8.0: 24/11/2022
C*    - Suppression de la génération des fichiers résultats par la routine SOS_TRPHI_OPTION
C*      en faveur de la fourniture de tableaux de résultats à la routine SOS_PROC permettant un appel du code par Python.
C*      (les fichiers résultats sont générés par SOS_ABS_MAIN).
C*   
C*  MOD:VERSION:8.1: 19/04/2023
C*    - Suppression de la variable IND_ANGOUT: gestion de la restriction aux angles utilisateur rebasculée 
C*      dans SOS_ABS_MAIN suite au passage au binding Python
C*
C*  MOD:VERSION:8.2: 14/11/2023
C*    - Correctif de l'information de dimension de I, Q et U dans les fichiers résultats
C******************************************************************************


C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisées :
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_OS_NBMU_MAX : Nombre maximal d'angles positifs pour la résolution de l'ETR.
C    CTE_THRESHOLD_Q_U_NULL : Valeur seuil sous laquelle  Q ou U sont forcées à zéro.
C    CTE_SEUIL_Z : Valeur seuil pour le calcul des angles de rotation.
C    CTE_SEUIL_X : Valeur seuil pour le calcul des éléments M21 et M31 de la 
C                  matrice de réflexion exprimee dans un repère lié au plan méridien.
C    CTE_SOLAR_DISC_SOLID_ANGLE : Angle solide du disque solaire (sr) pour 
C                                 la distance Terre - Soleil moyenne.
C
C Constantes spécifiques :   
C    INCTE_PI : valeur de PI                           
C    INCTE_VALEUR_INDEF : Valeur par défaut pour les grandeurs indéfinies
C                         (angle de polarisation si polarisation nulle,
C                          taux de polarisation si intensité nulle).
C
C Procédures:
C      - SOS_TRPHI_OPTION
C      - SOS_TRPHI
C      - SOS_POLAR
C      - SOS_GLITTE
C      - SOS_ANGLE
C      - SOS_REFLE
C      - SOS_MATRIC
C      - SOS_OUTPUT_HEADER
C      - SOS_OUTPUT_HEADER_POLAR_DIAG
C
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define INCTE_VALEUR_INDEF -999.


C============================================================================
C PROCEDURE: SOS_TRPHI_OPTION
C ==========
C
C      Cette procédure recombine les séries de Fourier des paramètres I, Q et U
C      pour une valeur d'azimut relatif ou un ensemble de valeurs d'azimuts relatifs,
C      pour le champ montant et descendant.
C
C      L'état de la polarisation est calculé : taux de polarisation, angle de polarisation.
C      
C      Si ITRPHI = 1, on obtient le plan de visées pour PHIos 
C                     (azimut relatif par rapport a l'azimut solaire)
C         Convention : 
C             L'angle de visée est positif pour le demi plan PHI = PHIos.
C	       L'angle de visée est negatif pour le demi plan PHI = PHIos + 180.
C
C
C      Si ITRPHI = 2, on obtient un diagramme polaire en fonction de l'angle 
C      de visée (teta) et de l'angle azimutal (phi) :
C
C         X =  phi    ,   Y = teta    ,   I  ,  Q  ,  U
C
C	   L'angle de visée croit du centre au bord de l'image.
C
C         L'utilisateur choisit le pas d'échantillonnage en azimut (PAS_PHI).
C
C
C      Pour le cas mer plate, le code ajoute la réflexion de la luminance solaire
C      directe, descendante du TOA à la surface, atténuée par l'atmosphère.
C      Réflexion symétrique à la direction d'incidence solaire par rapport au zénith.
C      NB : La luminance solaire est calculée pour un éclairement solaire au TOA 
C           valant PI, pour une distance Terre-Soleil moyenne,  
C           la constante CTE_SOLAR_DISC_SOLID_ANGLE donnant l'angle solide du disque
C           solaire apparent pour la distance Terre-Soleil moyenne. 
C
C
C Description des paramètres
C --------------------------
C     NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C
C     RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (E) Cosinus des angles.
C     GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE)  : (E) Poids des angles.
C
C     FICHOS   (CHARACTER*CTE_LENFIC2)  (E) :  Nom du fichier résultat des OS 
C                                     (répertoire + nom fichier + extension)
C     TAU (double)    (E) : Epaisseur optique totale de l'atmosphère 
C                          (éventuellement tronquée).
C     TAUOUT (double) (E) : Epaisseur optique au niveau de sortie des résultats pour le
C                           champ montant.
C     ZOUT (DOUBLE)   (E) : Altitude de sortie du champ de rayonnement 
C			         -1 : 300 km pour champ montant, 0 km pour champ descendant.
C                             ZOUT entre 0 et 300 km : altitude pour le champ montant et descendant.  
C     IGLI  (I4)      (E) : Option de Glitter :  0 : pas de Glitter.
C                                                1 : calcul du Glitter.
C     N0  (I4)        (E) : Numéro de mu de Gauss pour l'angle zénithal solaire. 
C                          (angle complémentaire : mus = -rmu(N0)) 
C                          (N0 > 0 => pour en tenir compte).
C     WIND (double)   (E) : Vitesse du vent en m/s.
C     IND_SURF (double) (E) : Indice de réfraction relatif surface / air.
C     IFRESNEL  (I4)  (E) : Option de calcul de réflexion de Fresnel sur mer plate :
C                             0 : pas de réflexion sur mer plate.
C                             1 : calcul de la réflexion de Fresnel sur mer plate.
C     IROUJEAN (I4)   (E) : Option de calcul de la réflexion de Roujean
C                             0 : pas de réflexion par la BRDF de Roujean.
C                             1 : calcul de la réflexion par la BRDF de Roujean.
C     K0_ROUJEAN (double) (E) : Paramètre k0 de la BRDF de Roujean.
C     K1_ROUJEAN (double) (E) : Paramètre k1 de la BRDF de Roujean.
C     K2_ROUJEAN (double) (E) : Paramètre k2 de la BRDF de Roujean.
C     IRONDEAUX (I4)  (E) : Option de calcul de la réflexion de Rondeaux
C                             0 : pas de réflexion par la BPDF de Rondeaux.
C                             1 : calcul de la réflexion par la BPDF de Rondeaux.
C     IBREON (I4)     (E) : Option de calcul de la réflexion de Bréon
C                             0 : pas de réflexion par la BPDF de Bréon.
C                             1 : calcul de la réflexion par la BPDF de Bréon.
C     INADAL (I4)     (E) : Option de calcul de la réflexion de Nadal
C                             0 : pas de réflexion par la BPDF de Nadal.
C                             1 : calcul de la réflexion par la BPDF de Nadal.
C     ALPHA_NADAL (double)  (E) : Paramètre alpha du modèle de Nadal.         	
C     BETA_NADAL (double)  (E) : Paramètre beta du modèle de Nadal.  
C     IMAIGNAN (I4)     (E) : Option de calcul de la réflexion de Maignan
C                             0 : pas de réflexion par la BPDF de Maignan.
C                             1 : calcul de la réflexion par la BPDF de Maignan. 
C     COEF_C_MAIGNAN  (double)  (E) : Paramètre C*exp(-NDVI) du modèle de Maignan.        
C     ITRPHI (I4)     (E) : Indice de type de sortie.
C     PHIOS (double)  (E) : Azimut relatif en degrés (pour ITRPHI=1).
C     PAS_PHI (I4)    (E) : Pas sur l'azimut en degrés (pour ITRPHI=2)
C                           Valeur entière.
C     IPOLAR  (I4)   (E) : Option de simulation sans polarisation.
C      			   0 : calculs sans polarisation.
C			          1 : calculs normaux avec polarisation.		 
C
C     PHI_FIN(0:360)	 (double) : (S) Table des azimuts relatifs du champ de sortie
C     THETA_FIN(0:CTE_OS_NBMU_MAX) (double) : (S) Table des angles zénithaux de visée du champ de sortie
C     SCA_UP_FIN(0:360,0:CTE_OS_NBMU_MAX) (double) : (S) Table des angles de diffusion pour le champ montant de sortie
C                                                        en fonction de Phi et Thetav
C     I_UP_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)   : (S) Table des valeurs I pour le champ montant de sortie  
C     Q_UP_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)   : (S) Table des valeurs Q pour le champ montant de sortie 
C     U_UP_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)   : (S) Table des valeurs U pour le champ montant de sortie 
C     POL_ANG_UP_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)  : (S) Table des angles de polarisation pour le champ montant de sortie
C     POL_RATE_UP_FIN(0:360,0:CTE_OS_NBMU_MAX) (double) : (S) Table des taux de polarisation pour le champ montant de sortie  
C     L_POL_UP_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)    : (S) Table des valeurs de luminance polarisée pour le champ montant de sortie 
C
C     SCA_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) (double) : (S) Table des angles de diffusion pour le champ descendant de sortie
C                                                          en fonction de Phi et Thetav
C     I_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)   : (S) Table des valeurs I pour le champ descendant de sortie  
C     Q_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)   : (S) Table des valeurs Q pour le champ descendant de sortie 
C     U_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)   : (S) Table des valeurs U pour le champ descendant de sortie 
C     POL_ANG_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)  : (S) Table des angles de polarisation pour le champ descendant de sortie
C     POL_RATE_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) (double) : (S) Table des taux de polarisation pour le champ descendant de sortie  
C     L_POL_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) (double)    : (S) Table des valeurs de luminance polarisée pour le champ descendant de sortie 
C
C     IER  (I4)       (S) : Code d'erreur =0 si pas d'erreur, =-1 sinon 	
C
C
C Variable d'environnement
C ------------------------
C   Aucune
C
C
C Description des fichiers utilisés
C ---------------------------------
C
C   -->Contenu du fichier résultat des OS
C       Fichier binaire non formaté, cree par SOS_OS
C       Il contient les paramètres de Stokes Q , U et I pour un ordre S
C       du développement en séries de Fourier, en fonction de l'angle mu
C       (valeurs en double précision).
C
C Common utilisé:
C --------------
C     Aucun   
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C       Erreur si le fichier résultat des OS n'existe pas
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C     
C  Routines utilisées par la procédure:
C ------------------------------------
C  Ce programme fait appel aux routines:
C      - SOS_TRPHI
C      - SOS_POLAR
C
C==============================================================================
  
      SUBROUTINE SOS_TRPHI_OPTION(NBMU,RMU,GA,
     &     FICHOS,TAU,TAUOUT,ZOUT,
     &     IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &     IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &     IRONDEAUX,IBREON,
     &     INADAL,ALPHA_NADAL,BETA_NADAL,
     &     IMAIGNAN,COEF_C_MAIGNAN,
     &     ITRPHI,PHIOS,PAS_PHI,IPOLAR,
     &     PHI_FIN, THETA_FIN, 
     &     SCA_UP_FIN, I_UP_FIN, Q_UP_FIN, U_UP_FIN,
     &     POL_ANG_UP_FIN, POL_RATE_UP_FIN, L_POL_UP_FIN,
     &     SCA_DOWN_FIN, I_DOWN_FIN, Q_DOWN_FIN, U_DOWN_FIN,
     &     POL_ANG_DOWN_FIN, POL_RATE_DOWN_FIN, L_POL_DOWN_FIN,
     &     IER)

   	 	
      IMPLICIT NONE

C* Définition des variables                 
C*-----------------------------------------
      DOUBLE PRECISION TAU ! Epaisseur optique totale de l'atmosphère (éventuellement tronquée).

      DOUBLE PRECISION TAUOUT ! Epaisseur optique au niveau de sortie des résultats pour le champ montant.

      DOUBLE PRECISION ZOUT ! Altitude de sortie du champ de rayonnement 
                            ! Si -1 : sortie à 300 km pour champ montant, 0 km pour champ descendant.
                            ! Si entre 0 et 300 km : sortie à ZOUT pour le champ montant et descendant.  
						
      DOUBLE PRECISION WIND ! Vitesse du vent en m/s. 			
      DOUBLE PRECISION IND_SURF ! Indice de réfraction relatif surface / air.
      DOUBLE PRECISION K0_ROUJEAN ! Paramètre k0 du modèle de BRDF de Roujean.
      DOUBLE PRECISION K1_ROUJEAN ! Paramètre k1 du modèle de BRDF de Roujean.
      DOUBLE PRECISION K2_ROUJEAN ! Paramètre k2 du modèle de BRDF de Roujean.
      DOUBLE PRECISION ALPHA_NADAL ! Paramètre alpha du modèle de BPDF de Nadal.
      DOUBLE PRECISION BETA_NADAL ! Paramètre beta du modèle de BPDF de Nadal.
      DOUBLE PRECISION COEF_C_MAIGNAN ! Paramètre C*exp(-NDVI) du modèle de Maignan. 
      DOUBLE PRECISION PHIOS ! Azimut relatif en degrés (pour ITRPHI=1).
      
      DOUBLE PRECISION TETA ! Angle de visée (en degrés).
      DOUBLE PRECISION PHI ! Azimut relatif (en radians).
      DOUBLE PRECISION PHI1, PHI2 ! Azimut relatif (en radians).
      
      DOUBLE PRECISION XIT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Paramètre de Stokes I 
                                                             ! en fonction de mu pour l'azimut phi.
                                                             !    mu < 0 : champ descendant.
                                                             !    mu > 0 : champ montant.
      DOUBLE PRECISION XQT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Paramètre de Stokes Q.
      DOUBLE PRECISION XUT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Paramètre de Stokes U.
   
      DOUBLE PRECISION XI !  Valeur du paramètre de Stokes I pour un angle
      DOUBLE PRECISION XQ
      DOUBLE PRECISION XU   	  
      
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) !Tableau des angles de Gauss (et utilisateur). 
      DOUBLE PRECISION GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  !Poids des angles de Gauss (et utilisateur).

      DOUBLE PRECISION ANGDIFF(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) !Angles de diffusion


      DOUBLE PRECISION XAN ! Angle de la polarisation (degrés)	
      DOUBLE PRECISION LPOL ! Intensité polarisée normalisée à l'éclairement solaire au TOA 	       
      DOUBLE PRECISION TPOL ! Taux de polarisation en %.

      CHARACTER*CTE_LENFIC2 FICHOS !Nom du fichier résultat des OS	

      INTEGER*4 IGLI ! Option de Glitter : 0 : pas de Glitter.
                     !                     1 : calcul du Glitter.

      INTEGER*4 IFRESNEL ! Indice précisant si les simulations introduisent la réflexion de Fresnel sur mer plate :
                         !    1 => pour un calcul de la réflexion de Fresnel.
                         !    Sinon => pas de prise en compte.      
	
      INTEGER*4 IROUJEAN ! Indice précisant si les simulations introduisent la réflexion par la BRDF de Roujean :
                         !    1 => pour un calcul avec la BRDF de Roujean.
                         !    Sinon => pas de prise en compte.      
      
      INTEGER*4 IRONDEAUX ! Indice précisant si les simulations introduisent la réflexion par la BPDF de Rondeaux :
                          !    1 => pour un calcul avec la BPDF de Rondeaux.
                          !    Sinon => pas de prise en compte.   
			  
      INTEGER*4 IBREON ! Indice précisant si les simulations introduisent la réflexion par la BPDF de Bréon :
                       !    1 => pour un calcul avec la BPDF de Bréon.
                       !    Sinon => pas de prise en compte.   
			  
      INTEGER*4 INADAL ! Indice précisant si les simulations introduisent la réflexion par la BPDF de Nadal :
                       !    1 => pour un calcul avec la BPDF de Nadal.
                       !    Sinon => pas de prise en compte.   
			  
      INTEGER*4 IMAIGNAN ! Indice précisant si les simulations introduisent la réflexion par la BPDF de Maignan :
                         !    1 => pour un calcul avec la BPDF de Maignan.
                         !    Sinon => pas de prise en compte.   
			             
      INTEGER*4 N0 ! Numéro de mu de Gauss pour l'angle zénithal solaire. 
                   !  (angle complémentaire : mus = -rmu(N0)) 
                   !   (N0 > 0 => pour en tenir compte).

      INTEGER*4 ITRPHI ! Indice de type de sortie.
                       !   1 =>  Plan de visée pour une valeur d'azimut.		           
                       !   2 =>  Diagramme polaire en fonction de la direction d'azimut et de visée.
			      
      INTEGER*4 PAS_PHI ! Pas sur l'azimut en degrés (pour ITRPHI=2, valeur entière)
      
      INTEGER*4 IPOLAR ! Option de simulation sans polarisation.
                       ! 0 : calculs sans polarisation.
                       ! 1 : calculs normaux avec polarisation.
			      
      INTEGER*4 IER ! code d'erreur =0 si pas d'erreur, =-1 sinon
      INTEGER*4 J ! Indice du tableau des angles de Gauss.
      INTEGER*4 IPHI ! Indice de boucle sur l'azimut relatif.
      INTEGER*4 NBMU ! nombre d'angles utiles					     
			
      INTEGER*4 I  ! indice de ligne des tableaux SCA_FIN, I_FIN, Q_FIN, U_FIN, POL_ANG_FIN, POL_RATE_FIN, L_POL_FIN
      INTEGER*4 JJ ! indice de colonne des tableaux THETA_FIN, SCA_FIN, I_FIN, Q_FIN, U_FIN, POL_ANG_FIN, POL_RATE_FIN, L_POL_FIN
      INTEGER*4 IP ! indice de ligne des tableaux SCA_FIN, I_FIN, Q_FIN, U_FIN, POL_ANG_FIN, POL_RATE_FIN, L_POL_FIN
       
   
      DOUBLE PRECISION PHI_FIN(0:360)		   !Table des azimuts relatifs du champ de sortie
      DOUBLE PRECISION THETA_FIN(0:CTE_OS_NBMU_MAX) !Table des angles zénithaux de visée du champ de sortie
            
      DOUBLE PRECISION SCA_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table des angles de diffusion pour le champ montant de sortie
      DOUBLE PRECISION SCA_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) !Table des angles de diffusion pour le champ descendant de sortie     

      DOUBLE PRECISION I_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)     !Table des valeurs I pour le champ montant de sortie    
      DOUBLE PRECISION Q_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)     !Table des valeurs Q pour le champ montant de sortie    
      DOUBLE PRECISION U_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)     !Table des valeurs U pour le champ montant de sortie  
      DOUBLE PRECISION I_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table des valeurs I pour le champ descendant de sortie      
      DOUBLE PRECISION Q_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table des valeurs Q pour le champ descendant de sortie          
      DOUBLE PRECISION U_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table des valeurs U pour le champ descendant de sortie   
      
      DOUBLE PRECISION POL_ANG_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table des angles de polarisation pour le champ montant de sortie
      DOUBLE PRECISION POL_ANG_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) !Table des angles de polarisation pour le champ descendant de sortie
     
      DOUBLE PRECISION POL_RATE_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table des taux de polarisation pour le champ montant de sortie 
      DOUBLE PRECISION POL_RATE_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) !Table des taux de polarisation pour le champ descendant de sortie        
  
      DOUBLE PRECISION L_POL_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table des valeurs de luminance polarisée pour le champ montant de sortie          
      DOUBLE PRECISION L_POL_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) !Table des valeurs de luminance polarisée pour le champ descendant de sortie    

      	
C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0
     
      
C* Sortie d'un plan azimutal
C------------------------------------------------------
      IF(ITRPHI.EQ.1) THEN

         !Demi plan PHIos + 180 degrés : Angles de visée négatifs par convention
         !----------------------------------------------------------------------
         PHI=INCTE_PI+PHIOS*INCTE_PI/180.D+00	
         CALL SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &        IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &        IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &        IRONDEAUX,IBREON,
     &        INADAL,ALPHA_NADAL,BETA_NADAL,
     &        IMAIGNAN,COEF_C_MAIGNAN,
     &        IPOLAR,XIT,XQT,XUT,ANGDIFF,IER)
         IF (IER.NE.0) GOTO 9999   ! cas d'erreur

         PHI_FIN(0) = PHIOS + 180

         DO J=NBMU,1,-1
            
            TETA=DACOS(RMU(J))*180.D+00/INCTE_PI
            
            !Enregistrement pour le champ montant
            !------------------------------------
            XI = XIT(J)
            XQ = XQT(J)
            XU = XUT(J)
            CALL SOS_POLAR(XI,XQ,XU,XAN,TPOL,LPOL)

C            JJ=NBMU-J
	    JJ=J-1
	    
            THETA_FIN(JJ) = -TETA
            SCA_UP_FIN(0, JJ) = ANGDIFF(J)

            I_UP_FIN(0, JJ) = XI
            Q_UP_FIN(0, JJ) = XQ
            U_UP_FIN(0, JJ) = XU
            POL_ANG_UP_FIN(0, JJ) = XAN
            POL_RATE_UP_FIN(0, JJ) = TPOL
            L_POL_UP_FIN(0, JJ) = LPOL
           
            !Enregistrement pour le champ descendant
            !---------------------------------------
            XI = XIT(-J)
            XQ = XQT(-J)
            XU = XUT(-J)
            CALL SOS_POLAR(XI,XQ,XU,XAN,TPOL,LPOL)

            THETA_FIN(JJ) = -TETA
            SCA_DOWN_FIN(0, JJ) = ANGDIFF(-J)
            
            I_DOWN_FIN(0, JJ) = XI
            Q_DOWN_FIN(0, JJ) = XQ
            U_DOWN_FIN(0, JJ) = XU
            POL_ANG_DOWN_FIN(0, JJ) = XAN
            POL_RATE_DOWN_FIN(0, JJ) = TPOL
            L_POL_DOWN_FIN(0, JJ) = LPOL
	     	    
         ENDDO
 

	    
         !Demi plan PHIos : Angles de visée positifs par convention
         !---------------------------------------------------------
         PHI=PHIOS*INCTE_PI/180.D+00
         CALL SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &        IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &        IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &        IRONDEAUX,IBREON,
     &        INADAL,ALPHA_NADAL,BETA_NADAL,
     &        IMAIGNAN,COEF_C_MAIGNAN,
     &        IPOLAR,XIT,XQT,XUT,ANGDIFF,IER)
         IF (IER.NE.0) GOTO 9999   ! cas d'erreur

         PHI_FIN(0) = PHIOS 
	 
         DO J=1,NBMU
            
            TETA=DACOS(RMU(J))*180.D+00/INCTE_PI

            !Enregistrement pour le champ montant
            !------------------------------------
            XI = XIT(J)
            XQ = XQT(J)
            XU = XUT(J)
            CALL SOS_POLAR(XI,XQ,XU,XAN,TPOL,LPOL)

            JJ=J-1
	    
            THETA_FIN(JJ) = TETA
            SCA_UP_FIN(1, JJ) = ANGDIFF(J)
            
            I_UP_FIN(1, JJ) = XI
            Q_UP_FIN(1, JJ) = XQ
            U_UP_FIN(1, JJ) = XU
            POL_ANG_UP_FIN(1, JJ) = XAN
            POL_RATE_UP_FIN(1, JJ) = TPOL
            L_POL_UP_FIN(1, JJ) = LPOL
        
            !Enregistrement pour le champ descendant
            !---------------------------------------
            XI = XIT(-J)
            XQ = XQT(-J)
            XU = XUT(-J)
            CALL SOS_POLAR(XI,XQ,XU,XAN,TPOL,LPOL)

            THETA_FIN(JJ) = TETA
            SCA_DOWN_FIN(1, JJ) = ANGDIFF(-J)
            
            I_DOWN_FIN(1, JJ) = XI
            Q_DOWN_FIN(1, JJ) = XQ
            U_DOWN_FIN(1, JJ) = XU
            POL_ANG_DOWN_FIN(1, JJ) = XAN
            POL_RATE_DOWN_FIN(1, JJ) = TPOL
            L_POL_DOWN_FIN(1, JJ) = LPOL	    
         
         ENDDO	    

      ENDIF  !ITRPHI=1
      


C* Sortie d'un diagramme polaire
C---------------------------------
      IF(ITRPHI.EQ.2) THEN        

         IP=0
	 
         DO IPHI=0,360,PAS_PHI
            
	    
            PHI = INCTE_PI * IPHI / 180.D+00
            CALL SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &           IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &           IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &           IRONDEAUX,IBREON,
     &           INADAL,ALPHA_NADAL,BETA_NADAL,
     &           IMAIGNAN,COEF_C_MAIGNAN,
     &           IPOLAR,XIT,XQT,XUT,ANGDIFF,IER)
            IF (IER.NE.0) GOTO 9999 ! cas d'erreur
            
            PHI = dble(IPHI)
            DO J=1,NBMU
               TETA=DACOS(RMU(J))*180.D+00/INCTE_PI
               !Enregistrement pour le champ montant
               !------------------------------------
               XI = XIT(J)
               XQ = XQT(J)
               XU = XUT(J)               
               CALL SOS_POLAR(XI,XQ,XU,XAN,TPOL,LPOL)

               PHI_FIN(IP) = PHI
               THETA_FIN(J-1) = TETA

               SCA_UP_FIN(IP, J-1) = ANGDIFF(J)
               
               I_UP_FIN(IP, J-1) = XI
               Q_UP_FIN(IP, J-1) = XQ
               U_UP_FIN(IP, J-1) = XU
               POL_ANG_UP_FIN(IP, J-1) = XAN
               POL_RATE_UP_FIN(IP, J-1) = TPOL
               L_POL_UP_FIN(IP, J-1) = LPOL

               !Cas nominal : champ descendant pour tous les angles  
               XI = XIT(-J)
               XQ = XQT(-J)
               XU = XUT(-J)
               CALL SOS_POLAR(XI,XQ,XU,XAN,TPOL,LPOL) 	 

               SCA_DOWN_FIN(IP, J-1) = ANGDIFF(-J)
               
               I_DOWN_FIN(IP, J-1) = XI
               Q_DOWN_FIN(IP, J-1) = XQ
               U_DOWN_FIN(IP, J-1) = XU
               POL_ANG_DOWN_FIN(IP, J-1) = XAN
               POL_RATE_DOWN_FIN(IP, J-1) = TPOL
               L_POL_DOWN_FIN(IP, J-1) = LPOL
      
             ENDDO              !MU
	     
	     !Incrementation de l'indice PHI
	     IP=IP+1
	     
	  ENDDO                 !PHI
          
      ENDIF !ITRPHI=2
      
C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur  retour du status -1 au programme appelant
C----------------------------------------------------------
  991 WRITE(6,*) "SOS_TRPHI_OPTION : Error while opening a file"
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 
     & "SOS_TRPHI_OPTION : Error while reading or writing on a file"
      IER=-1
      GOTO 9999
      
      
 9999 CONTINUE


      RETURN  	
      END  	!FIN DE LA PROCEDURE SOS_OPTION_TRPHI




C==============================================================================
C PROCEDURE: SOS_TRPHI
C ==========
C      Cette procédure recombine la décomposition en séries de Fourier des
C      paramètres de Stokes I, Q et U pour un azimut relatif fixe.
C      
C      Elle permet de recalculer analytiquement la réflexion directe sur mer
C      (glitter) supprimee des calculs des OS pour reduire l'ordre maximale
C      utile pour la décomposition.
C
C
C Description des paramètres
C --------------------------
C     FICHOS   (CHARACTER*CTE_LENFIC2)  (E) : Nom du fichier résultat des OS 
C                                   (répertoire + nom fichier + extension)     
C     NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C     RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)   (E) : Tableau des cosinus des angles.
C     TAU  (double)      (E) : Epaisseur optique totale de l'atmosphère 
C                              (éventuellement tronquée).
C     TAUOUT  (double)   (E) : Epaisseur optique au niveau de sortie
C                              des résultats pour le champ montant.
C     PHI  (double)      (E) : Azimut relatif en radians.     
C     IGLI  (I4)      (E) : Option de Glitter :  0 : pas de Glitter.
C                                                1 : calcul du Glitter.
C     N0  (I4)        (E) : Numéro de mu de Gauss pour l'angle zénithal solaire. 
C                          (angle complémentaire : mus = -rmu(N0)) 
C                          (N0 > 0 => pour en tenir compte).
C     WIND (double)   (E) : Vitesse du vent en m/s.
C     IND_SURF (double) (E) : Indice de réfraction relatif surface / air.
C     IFRESNEL  (I4)  (E) : Option de calcul de réflexion de Fresnel 
C                           sur mer plate :
C                             0 : pas de réflexion sur mer plate.
C                             1 : calcul de la réflexion de Fresnel sur mer plate.
C     IROUJEAN (I4)   (E) : Option de calcul de la réflexion de Roujean
C                             0 : pas de réflexion par la BRDF de Roujean.
C                             1 : calcul de la réflexion par la BRDF de Roujean.
C     K0_ROUJEAN (double) (E) : Paramètre k0 de la BRDF de Roujean.
C     K1_ROUJEAN (double) (E) : Paramètre k1 de la BRDF de Roujean.
C     K2_ROUJEAN (double) (E) : Paramètre k2 de la BRDF de Roujean.
C     IRONDEAUX (I4)  (E) : Option de calcul de la réflexion de Rondeaux
C                             0 : pas de réflexion par la BPDF de Rondeaux.
C                             1 : calcul de la réflexion par la BPDF de Rondeaux.
C     IBREON (I4)     (E) : Option de calcul de la réflexion de Bréon
C                             0 : pas de réflexion par la BPDF de Bréon.
C                             1 : calcul de la réflexion par la BPDF de Bréon.
C     INADAL (I4)     (E) : Option de calcul de la réflexion de Nadal
C                             0 : pas de réflexion par la BPDF de Nadal.
C                             1 : calcul de la réflexion par la BPDF de Nadal.
C     ALPHA_NADAL (double) (E) : Paramètre alpha du modèle de Nadal.         	
C     BETA_NADAL (double)  (E) : Paramètre beta du modèle de Nadal.      
C     IMAIGNAN (I4) : (E) Indice précisant si les simulations
C      	            introduisent le modèle de réflexion de Maignan :
C			      1 : pour un calcul de BPDF Maignan (il faut utiliser N0 > 0).
C			      Sinon : pas de prise en compte. 
C     COEF_C_MAIGNAN  (double)  (E) : Paramètre C*exp(-NDVI) du modèle de Maignan.      
C     IPOLAR  (I4)    (E) : Option de simulation sans polarisation.
C      			   0 : calculs sans polarisation.
C			          1 : calculs normaux avec polarisation.
C     XIT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)   (S) : 
C                              Paramètre de Stokes I en fonction de l'angle mu, 
C                              pour l'azimut phi.
C                                 mu < 0 : champ descendant.
C                                 mu > 0 : champ montant.
C     XQT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)   (S) : Paramètre de Stokes Q.
C     XUT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)   (S) : Paramètre de Stokes U.
C 
C     ANGDIFF(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)   (S) : Angle de diffusion pour les 
C                                                                 angles de visée.
C
C     IER  (I4)   (S) : code d'erreur =0 si pas d'erreur, =-1 sinon 
C
C 
C Description des fichiers utilisés
C ---------------------------------
C   -->Contenu du fichier résultat des OS
C       Fichier binaire non formaté, cree par SOS_OS
C       Il contient les paramètres de Stokes Q , U et I pour un ordre S
C       du développement en séries de Fourier, en fonction de l'angle mu
C       (valeurs en double précision).
C
C
C Common utilisé:
C --------------
C     Aucun
C     
C
C Cas d'erreur :
C ------------
C     -Ouverture d'un fichier
C       Erreur si le fichier des résultat des OS n'existe pas
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C     
C  Routines utilisées par la procédure:
C ------------------------------------
C  Ce programme fait appel aux routines:
C      - SOS_GLITTE
C      - SOS_ANGLE
C      - SOS_REFLEX
C      - SOS_MATRIC
C      - SOS_CALC_F_ROUJEAN
C      - SOS_CALC_MAIGNAN
C
C==============================================================================

      SUBROUTINE SOS_TRPHI(FICHOS,NBMU,RMU,TAU,TAUOUT,PHI,
     &                     IGLI,N0,WIND,IND_SURF,IFRESNEL,
     &                     IROUJEAN,K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                     IRONDEAUX,IBREON,
     &                     INADAL,ALPHA_NADAL,BETA_NADAL,
     &                     IMAIGNAN,COEF_C_MAIGNAN,
     &                     IPOLAR,XIT,XQT,XUT,ANGDIFF,IER)
           
      IMPLICIT NONE

C* Définition des variables       
C*----------------------------------------------------------------------     
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Tableau des mu de Gauss 

      DOUBLE PRECISION ANGDIFF(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) !Angles de diffusion

      DOUBLE PRECISION TAU ! Epaisseur optique totale de l'atmosphère (éventuellement tronquée).
				
      DOUBLE PRECISION TAUOUT ! Epaisseur optique au niveau de sortie des résultats pour le champ montant.	
				
      DOUBLE PRECISION PHI ! Azimut relatif en radians.			
      DOUBLE PRECISION WIND ! Vitesse du vent en m/s.		
      DOUBLE PRECISION IND_SURF ! Indice de réfraction relatif surface / air.
      DOUBLE PRECISION K0_ROUJEAN ! Paramètre k0 du modèle de BRDF de Roujean.
      DOUBLE PRECISION K1_ROUJEAN ! Paramètre k1 du modèle de BRDF de Roujean.
      DOUBLE PRECISION K2_ROUJEAN ! Paramètre k2 du modèle de BRDF de Roujean.
      DOUBLE PRECISION PHIRJ ! Azimut relatif en radians selon la convention du modèle de Roujean.	
      DOUBLE PRECISION F ! BRDF de Roujean * cos(angle incidence).
      	
      DOUBLE PRECISION ALPHA_NADAL ! Paramètre alpha du modèle de BPDF de Nadal.
      DOUBLE PRECISION BETA_NADAL ! Paramètre beta du modèle de BPDF de Nadal.	
      DOUBLE PRECISION F21FRESNEL ! Valeur de la fonction F21 de Fresnel.
      DOUBLE PRECISION F21NADAL   ! Valeur de la fonction F21 de Nadal.

      DOUBLE PRECISION COEF_C_MAIGNAN ! Paramètre C*exp(-NDVI) du modèle de Maignan.  
            				
      DOUBLE PRECISION XIT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) !Paramètre de Stokes I en fonction de l'angle mu, pour l'azimut phi.
                                                             !    mu < 0 : champ descendant.
                                                             !    mu > 0 : champ montant.
      DOUBLE PRECISION XQT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) !Paramètre de Stokes Q.
      DOUBLE PRECISION XUT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) !Paramètre de Stokes U.
     					     
      DOUBLE PRECISION I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Paramètre de Stokes I en fonction de l'angle mu,
                                                            ! pour un ordre IS du développement en séries de Fourier sur l'azimut
      DOUBLE PRECISION Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Paramètre de Stokes Q.
      DOUBLE PRECISION U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Paramètre de Stokes U.
            	
      DOUBLE PRECISION XPHI ! Produit de PHI par IS.
      DOUBLE PRECISION C0 !Cosinus de l'angle d'élévation solaire rmu(N0)
      DOUBLE PRECISION C1 !Cosinus de l'angle de visée rmu(J) > 0.
      DOUBLE PRECISION S0 !Sinus de l'angle d'elevation solaire rmu(N0)
      DOUBLE PRECISION S1 !Sinus de l'angle de visée rmu(J) > 0.

      DOUBLE PRECISION S12 ! Sin(Teta0) * Sin(Teta1) 

      DOUBLE PRECISION AT0 ! Transmission du TOA au sol pour la direction d'incidence.
					
      DOUBLE PRECISION AT(1:CTE_OS_NBMU_MAX) ! Transmission totale : 
                                             ! direction d'incidence (du TOA au sol) +
                                             ! direction d'émergence (du sol au niveau de sortie pour le champ montant). 
      
      DOUBLE PRECISION ATJ ! Valeur du tableau AT pour la direction J. 
      
      DOUBLE PRECISION SIGMA2 ! Carre du paramètre Sigma de la théorie de Cox et Munck.
      
      DOUBLE PRECISION P ! Probabilite de réflexion vers la direction (mu, phi).
      DOUBLE PRECISION COSDIF ! Cosinus de l'angle de diffusion.
      DOUBLE PRECISION COSKIP ! Cosinus de l'angle de rotation ki'.
      
      DOUBLE PRECISION R11 ! Eléments de la matrice de réflexion de Fresnel pour un repère lié au plan de réflexion.
      DOUBLE PRECISION R12 
      DOUBLE PRECISION R33
        
      DOUBLE PRECISION M11 ! Eléments de la matrice de réflexion dans le repère lié au plan méridien.
      DOUBLE PRECISION M21  
      DOUBLE PRECISION M31

      DOUBLE PRECISION COEF_SUN ! Valeur du rapport de l'éclairement solaire au TOA 
                                ! (PI) sur l'angle solide du disque solaire.   
       
      CHARACTER*CTE_LENFIC2 FICHOS ! Nom du fichier résultat des OS.
 
      INTEGER*4 IGLI ! Option de Glitter :
                     ! 0 : pas de Glitter / 1 : calcul du Glitter.

      INTEGER*4 IFRESNEL ! Indice précisant si les simulations introduisent la réflexion de Fresnel sur mer plate :
                         !    1 => pour un calcul de la réflexion de Fresnel.
                         !    Sinon => pas de prise en compte. 
      
      INTEGER*4 IROUJEAN ! Indice précisant si les simulations introduisent la réflexion par la BRDF de Roujean :
                         !    1 => pour un calcul avec la BRDF de Roujean.
                         !    Sinon => pas de prise en compte.      
      
      INTEGER*4 IRONDEAUX ! Indice précisant si les simulations introduisent la réflexion par la BPDF de Rondeaux :
                          !    1 => pour un calcul avec la BPDF de Rondeaux.
                          !    Sinon => pas de prise en compte.   
			  
      INTEGER*4 IBREON ! Indice précisant si les simulations introduisent la réflexion par la BPDF de Bréon :
                       !    1 => pour un calcul avec la BPDF de Bréon.
                       !    Sinon => pas de prise en compte.   
			  
      INTEGER*4 INADAL ! Indice précisant si les simulations introduisent la réflexion par la BPDF de Nadal :
                       !    1 => pour un calcul avec la BPDF de Nadal.
                       !    Sinon => pas de prise en compte.   

      INTEGER*4 IMAIGNAN ! Indice précisant si les simulations introduisent la réflexion par la BPDF de Maignan :
                         !    1 => pour un calcul avec la BPDF de Maignan.
                         !    Sinon => pas de prise en compte.   
	
      INTEGER*4 IPOLAR ! Option de simulation sans polarisation.
                       ! 0 : calculs sans polarisation.
                       ! 1 : calculs normaux avec polarisation.
								            			     
      INTEGER*4 N0 ! Numéro de mu de Gauss pour l'angle zénithal solaire. 
                   ! (angle complémentaire : mus = -rmu(N0)) 
                   ! (N0 > 0 => pour en tenir compte).

      INTEGER*4 IS ! Ordre de la décomposition en séries de Fourier.      
      INTEGER*4 IER ! code d'erreur =0 si pas d'erreur, =-1 sinon  
          
      INTEGER*4 J ! Indice du tableau des angles de Gauss.
      
      INTEGER*4 NBMU ! nombre d'angles utiles		
	
      DOUBLE PRECISION X ! Intermédiaire de calcul de Tan(TetaN) 
                         ! pour TetaN l'angle d'incidence à la normale de la facette.

C* Traitement                
C--------------    


C* Calcul de l'angle de diffusion pour les différentes directions d'observation
C-------------------------------------------------------------------------------
      C0 = RMU(N0)

      DO J= -NBMU, NBMU

         COSDIF = -C0*RMU(J)
     &	           +SIN(ACOS(C0))*SIN(ACOS(RMU(J)))*DCOS(PHI)

         ANGDIFF(J) = ACOS(COSDIF)*180.D+00/INCTE_PI

      ENDDO


C* Ouverture des fichiers 
C------------------------------------------------     
C     Ouverture du fichier résultat des OS
      OPEN(UNIT=16,FILE=FICHOS,FORM='UNFORMATTED',
     &      STATUS='OLD',ERR=991) 
     
C* Recombinaison pour l'azimut PHI des paramètres de Stokes I, Q et U 
C* à partir de la décomposition en séries de Fourier.  
C------------------------------------------------------------------
C*    Initialisation pour IS = 0    

C*    Lecture du fichier qui a été enregistré avec des 
C*    bornes de tableaux limitées au dimensionnement 
C*    utile -NBMU:NBMU
      READ(16,err=992) (Q3(J),J=-NBMU,NBMU),
     & 		  (U3(J),J=-NBMU,NBMU),
     &		         (I3(J),J=-NBMU,NBMU)
             
      
      DO 2 J= -NBMU, NBMU
         IF(J.EQ.0) GOTO 2
         XQT(J) = Q3(J)
         XUT(J) = U3(J)
         XIT(J) = I3(J)   
    2 CONTINUE

C*   Boucle sur l'ordre IS.    
      IS=1
      
C*      Lecture du fichier qui a été enregistré avec des 
C*      bornes de tableaux limitées au dimensionnement 
C*      utile -NBMU:NBMU
   10 READ(16,END=9,err=992) (Q3(J),J=-NBMU,NBMU),
     & 		        (U3(J),J=-NBMU,NBMU),
     &		               (I3(J),J=-NBMU,NBMU)
   
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

C* Ajout de la réflexion directe du faisceau solaire sur la mer. 
C-----------------------------------------------------------------     
      IF (IGLI.EQ.1) THEN

C* Transmission du TOA au sol pour la direction d'incidence
C------------------------------------------------------------ 
         C0 = RMU(N0)
         AT0=DEXP(-TAU/C0)

C* Transmission pour la direction d'incidence et pour chaque direction
C* d'émergence J du sol au niveau de sortie des résultats. 
C------------------------------------------------------------------      
         DO 11 J=1, NBMU
            AT(J) = AT0*DEXP(-(TAU-TAUOUT)/RMU(J))
  11     CONTINUE


C* Carre du paramètre Sigma de la théorie de Cox et Munck.
C---------------------------------------------------------- 
         SIGMA2 = .003 + .00512*WIND
         
         DO 101 J=1,NBMU
            C1=RMU(J)

C* Calcul de la probabilite P de réflexion dans la direction J 
C* pour la statistique de pentes de vague  
C----------------------------------------------------------    
            CALL SOS_GLITTE(SIGMA2,C0,C1,PHI,P)
      
C* Calcul des angles de diffusion et de rotation. 
C----------------------------------------------------------     
            CALL SOS_ANGLE(C0,C1,PHI,COSKIP,COSDIF)
      
C* Calcul des éléments de matrice de la réflexion de Fresnel 
C* pour un repère lié au plan de réflexion.
C----------------------------------------------------------   
            CALL SOS_REFLEX(COSDIF,IND_SURF,R11,R12,R33)
      
C* Calcul des éléments de matrice de la réflexion de Fresnel 
C* pour un repère lié au plan méridien (rotation). 
C----------------------------------------------------------        
            CALL SOS_MATRIC(COSKIP,R11,R12,M11,M21,M31)

C* Ajout de la réflexion directe pour un éclairement solaire 
C* normalisée (PI) attenue (ATJ) avec la probabilite P liée 
C* aux pentes de vagues.
C----------------------------------------------------------   
            ATJ =AT(J)
            XIT(J)= XIT(J)+M11*ATJ*P
	    
	    IF(IPOLAR.EQ.1) THEN
               XQT(J)= XQT(J)+M21*ATJ*P
               XUT(J)= XUT(J)+M31*ATJ*P
	    ENDIF
	
  101    CONTINUE
  
      ENDIF	!Test igli=1



C* Ajout de la réflexion du faisceau solaire direct sur mer plate  
C----------------------------------------------------------   

      IF (IFRESNEL.EQ.1) THEN
		
C* Test si le plan d'incidence correspond au plan principal solaire
C* et l'angle zénithal solaire fait parti des angles de sortie
C* (c'est a dire, s'il s'agit d'un angle de Gauss).     
          IF ((DCOS(PHI).EQ.1.D+00).AND.(N0.GT.0)) THEN

C* Transmission du TOA au sol pour la direction d'incidence 
              C0  = RMU(N0)
              AT0 = DEXP(-TAU/C0)

C* Transmission pour la direction d'incidence et de réflexion
C* au niveau de sortie des résultats     
              ATJ = AT0*DEXP(-(TAU-TAUOUT)/C0)	     

C* Coefficients de la matrice de réflexion de Fresnel pour
C* l'incidence solaire

              COSDIF = 1. - 2.*C0*C0
              CALL SOS_REFLEX(COSDIF,IND_SURF,R11,R12,R33) 

C* Valeur de la luminance solaire au TOA (normalisée à Es=PI)
C* == > Ls = Es/angle_solide_disque_solaire
              COEF_SUN = INCTE_PI / CTE_SOLAR_DISC_SOLID_ANGLE

C* Réflexion du faisceau solaire direct (Es = PI)	      
	       XIT(N0) = XIT(N0) + R11*COEF_SUN*ATJ 
	       IF(IPOLAR.EQ.1) XQT(N0) = XQT(N0) + R12*COEF_SUN*ATJ  
	      
	  ENDIF
	  
      ENDIF     !Test ifresnel=1



C* Ajout de la réflexion du faisceau solaire direct pour une BRDF
C* du modèle de Roujean  
C----------------------------------------------------------   

      IF (IROUJEAN.EQ.1) THEN
 
C* Transmission du TOA au sol pour la direction d'incidence 
         C0 = RMU(N0)
         AT0=DEXP(-TAU/C0)
	 
	  S0 = DSQRT(1.D+00-C0*C0)

C* Boucle sur les directions d'émergence  
	  DO 103 J=1,NBMU
	  
C*     Transmission pour la direction d'incidence et pour la direction
C*     d'émergence J du sol au niveau de sortie des résultats. 
	     C1 = RMU(J)
	     ATJ = AT0*DEXP(-(TAU-TAUOUT)/C1) 
	     
	     S1 = DSQRT(1.D+00-C1*C1)

C*     Calcul de la BRDF de Roujean (luminance normalisée*cos(tetav))	     
	     PHIRJ = INCTE_PI - PHI
            CALL SOS_CALC_F_ROUJEAN(K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &	                             C0,S0,C1,S1,PHIRJ,F)   
     
C*     Ajout de la réflexion directe pour un éclairement solaire 
C*     normalisée (PI) attenue (ATJ) de la réflexion de Roujean.
            XIT(J)= XIT(J) + ATJ*F/C1
            
103       CONTINUE

      ENDIF     !Test iroujean=1



C* Ajout de la réflexion du faisceau solaire direct pour une BPDF
C* du modèle de Rondeaux ou de Bréon
C----------------------------------------------------------   

      IF( (IRONDEAUX.EQ.1).OR.(IBREON.EQ.1).OR.(IMAIGNAN.EQ.1) )  THEN

C* Transmission du TOA au sol pour la direction d'incidence 
         C0 = RMU(N0)
         AT0=DEXP(-TAU/C0)

	  S0 = DSQRT(1.D+00-C0*C0)
	  
C* Boucle sur les directions d'émergence  
	  DO 105 J=1,NBMU
	  
C*     Transmission pour la direction d'incidence et pour la direction
C*     d'émergence J du sol au niveau de sortie des résultats. 
	     C1 = RMU(J)
	     ATJ = AT0*DEXP(-(TAU-TAUOUT)/C1) 	     

C***** Calcul de la BPDF de Rondeaux ou Bréon

C*        Calcul des angles de diffusion et de rotation.     
             CALL SOS_ANGLE(C0,C1,PHI,COSKIP,COSDIF)
      
C*       Calcul des éléments de matrice de la réflexion de Fresnel 
C*       pour un repère lié au plan de réflexion.
             CALL SOS_REFLEX(COSDIF,IND_SURF,R11,R12,R33)
      
C*       Calcul des éléments de matrice de la réflexion de Fresnel 
C*       pour un repère lié au plan méridien (rotation).    
             CALL SOS_MATRIC(COSKIP,R11,R12,M11,M21,M31)

C*       Ajout de la réflexion directe pour un éclairement solaire 
C*       normalisée (PI) attenue (ATJ) de la réflexion de Rondeaux 
C*       ou de Bréon.
	      IF (IRONDEAUX.EQ.1) P = 1. / ( 4.*(1+C1/C0) ) ! Modèle de Rondeaux
	      IF (IBREON.EQ.1) P = 1. / ( 4.*C1 )		   ! Modèle de Bréon
             IF (IMAIGNAN.EQ.1) THEN	  
  
	          S1 = DSQRT(1.D+00-C1*C1)
                 S12 = S0*S1
                 CALL SOS_CALCG_MAIGNAN(C0,C1,S12,PHI,COEF_C_MAIGNAN,P)
                 P = P / ( 4.*C1 )	   ! Modèle de Maignan

           ENDIF 

	    XIT(J)= XIT(J)+M11*ATJ*P
	    IF(IPOLAR.EQ.1) THEN
	       XQT(J)= XQT(J)+M21*ATJ*P
              XUT(J)= XUT(J)+M31*ATJ*P
	    ENDIF         
	           
105       CONTINUE


      ENDIF     !Test irondeaux ou ibreon ou imaignan = 1
      



C* Ajout de la réflexion du faisceau solaire direct pour une BPDF
C* du modèle de Nadal  
C----------------------------------------------------------   

      IF (INADAL.EQ.1) THEN

C* Transmission du TOA au sol pour la direction d'incidence 
          C0 = RMU(N0)
          AT0=DEXP(-TAU/C0)
	  
C* Boucle sur les directions d'émergence  
	  DO 107 J=1,NBMU
	  
C*     Transmission pour la direction d'incidence et pour la direction
C*     d'émergence J du sol au niveau de sortie des résultats. 
	      C1 = RMU(J)
	      ATJ = AT0*DEXP(-(TAU-TAUOUT)/C1) 	     

C***** Calcul de la BPDF de Nadal

C*        Calcul des angles de diffusion et de rotation.     
             CALL SOS_ANGLE(C0,C1,PHI,COSKIP,COSDIF)
      
C*       Calcul des éléments de matrice de la réflexion de Fresnel 
C*       pour un repère lié au plan de réflexion.
             CALL SOS_REFLEX(COSDIF,IND_SURF,R11,R12,R33)
      
C*       Calcul des éléments de matrice de la réflexion de Fresnel 
C*       pour un repère lié au plan méridien (rotation).    
             CALL SOS_MATRIC(COSKIP,R11,R12,M11,M21,M31)

C*       Ajout de la réflexion directe pour un éclairement solaire 
C*       normalisée (PI) attenue (ATJ) de la réflexion de Nadal

C*         F21 de Fresnel (convention de signe de Nadal)                
             F21FRESNEL = -R12     
      
C*         Calcul de F21 de Nadal
             F21NADAL = -BETA_NADAL * F21FRESNEL / (C0 + C1)     
	      F21NADAL = ALPHA_NADAL * ( 1.D+00 - DEXP(F21NADAL) )

C*         Calcul de P               
             IF(F21FRESNEL.LT.1.0D-10) THEN	!Développement limité au voisinage 
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



C******************************************************************
C* Mise à zéro des termes de polarisation trop faible
C*
C* --> Evite des erreurs d'arrondi (U non nul) 
C*     pour le passage de PHI=0 à PHI=PI  
C*
C* Mise à zéro du terme I si inférieur à 1.E-99.
C****************************************************************** 
      DO J= -NBMU,NBMU

         IF (XIT(J).LE.1.D-99) XIT(J)=0.D+00
	  IF (DABS(XQT(J)).LT.CTE_THRESHOLD_Q_U_NULL) XQT(J)=0.D+00
         IF (DABS(XUT(J)).LT.CTE_THRESHOLD_Q_U_NULL) XUT(J)=0.D+00

      ENDDO
      
 

C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur : retour du status 1 au programme appelant
C------------------------------------------------------
  991 WRITE(6,*) "SOS_TRPHI : Error while opening a file"
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 
     &  "SOS_TRPHI : Error while reading or writing on a file"
      IER=-1 
      GOTO 9999

 9999 CONTINUE

C* Fermeture fichier
C-------------------
      CLOSE(16)
            
      RETURN  	!FIN DE LA PROCEDURE SOS_TRPHI
      END



C==============================================================================
C PROCEDURE: SOS_GLITTE
C ==========
C      Cette procédure calcule pour une incidence (mup, phip) sur des vagues
C      la probabilite de réflexion vers la direction (mu, phi).
C      L'agitation de la mer depend de la vitesse du vent. 
C
C      Calculs bases sur la théorie de Cox et Munck (1954).
C
C
C Description des paramètres
C --------------------------
C     SIG	  (double)  (E)  Carre du paramètre Sigma de la théorie de 
C      			   Cox et Munck.
C     C0	  (double)  (E)  Valeur absolue du cosinus de l'angle zénithal
C      			   d'incidence (mup).
C     C1	  (double)  (E)  Valeur absolue du cosinus de l'angle zénithal 
C     			          d'émergence (mu).
C     PHI	  (double)  (E)  Azimut relatif entre la direction d'émergence
C                              et celle d'incidence phi - phip (en radians).
C     P	  (double)  (S)  Probabilite de réflexion vers la direction (mu, phi).
C
C Common utilisé:
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
      
C* Définition des variables   
C*-----------------------------------------     
      DOUBLE PRECISION SIG !Carre du paramètre Sigma de la théorie de Cox et Munck.

      DOUBLE PRECISION C0 !Valeur absolue du cosinus de l'angle zénithal d'incidence (mup).
      DOUBLE PRECISION C1 !Valeur absolue du cosinus de l'angle zénithal d'émergence (mu).
      DOUBLE PRECISION PHI !Azimut relatif entre la direction d'émergence et celle d'incidence phi - phip (en radians).
      
      DOUBLE PRECISION P !Probabilite de réflexion vers la direction (mu, phi).
      
      DOUBLE PRECISION X1
      DOUBLE PRECISION X2
      DOUBLE PRECISION X3
      DOUBLE PRECISION C0N ! Cosinus de l'angle zénithal de la normale à la facette des vagues.
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
           PP=(1/SIG)*DEXP(XXX)
           P=PP/(4*C1*C0N**4)
       ENDIF     
 
       RETURN
       END  	!FIN DE LA PROCEDURE SOS_GLITTE



C==============================================================================
C PROCEDURE: SOS_ANGLE
C ==========
C      Cette procédure calcule :
C         _ le cosinus de l'angle de diffusion entre la direction d'incidence 
C	     solaire et la direction de visée.
C	   _ le cosinus de l'angle de rotation KI'. 
C
C Description des paramètres
C --------------------------
C       C0      (double)  (E) Valeur absolue du cosinus de l'angle zénithal solaire.
C       C1	  (double)  (E) Valeur absolue du cosinus de l'angle zénithal de visée.
C       PHI	  (double)  (S) Azimut relatif entre la direction de visée
C                             et celle du Soleil (en radians).
C       COSKIP  (double)  (S)  Cosinus de l'angle de rotation Ki'.
C       COSDIF  (double)  (S)  Cosinus de l'angle de diffusion.
C
C Common utilisé:
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

C* Définition des variables   
C*----------------------------------------
      DOUBLE PRECISION C0 !Valeur absolue du cosinus de l'angle zénithal solaire.
      DOUBLE PRECISION C1 !Valeur absolue du cosinus de l'angle zénithal de visée.
      DOUBLE PRECISION PHI !Azimut relatif entre la direction de visée et celle du Soleil (en radians).
      
      DOUBLE PRECISION COSKIP !(S) Cosinus de l'angle de rotation Ki'.
      DOUBLE PRECISION COSDIF !(S) Cosinus de l'angle de diffusion.
            				 
      DOUBLE PRECISION S ! 1 ou -1.
      DOUBLE PRECISION Z	

C* Traitement                
C--------------         
      S=1.
      IF (DSIN(PHI).GT.0.D+00)  S=-1.
      
      COSDIF=-C0*C1+DSQRT(1-C0**2)*DSQRT(1-C1**2)*DCOS(PHI)
      Z=S*(DSQRT(1-COSDIF**2))*(DSQRT(1-C1*C1))
      COSKIP=0.
      
      IF (DABS(Z).GT. CTE_SEUIL_Z)  COSKIP=(C1*COSDIF+C0)/Z
      
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_ANGLE





C==============================================================================
C PROCEDURE: SOS_REFLEX
C ==========
C      Cette procédure calcule les éléments R11, R12 et R33 de la matrice
C      de réflexion de Fresnel exprimee dans un repère lié au plan de réflexion.
C
C      La matrice de Fresnel s'ecrit :
C                  | R11(teta)  R12(teta)    0      |
C         F(teta)= | R12(teta)  R11(teta)    0      |
C                  |    0          0      R33(teta) |
C
C         avec :
C           _ teta : angle de réflexion par rapport à la normale au dioptre,
C           _ R11(teta) = 0.5 * (RL*RL + RR*RR),
C           _ R12(teta) = 0.5 * (RL*RL - RR*RR),
C           _ R33(teta) = RL*RR,
C
C         L'indice de réfraction relatif du dioptre, entre le milieu de 
C         réflexion et le milieu de transmission, est IND.
C
C         On calcule RL et RR en fonction de IND et teta, par :
C           _ RL : le coefficient de réflexion parallele
C
C               RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C           _ RR : le coefficient de réflexion perpendiculaire
C
C               RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C
C
C Description des paramètres
C --------------------------
C       COSDIF	  (double)  (E)  Cosinus de l'angle de diffusion.
C       IND	  (double)  (E)  indice de réfraction relatif 
C     				 (eau / air pour une réflexion sur mer).
C       R11	  (double)  (S) Elément de la matrice de réflexion de Fresnel.
C       R12	  (double)  (S) idem R11
C       R33	  (double)  (S) idem R11
C
C
C Common utilisé:
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

C* Définition des variables                
C*---------------------------    
      DOUBLE PRECISION COSDIF ! Cosinus de l'angle de diffusion.
      
      DOUBLE PRECISION IND !indice de réfraction relatif (eau / air pour une réflexion sur mer).
      
      DOUBLE PRECISION R11 !Eléments de la matrice de réflexion de Fresnel.
      DOUBLE PRECISION R12
      DOUBLE PRECISION R33
        
      DOUBLE PRECISION IND2 ! Carre de l'indice de réfraction relatif.
				
      DOUBLE PRECISION COSW ! Cosinus de l'angle de réflexion W 
                            ! (angle d'incidence par rapport à la normale du dioptre)
      
      DOUBLE PRECISION V ! Sinus carre de W.
      DOUBLE PRECISION X	
      
      DOUBLE PRECISION RL ! Coefficient de réflexion parallele.
      DOUBLE PRECISION RR ! Coefficient de réflexion perpendiculaire.
      

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
C      Cette procédure calcule les éléments M11, M21 et M31 de la matrice
C      de réflexion exprimee dans un repère lié au plan méridien (1ere colonne).
C      
C
C Description des paramètres
C --------------------------
C       COSKIP  (double)  (E) Cosinus de l'angle de rotation KI'      
C       R11	  (double)  (E) Elément de la matrice de réflexion de Fresnel
C     				  dans le repère lié au plan de réflexion.
C       R12	  (double)  (E) idem R11
C       M11	  (double)  (S) Elément de la matrice de réflexion dans le 
C                             repère lié au plan méridien.
C       M21	  (double)  (S) idem M11
C       M31	  (double)  (S) idem M11
C  
C
C Common utilisé:
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

C* Définition des variables   
C*------------------------------------------------   
      DOUBLE PRECISION COSKIP !Cosinus de l'angle de rotation KI'.
      
      DOUBLE PRECISION R11 ! Eléments de la matrice de réflexion de Fresnel dans le repère lié au plan de réflexion.
      DOUBLE PRECISION R12 
      
      DOUBLE PRECISION M11 !Elément de la matrice de réflexion dans le repère lié au plan méridien.
      DOUBLE PRECISION M21
      DOUBLE PRECISION M31
   
      DOUBLE PRECISION C2 ! Cosinus de 2*KI'.
      DOUBLE PRECISION S2 ! Sinus de 2*KI'.
      DOUBLE PRECISION X

C* Traitement                
C--------------               
      X=1.-DABS(COSKIP)
      C2=1.
      S2=0.
      IF (X.GE.CTE_SEUIL_X) THEN
         C2=2.*COSKIP*COSKIP-1.
         S2=2.*COSKIP*DSQRT(1.-COSKIP*COSKIP)
      ENDIF
      
      IF (COSKIP.EQ.0.D+00) R12=0.
      
      M11=R11
      M21=C2*R12
      M31=S2*R12
      
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_MATRIC




C==============================================================================
C PROCEDURE: SOS_OUTPUT_HEADER
C ==========
C      Cette procédure écrit des informations d'entête pour les fichiers
C      résultats de simulation du champ de luminance.
C      
C
C Description des paramètres
C --------------------------
C       FICID	  (I4)     (E) : Indice du fichier ouvert en écriture   
C       UPDOWN  (I4)     (E) : Indice pour champ montant (1) ou descendant
C       PHI1    (double) (E) : Direction azimutale définie par l'utilisateur
C       PHI2    (double) (E) : Direction PHI1 + PI
C       ZALT    (double) (E) : Altitude of the level of the output radiance
C
C Common utilisé:
C --------------
C     Aucun      
C
C Cas d'erreur :
C ------------
C     Aucun      
C
C==============================================================================
      SUBROUTINE SOS_OUTPUT_HEADER(FICID, UPDOWN, PHI1, PHI2, ZALT)
      
      IMPLICIT NONE

C* Définition des variables   

      INTEGER FICID	
      INTEGER UPDOWN ! Indice précisant si champ montant (1) ou descendant (autre)	
      DOUBLE PRECISION PHI1 ! Azimut relatif (radians)
      DOUBLE PRECISION PHI2 ! Azimut relatif PHI1 + Pi (radians)
      DOUBLE PRECISION ZALT ! Altitude of the radiance (km)

C* Traitement                
C--------------    

      IF (UPDOWN.EQ.1) THEN
          WRITE(FICID,'(A55)') 
     &    '# UPWARD RADIANCE FIELD VERSUS THE VIEWING ZENITH ANGLE'
      ELSE
          WRITE(FICID,'(A57)') 
     &    '# DOWNWARD RADIANCE FIELD VERSUS THE VIEWING ZENITH ANGLE'
      ENDIF

      WRITE(FICID,'(A43)') 
     &    '# (RELATIVE AZIMUTH AND ALTITUDE ARE FIXED)'
      WRITE(FICID,'(A102)') '#--------------------------------------'//
     &'---------------------------------------------------------------'
	  
      WRITE(FICID,'(A30)')  '# Relative azimuth (degrees) :'
      WRITE(FICID,'(A1)')   '#'
      WRITE(FICID,'(A36)')  '#      Relative azimuth convention :'

      IF (UPDOWN.EQ.1) THEN
         WRITE(FICID,'(A58)')  
     &   '#        180° <-> Satellite and Sun in the same half-plane'
         WRITE(FICID,'(A96)')  
     &     '#          0° <-> Satellite and Sun in opposite'//
     &     ' half-planes with respect to the zenith direction'          
      ELSE
          WRITE(FICID,'(A66)')  
     &     '#        180° <-> Viewing direction and Sun'//
     &     ' in the same half-plane'
         WRITE(FICID,'(A104)')  
     &     '#          0° <-> Viewing direction and Sun in opposite'//
     &     ' half-planes with respect to the zenith direction'              
      ENDIF

      WRITE(FICID,'(A1)')   '#'
      WRITE(FICID,'(A45)')
     &	    '#      Simulated relative azimuth (degrees) :'
      WRITE(FICID,10)  '#          for VZA < 0 (sign convention):',
     &	                  PHI2
      WRITE(FICID,10)  '#          for VZA > 0 (sign convention):',
     &	                  PHI1

      WRITE(FICID,'(A1)')   '#'
      WRITE(FICID,20)  
     &'# Value of the selected altitude for the output (km) :', ZALT
     
      WRITE(FICID,'(A1)')   '#'
      WRITE(FICID,'(A22)') '# Columns parameters :'
      WRITE(FICID,'(A48)')  
     &'#   VZA     :  Viewing Zenith Angle (in degrees)'
      WRITE(FICID,'(A44)')  
     &'#   SCA_ANG :  Scattering angle (in degrees)'
      WRITE(FICID,'(A64)')  
     &'#   I       :  Stokes parameter I at output altitude z'// 
     &' (in sr-1)'
      WRITE(FICID,'(A85)')  
     &'#              normalised to the extraterrestrial solar'// 
     &' irradiance (PI * L(z) / Esun)'
      WRITE(FICID,'(A64)')  
     &'#   Q       :  Stokes parameter Q at output altitude z'// 
     &' (in sr-1)' 
      WRITE(FICID,'(A66)')  
     &'#              normalised to the extraterrestrial solar'// 
     &' irradiance'
      WRITE(FICID,'(A64)')  
     &'#   U       :  Stokes parameter U at output altitude z'// 
     &' (in sr-1)'
      WRITE(FICID,'(A66)')  
     &'#              normalised to the extraterrestrial solar'// 
     &' irradiance '
      WRITE(FICID,'(A88)')   
     &	    '#   POL_ANG :  Polarization angle (in degrees).'// 
     &'  Note: if undefined the value is -999.00'
      WRITE(FICID,'(A44)')   
     &	    '#   POL_RATE:  Degree of polarization (in %)'
      WRITE(FICID,'(A55)')  
     &	    '#   IPOL    :  Polarized intensity at level z (in sr-1)'
      WRITE(FICID,'(A88)') 
     &	    '#              normalised to the extraterrestrial solar'// 
     &      ' irradiance (PI * Lpol(z) / Esun)'
      WRITE(FICID,'(A102)') '#--------------------------------------'//
     &'---------------------------------------------------------------'
	  

      WRITE(FICID,'(A91)') '#   VZA     SCA_ANG'                //
     &               '       I              Q              U'   //
     &               '         POL_ANG  POL_RATE    IPOL'

      WRITE(FICID,'(A92)') '#(degrees) (degrees)    (sr-1)'    //
     &'         (sr-1)         (sr-1)      (degrees)  (%)      (sr-1)'

 10   FORMAT(A41,1X,F7.3)
 20   FORMAT(A54,1X,F7.3)
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_OUTPUT_HEADER
	  




C==============================================================================
C PROCEDURE: SOS_OUTPUT_HEADER_POLAR_DIAG
C ==========
C      Cette procédure écrit des informations d'entête pour les fichiers
C      résultats de simulation du champ de luminance (cas d'un diagramme polaire)
C      
C
C Description des paramètres
C --------------------------
C       FICID	  (I4)     (E) : Indice du fichier ouvert en écriture   
C       UPDOWN  (I4)     (E) : Indice pour champ montant (1) ou descendant
C       ZALT    (double) (E) : Altitude of the level of the output radiance
C
C Common utilisé:
C --------------
C     Aucun      
C
C Cas d'erreur :
C ------------
C     Aucun      
C
C==============================================================================
      SUBROUTINE SOS_OUTPUT_HEADER_POLAR_DIAG(FICID, UPDOWN, ZALT)
      
      IMPLICIT NONE

C* Définition des variables   

      INTEGER FICID	
      INTEGER UPDOWN ! Indice précisant si champ montant (1) ou descendant (autre)	
      DOUBLE PRECISION ZALT ! Altitude of the radiance (km)

C* Traitement                
C--------------    

      IF (UPDOWN.EQ.1) THEN
          WRITE(FICID,'(A72)') 
     &    '#UPWARD RADIANCE FIELD VERSUS THE AZIMUTH ANGLE'// 
     &    ' AND VIEWING ZENITH ANGLE'
      ELSE
          WRITE(FICID,'(A74)') 
     &    '#DOWNWARD RADIANCE FIELD VERSUS THE AZIMUTH ANGLE'// 
     &    ' AND VIEWING ZENITH ANGLE'
      ENDIF

      WRITE(FICID,'(A17)') 
     &    '#(ALTITUDE FIXED)'
      
      WRITE(FICID,'(A102)') '#--------------------------------------'//
     &'---------------------------------------------------------------'
	  
      WRITE(FICID,'(A31)')  '# Relative azimuth convention :'
      WRITE(FICID,'(A58)')  
     &   '#        180° <-> Satellite and Sun in the same half-plane'
      WRITE(FICID,'(A96)')  
     &     '#          0° <-> Satellite and Sun in opposite'//
     &     ' half-planes with respect to the zenith direction'
      WRITE(FICID,'(A1)')   '#'
      WRITE(FICID,20)  
     &'# Value of the selected altitude for the output (km) :', ZALT
     
      WRITE(FICID,'(A1)')   '#'
      WRITE(FICID,'(A22)') '# Columns parameters :'
      WRITE(FICID,'(A50)')  
     &'#   PHI     :  Relative azimuth Angle (in degrees)'
      WRITE(FICID,'(A48)')  
     &'#   VZA     :  Viewing Zenith Angle (in degrees)'
      WRITE(FICID,'(A44)')  
     &'#   SCA_ANG :  Scattering angle (in degrees)'
      WRITE(FICID,'(A64)')  
     &'#   I       :  Stokes parameter I at output altitude z'// 
     &' (in sr-1)'
      WRITE(FICID,'(A85)')  
     &'#              normalised to the extraterrestrial solar'// 
     &' irradiance (PI * L(z) / Esun)'
      WRITE(FICID,'(A64)')  
     &'#   Q       :  Stokes parameter Q at output altitude z'// 
     &' (in sr-1)' 
      WRITE(FICID,'(A66)')  
     &'#              normalised to the extraterrestrial solar'// 
     &' irradiance'
      WRITE(FICID,'(A64)')  
     &'#   U       :  Stokes parameter U at output altitude z'// 
     &' (in sr-1)'
      WRITE(FICID,'(A66)')  
     &'#              normalised to the extraterrestrial solar'// 
     &' irradiance '
      WRITE(FICID,'(A88)')   
     &	    '#   POL_ANG :  Polarization angle (in degrees).'// 
     &'  Note: if undefined the value is -999.00'
      WRITE(FICID,'(A44)')   
     &	    '#   POL_RATE:  Degree of polarization (in %)'
      WRITE(FICID,'(A55)')  
     &	    '#   IPOL    :  Polarized intensity at level z (in sr-1)'
      WRITE(FICID,'(A88)') 
     &	    '#              normalised to the extraterrestrial solar'// 
     &      ' irradiance (PI * Lpol(z) / Esun)'

      WRITE(FICID,'(A102)') '#--------------------------------------'//
     &'---------------------------------------------------------------'
	  

      WRITE(FICID,'(A99)') '#   PHI      VZA     SCA_ANG'         //
     &               '        I              Q              U'   //
     &               '       POL_ANG  POL_RATE    IPOL'

      WRITE(FICID,'(A100)')
     &'#(degrees) (degrees) (degrees)    (sr-1) '    //
     &'        (sr-1)         (sr-1)    (degrees)  (%)      (sr-1)'

 10   FORMAT(A41,1X,F7.3)
 20   FORMAT(A54,1X,F7.3)
      RETURN
      END  	!FIN DE LA PROCEDURE SOS_OUTPUT_HEADER_POLAR_DIAG
	  


C==============================================================================
C PROCEDURE: SOS_POLAR
C ==========
C      Routine calculant l'angle de la polarisation, le taux de polarisation
C      et l'intensité polarisée pour un vecteur (I,Q,U) en entrée
C
C      L'angle de polarisation XAN est définie de -PI/2 à +PI/2.
C
C      Les cas Q>0 et Q<0 sont établis à partir de l'analyse des valeurs de
C      cos(2*XAN) et sin(2*XAN) couplés au domaine de définition de XAN.
C
C      Le cas Q = 0 est établi à partir des valeurs possibles de XAN pour Q = 0.
C      Q = U = 0 correspond à une absence de polarisation ==> Angle indéterminé.
C
C Description des paramètres 			
C ------------------------------
C
C     XI (double)     (E) : Valeur du paramètre de Stokes I
C     XQ (double)     (E) : Valeur du paramètre de Stokes Q
C     XU (double)     (E) : Valeur du paramètre de Stokes U
C     XAN (double)    (S) : Valeur du l'angle de la polarisation (degrés)
C     TPOL (double)   (S) : Valeur du taux de polarisation (%)
C     LPOL (double)   (S) : Intensité polarisée normalisée à l'éclairement solaire au TOA 
C                           (PI * Lpol(z) / Esun)	       
C
C Common utilisé:
C --------------
C  Aucun
C
C Cas d'erreur :
C --------------
C  Aucun
C
C  Sous programmes utilisés:
C -------------------------
C  Aucun
C
C
C MOD:VERSION:1.0: 22/08/2022: Initial version of the SOS-ABS code 
C    Imported from OSOAA V1.6 (routine version 1.1: 04/06/2018)
C
C==============================================================================

      SUBROUTINE SOS_POLAR(XI,XQ,XU,XAN,TPOL,LPOL)

      IMPLICIT NONE     
                       	    
C* Définition des variables                         
C*-----------------------------------------------------------------
      
      DOUBLE PRECISION XI ! (E) Paramètre de Stokes I 
      DOUBLE PRECISION XQ ! (E) Paramètre de Stokes Q.
      DOUBLE PRECISION XU ! (E) Paramètre de Stokes U.
      
      DOUBLE PRECISION XT ! Rapport U / Q
     					     
      DOUBLE PRECISION XAN ! (S) Angle de la polarisation (degrés)
	
      DOUBLE PRECISION LPOL ! (S) Intensité polarisée normalisée à l'éclairement solaire au TOA 	       
      DOUBLE PRECISION TPOL ! (S) Taux de polarisation en %.
	

C* Traitements                       
C*------------      
      
      IF (XQ.NE.0.) THEN !Cas XQ!=0
 
         XT= XU/XQ
		 
         IF (XQ.GT.0.) THEN ! Cas si XQ > 0

            XAN= 90.*ATAN(XT)/INCTE_PI    
   
         ELSE ! Cas si XQ < 0

	     IF (XU.GT.0.) THEN
               XAN=90.+90.*DATAN(XT)/INCTE_PI
            ELSE
              XAN= -90.+90.*DATAN(XT)/INCTE_PI
            ENDIF

         ENDIF             

      ELSE !Cas XQ=0

          IF (XU.GT.0.) THEN
              XAN=  45.
          ELSE
             IF (XU.LT.0) THEN
                XAN= -45.
             ELSE !Cas XQ=XU indéfinie
                XAN= INCTE_VALEUR_INDEF
             ENDIF
          ENDIF

      ENDIF
	
      LPOL = SQRT(XQ*XQ+XU*XU)
      
      IF (XI.NE.0.D+00) THEN
         TPOL = 100.*LPOL/XI
      ELSE
         TPOL = INCTE_VALEUR_INDEF
      ENDIF
            
9999  RETURN   
            
      END      !FIN DE LA PROCEDURE SOS_POLAR   
      
