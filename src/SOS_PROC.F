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
C* FICHIER: SOS_PROC.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: 
C*      Cette procédure réalise le contrôle des paramètres de simulation.
C*      Elle exécute successivement les différents modules du code:
C*         - Définition des angles utilisées : appel de SOS_ANGLES
C*         - Calcul des propriétés radiatives des aérosols : appel de SOS_AEROSOLS
C*         - Estimation de l'AOT à la longueur d'onde de référence
C*         - Calcul des propriétés radiatives de la surface : appel de SOS_SURFACE
C*         - Préparation des données du profil d'absorption gazeuse : appel de SOS_PREPA_ABSPROFILE
C*         - Préparation des données radiatives pour le calcul des OS : appel de SOS_PREPA_OS
C*
C*         Pour le mode CKD fin:
C*              Pour chaque exponentielle de la décomposition CKD:
C*              - Génération du profil en épaisseur optique d'absorption gazeuse pour une  exponentielle : 
C*                appel de SOS_ABSPROFILE
C*              - Génération du profil en épaisseur optique (diffusion + absorption) : appel de SOS_PROFILE
C*              - Calcul du champ de luminance : appel de SOS
C*              - Aggrégation des luminances (et transmissions) selon la décomposition CKD : appel de SOS_AGGREGATE
C*
C*         Pour le mode CKD simplifié: Calcul de l'épaisseur optique d'absorption par CKD 
C*              - Génération du profil en épaisseur optique d'absorption gazeuse : appel de SOS_ABSPROFILE
C*              - Génération du profil en épaisseur optique (diffusion + absorption) : appel de SOS_PROFILE
C*              - Calcul du champ de luminance : appel de SOS
C*
C*         Recomposition de la luminance à partir des résultats de la décomposition en série de Fourier selon l'azimut : 
C*         appel de SOS_TRPHI_OPTION
C*
C*      Procédure pouvant être appelée depuis un programme Python (binding Python).
C*
C* AUTEURS: 
C*   CS GROUP France, X.Lenot, T.Erudel, B. Lafrance
C*
C* DATE: 23/11/2023
C*
C* MOD:VERSION:1.0: Création du fichier
C*
C* MOD:VERSION:1.1: 05/04/2023
C*    - Introduction d'un controle de la somme des coefficients CKD AIK et normalisation pour obtenir strictement 1.
C*   
C*  MOD:VERSION:1.2: 19/04/2023
C*    - Suppression de la variable IND_ANGOUT à l'appel de SOS_FIN_OS
C*    - Ajout de la variable IND_ANGOUT_FIN aux paramètres de sortie de SOS_PROC
C*   
C*  MOD:VERSION:1.3: 10/08/2023
C*    - Récupération du coeffficient de troncature de la fonction de phase aérosols en sortie de la routine SOS_AEROSOLS
C*    - Ajout du coeffficient de troncature de la fonction de phase aérosols en paramètre de sortie 
C*      (binding Python et routine vers programme principal Fortran)
C*
C*  MOD:VERSION:1.4: 15/09/2023
C*    - Correction des boucles sur indices pour l'initialisation à 1 de la table KDIS_AI(I,J,K)
C*    - Suppression d'une virgule après la commande WRITE(...) en un endroit (sortie fichier log).
C*      Impact nul sur les résultats de simulations, mais levant un warning à la compilation.
C*    - Suppression de l'appel à la routine SOS_FIN_OS en faveur d'une appel direct de SOS_TRPHI_OPTION
C*
C*  MOD:VERSION:1.5: 14/11/2023
C*    - Suppression des tests épaisseurs optiques moléculaires et aérosols > 0.0001 (O. à la place) pour vérifier les valeurs HR et HA fournies.
C*    - Correcte dénomination de la formule de calcul de l'épaisseur optique moléculaire (CNES, Perbos, 1982) au lieu de Hansen & Travis
C*
C----------------------------------------------------------------------------
C Définition des constantes   
C---------------------------------------------------------------------------- 
C Constantes utilisées :
C    CTE_LENDIR : Longueur des noms de répertoires.
C    CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_LENCOM : Longueur max d'une ligne de commandes
C    CTE_WAMIN: Longueur d'onde (en micron) minimum
C    CTE_WAMAX: Longueur d'onde (en micron) maximum
C    CTE_MDF : Facteur de dépolarisation moléculaire
C    CTE_CKD_NWVL_MAX: Nombre max de nombres d'ondes pour le calcul de l'absorption
C    CTE_CKD_NAI_MAX: Nombre max d'exponientielles par gaz absorbant (méthode CKD)
C    CTE_CKD_NT_MAX   : Nombre maximal de températures des données CKD
C    CTE_CKD_NP_MAX   : Nombre maximal de pressions des données CKD
C    CTE_HT_STD_PSURF: Pression atmosphèrique standard
C    CTE_NBABS: Nombre de gaz absorbants
C    CTE_ABS_NBLEV : Nombre de couches dans le profil de gaz absorbants
C    CTE_ABS_NBCOL : Nombre de colonnes des fichiers de profils des gaz d'absorption 
C    CTE_NOT_DEFINED_VALUE_INT: valeur par défaut pour les entiers
C    CTE_NOT_DEFINED_VALUE_DBLE: valeur par défaut pour les réels
C    CTE_DEFAULT_FICGRANU: nom par défaut du fichier granulométrie 
C    CTE_REMOVE_PREVIOUS_BIN : Option de suppression ou non d'un fichier résultat 
C    CTE_OS_NBMU_MAX  : Nombre maximal d'angles de Gauss positifs pour les calculs de luminance
C    CTE_OS_NB_MAX :    Valeur maximale pour l'ordre limite du développement en  polynomes de Legendre 
C    CTE_TOA_ALT : Altitude du niveau TOA (km).
C    CTE_DEFAULT_IGMAX : Nombre maximal d'ordres de diffusion par défaut
C
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define INCTE_IDLOG_MAIN 77
#define INCTE_IDLOG_PROFIL 88
#define INCTE_IDLOG_OS 99     
#define INCTE_NB_Modeles 4
#define INCTE_NB_PROFIL 2


C
C Description des paramètres
C ---------------------------
C
C   Entrées:
C			   
C     RESROOT (CHARACTER) : (E) Répertoire de travail et de sortie des résultats
C				    (chemin complet)
C
C     FICMAIN_LOG (CHARACTER) : (E) Nom du fichier Trace du programme principal.
C                                   (nom fichier + extension) 
C
C     WA_SIMU (DOUBLE) : (E) Longueur d'onde de simulation de la luminance en microns.
C
C     NBMU_GAUSS_LUM (I4) : (E) Nombre d'angles de Gauss pour les calculs de luminance.
C
C     FICANGLES_USER_LUM (CHARACTER) : (E) nom complet du fichier des angles "utilisateurs" pour les calculs de luminances.
C                 			        (répertoire + nom fichier + extension)
C
C     TETAS (DOUBLE) : (E) Angle zenithal solaire en degrés (0 < tetas < 90).
C      		      (angle complémentaire : mus= -cos(tetas))
C
C     FICANGLES_RES_LUM (CHARACTER) : (E) Fichier angles pour le calcul des luminances (dans RESROOT/SOS)
C
C     NBMU_GAUSS_MIE (I4) : (E) Nombre d'angles de Gauss pour les calculs de fonctions de phase
C
C     FICANGLES_USER_MIE (CHARACTER) : (E) Nom complet du fichier des angles "utilisateurs" pour les calculs de Mie.
C                 			        (répertoire + nom fichier + extension)
C
C     FICANGLES_RES_MIE : (E) Fichier angles pour les calcul de Mie (dans RESROOT/SOS)
C
C     FICANGLOG (E) : Nom complet du fichier trace (repertoire + nom fichier + extension)

C     WAREF_AOT (DOUBLE) : (E) Longueur d'onde de référence pour l'AOT.
C
C     AOT_REF (DOUBLE) : (E) Epaisseur optique des aérosols à la longueur d'onde de référence (non tronquée).
C
C     ITRONC_AER (I4) : (E) Option de troncature de la fonction de phase des aérosols 
C
C     FICGRANU_LOG (CHARACTER) : (E) Nom du fichier TRACE des calculs d'IOP aérosols (dans RESROOT/LOG)
C
C     FICMIE_LOG (CHARACTER) : (E) Nom du fichier TRACE des calculs de Mie (dans RESROOT/LOG)
C
C     DIR_MIE (CHARACTER) : (E) Répertoire de sortie des calculs de Mie
C
C     FICGRANU (CHARACTER) : (E) Nom du fichier résultats des calculs d'IOP aérosols (dans RESROOT/SOS)
C
C     IMOD_AER (I2) : (E) Indice de type de granulometrie mono ou multi-modale
C			        0 : Mono-modal models.
C			        1 : WMO models.
C			        2 : Shettle & Fenn models.
C		               3 : LND bi-modal models.
C                            4 : Use of external phase functions
C
C     RN_WA, IN_WA, RN_WAREF, IN_WAREF (DOUBLE) : (E) Indice de réfraction des aérosols dans le cas mono-mode 
C                    Partie réelle et partie imaginaire
C                    Pour la longueur d'onde de simulation de la luminance et pour celle de référence de l'AOT
C
C     IGRANU (I2) : (E) Indice de modèle mono-mode
C			        1 : LND
C			        2 : Junge Law
C
C     LND_RADIUS_MMD_AER, LND_LNVAR_MMD_AER (DOUBLE) : (E) Rayon modal et variance de la loi LND d'un cas mono-modal.
C
C     JD_SLOPE_MMD_AER(DOUBLE) : (E) Pente de la loi de Junge d'un cas mono-modal
C
C     JD_RMIN_MMD_AER, JD_RMAX_MMD_AER (DOUBLE) : (E) Rayon minimal et maximal de la loi de Junge.
C
C     IModele_WMO (I2) : (E) Type du modèle WMO
C      		        1 : Continental WMO model.
C			        2 : Maritime WMO model.
C			        3 : Urban WMO model.
C		               4 : WMO model by user définition.
C
C     C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO (DOUBLE) : (E) Valeurs utilisateur du mélange WMO de l'utilisateur
C                 Proportion pour particules Dust-Like, Water Soluble, Oceanic et Soot 
C
C     IModele_SF (I2) : (E) Type du modèle Shettle & Fenn
C      		1 : Tropospheric S&F model.
C			2 : Urban S&F model.
C			3 : Maritime S&F model.
C		       4 : Coastal S&F model.
C
C     RH (DOUBLE) : (E) Humidité relative (%) du modèle Shettle & Fenn
C
C     MODE_PARAM_BILND (I2) : (E) Indice du choix de description d'un mélange de 2 LND.
C      		1 : Use of user volumetric concentrations.'
C      		2 : Use of user ratio AOT_coarse / AOT_total.'
C     
C     USER_CV_COARSE, USER_CV_FINE (DOUBLE) : (E) Valeur utilisateur de la concentration volumique 
C                                                 de la "LND coarse" et "LND fine"
C
C     RTAUct_WAREF (DOUBLE) : (E) Rapport des épaisseurs optique AOT_coarse / AOT_total pour la longueur d'onde de référence
C
C     BMD_CM_MRWA, BMD_CM_MIWA, BMD_CM_MRWAREF, BMD_CM_MIWAREF (DOUBLE) : (E) 
C                    Indice de réfraction pour le mode Coarse d'un modèle bimodal d'aérosols
C                    Partie réelle et partie imaginaire
C                    Pour la longueur d'onde de simulation de la luminance et pour celle de référence de l'AOT
C
C     BMD_CM_RMODAL, BMD_CM_VAR (DOUBLE) : (E) Rayon modal et variance de la loi LND 
C                                              pour le mode Coarse d'un modèle bimodal d'aérosols.

C     BMD_FM_MRWA, BMD_FM_MIWA, BMD_FM_MRWAREF, BMD_FM_MIWAREF(DOUBLE) : (E) 
C                    Indice de réfraction pour le mode Fine d'un modèle bimodal d'aérosols
C                    Partie réelle et partie imaginaire
C                    Pour la longueur d'onde de simulation de la luminance et pour celle de référence de l'AOT

C     BMD_FM_RMODAL, BMD_FM_VAR (DOUBLE) : (E) Rayon modal et variance de la loi LND 
C                                              pour le mode Fine d'un modèle bimodal d'aérosols.
C
C     FICEXTDATA_AER (CHARACTER) : (E) Nom du fichier utilisateur donnant la fonction de phase aérosols
C                                      et les coefficients de diffusion et d'extinction.
C                 			    (chemin complet : répertoire + nom fichier )
C
C     FICMIXTURE_AER (CHARACTER) : (E) Nom du fichier de définition d'un mélange de modèles d'aérosols par l'utilisateur
C                 			    (chemin complet : répertoire + nom fichier )
C
C     FICUSER_AER (CHARACTER) : (E) Nom du fichier aérosols fourni par l'utilisateur
C                 			    (chemin complet : répertoire + nom fichier )
C
C     FICPROFIL_LOG (CHARACTER) : (E) Nom du fichier Trace des profils par ordre CKD (dans RESROOT/LOG)
C
C     TR (DOUBLE) : (E) Epaisseur optique moléculaire à la longueur de simulation de la luminance
C
C     HR, HA (DOUBLE) : (E) Echelle de hauteur respectivement des molécules et des aérosols (km)
C
C     IPROFIL (I2) : (E) Type du profil d'aérosols
C      		1 : Aerosol profile described by an exponential decrease using a heigth scale
C			2 : Homogeneous aerosol layer between two altitudes.
C
C     ZMIN, ZMAX (DOUBLE) : (E) Altitudes minimale et maximale de la couche molécules + aérosols (km)
C                               dans le cas d'un mélange d'aérosols entre deux altitudes.
C
C     PSURF (DOUBLE) : (E) Pression de surface (mbar)
C
C     H2O (DOUBLE) : (E) Contenu intégré en vapeur d'eau (g/cm2)
C
C     O3 (DOUBLE) : (E) Contenu intégré en ozone (Dobson)
C
C     CO2 (DOUBLE) : (E) Concentration du dioxyde de carbone en surface (ppmv)
C
C     CH4 (DOUBLE) : (E) Concentration du méthane en surface (ppmv)
C
C     ABSPROFIL (I2) : (E) Type du profil en absorption
C			      0:User
C			      1:TROPICAL
C			      2:MLS
C			      3:MLW
C			      4:SAS
C			      5:SAW
C			      6:USST62
C			      7:No absorption
C
C     FICABSPROFIL (CHARACTER) : (E) Fichier profil en absorption fourni par utilisateur (cas profil utilisateur)
C                 			    (chemin complet : répertoire + nom fichier )
C
C     NUSTEP (DOUBLE) : (E) Pas du nombre d'onde (cm-1) dans le fichier de coeff CKD
C                           Doit valoir 1, 5 ou 10 cm-1
C
C     ISURF (I2) : (E) Type de surface
C      		0 : lambertian surface for an albedo RHO.
C			1 : lambertian surface + glitter (roughness water).
C			2 : lambertian surface + Fresnel's reflexion (plane water).
C		       3 : lambertian surface + Roujean's BRDF.
C			4 : lambertian surface + Roujean's BRDF + Rondeaux's BPDF.
C			5 : lambertian surface + Roujean's BRDF + Breon's BPDF .
C			6 : lambertian surface + Roujean's BRDF + Nadal's BPDF.
C
C     DIR_SURF (CHARACTER) : (E) Répertoire de stokage des fichiers de BRDF/BPDF des modèles de surface
C                 		     (chemin complet du répertoire )
C
C     FICSURF_LOG (CHARACTER) : (E) Nom du fichier Trace des calculs de BRDF/BPDF (dans RESROOT/LOG)
C
C     SURF_IND (DOUBLE) : (E) Indice de réfraction de la surface / air.	
C
C     WIND (DOUBLE) : (E) Vitesse du vent en surface (m/s)

C     K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN (DOUBLE) : (E) Paramètres du modèle de BRDF de Roujean
C
C     ALPHA_NADAL,  BETA_NADAL (DOUBLE) : (E) Paramètres du modèle de BPDF de Nadal
C
C     COEF_C_MAIGNAN (DOUBLE) : (E) Paramètre C*exp(-NDVI) du modèle de BPDF de Maignan
C
C     RHO (DOUBLE) : (E) Albédo de surface
C
C     FICSURF (CHARACTER) : (E) Nom d'un fichier de matrice de réflexion fourni par l'utilisateur
C                 		    (chemin complet : répertoire + nom fichier )
C
C     FICSOS_LOG (CHARACTER) : (E) Nom du fichier Trace des calculs de résolution de l'ETR par SOS_OS (dans RESROOT/LOG)
C
C     FICSOS_RES_BIN (CHARACTER) : (E) Nom du fichier résultat des calculs d'ordre successifs (dans RESROOT/SOS)
C                                      pour un traitement unitaire (i.e une exponentielle CKD si calcul d'absorption gazeuse)
C
C
C   Sorties:
C
C     FICTRANS  (CHARACTER) : (S) Nom du fichier optionnel de sortie des transmissions (dans RESROOT/SOS)
C
C     FICFLUX  (CHARACTER) : (S) Nom du fichier optionnel de sortie des flux up TOA et down BOA (dans RESROOT/SOS)
C
C     ZOUT (DOUBLE) : (E) Altitude du niveau de sortie du champ de luminance
C			          -1 =>  Sortie standard: Luminance émergeante au TOA + descendante au niveau du sol.
C			      z (km) =>  Luminance montante et descendante pour l'altitude z (0 <= z <= CTE_TOA_ALT).
C
C     IGMAX (I4) : (E) Ordre maximal d'interaction
C
C     IPOLAR (I4) : (E) Option de coupure de la polarisation.
C                          0 : simulation sans polarisation.
C			      1 : simulation normale avec polarisation.
C
C     ITRPHI (I4) : (E) Option de type de sortie : 
C                          1 : Plan de visée pour une valeur d'azimut.
C                          2 : Diagramme polaire en fonction de la direction d'azimut et de visée.
C
C     PHIOS (DOUBLE) : (E) Azimut relatif en degrés pour le cas d'une sortie dans un plan azimutal.
C
C     PAS_PHI (I4) : (E) Pas sur l'azimut en degrés (valeur entière) dans le cas d'un diagramme polaire.
C
C     IMODE_CKD_CALCUL (I4) : (E) Option de mode d'application de la décomposition CKD
C          1 : Calcul le plus fin avec un calcul de luminance par ordre de la décomposition
C              puis recomposition selon la pondération CKD.
C          2 : Calcul pour une épaisseur optique d'absorption estimée par l'approche CKD
C              --> Un unique calcul de luminance
C
C     TRACE (BOOLEAN) : (E) Vrai si génération d'un fichier Trace du programme principal
C
C
C     IER (I4) : (S) Code de retour d'une erreur (1 si erreur, 0 sinon)
C
C     LUM_NBMU (I4) : (S) Nombre d'angles effectifs pour les calculs de luminances   
C   
C     IND_ANGOUT_FIN (DOUBLE) : (S) Indices des tableaux d'angles pour lesquels il doit y avoir  
C                                   une sortie dans le fichier résultat des angles utilisateur
C      			        (1 : pour une sortie requise)  
C
C     PHI_FIN (DOUBLE) : (S) Table des valeurs d'azimuts relatifs du champ de sortie
C  
C     THETA_FIN (DOUBLE) : (S) Table des valeurs d'angles zénithaux de visée du champ de sortie
C  
C     SCA_UP_FIN, SCA_DOWN_FIN (DOUBLE) : (S) Tables des angles de diffusion pour les champs montant et descendant de sortie 
C
C     I_UP_FIN, Q_UP_FIN, U_UP_FIN (DOUBLE) : (S) Tables des valeurs I,Q et U pour les champs montant de sortie 
C  
C     I_DOWN_FIN, Q_DOWN_FIN, U_DOWN_FINC (DOUBLE) : (S) Tables des valeurs I,Q et U pour les champs descendant de sortie   
C   
C     POL_ANG_UP_FIN, POL_ANG_DOWN_FIN (DOUBLE) : (S) Tables des angles de polarisation 
C                                                     pour les champs montant et descendant de sortie 
C
C     POL_RATE_UP_FIN, POL_RATE_DOWN_FIN (DOUBLE) : (S) Tables des taux de polarisation 
C                                                       pour les champs montant et descendant de sortie 
C
C     L_POL_UP_FIN, L_POL_DOWN_FIN (DOUBLE) : (S) Table des valeurs de luminance polarisée
C                                                 pour les champs montant et descendant de sortie   
C     
C     TDIR_VRAI (DOUBLE) : (S) Transmission directe pour la direction d'incidence solaire de l'atmosphère vraie
C
C     FLUX_DIFF_DOWN (DOUBLE) : (S) Flux diffus descendant éclairant la surface, normalisé par le flux solaire incident au TOA 
C      
C     FLUX_DOWN (DOUBLE) : (S) Flux total descendant au niveau de la surface, normalisé par le flux solaire incident au TOA 
C     
C     EPLUS (DOUBLE) : (S) Flux diffus émergeant au TOA, normalisé par le flux solaire incident au TOA  
C
C 
C Description des fichiers utilisés
C ---------------------------------
C
C     Génération du fichier Trace principal en format texte.
C
C     Fichier FICTRANS : 
C        --> Fournit la tranmission directe et diffuse du TOA au sol pour la direction d'incidence solaire.
C            Fournit également la transmission diffuse montante de la surface au TOA pour les différents angles de visée.
C
C     Fichier FICFLUX : 
C        --> Fournit le flux direct et diffus descendant au niveau de la surface.
C            Fournit aussi le flux diffus montant au TOA.
C               NB: les valeurs de flux sont normalisées par l'éclairement solaire incident au TOA (Esun * cos SZA).
C            Le profil atmosphérique est rappelé : 
C                  Epaisseur optique moléculaire (MOT), 
C                  Epaisseur optique des aérosols (AOT),  
C                  Epaisseur optique d'absorption gazeuse (GOT).
C
C Common utilisé:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C
C    Divers contrôle de cohérence du paramètrage
C    
C    Contrôle de retour en erreur d'une routine appelée:
C          SOS_ANGLES, SOS_AEROSOLS, SOS_SURFACE
C          SOS_PREPA_ABSPROFILE, SOS_PREPA_OS, SOS_ABSPROFILE, SOS_PROFILE
C          SOS, SOS_AGGREGATE et SOS_TRPHI_OPTION
C
C--------------------------------------------------------------------------------------------------------------------
      SUBROUTINE SOS_PROC(RESROOT,
     &     FICMAIN_LOG,
     &     WA_SIMU,
     &     NBMU_GAUSS_LUM,
     &     FICANGLES_USER_LUM,
     &     TETAS,
     &     FICANGLES_RES_LUM,
     &     NBMU_GAUSS_MIE,
     &     FICANGLES_USER_MIE,
     &     FICANGLES_RES_MIE,
     &     FICANGLOG,
     &     WAREF_AOT, AOT_REF,
     &     ITRONC_AER,
     &     FICGRANU_LOG, FICMIE_LOG,
     &     DIR_MIE,
     &     FICGRANU,
     &     IMOD_AER,
     &     RN_WA, IN_WA, RN_WAREF, IN_WAREF,
     &     IGRANU,
     &     LND_RADIUS_MMD_AER, LND_LNVAR_MMD_AER,
     &     JD_SLOPE_MMD_AER, JD_RMIN_MMD_AER, JD_RMAX_MMD_AER,
     &     IModele_WMO, C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &     IModele_SF, RH,
     &     MODE_PARAM_BILND, USER_CV_COARSE, USER_CV_FINE,
     &     RTAUct_WAREF,
     &     BMD_CM_MRWA, BMD_CM_MIWA, BMD_CM_MRWAREF, BMD_CM_MIWAREF,
     &     BMD_CM_RMODAL, BMD_CM_VAR,
     &     BMD_FM_MRWA, BMD_FM_MIWA, BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &     BMD_FM_RMODAL, BMD_FM_VAR,
     &     FICEXTDATA_AER,
     &     FICMIXTURE_AER,
     &     FICUSER_AER,
     &     FICPROFIL_LOG,
     &     TR, HR, HA,
     &     IPROFIL,
     &     ZMIN, ZMAX,
     &     PSURF, H2O, O3, CO2, CH4,
     &     ABSPROFIL, FICABSPROFIL,
     &     NUSTEP,
     &     ISURF,
     &     DIR_SURF,
     &     FICSURF_LOG,
     &     SURF_IND,
     &     WIND,
     &     K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN,
     &     ALPHA_NADAL,  BETA_NADAL,
     &     COEF_C_MAIGNAN,
     &     RHO,
     &     FICSURF,
     &     FICSOS_LOG,
     &     FICSOS_RES_BIN,
     &     FICTRANS, FICFLUX,
     &     ZOUT,
     &     IGMAX,
     &     IPOLAR,
     &     ITRPHI, PHIOS, PAS_PHI,
     &     IMODE_CKD_CALCUL,
     &     IER, TRACE,
     &     LUM_NBMU,
     &     IND_ANGOUT_FIN,
     &     PHI_FIN, THETA_FIN, 
     &     SCA_UP_FIN, I_UP_FIN, Q_UP_FIN, U_UP_FIN,
     &     POL_ANG_UP_FIN, POL_RATE_UP_FIN, L_POL_UP_FIN,
     &     SCA_DOWN_FIN, I_DOWN_FIN, Q_DOWN_FIN, U_DOWN_FIN,
     &     POL_ANG_DOWN_FIN, POL_RATE_DOWN_FIN, L_POL_DOWN_FIN,
     &     TDIR_VRAI, FLUX_DIFF_DOWN, FLUX_DOWN, EPLUS,
     &     COEF_TRONCA)
      
      IMPLICIT NONE
      
C Constantes
C -----------------            
      
! Numéro identifiant du fichier Trace pour PROC.F
      INTEGER*2 IDLOG	    
      PARAMETER(IDLOG=INCTE_IDLOG_MAIN)

! Numéro identifiant du fichier Trace pour PROC.F
      INTEGER*4 IDLOGPROF
      PARAMETER(IDLOGPROF=INCTE_IDLOG_PROFIL)

! Numéro identifiant du fichier Trace pour PROC.F   
      INTEGER*4 IDLOGOS
      PARAMETER(IDLOGOS=INCTE_IDLOG_OS)
    
C Constantes AEROSOLS
C -------------------

!Nombre de modèles WMO ou Shettle & Fenn possibles.
      INTEGER CTE_NB_Modeles
      PARAMETER (CTE_NB_Modeles=INCTE_NB_Modeles)

C Constantes PROFIL
C -----------------

!Valeur max du paramètre profil 
      INTEGER CTE_NB_PROFIL
      PARAMETER (CTE_NB_PROFIL=INCTE_NB_PROFIL)


     	            
C* Définition des variables                         
C*-----------------------------------------------------------------                                    
      
C-- Variables liées à l'arborescence de travail

C     Longueur de la chaine de caractères RESROOT
      INTEGER*4  LDIR

C     Longueur de la chaine de caractères DIR_SURF
      INTEGER*4  LDIRSURF

C     Longueur de la chaine de caractères DIR_SURF/GLITTER
      INTEGER*4  LDIRGLIT

C     Longueur de la chaine de caractères DIR_SURF/ROUJEAN
      INTEGER*4  LDIRROUJ

C     Longueur de la chaine de caractères DIR_SURF/RH
      INTEGER*4  LDIRRH

C     Longueur de la chaine de caractères DIR_SURF/BREON
      INTEGER*4  LDIRBREON

C     Longueur de la chaine de caractères DIR_SURF/BREON
      INTEGER*4  LDIRNADAL

C     Longueur de la chaine de caractères DIR_SURF/MAIGNAN
      INTEGER*4  LDIRMAIGNAN

C     Longueur de la chaine de caractères DIR_MIE
      INTEGER*4  LDIRMIE

C     Longueur de la chaine de caractères DIR_SOS
      INTEGER*4  LDIRSOS

C     Longueur de la chaine de caractères DIR_LOG
      INTEGER*4  LDIRLOG

C     Longueur de la chaine de caractères DIR_TMP
      INTEGER*4  LDIRTMP

C     Longueur de la chaine de caractères d'un nom de fichier
      INTEGER*4 LFIC

C     Répertoire de sortie       
      CHARACTER*CTE_LENDIR RESROOT

C     Répertoire de sortie des calculs de surface
      CHARACTER*CTE_LENDIR DIR_SURF

C     Répertoire de sortie des calculs de glitter
      CHARACTER*CTE_LENDIR DIR_GLIT

C     Répertoire de sortie des calculs de la BRDF Roujean
      CHARACTER*CTE_LENDIR DIR_ROUJ

C     Répertoire de sortie des calculs de la BRDF Rondeaux & Herman
      CHARACTER*CTE_LENDIR DIR_RH

C     Répertoire de sortie des calculs de la BRDF Breon
      CHARACTER*CTE_LENDIR DIR_BREON

C     Répertoire de sortie des calculs de la BRDF Nadal
      CHARACTER*CTE_LENDIR DIR_NADAL

C     Répertoire de sortie des calculs de la BRDF Maignan
      CHARACTER*CTE_LENDIR DIR_MAIGNAN

C     Répertoire de sortie des calculs de Mie
      CHARACTER*CTE_LENDIR DIR_MIE

C     Répertoire de sortie des ordres successifs
      CHARACTER*CTE_LENDIR DIR_SOS

C     Répertoire de sortie des fichiers traces
      CHARACTER*CTE_LENDIR DIR_LOG

C     Répertoire de sortie des fichiers temporaires
      CHARACTER*CTE_LENDIR DIR_TMP


C-- Variables globales

C     Longueur d'onde de simulation des luminances (en microns)
      DOUBLE PRECISION WA_SIMU

C     Longueur de la chaine de caractères COMMAND
      INTEGER*4 LCOM   
    
C     Nom du fichier trace principal 
      CHARACTER*CTE_LENFIC2 FICMAIN_LOG

C     Chaine de commande pour appel système.
      CHARACTER*CTE_LENCOM  COMMAND

C     Indices de tableau
      INTEGER*4 I,J,K

C     Code d'erreur =0 si pas d'erreur, =-1 sinon
      INTEGER*4 IER

C     Vrai si écriture dans le fichier trace
      LOGICAL TRACE

C     Vrai si écriture dans le fichier trace du profil
      LOGICAL TRACEPROF

C     Vrai si écriture dans le fichier trace des OS
      LOGICAL TRACEOS

C     Vérifie l'existence d'un fichier
      LOGICAL EX

                  
C-- Variables liées à la définition des angles  

C     Angle zénithal solaire en degrés (0 < tetas < 90). 
C     (angle complémentaire : mus= -cos(tetas))    		
      DOUBLE PRECISION TETAS

C     Nombre d'angles de Gauss utilisés pour les calculs de luminance   
      INTEGER*4 NBMU_GAUSS_LUM

C     Nombre d'angles de Gauss utilisés pour les calculs de Mie    
      INTEGER*4 NBMU_GAUSS_MIE  


      !Calculs de luminance :
C     Chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_LUM

C     Chemin et nom du fichier des angles utilisés
      CHARACTER*CTE_LENFIC2 FICANGLES_RES_LUM

      !Calculs de fonctions de phase :
C     Chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_MIE

C     Chemin et nom du fichier des angles utilisés
      CHARACTER*CTE_LENFIC2 FICANGLES_RES_MIE
  
C     Nom du fichier trace
      CHARACTER*CTE_LENFIC2 FICANGLOG
 
      
C---- Variables liées au modèle d'aérosols   

C     Epaisseur optique des aérosols à la longueur d'onde de référence (non tronquée).
      DOUBLE PRECISION AOT_REF

C     Longueur d'onde de référence
      DOUBLE PRECISION WAREF_AOT


      !Cas mono-mode aérosols : 
      !----------------------       
C     Partie réelle de l'indice de réfraction 		      		
      DOUBLE PRECISION RN_WA

C     Partie imaginaire de l'indice de réfraction
      DOUBLE PRECISION IN_WA

C     Partie réelle de l'indice de réfraction à la longueur d'onde de reference
      DOUBLE PRECISION RN_WAREF

C     Partie imaginaire de l'indice de réfraction à la longueur d'onde de référence    
      DOUBLE PRECISION IN_WAREF

C     Rayon modal de la LND
      DOUBLE PRECISION LND_RADIUS_MMD_AER
	
C     Variance de la LND
      DOUBLE PRECISION LND_LNVAR_MMD_AER

C     Pente de la loi de Junge.	
      DOUBLE PRECISION JD_SLOPE_MMD_AER

C     Rayon minimal de la loi de Junge.
      DOUBLE PRECISION JD_RMIN_MMD_AER

C     Rayon maximal de la loi de Junge.
      DOUBLE PRECISION JD_RMAX_MMD_AER
      
C     Noms généraux des paramètres de granulométrie pour modèle mono-mode   
      DOUBLE PRECISION VARGRANU1_MMD_AER
      DOUBLE PRECISION VARGRANU2_MMD_AER
      DOUBLE PRECISION VARGRANU3_MMD_AER


      !Proportion volumique des constituants WMO
      !-----------------------------------------     
C     Valeur utilisateur pour particule DL du modèle WMO
      DOUBLE PRECISION C_WMO_DL

C     Valeur utilisateur pour particule WS du modèle WMO
      DOUBLE PRECISION C_WMO_WS

C     Valeur utilisateur pour particule OC du modèle WMO
      DOUBLE PRECISION C_WMO_OC

C     Valeur utilisateur pour particule SO du modèle WMO
      DOUBLE PRECISION C_WMO_SO

C     Pourcentage d'humidité relative des modèles de Shettle & Fenn
      DOUBLE PRECISION RH

    
      !Paramètres pour un modèle bimodal de LND
      !-----------------------------------------    
C     Valeur utilisateur de la concentration volumique de la "LND coarse"
      DOUBLE PRECISION USER_CV_COARSE
   
C     Valeur utilisateur de la concentration volumique de la "LND fine"
      DOUBLE PRECISION USER_CV_FINE

C     Rapport des épaisseurs optique AOT_coarse / AOT_total pour la longueur d'onde de référence
      DOUBLE PRECISION RTAUct_WAREF
 
C     Valeur lue pour RN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MRWA

C     Valeur lue pour IN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MIWA

C     Valeur lue pour RN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MRWAREF

C     Valeur lue pour RN_LND_WAREF(2) : fine mode 
      DOUBLE PRECISION BMD_CM_MIWAREF

C     Valeur lue pour RMODAL_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_RMODAL

C     Valeur lue pour VAR_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_VAR

C     Valeur lue pour RN_LND(2) : fine mode  
      DOUBLE PRECISION BMD_FM_MRWA

C     Valeur lue pour IN_LND(2) : fine mode   					 
      DOUBLE PRECISION BMD_FM_MIWA

C     Valeur lue pour IN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MRWAREF

C     Valeur lue pour IN_LND_WAREF(2) : fine mode 
      DOUBLE PRECISION BMD_FM_MIWAREF 
   
C     Valeur lue pour RMODAL_LND(2) : fine mode   
      DOUBLE PRECISION BMD_FM_RMODAL

C     Valeur lue pour VAR_LND(2) : fine mode 
      DOUBLE PRECISION BMD_FM_VAR

C     Option de troncature de la fonction de phase des aérosols (1 pour la réaliser)                            
      INTEGER*4 ITRONC_AER

C     Option de troncature (=ITRONC_AER en entrée d'appel de SOS_AEROSOLS 
C     mais peut basculer de 1 à 0 en sortie si la troncature est automatiquement annulée)
      INTEGER*4 ITRONC_AER_IN

C     Granulometrie mono ou multi-modale
      INTEGER*2 IMOD_AER

C     Numéro du modèle de granulométrie
      INTEGER*4 IGRANU

C     Type du modèle WMO
      INTEGER*2 IModele_WMO

C     Type du modèle Shettle & Fenn
      INTEGER*2 IModele_SF

C     Indice du choix de description du mélange de 2 LND
      INTEGER*2 MODE_PARAM_BILND

C     Nom du fichier TRACE des calculs de MIE       
      CHARACTER*CTE_LENFIC2 FICMIE_LOG

C     Nom du fichier TRACE des calculs de propriétés aérosols GRANU 
      CHARACTER*CTE_LENFIC2 FICGRANU_LOG

C     Nom du fichier résultat GRANU 
      CHARACTER*CTE_LENFIC2 FICGRANU

C     Chemin et nom du fichier utilisateur donnant la fonction de phase 
C     et les coefficients de diffusion et d'extinction
      CHARACTER*CTE_LENFIC2 FICEXTDATA_AER

C     Chemin et nom du fichier aérosols dénommé par l'utilisateur
      CHARACTER*CTE_LENFIC2 FICUSER_AER

C     Chemin et nom du fichier de définition d'un mélange de modèles d'aérosols par l'utilisateur  
      CHARACTER*CTE_LENFIC2 FICMIXTURE_AER  
	
 
			
C --- Variables liées à la définition du profil atmosphérique    

C     Epaisseur optique rayleigh
      DOUBLE PRECISION TR

C     Echelle de hauteur du profil moléculaire (km)
      DOUBLE PRECISION HR

C     Echelle de hauteur du profil d'aérosols (km)
      DOUBLE PRECISION HA

C     Altitudes minimale et maximale de la couche molécules + aérosols (km)
      DOUBLE PRECISION ZMIN
      DOUBLE PRECISION ZMAX
  
C     Pression atmosphérique de surface (mbar)
      DOUBLE PRECISION PSURF

C     Contenu intégré en vapeur d'eau (g/cm2)
      DOUBLE PRECISION H2O  

C     Contenu intégré en ozone (Dobson)
      DOUBLE PRECISION O3

C     Concentration du dioxyde de carbone en surface (ppmv)
      DOUBLE PRECISION CO2

C     Concentration du méthane en surface (ppmv)
      DOUBLE PRECISION CH4

C     Transmission directe pour la direction d'incidence solaire de l'atmosphère équivalente 
      DOUBLE PRECISION TDIR_TRONC

C     Transmission directe pour la direction d'incidence solaire de l'atmosphere vraie
      DOUBLE PRECISION TDIR_VRAI

C     Altitude (km)
      DOUBLE PRECISION ZALT

C     Epaisseur optique moléculaire, aérosols et d'absorption gazeuse à ZALT
      DOUBLE PRECISION TR_ZALT
      DOUBLE PRECISION TA_ZALT
      DOUBLE PRECISION TG_ZALT

C     Type du profil aérosol
      INTEGER*2 IPROFIL

C     Type du profil en absorption
      INTEGER*2 ABSPROFIL

C     Nombre de couches atmosphériques
      INTEGER*4 NT

C     Fichier profil en absorption fourni par utilisateur     
      CHARACTER*CTE_LENFIC2 FICABSPROFIL

C     Nom du fichier résultat PROFIL temporaire pour un ordre de décomposition CKD
      CHARACTER*CTE_LENFIC2 FICPROFIL_TMP

C     Chemin et nom du fichier Trace des profils par ordre CKD
      CHARACTER*CTE_LENFIC2 FICPROFIL_LOG

 
C --- Variables liées à la définition du choix des données de coefficients CKD
C     Pas du nombre d'onde (cm-1) dans le fichier de coeff CKD
      DOUBLE PRECISION  NUSTEP

            
C --- Variables liées à la définition du modèle de surface 

C     Type de surface
      INTEGER*2 ISURF 

C     Albédo de surface
      DOUBLE PRECISION RHO

C     Vitesse du vent (m/s
      DOUBLE PRECISION WIND 

C     Indice de réfraction du dioptre air/surface
      DOUBLE PRECISION SURF_IND

C     Paramètres du modèle de Roujean	   
      DOUBLE PRECISION K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN

C     Paramètres du modèle de  Nadal 
      DOUBLE PRECISION ALPHA_NADAL, BETA_NADAL

C     Paramètre C*exp(-NDVI) du modèle de Maignan
      DOUBLE PRECISION COEF_C_MAIGNAN

C     Nom du fichier trace SURFACE
      CHARACTER*CTE_LENFIC2 FICSURF_LOG

C     Nom d'un fichier contenant une matrice de réflexion existante
      CHARACTER*CTE_LENFIC2 FICSURF

  
C --- Variables liées aux calculs SOS

C     Type de sortie : -1 =>  Sortie standard: Luminance émergeante au TOA + descendante au niveau du sol.
C                       z =>  Luminance montante et descendante pour l'altitude z (0 <= z <= CTE_TOA_ALT).  
      DOUBLE PRECISION ZOUT
 
C     Ordre maximal de diffusion		
      INTEGER*4 IGMAX

C     Option de simulation sans polarisation : 0 : calculs sans polarisation
C                                              1 : calculs normaux avec polarisation
      INTEGER*4 IPOLAR

C     Option de mode d'application de la décomposition CKD
C          1 : Calcul le plus fin avec un calcul de luminance par ordre de la décomposition
C              puis recomposition selon la pondération CKD.
C          2 : Calcul pour une épaisseur optique d'absorption estimée par l'approche CKD
C              --> Un unique calcul de luminance
      INTEGER*4 IMODE_CKD_CALCUL

C     Indice de type de sortie : 1 =>  Plan de visée pour une valeur d'azimut.
C                                2 =>  Diagramme polaire en fonction de la direction d'azimut et de visée.
      INTEGER*4 ITRPHI        

C     Pas sur l'azimut en degrés (pour ITRPHI=2)			
      INTEGER*4 PAS_PHI

C     Azimut relatif en degrés (pour ITRPHI=1)     			
      DOUBLE PRECISION PHIOS

C     Facteur de dépolarisation moléculaire 
      DOUBLE PRECISION RON

C     Chemin et nom du fichier Trace des calculs de SOS_OS    
      CHARACTER*CTE_LENFIC2 FICSOS_LOG

C     Nom du fichier optionnel de sortie des transmissions
      CHARACTER*CTE_LENFIC2 FICTRANS

C     Nom du fichier optionnel de sortie des flux up TOA et down BOA
      CHARACTER*CTE_LENFIC2 FICFLUX

C     Nom du fichier résultat des OS (pour un traitement unitaire)
      CHARACTER*CTE_LENFIC2 FICSOS_RES_BIN

C     Nom du fichier temporaire resultat des OS (calcul aggrégé)
      CHARACTER*CTE_LENFIC2 FICOS_AGG_TMP



C --- Variables internes

C -   Pour les aérosols

C     Longueur d'onde de calcul des paramètres aérosols (WA_SIMU ou WAREF)
      DOUBLE PRECISION WA

C     Section efficace d'extinction des aérosols pour la longueur d'onde calculée (simu ou ref)
      DOUBLE PRECISION KMAT1

C     Section efficace d'extinction des aérosols pour la longueur d'onde de référence de l'ep. opt. aérosols				     
      DOUBLE PRECISION KMAT1_WAREF
	
C     Epaisseur optique des aérosols pour la longueur d'onde de simulation			     
      DOUBLE PRECISION TA

C     Albédo de simple diffusion des aérosols pour la longueur d'onde de simulation				     
      DOUBLE PRECISION PIZAER

C     Coefficient de troncature de la fonction de phase des aérosols				     
      DOUBLE PRECISION COEF_TRONCA


C -   Pour les surfaces

C     Répertoire résultats des calculs de BRDF/BPDF de surface
      CHARACTER*CTE_LENDIR  SURF_OUTDIR
      
C -   Pour les profils

C     Nombre d'onde (cm-1)
      DOUBLE PRECISION NU

C     Profil en absorption
      DOUBLE PRECISION ABS_USERPROFIL(CTE_ABS_NBLEV,CTE_ABS_NBCOL)

C     Altitudes des niveaux du profil d'absorption (ordre décroissant)
      DOUBLE PRECISION ALTABS(1:CTE_ABS_NBLEV)

C     Profil d'épaisseur optique d'absorption gazeuse
      DOUBLE PRECISION TAUABS(1:CTE_ABS_NBLEV)

C     Profil de transmission verticale (masse d'air fictive de 1)        
      DOUBLE PRECISION TRSCKD(1:CTE_ABS_NBLEV)

C     Nombre de particules/cm2 par espèce gazeuse (H2O CO2 O3 N2O CO CH4 O2 NO2) et par couche du profil atmosphérique
      DOUBLE PRECISION RO(CTE_NBABS,CTE_ABS_NBLEV)
 
C     Tableau identifiant les gaz absorbant pris en compte (H2O CO2 O3 N2O CO CH4 O2 NO2)   
C     Inclue les gaz avec uniquement une absorption par continuum (NO2).
      INTEGER*2 IABS(1:CTE_NBABS)

C     Indice des tableaux CKD (dépend de NU)
      INTEGER*4 LAMB1

C -   Coefficients CKD

C     Nombre d'exponentielles par longueur d'onde et espèce       
      INTEGER*4  NEXP(CTE_NBABS,CTE_CKD_NWVL_MAX)

C     Valeur des coefficients ai par longueur d'onde, espèce et ordre i de la décomposition en exponentielles     
      DOUBLE PRECISION  KDIS_AI(CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX) 
     
C     Valeur des coefficients ki par longueur d'onde, espèce, ordre i de la décomposition en exponentielles, 
C     pression et température                                    
      DOUBLE PRECISION  KDIS_KI(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                          CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX )
 
C     Valeur des coefficients ki de la vapeur d'eau par longueur d'onde,  ordre i de la décomposition en exponentielles,
C     concentration, pression et température
      DOUBLE PRECISION  KDIS_KI_H2O(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                              CTE_CKD_NC_MAX, CTE_CKD_NAI_MAX, 
     &                              CTE_CKD_NWVL_MAX )

   
C     Nombre de valeurs de température, pression dans les fichiers CKD   
      INTEGER NB_TEMP, NB_PRES

C     Tableau des valeurs de température lues
      DOUBLE PRECISION TAB_TEMP(CTE_CKD_NT_MAX)

C     Tableau des valeurs de pression lues 
      DOUBLE PRECISION TAB_PRES(CTE_CKD_NP_MAX)

C     Nombre de valeurs de concentrations H2O dans le fichier CKD H2O  
      INTEGER NB_CONC_H2O 

C     Tableau des valeurs de concentrations H2O lues dans le fichier CKD H2O (ppmv) 
      DOUBLE PRECISION TAB_CONC_H2O(CTE_CKD_NC_MAX) 

C     Indice d'exponentielles      
      INTEGER*4 IK1,IK2,IK3,IK4,IK5,IK6,IK7,IK8

C     Coefficient de pondération des luminances
      DOUBLE PRECISION  AIK

C     Sommation des coefficients de pondération des luminances
      DOUBLE PRECISION sumAIK
   
C -   Pour les OS

C     Ordre des développements des fonctions de phase en Polynomes de Legendre. 
      INTEGER*4 OS_NB

C     Albédo de simple diffusion 
      DOUBLE PRECISION PIZ

C     Albédo de simple diffusion après troncature
      DOUBLE PRECISION PIZTR

C     Coefficient de troncature de la fonction de phase aérosols
      DOUBLE PRECISION COEFTRONCA

C     Nombre d'angles effectifs des simulations
      INTEGER*4 LUM_NBMU

C     Cosinus des angles de Gauss (et utilisateur)
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  

C     Poids des angles de Gauss (et utilisateur)
      DOUBLE PRECISION GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
				
C     Paramètres alpha(k), beta(k), gamma(k) et zeta(k) de la matrice de phase aérosols
      DOUBLE PRECISION AER_ALPHA(0:CTE_OS_NB_MAX)
      DOUBLE PRECISION AER_BETA(0:CTE_OS_NB_MAX)
      DOUBLE PRECISION AER_GAMMA(0:CTE_OS_NB_MAX)
      DOUBLE PRECISION AER_ZETA(0:CTE_OS_NB_MAX)

C     Flux diffus éclairant le sol et émergeant au TOA, normalisé par le flux solaire incident au TOA
      DOUBLE PRECISION EMOINS
      DOUBLE PRECISION EPLUS

C     Transmission diffuse pour la direction d'incidence solaire de l'atmosphère équivalente
      DOUBLE PRECISION TDIFMUS

C     Transmission diffuse pour les angles de Gauss de l'atmosphère équivalente				
      DOUBLE PRECISION TDIFMUG(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  

C     Epaisseur optique totale au sol (tronquée)
      DOUBLE PRECISION TTOT_TRONC

C     Epaisseur optique totale au sol (non tronquée).
      DOUBLE PRECISION TTOT_VRAI

C     Epaisseur optique (tronquée) au niveau de sortie des résultats pour le champ montant 
      DOUBLE PRECISION TAUOUT

C     Variables précédentes de flux et épaisseurs optiques calculées 
C     pour un octuor ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8 de décomposition CKD)
      DOUBLE PRECISION EMOINS_TMP
      DOUBLE PRECISION EPLUS_TMP
      DOUBLE PRECISION TDIFMUS_TMP		
      DOUBLE PRECISION TDIFMUG_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      DOUBLE PRECISION TTOT_TRONC_TMP
      DOUBLE PRECISION TTOT_VRAI_TMP
      DOUBLE PRECISION TAUOUT_TMP

C     Numéro de mu de Gauss pour l'angle zénithal solaire
C     angle complémentaire : mus = -rmu(N0)
C     N0 > 0 => pour en tenir compte			
      INTEGER*4 N0
   
C     Indice précisant si les simulations introduisent la réflexion de Fresnel sur mer plate :
C         1 => pour un calcul de la réflexion de Fresnel.
C         Sinon => pas de prise en compte.    
      INTEGER*4 IFRESNEL

C     Indice précisant si les simulations utilisent un modèle de réflexion BRDF ou BPDF :
C         1 => pour calcul de BRDF ou BPDF (il faut utiliser N0 > 0).
C         Sinon => pas de simulation de BRDF ou BPDF. 	
      INTEGER*4 IMAT_SURF

C     Indice précisant si les simulations introduisent le Glitter :
C         1 => pour calcul du Glitter (il faut alors utiliser N0 > 0).
C         Sinon => pas de simulation du Glitter.  
      INTEGER*4 IGLI
          
C     Indice précisant si les simulations utilisent le modèle de réflexion de Roujean :
C         1 => pour calcul de BRDF Roujean (il faut utiliser N0 > 0).
C         Sinon => pas de simulation de BRDF Roujean. 
      INTEGER*4 IROUJEAN

C     Indice précisant si les simulations utilisent le modèle de réflexion de Rondeaux - Herman :
C         1 => pour calcul de BPDF RH (il faut utiliser N0 > 0).
C         Sinon => pas de simulation de BPDF RH. 			 
      INTEGER*4 IRONDEAUX
			 
C     Indice précisant si les simulations utilisent le modèle de réflexion de Breon :
C         1 => pour calcul de BPDF Breon (il faut utiliser N0 > 0).
C         Sinon => pas de simulation de BPDF Breon. 			
      INTEGER*4 IBREON

C     Indice précisant si les simulations utilisent le modèle de réflexion de Nadal :
C         1 => pour calcul de BPDF Nadal (il faut utiliser N0 > 0).
C         Sinon => pas de simulation de BPDF Nadal. 
      INTEGER*4 INADAL

C     Indice précisant si les simulations utilisent le modèle de réflexion de Maignan :
C         1 => pour calcul de BPDF Nadal (il faut utiliser N0 > 0).
C         Sinon => pas de simulation de BPDF Maignan. 
      INTEGER*4 IMAIGNAN

C     Indices des tableaux d'angles pour lesquels il doit y avoir une sortie dans le fichier résultat des angles utilisateur 
C     ==> 1 : pour une sortie requise
C     Indices à partir de 0 pour passage à Python.
      INTEGER*4 IND_ANGOUT_FIN(0:CTE_OS_NBMU_MAX)  
				
C     Longueur de la chaine de caractères FICSOS_RES_BIN
      INTEGER*4 LFICOS

C     Nom du fichier temporaire résultat des OS		 
      CHARACTER*CTE_LENFIC2 FICOS_TMP
	
C     Table des valeurs d'azimuts relatifs du champ de sortie
      DOUBLE PRECISION PHI_FIN(0:360)

C     Table des valeurs d'angles zénithaux de visée du champ de sortie
      DOUBLE PRECISION THETA_FIN(0:CTE_OS_NBMU_MAX)

C     Tables des angles de diffusion pour les champs montant et descendant de sortie 
      DOUBLE PRECISION SCA_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)
      DOUBLE PRECISION SCA_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)     

C     Tables des valeurs I,Q et U pour les champs montant et descendant de sortie    
      DOUBLE PRECISION I_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)        
      DOUBLE PRECISION Q_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)       
      DOUBLE PRECISION U_UP_FIN(0:360,0:CTE_OS_NBMU_MAX) 
      DOUBLE PRECISION I_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)     
      DOUBLE PRECISION Q_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)      
      DOUBLE PRECISION U_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) 

C     Tables des angles de polarisation pour les champs montant et descendant de sortie    
      DOUBLE PRECISION POL_ANG_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)       
      DOUBLE PRECISION POL_ANG_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)  

C     Tables des taux de polarisation pour les champs montant et descendant de sortie      
      DOUBLE PRECISION POL_RATE_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)  
      DOUBLE PRECISION POL_RATE_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)  

C     Table des valeurs de luminance polarisée pour les champs montant et descendant de sortie         
      DOUBLE PRECISION L_POL_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)       
      DOUBLE PRECISION L_POL_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)

C     Flux diffus descendant au niveau de la surface      
      DOUBLE PRECISION FLUX_DIFF_DOWN

C     Flux total descendant au niveau de la surface           
      DOUBLE PRECISION FLUX_DOWN
  

    
C =====================================================================
C --- Définition des variables I/O pour usage du binding Python
C =====================================================================

Cf2py intent(in) resroot, ficmain_log,wa_simu, nbmu_gauss_lum, ficangles_user_lum, tetas
Cf2py intent(in)  ficangles_res_lum, nbmu_gauss_mie,ficangles_user_mie, ficangles_res_mie 
Cf2py intent(in) ficanglog, waref_aot, aot_ref, itronc_aer, ficgranu_log, ficmie_log, dir_mie
Cf2py intent(in)  ficgranu,imod_aer,rn_wa, in_wa, rn_waref, in_waref, igranu, lnd_radius_mmd_aer
Cf2py intent(in)  lnd_lnvar_mmd_aer, jd_slope_mmd_aer, jd_rmin_mmd_aer, jd_rmax_mmd_aer
Cf2py intent(in) imodel_wmo, c_wmo_dl, c_wmo_ws, c_wmo_oc, c_wmo_so 
Cf2py intent(in) imodele_sf, rh, mode_param_bilnd, user_cv_coarse, user_cv_fine 
Cf2py intent(in) rtauct_waref,bmd_cm_mrwa, bmd_cm_miwa, bmd_cm_mrwaref, bmd_cm_miwaref 
Cf2py intent(in) bmd_cm_rmodal, bmd_cm_var, bmd_fm_mrwa, bmd_fm_miwa, bmd_fm_mrwaref
Cf2py intent(in) bmd_fm_miwaref, bmd_fm_rmodal, bmd_fm_var, ficextdata_aer, ficmixture_aer 
Cf2py intent(in) ficuser_aer, ficprofil_log, tr, hr, ha, iprofil, zmin, zmax,psurf, h2o
Cf2py intent(in)  o3, co2, ch4, absprofil, ficabsprofil, nustep, isurf, dir_surf, ficsurf_log
Cf2py intent(in)  surf_ind, wind, k0_roujean, k1_roujean, k2_roujean, alpha_nada, beta_nada
Cf2py intent(in) coef_c_maignan, rho, ficsurf, ficsos_log, ficsos_res_bin, fic_trans 
Cf2py intent(in) fic_flux, zout, igmax, ipolar, itrphi, phios, pas_phi, imode_ckd_calcul, ier, trace

Cf2py intent(out) lum_nbmu, ind_angout_fin, phi_fin, theta_fin, sca_up_fin, i_up_fin, q_up_fin, u_up_fin, pol_ang_up_fin, pol_rate_up_fin, l_pol_up_fin
Cf2py intent(out) sca_down_fin, i_down_fin, q_down_fin, u_down_fin, pol_ang_down_fin, pol_rate_down_fin, l_pol_down_fin
Cf2py intent(out) tdir_vrai, flux_diff_down, flux_down, eplus
Cf2py intent(out) coef_tronca

                        
      
C =========================================
C --- Initialisations 
C =========================================

C --- Initialisation de l'indicateur d'erreur de traitement
C----------------------------------------------------------
      IER=0
 
C     Répertoire de sortie
      DIR_ROUJ="NO_ROUJEAN_DIR"

C     Répertoire de sortie
      SURF_OUTDIR="NO_OUTPUT_DIR"


C     Profil absorption en fonction de l'altitude
      DO I=1,CTE_ABS_NBLEV
         ALTABS(I)=CTE_NOT_DEFINED_VALUE_DBLE
         TAUABS(I)=CTE_NOT_DEFINED_VALUE_DBLE
      ENDDO
 
C     Nombre de coefficients pour chaque gaz absorbant et longueur d'onde
C     Par défaut, on le fixe à 1 => cas sans absorption ou continuum seul
      DO J=1,CTE_CKD_NWVL_MAX
         DO I=1,CTE_NBABS
            NEXP(I,J)=1
         ENDDO
      ENDDO

C     Poids de chaque exponentielle pour chaque gaz absorbant et longueur d'onde
C     Par défaut, on le fixe à 1 => cas sans absorption ou continuum seul
      DO K=1,CTE_CKD_NWVL_MAX	
         DO J=1,CTE_NBABS
            DO I=1,CTE_CKD_NAI_MAX
               KDIS_AI(I,J,K)=1.		
            ENDDO
         ENDDO
      ENDDO

C     Somme des coefficients AIK par défaut
      sumAIK=1.D+00

C     Indice correspondant au nombre d'onde / 10
      LAMB1=1
            
C --- Initialisation des valeurs par défaut des paramètres SOS

C     Facteur de dépolarisation moléculaire.
      RON=CTE_MDF

C     Init des épaisseurs optiques pour aggrégation
      TTOT_VRAI=0.
      TTOT_TRONC=0.
      TAUOUT=0.
      
C     Init des transmissions diffuses     
      TDIFMUS=0.
      DO I=-CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX
         TDIFMUG(I)=0.
      ENDDO	  
      EMOINS=0.
      EPLUS=0.


C =========================================================
C --- Définition de l'arborescence de travail
C =========================================================           				

C ---- Controle de la définition de l'espace de travail (obligatoire)  
      IF (RESROOT.EQ."UNDEFINED_REPERTORY") GOTO 1000
   
      LDIR = INDEX(RESROOT,' ') - 1
      IF((LDIR.LE.0).OR.(LDIR.GT.CTE_LENDIR)) GOTO 1100

C ---- Controle de la définition du répertoire de résultats MIE 
C     (obligatoire si AOT non nulle et si on n'utilise pas un fichier de 
C      données externes : fonctions de phase ou directement les coefficients Greek)

C     Définition de l'ep opt aérosols à la longueur d'onde de référence : obligatoire  
      IF (AOT_REF.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2301  	

      IF (AOT_REF.GT.0.D+00) THEN
		
C        Définition du type de modèle d'aérosols : obligatoire  
         IF ((IMOD_AER.EQ.CTE_NOT_DEFINED_VALUE_INT)
     &        .AND.(FICUSER_AER.EQ."NO_USER_AEROSOLS" )) GOTO 2304

C        Définition de la longueur d'onde de référence pour l'ep opt aérosols : obligatoire  
         IF (WAREF_AOT.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2302  

         IF ((IMOD_AER.NE.4.)
     &        .AND.(FICUSER_AER.EQ.'NO_USER_AEROSOLS')) THEN

            IF (DIR_MIE.EQ."UNDEFINED_REPERTORY") GOTO 1001

            !Création de l'arborescence de résultats aérosols
            LDIRMIE = INDEX(DIR_MIE,' ') - 1
            IF((LDIRMIE.LE.0).OR.(LDIRMIE.GT.CTE_LENDIR)) GOTO 1101
            COMMAND='mkdir -p '//DIR_MIE(1:LDIRMIE)
            LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1202
            CALL SYSTEM(COMMAND)
         ENDIF    
         
      ENDIF


C ---- Définition de l'arborescence de résultats de surface
      IF((ISURF.NE.0).AND.(ISURF.NE.2)) THEN

         !Controle de la définition du répertoire de résultats DIR_SURF (requis si BRDF ou BPDF calculée)
         IF (DIR_SURF.EQ."UNDEFINED_REPERTORY")  GOTO 1002

         LDIRSURF = INDEX(DIR_SURF,' ') - 1
         IF((LDIRSURF.LE.0).OR.(LDIRSURF.GT.CTE_LENDIR)) GOTO 1102
         
         !Définition des répertoires 
         !--Répertoire GLITTER
         LDIRGLIT= LDIRSURF + 8  	
         IF(LDIRGLIT.GT.CTE_LENDIR) GOTO 1103
         DIR_GLIT(1:LDIRGLIT) = DIR_SURF(1:LDIRSURF)//'/'//'GLITTER'
      
         !--Répertoire ROUJEAN
         LDIRROUJ= LDIRSURF + 8  	
         IF(LDIRROUJ.GT.CTE_LENDIR) GOTO 1104
         DIR_ROUJ(1:LDIRROUJ) = DIR_SURF(1:LDIRSURF)//'/'//'ROUJEAN'      	
      
         !--Répertoire Rondeaux & Herman
         LDIRRH= LDIRSURF + 3  	
         IF(LDIRRH.GT.CTE_LENDIR) GOTO 1105
         DIR_RH(1:LDIRRH) = DIR_SURF(1:LDIRSURF)//'/'//'RH'     

         !--Répertoire BREON
         LDIRBREON= LDIRSURF + 6  	
         IF(LDIRBREON.GT.CTE_LENDIR) GOTO 1106
         DIR_BREON(1:LDIRBREON) = DIR_SURF(1:LDIRSURF)//'/'//'BREON'     
  
         !--Répertoire NADAL
         LDIRNADAL= LDIRSURF + 6  	
         IF(LDIRNADAL.GT.CTE_LENDIR) GOTO 1107
         DIR_NADAL(1:LDIRNADAL) = DIR_SURF(1:LDIRSURF)//'/'//'NADAL'     

         !--Répertoire MAIGNAN
         LDIRMAIGNAN= LDIRSURF + 8  	
         IF(LDIRMAIGNAN.GT.CTE_LENDIR) GOTO 11071
         DIR_MAIGNAN(1:LDIRMAIGNAN) 
     &   = DIR_SURF(1:LDIRSURF)//'/'//'MAIGNAN'

         !Création de l'arborescence de résultats de surface
         !--Répertoire SURFACE
         COMMAND='mkdir -p '//DIR_SURF(1:LDIRSURF)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1203
         CALL SYSTEM(COMMAND)

         !--Répertoire GLITTER
         COMMAND='mkdir -p '//DIR_GLIT(1:LDIRGLIT)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1204
         CALL SYSTEM(COMMAND)
     
         !--Répertoire ROUJEAN
         COMMAND='mkdir -p '//DIR_ROUJ(1:LDIRROUJ)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1205
         CALL SYSTEM(COMMAND)
     
         !--Répertoire Rondeaux & Herman
         COMMAND='mkdir -p '//DIR_RH(1:LDIRRH)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1206
         CALL SYSTEM(COMMAND)

         !--Répertoire BREON
         COMMAND='mkdir -p '//DIR_BREON(1:LDIRBREON)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1207
         CALL SYSTEM(COMMAND)
  
         !--Répertoire NADAL
         COMMAND='mkdir -p '//DIR_NADAL(1:LDIRNADAL)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1208
         CALL SYSTEM(COMMAND)

         !--Répertoire MAIGNAN
         COMMAND='mkdir -p '//DIR_MAIGNAN(1:LDIRMAIGNAN)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 12081
         CALL SYSTEM(COMMAND)

      ENDIF
		
C ---- Définition de l'arborescence de résultats SOS

      LDIRSOS= LDIR + 4  	
      IF(LDIRSOS.GT.CTE_LENDIR) GOTO 1108
      DIR_SOS(1:LDIRSOS) = RESROOT(1:LDIR)//'/'//'SOS'	

C ---- Définition du répertoire de log

      LDIRLOG= LDIR + 4  	
      IF(LDIRLOG.GT.CTE_LENDIR) GOTO 1109
      DIR_LOG(1:LDIRLOG) = RESROOT(1:LDIR)//'/'//'LOG'	


C ---- Définition du répertoire temporaire des traitements

      LDIRTMP= LDIR + 4  	
      IF(LDIRTMP.GT.CTE_LENDIR) GOTO 1110
      DIR_TMP(1:LDIRTMP) = RESROOT(1:LDIR)//'/'//'TMP'



C =========================================================
C --- Définition des fichiers temporaires par ordre CKD
C =========================================================           				

C     Fichier temporaire de profil
      FICPROFIL_TMP = DIR_TMP(1:LDIRTMP)//'/'//"PROFIL_TMP"  

C     Fichier temporaire de résultats des OS (1 simulation)
      FICOS_TMP = DIR_TMP(1:LDIRTMP)//'/'//"FICOS_TMP"  

C     Fichier temporaire de résultats des OS (simulations aggrégées)
      FICOS_AGG_TMP = DIR_TMP(1:LDIRTMP)//'/'//"FICOS_AGG_TMP.bin"


C =========================================================
C --- Création de l'arborescence de travail
C =========================================================           				

C ---- Espace de travail

      COMMAND='mkdir -p '//RESROOT(1:LDIR)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1201
      CALL SYSTEM(COMMAND)
	
C ---- Création de l'arborescence de résultats SOS

      COMMAND='mkdir -p '//DIR_SOS(1:LDIRSOS)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1209
      CALL SYSTEM(COMMAND)

C ---- Création du répertoire de log

      COMMAND='mkdir -p '//DIR_LOG(1:LDIRLOG)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1210
      CALL SYSTEM(COMMAND)


C ---- Création du répertoire des fichiers temporaires

      COMMAND='mkdir -p '//DIR_TMP(1:LDIRTMP)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1211
      CALL SYSTEM(COMMAND)

      
C =========================================================
C --- Ouverture du fichier Trace
C =========================================================        

C     ---- Définition du fichier Trace principal
      IF (TRACE) THEN

         LFIC = INDEX(FICMAIN_LOG,' ') - 1
         IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 1301
         FICMAIN_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICMAIN_LOG(1:LFIC)
       
C ---- Création du fichier
         OPEN (IDLOG,FILE=FICMAIN_LOG,ERR=1302)

         WRITE (IDLOG,*,ERR=1303) "WORKING FOLDER DEFINITION"
         WRITE (IDLOG,*,ERR=1303) "-------------------------"
         WRITE (IDLOG,*,ERR=1303) "Working folder :", RESROOT(1:LDIR)
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &        "Surface result folder :", DIR_SURF(1:LDIRSURF)
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &        "Mie result folder :", DIR_MIE(1:LDIRMIE)
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &        "SOS result folder :", DIR_SOS(1:LDIRSOS)     
         WRITE (IDLOG,*,ERR=1303) " "
      ENDIF  ! Fin Trace 
  
	    
C =========================================================
C --- Tests d'existence et de cohérence des paramètres
C =========================================================     

C --------------------------------------------------------------        
C --- Controle des paramètres de simulation
C --------------------------------------------------------------  
C 
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &         "--> Starting the check of user parameters \n"
      ENDIF

C --- Controle de la définition de la longueur d'onde de simulation
	     
      IF (WA_SIMU.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2100           
      IF ((WA_SIMU.LT.CTE_WAMIN).OR.(WA_SIMU.GT.CTE_WAMAX))  GOTO 2101
      
      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "************** ",
     &         "Simulation parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "
          
C        Longueur d'onde de simulation
	  WRITE (IDLOG,*,ERR=1303) 
     &         "Simulation wavelength (microns): ",WA_SIMU
	  WRITE (IDLOG,*,ERR=1303) " "
	  
      ENDIF 
      
C --------------------------------------------------------------        
C --- Controle du passage des arguments de définition des angles
C --------------------------------------------------------------  
C  
C ---- Test d'existence et de validité des paramètres obligatoires

C     Angle solaire : obligatoire           
      IF (TETAS.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2200
      IF ((TETAS.LT.0.D-00).OR.(TETAS.GE.90.D-00)) GOTO 2201
         
	              
C ---- Affichage des paramètres utilisés dans le fichier de trace
	     
      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "************** ",
     &         "Angles parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "	  
	  WRITE (IDLOG,*,ERR=1303) "***** Mandatory Parameters "
	  
C        Angle solaire
	  WRITE (IDLOG,*,ERR=1303) "Solar Zenith Angle: ",TETAS
          
C        Fichier résultat pour les luminances       
	  WRITE (IDLOG,*,ERR=1303) "User's filename ",
     &         "for the radiance angles définition :"
         WRITE (IDLOG,*,ERR=1303) FICANGLES_RES_LUM
               
C        Fichier résultat pour les fonctions de phase
	  WRITE (IDLOG,*,ERR=1303) "User's filename ",
     &         "for the Mie angles définition :"
         WRITE (IDLOG,*,ERR=1303) FICANGLES_RES_MIE
	  	   
C        Fichier de log pour le calcul des angles
	  WRITE (IDLOG,*,ERR=1303) "User's log file: "
         WRITE (IDLOG,*,ERR=1303) FICANGLES_RES_MIE
               
	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "***** Optional Parameters "

C        Nombre d'angles de Gauss pour les simulations de luminances
	  IF (NBMU_GAUSS_LUM.NE.CTE_NOT_DEFINED_VALUE_INT) THEN 
	     WRITE (IDLOG,*,ERR=1303) "Nb Gauss Angles for ",
     &            "luminance calculation (user value):",NBMU_GAUSS_LUM
	  ELSE
	     WRITE (IDLOG,*,ERR=1303) "Nb Gauss Angles for ",
     &            "luminance calculation (default value)"
	  ENDIF
	 
C        Nombre d'angles de Gauss pour les calculs de mie
	  IF (NBMU_GAUSS_MIE.NE.CTE_NOT_DEFINED_VALUE_INT) THEN 
	     WRITE (IDLOG,*,ERR=1303) "Nb Gauss Angles for ",
     &            "Mie calculation (user value):",NBMU_GAUSS_MIE
	  ELSE
	     WRITE (IDLOG,*,ERR=1303) "Nb Gauss Angles for ",
     &            "Mie calculation (default value)"
	  ENDIF
          
C        Fichier des angles utilisateurs pour les calculs de luminance
	  IF (FICANGLES_USER_LUM.NE.'NO_USER_ANGLES') THEN 
	     WRITE (IDLOG,*,ERR=1303) "User's file of specific angles",
     &            " for luminance calculations (user value):",
     &            FICANGLES_USER_LUM
	  ENDIF		 
		 	 
C	  Fichier des angles utilisateurs pour les calculs de Mie
	  IF (FICANGLES_USER_MIE.NE.'NO_USER_ANGLES') THEN 
	     WRITE (IDLOG,*,ERR=1303) "User's file of specific angles",
     &            " for MIE calculations (user value):",
     &            FICANGLES_USER_MIE
	  ENDIF		  	 
          
	  WRITE (IDLOG,*,ERR=1303) " "
          
      ENDIF       ! Fin Trace

      
C --------------------------------------------------------------  
C --- Controle du passage des arguments de modèle d'aérosols
C --------------------------------------------------------------  
C
C ---- Test d'existence et de validité des paramètres obligatoires
   	
C     Cas AOT_REF > 0 :
      IF (( AOT_REF.GT.0.D+00).AND.
     & ( FICUSER_AER.EQ."NO_USER_AEROSOLS" )) THEN
      
C        Controle de la définition des paramètres aérosols
C        Modèle aérosol
	  IF ((IMOD_AER.LT.0).OR.(IMOD_AER.GT.5)) GOTO 2305

C        Option de suppression de la troncature
         IF ((ITRONC_AER.LT.0).OR.(ITRONC_AER.GT.1)) GOTO 23141
	 
C        ----------------------------------------------------
C        Cas d'un modèle d'aérosols mono-modal : IMOD_AER = 0  
C        ----------------------------------------------------
	  IF (IMOD_AER.EQ.0) THEN 

C            Parties réelles et imaginaires de l'indice de réfraction
	     IF ((RN_WA.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &          (IN_WA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2309

C           Validite de la partie imaginaire
	     IF (IN_WA.GT.0.D+00) GOTO 2310

C           Validite du type de granulométrie
            IF (IGRANU.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2311

	     IF ((IGRANU.LT.1).OR.(IGRANU.GT.2)) GOTO 2312	 

	     ! Contrôle de la récupération des paramètres du modèle de LND
	     IF (IGRANU.EQ.1) THEN
	         IF (LND_RADIUS_MMD_AER.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 23131
		  IF (LND_LNVAR_MMD_AER.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 23131

                VARGRANU1_MMD_AER = LND_RADIUS_MMD_AER
	         VARGRANU2_MMD_AER = LND_LNVAR_MMD_AER
	         VARGRANU3_MMD_AER = CTE_NOT_DEFINED_VALUE_DBLE
	     ENDIF     

	     ! Contrôle de la récupération des paramètres obligatoires du modèle de Junge
	     IF (IGRANU.EQ.2) THEN
	         IF (JD_SLOPE_MMD_AER.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 23132
		  IF (JD_RMIN_MMD_AER.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 23132

                VARGRANU1_MMD_AER = JD_RMIN_MMD_AER
	         VARGRANU2_MMD_AER = JD_SLOPE_MMD_AER
	         VARGRANU3_MMD_AER = JD_RMAX_MMD_AER
	     ENDIF	    

C	     Cas où la longueur d'onde de référence est différente de celle de simu
C	     ==> Controle de la récupération de la valeur de l'indice de réfraction
C	         pour la longueur d'onde de référence.
	      IF (WA_SIMU.NE.WAREF_AOT) THEN
	         IF((RN_WAREF.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &	            (IN_WAREF.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2314
	      ELSE
                RN_WAREF = RN_WA
                IN_WAREF = IN_WA
	     ENDIF 
	 		 	
         ENDIF
	 
C         ----------------------------------------------------
C         Cas d'un modèle d'aérosols WMO : IMOD_AER = 1  
C         ----------------------------------------------------	 
 	  
	  IF (IMOD_AER.EQ.1) THEN 	!Cas WMO	

	     ! Controle la validité du paramètre IModele
	     ! -----------------------------------------	
	     IF (IModele_WMO.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2315	

	     ! Controle de la récupération du choix de modèle WMO   
	     IF ((IModele_WMO.LT.1).OR.(IModele_WMO.GT.CTE_NB_Modeles))
     &         GOTO 2316	

	     ! Si modèle libre de l'utilisateur, controle de la récupération 
	     ! des proportions des 4 composants
	     IF (IModele_WMO.EQ.4) THEN     
	         IF ((C_WMO_DL.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &		      (C_WMO_WS.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &		      (C_WMO_OC.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &		      (C_WMO_SO.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2317
	     ENDIF !Cas modèle utilisateur WMO
	  
	  ENDIF !Fin cas IMOD_AER = 1 (WMO)  
	  
C        ----------------------------------------------------
C        Cas d'un modèle d'aérosols S&F : IMOD_AER = 2  
C        ----------------------------------------------------
        
	  IF (IMOD_AER.EQ.2) THEN !     cas S&F
         
	     ! Controle de la récupération du choix de modèle S&F 
	     IF (IModele_SF.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2318
	     
	     ! Controle de la récupération de la valeur de l'humidité relative
	     IF (RH.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2319

            !   Controle de la validite des parametres IMODELE et RH
            IF ((IModele_SF.LT.1).OR.(IModele_SF.GT.CTE_NB_Modeles))
     &         GOTO 2320        !Erreur choix modèle
            IF ((RH.LT.0.D-00).OR.(RH.GT.99.D-00)) GOTO 2321 !Erreur sur % humidité

	  ENDIF
         !Fin cas IMOD_AER = 2 (S&F)         

C        ----------------------------------------------------
C        Cas d'un modèle d'aérosols bimode LND : IMOD_AER = 3  
C        ----------------------------------------------------	  
         !Cas bi LND
	  IF (IMOD_AER.EQ.3) THEN

             !Controle de la récupération des paramètres specifiques au calcul bimodal LND
             IF (MODE_PARAM_BILND.EQ.
     &           CTE_NOT_DEFINED_VALUE_INT)  GOTO 2322

             !controle de la validite du paramètre MODE_PARAM_BILND
             IF ((MODE_PARAM_BILND.LT.1).OR.(MODE_PARAM_BILND.GT.2)) 
     &         GOTO 2323
  



	     
             !Récupération de l'information de la description du mélange
	      !    Cas 1 : Utilisation des concentrations volumiques de l'utilisateur
             !    --> Controle que les paramètres utilisateurs sont renseignés pour l'option
             IF (MODE_PARAM_BILND.EQ.1) THEN          
	         IF (USER_CV_COARSE.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2324
	         IF (USER_CV_FINE.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2325  
             ENDIF
    
            !    Cas 2 : Utilisation du rapport AOT_coarse / AOT_totale
            !    --> Controle que les paramètres sont tous renseignes pour l'option
            IF (MODE_PARAM_BILND.EQ.2) THEN 
	         IF (RTAUct_WAREF.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2326
            ENDIF
   
            !Récupération de l'information d'indice et granulométrie des 2 LND
            !-----------------------------------------------------------------         	 
            IF (BMD_CM_MRWA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2327
            IF (BMD_CM_MIWA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2327 
            IF (BMD_CM_RMODAL.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2327
            IF (BMD_CM_VAR.EQ.CTE_NOT_DEFINED_VALUE_DBLE)    GOTO 2327 
            IF (BMD_FM_MRWA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2328
            IF (BMD_FM_MIWA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2328 
            IF (BMD_FM_RMODAL.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2328
            IF (BMD_FM_VAR.EQ.CTE_NOT_DEFINED_VALUE_DBLE)    GOTO 2328            
     
            !Récupération de l'information d'indice des 2 LND pour WAREF 
            !-----------------------------------------------------------------     
            IF (MODE_PARAM_BILND.EQ.2) THEN 
	     
	        IF (WA_SIMU.NE.WAREF_AOT) THEN
	           IF (BMD_CM_MRWAREF.EQ.
     &                CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2329
	           IF (BMD_CM_MIWAREF.EQ.
     &                CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2329
	           IF (BMD_FM_MRWAREF.EQ.
     &                CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2329
	           IF (BMD_FM_MIWAREF.EQ.
     &                CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2329   

	        ELSE !Si waref=wa ==> On impose de récupérer l'info d'indice pour wa 
		      !                et on l'affecte à waref
		    BMD_CM_MRWAREF = BMD_CM_MRWA
		    BMD_CM_MIWAREF = BMD_CM_MIWA
		    BMD_FM_MRWAREF = BMD_FM_MRWA
		    BMD_FM_MIWAREF = BMD_FM_MIWA  

	        ENDIF ! Fin test si Wa != WaRef
	        
            ENDIF !Fin cas description du mélange bi LND par rapport d'AOT  
	           
	  ENDIF


C         ----------------------------------------------------
C         Cas d'un modèle d'aérosols fourni par l'utilisateur : IMOD_AER = 4  
C         ----------------------------------------------------	  
	  IF (IMOD_AER.EQ.4) THEN 	!Cas fonctions externes

             ! Controle de la récupération du nom du fichier de données externes
	     IF (FICEXTDATA_AER.EQ.
     &          "NO_USER_AEROSOLS_PHAZE_FCT") GOTO 2330

	     ! Controle que la longueur d'onde de référence pour les aérosols et
	     ! identique à celle de simulation du champs de luminance 
	     IF (WA_SIMU.NE.WAREF_AOT) GOTO 2331

	  ENDIF !Fin cas IMOD_AER = 4 (fonctions phases externes)

	  
C         ---------------------------------------------------------------
C         Cas d'un mélange de modèles d'aérosols défini par l'utilisateur: IMOD_AER = 5  
C         ---------------------------------------------------------------	 
 	  
	  IF (IMOD_AER.EQ.5) THEN 	!Cas mélange de modèles par l'utilisateur	
	  
	     ! Controle de la récupération du nom du fichier de définition du mélange
	     IF (FICMIXTURE_AER.EQ."NO_USER_AEROSOLS_MIXTURE") GOTO 2340
  
	  ENDIF !Fin cas IMOD_AER = 5 (Mélange utilisateur)  

	  	  	  	  	             
      ENDIF ! FIN cas test sur AOT > 0 et pas de fichier aérosols utilisateur



       
C     Cas AOT_REF > 0 et DEF_FICUSER_AER = 1 
C     => Utilisation du fichier utilisateur :
C     Cas d'erreur : le fichier résultats ne doit pas être défini
C                    les longueurs d'onde de simu et de référence doivent être identiques

      IF (( AOT_REF.GT.0.D+00)
     &      .AND.(FICUSER_AER.NE."NO_USER_AEROSOLS" )) THEN
          
	  ! Controle que la longueur d'onde de référence pour les aérosols est
	  ! identique à celle de simulation du champs de luminance 
         IF (WA_SIMU.NE.WAREF_AOT) GOTO 2350

	  ! Controle de la non définition par l'utilisateur du fichier résultat des calculs d'IOP
	   IF (FICGRANU.NE.CTE_DEFAULT_FICGRANU) GOTO 2351
	  
      ENDIF
      
C ---- Affichage des paramètres utilisés dans le fichier de trace

      IF (TRACE) THEN

         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) "************** ",
     &         "Aerosols parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "
	  
         WRITE (IDLOG,*,ERR=1303) "***** Mandatory Parameters "
	  
C        AOT de référence
         WRITE (IDLOG,*,ERR=1303) 
     &         "AOT for the reference wavelength:",AOT_REF

C	  Cas AOT > 0
         IF (( AOT_REF.GT.0.D+00)
     &         .AND.(FICUSER_AER.EQ."NO_USER_AEROSOLS")) THEN

C           Longueur d'onde de référence	
            WRITE (IDLOG,*,ERR=1303) 
     &            "Reference Wavelength (microns): ",WAREF_AOT
      
C           Type de model d'aérosols
            WRITE (IDLOG,*,ERR=1303) "Aerosol Model Type: ",IMOD_AER


C	     Cas modèle mono-mode	
            IF (IMOD_AER.EQ.0) THEN

                WRITE (IDLOG,*,ERR=1303)
     &               "Mono-modal aerosol Model"
		   
C         	  Partie réelle de l'indice de réfraction à WA_SIMU
                WRITE (IDLOG,*,ERR=1303)
     &               "Real part of the refractive index: ",RN_WA
     
C         	  Partie imaginaire de l'indice de réfraction à WA_SIMU
                WRITE (IDLOG,*,ERR=1303)
     &               "Imaginary part of the refractive index: ",IN_WA
     
C		  Type de granulométrie
                WRITE (IDLOG,*,ERR=1303)
     &               "Particle size distribution model: ",IGRANU

C		  Paramètres de granulométrie
                IF (IGRANU.EQ.1) THEN
                   WRITE (IDLOG,*,ERR=1303) "LND parameters: "
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Modal radius (microns): ", VARGRANU1_MMD_AER 
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Variance (microns): ", VARGRANU2_MMD_AER 
                ELSE
                   WRITE (IDLOG,*,ERR=1303) "JUNGE's law parameters: "
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Slope: ", VARGRANU1_MMD_AER 
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Min radius (microns): ", VARGRANU2_MMD_AER 
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Max radius (microns): ", VARGRANU3_MMD_AER 
                ENDIF

	     	  IF (WA_SIMU.NE.WAREF_AOT) THEN
	     
C                  Partie réelle de l'indice de réfraction à WA_REF
                   WRITE (IDLOG,*,ERR=1303)
     &                    "Real part of the refractive index",
     &                    " (ref wvl): ",RN_WAREF
     
C         	     Partie imaginaire de l'indice de réfraction à WA_SIMU
                   WRITE (IDLOG,*,ERR=1303)
     &                    "Imaginary part of the refractive index ",
     &                    "(ref wvl): ",IN_WAREF
                ENDIF
                  
            ENDIF ! Fin IMOD_AER.EQ.0
		
		
C	     Cas modèle WMO		
            IF (IMOD_AER.EQ.1) THEN

	  	  WRITE (IDLOG,*,ERR=1303) "WMO aerosol Model"
		
C         	  Modèle WMO
                WRITE (IDLOG,*,ERR=1303)
     &                 "WMO Model type: ",IModele_WMO
                  
	     	  IF (IModele_WMO.EQ.4) THEN
	     
C         	     Pourcentage volumique des particules DL du modèle WMO
                   WRITE (IDLOG,*,ERR=1303)
     &                    "DL particules ratio",C_WMO_DL
     
C         	     Pourcentage volumique des particules WS du modèle WMO
                   WRITE (IDLOG,*,ERR=1303)
     &                    "WS particules ratio",C_WMO_WS

C         	     Pourcentage volumique des particules SO du modèle WMO
                   WRITE (IDLOG,*,ERR=1303)
     &                    "SO particules ratio",C_WMO_SO

C         	     Pourcentage volumique des particules OC du modèle WMO
                   WRITE (IDLOG,*,ERR=1303)
     &                    "OC particules ratio",C_WMO_OC

		  ENDIF
     
            ENDIF            ! Fin IMOD_AER.EQ.1

C	     Cas modèle SF		
            IF (IMOD_AER.EQ.2) THEN

	  	  WRITE (IDLOG,*,ERR=1303) "S&F aerosol Model"
					
C               Modèle S&F
                WRITE (IDLOG,*,ERR=1303)
     &                 "S&F Model type: ",IModele_SF
          
C         	  Humidite relative
                WRITE (IDLOG,*,ERR=1303)
     &                 "Relative Humidity: ",RH

            ENDIF  ! Fin IMOD_AER.EQ.2
		

C	     Cas model bi-modal

            IF (IMOD_AER.EQ.3) THEN

	  	  WRITE (IDLOG,*,ERR=1303)"Bi-modal aerosol Model"
		
C         	  Mode de calcul bimodal
                WRITE (IDLOG,*,ERR=1303)
     &                 "Calculation mode: ",MODE_PARAM_BILND
          
		  IF (MODE_PARAM_BILND.EQ.1) THEN        

C	            concentration volumique LND coarse
	  	     WRITE (IDLOG,*,ERR=1303)
     &                 "User's volumic concentration of LND coarse: ",
     &                    USER_CV_COARSE

C		     concentration volumique LND fine
	  	     WRITE (IDLOG,*,ERR=1303)
     &                    "User's volumic concentration of LND fine: ",
     &                    USER_CV_FINE
 			  
            	  ENDIF

		  IF (MODE_PARAM_BILND.EQ.2) THEN        

C		     Rapport AOT coarse / AOT totale
	  	     WRITE (IDLOG,*,ERR=1303)
     &                    "User's ratio (coarse AOT / total AOT) : ",
     &                    RTAUct_WAREF
 			  
            	  ENDIF

C         	  Partie réelle de l'indice de réfraction du mode Coarse
                WRITE (IDLOG,*,ERR=1303)
     &                 "Real part of the refractive index ",
     &                 "of Coarse model: ",BMD_CM_MRWA
     
C         	  Partie imaginaire de l'indice de réfraction du mode Coarse
                WRITE (IDLOG,*,ERR=1303)
     &                 "Imaginary part of the refractive index ",
     &                 "of Coarse model: ",BMD_CM_MIWA

C         	  Rayon modal du mode Coarse
                WRITE (IDLOG,*,ERR=1303)
     &          "Coarse modal radius of Coarse model: ",BMD_CM_RMODAL
     
C         	  Variance de LND pour le mode Coarse
                WRITE (IDLOG,*,ERR=1303)
     &          "variance LND of the Coarse Mode ",BMD_CM_VAR

C         	  Partie réelle de l'indice de réfraction du mode Fine
                 WRITE (IDLOG,*,ERR=1303)
     &                 "Real part of the refractive index",
     &                 " of Fine model: ",BMD_FM_MRWA
     
C         	  Partie imaginaire de l'indice de réfraction du mode Fine
                WRITE (IDLOG,*,ERR=1303)
     &                 "Imaginary part of the refractive index ",
     &                 "of Fine model: ",BMD_FM_MIWA

C         	  Rayon modal du mode Fine
                WRITE (IDLOG,*,ERR=1303)
     &          "Coarse modal radius of Fine model: ",BMD_FM_RMODAL
     
C         	  Variance de LND pour le mode Fine
                WRITE (IDLOG,*,ERR=1303)
     &          "variance LND of the Fine Mode ",BMD_FM_VAR
     
                IF (WA_SIMU.NE.WAREF_AOT) THEN

C                  Partie réelle de l'indice de réfraction du mode Coarse
                   WRITE (IDLOG,*,ERR=1303)
     &		     "Real part of the refractive index",
     &		     "of Coarse model (ref wvl): ",BMD_CM_MRWAREF
     
C         	     Partie imaginaire de l'indice de réfraction du mode Coarse
                   WRITE (IDLOG,*,ERR=1303)
     &		     "Imaginary part of the refractive index ",
     &		     "of Coarse model (ref wvl): ",BMD_CM_MIWAREF


C         	     Partie réelle de l'indice de réfraction du mode Fine
                   WRITE (IDLOG,*,ERR=1303)
     &		     "Real part of the refractive index",
     &		     "of Fine model (ref wvl): ",BMD_FM_MRWAREF
     
C         	     Partie imaginaire de l'indice de réfraction du mode Fine
                   WRITE (IDLOG,*,ERR=1303)
     &		     "Imaginary part of the refractive index ",
     &		     "of Fine model (ref wvl): ",BMD_FM_MIWAREF

		  ENDIF
			
	     ENDIF  ! Fin IMOD_AER.EQ.3
	 

C	     Cas données utilisateur		
	     IF (IMOD_AER.EQ.4) THEN

	  	  WRITE (IDLOG,*,ERR=1303)
     &		  "User gives aerosol radiave properties"
					
C         	  Fichier utilisateur
                WRITE (IDLOG,*,ERR=1303)
     &		  "User's filename: ",FICEXTDATA_AER
     
	     ENDIF    



C	     Cas modèle de mélange défini par l'utilisateur	
	     IF (IMOD_AER.EQ.5) THEN

	  	  WRITE (IDLOG,*,ERR=1303)
     &		  "User defines the mixture of aerosol models"
					
C         	  Fichier utilisateur
                WRITE (IDLOG,*,ERR=1303)
     &		  "User's filename : ",FICMIXTURE_AER
     
	     ENDIF   

		
	  ENDIF ! fin test AOT > 0 et DEF_FICUSER_AER = 0
	 

C	  Cas d'un fichier aérosol fourni par l'utilisateur
	  IF (( AOT_REF.GT.0.D+00)
     &       .AND.(FICUSER_AER.NE.'NO_USER_AEROSOLS')) THEN

	     WRITE (IDLOG,*,ERR=1303) " "	  
	     WRITE (IDLOG,*,ERR=1303) 
     &      "***** Use of an aerosol User File : ",FICUSER_AER
	     WRITE (IDLOG,*,ERR=1303) "=> No Aerosol Processing"
	  
	  ELSE	  

C         Fichier de résultats 	

	     WRITE (IDLOG,*,ERR=1303) " "	  
            WRITE (IDLOG,*,ERR=1303) "***** Optional Parameters "
	  
	     IF (FICGRANU.NE.CTE_DEFAULT_FICGRANU) THEN 
	     	  WRITE (IDLOG,*,ERR=1303) "User's filename of the ",
     &       	  "result file (user value):", FICGRANU
	     ELSE
	     	  WRITE (IDLOG,*,ERR=1303) "User's filename of the ",
     &       	  "result file (default value):", FICGRANU  
	     ENDIF		  	 
 
   	  ENDIF  ! fin test utilisation fichier utilisateur

	  WRITE (IDLOG,*,ERR=1303) " "	  
	  WRITE (IDLOG,*,ERR=1303) "***** Optional Parameters "	  
	  WRITE (IDLOG,*,ERR=1303) "Phase function troncature option:",
     &       	    ITRONC_AER
      ENDIF  ! fin écriture trace
          

      
C ------------------------------------------------------------------        
C --- Controle du passage des arguments pour les calculs de surface
C ------------------------------------------------------------------  
C  
C ---- Test d'existence et de validité des paramètres obligatoires

C     Argument RHO (albédo de surface) obligatoire, controle de sa valeur         
      IF (RHO.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2401        
      IF (RHO.LT.0.D-00) GOTO 2402
 
C     Argument ISURF (type de surface) obligatoire, controle de sa valeur  
      IF (ISURF.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2403
     
C     Argument ISURF obligatoire, controle de sa valeur
      IF((ISURF.NE.0).AND.(ISURF.NE.1).AND.(ISURF.NE.2)
     &     .AND.(ISURF.NE.3).AND.(ISURF.NE.4).AND.(ISURF.NE.5)
     &     .AND.(ISURF.NE.6).AND.(ISURF.NE.7)) GOTO 2404
            
C     Argument SURF_IND requis si ISURF = 1, 2, 4, 5, 6 ou 7      
      IF ( (SURF_IND.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &     .AND.((ISURF.EQ.1).OR.(ISURF.EQ.2)
     &     .OR.(ISURF.EQ.4).OR.(ISURF.EQ.5).OR.(ISURF.EQ.6)
     &     .OR.(ISURF.EQ.7)) ) GOTO 2405  

C     Cas du GLITTER => Argument WIND obligatoire              
      IF (ISURF.EQ.1) THEN
         IF (WIND.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2406
         IF (WIND.LT.0.D-00) GOTO 24061
      ENDIF	            
      
C     Cas du modèle BRDF de Roujean => paramètres K0, K1 et K2 obligatoires
      IF (ISURF.GE.3) THEN
         IF ( (K0_ROUJEAN.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &        .OR.(K1_ROUJEAN.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &        .OR.(K2_ROUJEAN.EQ.CTE_NOT_DEFINED_VALUE_DBLE) ) THEN
             GOTO 2407
         ENDIF 	 
      ENDIF	                     

C     Cas BPDF de NADAL => paramètres ALPHA_NADAL et BETA_NADAL doivent être définis    
      IF (ISURF.EQ.6) THEN

         !Suppression de la possibilité d'utiliser le modèle Nadal (non applicable au formalisme SOS)
         IF (TRACE) THEN
            WRITE (IDLOG,*,ERR=1303) " "
            WRITE (IDLOG,*,ERR=1303) "WARNING:"
            WRITE (IDLOG,*,ERR=1303) 
     &           "The Nadal's BPDF model is not supported " 
            WRITE (IDLOG,*,ERR=1303) "==> Select another surface model"
         ENDIF 

         WRITE(6,*) "  "
         WRITE(6,*) "  "
         WRITE(6,*) "  The Nadal's BPDF model is not supported"
         WRITE(6,*) "  ==> Select another surface model"
         GOTO 9999
        
         IF ((ALPHA_NADAL.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &   .OR.(BETA_NADAL.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2408

C        Validité des coefficients ALPHA_NADAL et BETA_NADAL		
         IF((ALPHA_NADAL.LT.0.).OR.(ALPHA_NADAL.GE.10.)) GOTO 2409
	  IF((BETA_NADAL.LT.0.).OR.(BETA_NADAL.GE.100.)) GOTO 2410
	 
      ENDIF 


C     Cas du modèle de BPDF de Maignan => Argument COEF_C_MAIGNAN obligatoire              
      IF (ISURF.EQ.7) THEN
         IF (COEF_C_MAIGNAN.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2411        
      ENDIF


C ---- Affichage des paramètres utilisés dans le fichier de trace

      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "************** ",
     &       "Surface parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "

C        Albédo de surface
	  WRITE (IDLOG,*,ERR=1303) "Surface Albedo: ",RHO

C        Modèle de surface
	  WRITE (IDLOG,*,ERR=1303) "Surface Type: ",ISURF
	  
C	  Indice de réfraction
         IF (SURF_IND.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
	      WRITE (IDLOG,*,ERR=1303) "Refractive Index ",SURF_IND
	  ENDIF
	  
C	  Cas GLITTER
      	  IF (ISURF.EQ.1) THEN
	      WRITE (IDLOG,*,ERR=1303) "GLITTER case"
	      WRITE (IDLOG,*,ERR=1303) "Wind: ",WIND
      	  ENDIF	            

      	  IF (ISURF.EQ.2) THEN
	      WRITE (IDLOG,*,ERR=1303)
     &       "Fresnel's reflexion on flat water"
      	  ENDIF
  
C	  Cas BRDF de Roujean
      	  IF (ISURF.GE.3) THEN
	      WRITE (IDLOG,*,ERR=1303) "Roujean BRDF parameters"
	      WRITE (IDLOG,*,ERR=1303) "K0: ",K0_ROUJEAN           
	      WRITE (IDLOG,*,ERR=1303) "K1: ",K1_ROUJEAN 
	      WRITE (IDLOG,*,ERR=1303) "K2: ",K2_ROUJEAN 
         ENDIF	                     

      	  IF (ISURF.EQ.4) THEN
	      WRITE (IDLOG,*,ERR=1303) "Rondeaux's BPDF model"
      	  ENDIF

      	  IF (ISURF.EQ.5) THEN
	      WRITE (IDLOG,*,ERR=1303) "Breon's BPDF model"
      	  ENDIF

C     	  Cas BPDF de NADAL => paramètres ALPHA_NADAL et BETA_NADAL    
      	  IF (ISURF.EQ.6) THEN
	      WRITE (IDLOG,*,ERR=1303) "Nadal's BPDF model"
             WRITE (IDLOG,*,ERR=1303) "Parameters:"		
             WRITE (IDLOG,*,ERR=1303) "  ALPHA_NADAL: ",ALPHA_NADAL
             WRITE (IDLOG,*,ERR=1303) "  BETA_NADAL: ",BETA_NADAL	 
      	  ENDIF 	 

      	  IF (ISURF.EQ.7) THEN
	      WRITE (IDLOG,*,ERR=1303) "Maignan's BPDF model"
             WRITE (IDLOG,*,ERR=1303) "  C coef: ",COEF_C_MAIGNAN 
      	  ENDIF
	  
      ENDIF  ! Fin de trace

C ------------------------------------------------------------------        
C --- Controle du passage des arguments pour les calculs du profil
C ------------------------------------------------------------------  
C  
C ---- Test d'existence et de validité des paramètres obligatoires

C     Epaisseur optique Rayleigh      
      IF ((TR.NE.CTE_NOT_DEFINED_VALUE_DBLE)
     &     .AND.(TR.LT.0.D-00)) GOTO 2502
      
C     Echelle de Hauteur Rayleigh    

      IF ((TR.NE.CTE_NOT_DEFINED_VALUE_DBLE)
     &     .AND.(TR.GT.0.)) THEN
         IF (HR.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2503  
         IF (HR.LE.0.D-00) GOTO 2504
      ENDIF
   
C     Type de profil pour diffusion aérosols     
      IF (IPROFIL.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2505   
      IF ((IPROFIL.LT.1).OR.(IPROFIL.GT.CTE_NB_PROFIL)) GOTO 2506

C     Echelle de Hauteur aérosols si IPROFIL=1
      IF (IPROFIL.EQ.1) THEN
         IF ((HA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &        .AND.(AOT_REF.GT.0.)) GOTO 2507
         IF ((HA.LE.0.D-00).AND.(AOT_REF.GE.0.))GOTO 2508	  
      ENDIF
       
C     Altitudes min et max du profil si IPROFIL=2      
      IF (IPROFIL.EQ.2) THEN
         IF ((ZMIN.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &      .OR.(ZMAX.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2509
      ENDIF

C     Type de profil d'absorption gazeuse     
            
      IF (ABSPROFIL.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2510      
      IF ((ABSPROFIL.LT.0).OR.(ABSPROFIL.GT.7)) GOTO 2511

C     Fichier utilisateur du profil d'absorption si ABSPROFIL = 0
	  
      IF (ABSPROFIL.EQ.0) THEN
      	  IF (FICABSPROFIL.EQ."NO_USER_ABS_PROFILE_FILE") GOTO 2512
      ENDIF
     
C     Pas d'absorption gazeuse si IPROFIL=2
      IF ((IPROFIL.EQ.2).AND.(ABSPROFIL.NE.7)) GOTO 2513
      
C     Résolution spectrale à préciser si ABSPROFIL != 7
C     Mode de calcul des luminances à préciser si ABSPROFIL != 7	
      IF (ABSPROFIL.NE.7) THEN  ! Cas avec absorption gazeuse

      	  IF (NUSTEP.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2514
		
         IF ((NUSTEP.NE.1).AND.(NUSTEP.NE.5)
     &      .AND.(NUSTEP.NE.10)) GOTO 25141

         IF (IMODE_CKD_CALCUL.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2515
                
      ELSE
          IMODE_CKD_CALCUL = 2	!Mode par défaut pour cas sans absorption gazeuse  
      ENDIF


      
C ---- Affichage des paramètres utilisés dans le fichier de trace
      
      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "************** ",
     &       "Profil parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "

C        Epaisseur optique Rayleigh
	  WRITE (IDLOG,*,ERR=1303) "Molecular Optical Thickness: ",TR

         IF (TR.EQ.CTE_NOT_DEFINED_VALUE_DBLE) THEN
             WRITE (IDLOG,*,ERR=1303)
     &       "--> Calculated by CNES formulation (Perbos, 1982)"

             WRITE (IDLOG,*,ERR=1303) "Atmospheric pressure: ",PSURF
             IF(PSURF.EQ.CTE_HT_STD_PSURF) THEN
		  WRITE (IDLOG,*,ERR=1303) "--> Default value "
             ELSE
		  WRITE (IDLOG,*,ERR=1303) "--> User value "
             ENDIF
  	  ENDIF


C        Echelle de Hauteur Rayleigh     
	  IF (HR.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
	      WRITE (IDLOG,*,ERR=1303) "Molecular Height scale ",HR
	  ENDIF
	  
C	  Type de profil
      	  WRITE (IDLOG,*,ERR=1303) "Aerosol profile type: ",IPROFIL           
  
C        Echelle de Hauteur aérosol     
	  IF (HR.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
	      WRITE (IDLOG,*,ERR=1303) "Aerosol Height Scale ",HA
	  ENDIF                   

C        Altitudes min et max du profil si IPROFIL=2        
      	  IF (IPROFIL.EQ.2) THEN
	      WRITE (IDLOG,*,ERR=1303) "Profil min altitude ",ZMIN
	      WRITE (IDLOG,*,ERR=1303) "Profil max altitude ",ZMAX
      	  ENDIF	  

C     	  Type de profil d'absorption gazeuse            
      	  WRITE (IDLOG,*,ERR=1303) "Absorption Profil type ",ABSPROFIL
          
C        Fichier utilisateur du profil d'absorption si ABSPROFIL = 0	  
          IF (ABSPROFIL.EQ.0) THEN
      	      WRITE (IDLOG,*,ERR=1303) "Absorption user Profile file ",
     &	      FICABSPROFIL
          ENDIF
	  
      ENDIF  ! Fin de trace
     
   	
C -------------------------------------------------------------------  
C --- Controle des paramètres propres aux calculs d'Ordres Successifs 
C     et à la production des fichiers résultats
C -------------------------------------------------------------------  
C  
C ---- Test d'existence et de validité des paramètres obligatoires

C     Controle ordre maximal d'interaction       
      IF (IGMAX.LT.1) GOTO 2604

C     Option de visée : obligatoire  
      IF (ITRPHI.EQ.CTE_NOT_DEFINED_VALUE_INT ) GOTO 2605

C     Option de visée : obligatoire
      IF ((ITRPHI.LT.1).OR.(ITRPHI.GT.2)) GOTO 2606

      !Si choix ITRPHI=1 : il faut définir l'azimut relatif	 
      IF (ITRPHI.EQ.1) THEN 
         IF (PHIOS.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2607
      ENDIF
          	
      !Si choix ITRPHI=2 : il faut définir le pas en azimut	 
      IF (ITRPHI.EQ.2) THEN 
         IF (PAS_PHI.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2608
         IF (PAS_PHI.LE.0) GOTO 2609
      ENDIF
      
      
C ---- Test d'existence et de validité des paramètres optionnels

C     Option de coupure de la polarisation   
      IF ((IPOLAR.LT.0).OR.(IPOLAR.GT.1)) GOTO 2610    
     
C     Niveau de sortie des simulations de luminance : optionnel  
C     Si non défini : la sortie est standard : 
C        Champ émergeant au TOA, champ descendant au SOL : ZOUT = -1     
      IF (((ZOUT.LT.0.D+00).AND.(ZOUT.NE.-1.D+00))
     &     .OR.(ZOUT.GT.CTE_TOA_ALT)) THEN
         GOTO 2611
      ENDIF      
      
      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &   "************** SOS parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "

	  WRITE (IDLOG,*,ERR=1303) "***** Mandatory Parameters "
	  
C        Facteur de dépolarisation moléculaire
	  WRITE (IDLOG,*,ERR=1303) "Molecular depolarization factor: ",
     &	  RON

C        Ordre maximal d'interaction :     
	  WRITE (IDLOG,*,ERR=1303) "Maximum Number of interactions: ",
     &	  IGMAX

C	  Option de visée
	  WRITE (IDLOG,*,ERR=1303) "Viewing mode: ",ITRPHI

C	  Cas ITRPHI=1 => Affichage azimut relatif		
	  IF (ITRPHI.EQ.1) THEN
	  	WRITE (IDLOG,*,ERR=1303) "Relative azimut: ",PHIOS
	  ENDIF	  

C	  Cas ITRPHI=2 => Affichage du pas d'échantillonnage en phi	  	
	  IF (ITRPHI.EQ.2) THEN
	      WRITE (IDLOG,*,ERR=1303) "Azimut sampling step: ",PAS_PHI
	  ENDIF		  

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "***** Optional Parameters "	  

C        Option de coupure de la polarisation : optionnelle
	  IF (IPOLAR.EQ.0) THEN
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "Simulation without polarization"
	  ELSE
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "Simulation including polarization"
	  ENDIF	  

C        Niveau de sortie des luminances : optionnelle
         IF (ZOUT.NE.-1.D+00) THEN
	     WRITE (IDLOG,*,ERR=1303) "Output Altitude (km): ",
     &	     ZOUT
	  ENDIF	  

      ENDIF  ! Fin trace
	 

C =========================================================
C --- Gestion de l'arborescence des fichiers de sortie
C =========================================================     
 
      
C --- Fichiers de trace

      !Fichier Trace Angles
      IF (FICANGLOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICANGLOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3001
          FICANGLOG = DIR_LOG(1:LDIRLOG)//'/'//FICANGLOG(1:LFIC)
      ENDIF 
       
      !Fichier Trace des calculs de Mie des aérosols
      IF (FICMIE_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICMIE_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3002
          FICMIE_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICMIE_LOG(1:LFIC)
      ENDIF 

      !Fichier Trace des calculs de phase des aérosols
      IF (FICGRANU_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICGRANU_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3003
          FICGRANU_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICGRANU_LOG(1:LFIC)
      ENDIF 
       
      !Fichier Trace des calculs de profils
      IF (FICPROFIL_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICPROFIL_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3004
          FICPROFIL_LOG = DIR_LOG(1:LDIRLOG)//'/'
     &	  		 //FICPROFIL_LOG(1:LFIC)
          TRACEPROF=.TRUE.
      ELSE
          TRACEPROF=.FALSE.	  
      ENDIF     
       
      !Fichier Trace des calculs de matrices de surface
      IF (FICSURF_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICSURF_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3005
          FICSURF_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICSURF_LOG(1:LFIC)
      ENDIF     
       
      !Fichier Trace des calculs SOS
      IF (FICSOS_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICSOS_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3006
          FICSOS_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICSOS_LOG(1:LFIC)
	   TRACEOS=.TRUE.
      ELSE
          TRACEOS=.FALSE.	   
      ENDIF           

      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &       "************** Trace filenames ***************"
         WRITE (IDLOG,*,ERR=1303) " "
	  
C        Fichier Trace Angles
         WRITE (IDLOG,*,ERR=1303) "Angle trace filename: ",FICANGLOG
	
C        Fichier Trace des calculs de Mie des aérosols
         WRITE (IDLOG,*,ERR=1303) "Mie trace filename: ",FICMIE_LOG
	
C        Fichier Trace des calculs de phase des aérosols
         WRITE (IDLOG,*,ERR=1303) "Aerosol trace filename: ",
     &	  FICGRANU_LOG

C        Fichier Trace des calculs de profils
         WRITE (IDLOG,*,ERR=1303) "Profil trace filename: ",
     &	  FICPROFIL_LOG
	  
C        Fichier Trace des calculs de matrices de surface
         WRITE (IDLOG,*,ERR=1303) "Surface trace filename: ",
     &	  FICSURF_LOG

C        Fichier Trace des calculs SOS
         WRITE (IDLOG,*,ERR=1303) "SOS trace filename: ",FICSOS_LOG
              	  
      ENDIF	! fin trace      
                                     
C ---- Fichiers de sortie
       
      !Fichier Angles pour luminances
      LFIC = INDEX(FICANGLES_RES_LUM,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3007
      FICANGLES_RES_LUM = DIR_SOS(1:LDIRSOS)//'/'
     &                    //FICANGLES_RES_LUM(1:LFIC)
         
      !Fichier Angles pour calculs de Mie
      LFIC = INDEX(FICANGLES_RES_MIE,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3008
      FICANGLES_RES_MIE = DIR_SOS(1:LDIRSOS)//'/'
     &                    //FICANGLES_RES_MIE(1:LFIC)

        
      !Fichier matrice de phase des aérosols
      LFIC = INDEX(FICGRANU,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3009
      FICGRANU = DIR_SOS(1:LDIRSOS)//'/'//FICGRANU(1:LFIC)
             
                
      !Fichier binaire de luminance en séries de Fourier
      IF (FICSOS_RES_BIN.NE.'NO_OUTPUT') THEN
         LFIC = INDEX(FICSOS_RES_BIN,' ') - 1
         IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3011
         FICSOS_RES_BIN = DIR_SOS(1:LDIRSOS)//'/'
     &        //FICSOS_RES_BIN(1:LFIC)
      ENDIF

      !Test d'existence du fichier
      LFICOS=INDEX(FICSOS_RES_BIN,' ')-1
      INQUIRE(FILE=FICSOS_RES_BIN(1:LFICOS),EXIST=EX)
      IF (EX.EQV..TRUE.) THEN
       
c         WRITE(6,*) 
c     &       'File -SOS.ResBin exists from previous calculations:',
c     &       FICSOS_RES_BIN(1:LFICOS)
          
	  IF (CTE_REMOVE_PREVIOUS_BIN_FILE.EQ.0) THEN
	     GOTO 30111 
	  ELSE
c	     WRITE(6,*)'  Cte CTE_REMOVE_PREVIOUS_BIN_FILE is set to 1:'
c	     WRITE(6,*)'  The file is deleted before new computations'
            OPEN(UNIT=17,FILE=FICSOS_RES_BIN,FORM='UNFORMATTED') 
            CLOSE(17,STATUS='DELETE')
	  ENDIF   
	  
      ENDIF    
       
        
      !Fichier de transmissions
      LFIC = INDEX(FICTRANS,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3012
      IF (FICTRANS.NE.'NO_OUTPUT') THEN
          FICTRANS = DIR_SOS(1:LDIRSOS)//'/'//FICTRANS(1:LFIC)
      ENDIF
 
      !Fichier de flux
      LFIC = INDEX(FICFLUX,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 30121
      IF (FICFLUX.NE.'NO_OUTPUT') THEN
          FICFLUX = DIR_SOS(1:LDIRSOS)//'/'//FICFLUX(1:LFIC)
      ENDIF

      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &       "************** Output filenames ***************"
         WRITE (IDLOG,*,ERR=1303) " "
	  

C        Fichier Angles pour luminances
         WRITE (IDLOG,*,ERR=1303) "Radiance Angles: ",FICANGLES_RES_LUM
       	
C        Fichier Angles pour calculs de Mie
         WRITE (IDLOG,*,ERR=1303) "Mie Angles: ",FICANGLES_RES_MIE
	
C        Fichier matrice de phase des aérosols
         WRITE (IDLOG,*,ERR=1303) "Aerosol phase matrix: ",
     &	  FICGRANU         
	  
C        Fichier binaire de luminance en séries de Fourier
         WRITE (IDLOG,*,ERR=1303) "SOS binary output: ",FICSOS_RES_BIN

C        Fichier de transmissions
         WRITE (IDLOG,*,ERR=1303) "Transmittances: ",FICTRANS

C        Fichier de flux
         WRITE (IDLOG,*,ERR=1303) "Flux: ",FICFLUX
	      	        	  
      ENDIF	! fin trace            

C =========================================================
C --- Préparation des angles pour les traitements
C =========================================================   

      WRITE(6,*) "DEFINING ANGLES FOR COMPUTATIONS"

      IF (FICANGLES_USER_LUM.EQ.'NO_USER_ANGLES') THEN
          LFIC=14
      ELSE
          LFIC=INDEX(FICANGLES_USER_LUM,' ')
      ENDIF
      FICANGLES_USER_LUM=FICANGLES_USER_LUM(1:LFIC)

      IF (FICANGLES_USER_MIE.EQ.'NO_USER_ANGLES') THEN
          LFIC=14
      ELSE
          LFIC=INDEX(FICANGLES_USER_MIE,' ')
      ENDIF
      FICANGLES_USER_MIE=FICANGLES_USER_MIE(1:LFIC)

      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "ANGLES CALCULATION"
	  WRITE (IDLOG,*,ERR=1303) "------------------"
         WRITE (IDLOG,*,ERR=1303) 
     &	  "--> SOS_ANGLES is called with the parameters :"
	  WRITE (IDLOG,*,ERR=1303) "NBMU_GAUSS_LUM ",NBMU_GAUSS_LUM
	  WRITE (IDLOG,*,ERR=1303) "TETAS ",TETAS
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_USER_LUM ",
     &	  FICANGLES_USER_LUM
	  WRITE (IDLOG,*,ERR=1303) "NBMU_GAUSS_MIE ",NBMU_GAUSS_MIE
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_USER_MIE ",
     &	  FICANGLES_USER_MIE
	  WRITE (IDLOG,*,ERR=1303) "FICANGLOG ",FICANGLOG
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_RES_LUM ",
     &	  FICANGLES_RES_LUM
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_RES_MIE ",
     &	  FICANGLES_RES_MIE	  
       ENDIF

       CALL SOS_ANGLES(NBMU_GAUSS_LUM, TETAS, FICANGLES_USER_LUM,                 
     &      NBMU_GAUSS_MIE, FICANGLES_USER_MIE, 
     &      FICANGLOG,
     &      FICANGLES_RES_LUM, FICANGLES_RES_MIE, IER)
     
      IF (IER.NE.0) THEN ! cas d'erreur
      		
         IF (TRACE) THEN
	    WRITE (IDLOG,*,ERR=1303) 
     &	    "   ==> Error in subroutine SOS_ANGLES"
 	  ENDIF
	  GOTO 4000   
	  
      ENDIF
     
C =========================================================
C --- Préparation de la matrice de phase des aérosols
C =========================================================    
      

      WRITE(6,*) "AEROSOLS RADIATIVE PROPERTIES COMPUTATION"
 
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "-----------------------------------------"
	  WRITE (IDLOG,*,ERR=1303)
     &	  "AEROSOLS RADIATIVE PROPERTIES COMPUTATION"
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "-----------------------------------------"	  
      ENDIF !Fin Trace 	     


C --- Si l'épaisseur optique est nulle ou si la longueur d'onde de 
C --- référence est la même que celle de simu, alors on lance la 
C --- routine SOS_aérosols une seule fois. Sinon, il faut l'exécuter deux 
C --- fois et déterminer l'épaisseur optique à la longueur d'onde de simu.

C     Pour le premier lancement de la routine, on affecte à TA 
C     (ep optique de simu) la valeur de AOT_REF
      TA = AOT_REF

C     Détermination de la longueur d'onde pour le lancement de la routine
      		            
      IF (AOT_REF.EQ.0.D+00) THEN  ! Cas AOT_REF nulle
      
         ! Dans ce cas, on effectue les calculs directement à la longueur d'onde de simu
	  WA = WA_SIMU 	
	  
      ELSE ! Cas AOT_REF non nulle

	  ! Les traitements sont faits à la longueur d'onde de référence  
         ! Un deuxieme traitement sera ensuite fait à la longueur d'onde
	  ! de simu si elle est différente de la référence
	  WA = WAREF_AOT 
			  
      ENDIF

C --- Gestion de la trace (affichage paramètres de calcul
      
      IF (TRACE) THEN
	 
         IF (( FICUSER_AER.EQ.'NO_USER_AEROSOLS' )
     &        .OR.( AOT_REF.EQ.0.D+00)) THEN
	 
           WRITE (IDLOG,*,ERR=1303) " "
	    WRITE (IDLOG,*,ERR=1303) 
     &	    "=> First Calculation of aerosol parameters"
           WRITE (IDLOG,*,ERR=1303) 
     &	    "--> SOS_AEROSOLS is called with the parameters :"
	    WRITE (IDLOG,807,ERR=1303) FICANGLES_RES_MIE
	    WRITE (IDLOG,841,ERR=1303) WA
	    WRITE (IDLOG,842,ERR=1303) TA
	    WRITE (IDLOG,843,ERR=1303) WAREF_AOT
	    WRITE (IDLOG,844,ERR=1303) AOT_REF
	    WRITE (IDLOG,845,ERR=1303) ITRONC_AER
	    WRITE (IDLOG,846,ERR=1303) IMOD_AER
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,*,ERR=1303) "      mono-modal model :"
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
 	    WRITE (IDLOG,847,ERR=1303) RN_WAREF, IN_WAREF 
	    WRITE (IDLOG,849,ERR=1303) IGRANU
	    WRITE (IDLOG,850,ERR=1303) VARGRANU1_MMD_AER, 
     &	     			          VARGRANU2_MMD_AER, 
     &	     			          VARGRANU3_MMD_AER
	    WRITE (IDLOG,*,ERR=1303) "      -----------"
	    WRITE (IDLOG,*,ERR=1303) "      WMO model :"
	    WRITE (IDLOG,*,ERR=1303) "      -----------"
           WRITE (IDLOG,851,ERR=1303) IModele_WMO
	    WRITE (IDLOG,852,ERR=1303) C_WMO_DL,C_WMO_WS,C_WMO_OC,
     &	    			          C_WMO_SO
	    WRITE (IDLOG,*,ERR=1303) "      ----------------------"
	    WRITE (IDLOG,*,ERR=1303) "      Shettle & Fenn model :"
	    WRITE (IDLOG,*,ERR=1303) "      ----------------------"
	    WRITE (IDLOG,853,ERR=1303) IModele_SF
	    WRITE (IDLOG,854,ERR=1303) RH
	    WRITE (IDLOG,*,ERR=1303) "      -----------------------"
	    WRITE (IDLOG,*,ERR=1303) "      Bi-modal distribution :"
	    WRITE (IDLOG,*,ERR=1303) "      -----------------------"
	    WRITE (IDLOG,855,ERR=1303) MODE_PARAM_BILND
	    WRITE (IDLOG,856,ERR=1303) USER_CV_COARSE, USER_CV_FINE
	    WRITE (IDLOG,857,ERR=1303) RTAUct_WAREF
	    WRITE (IDLOG,*,ERR=1303) "      Coarse mode :"
	    WRITE (IDLOG,858,ERR=1303) BMD_CM_MRWAREF, BMD_CM_MIWAREF
	    WRITE (IDLOG,859,ERR=1303) BMD_CM_MRWAREF, BMD_CM_MIWAREF
	    WRITE (IDLOG,860,ERR=1303) BMD_CM_RMODAL, BMD_CM_VAR
	    WRITE (IDLOG,*,ERR=1303) "      Fine mode :"
	    WRITE (IDLOG,858,ERR=1303) BMD_FM_MRWAREF, BMD_FM_MIWAREF
	    WRITE (IDLOG,859,ERR=1303) BMD_FM_MRWAREF, BMD_FM_MIWAREF
	    WRITE (IDLOG,860,ERR=1303) BMD_FM_RMODAL, BMD_FM_VAR
	    WRITE (IDLOG,*,ERR=1303) "      ---------------"
 	    WRITE (IDLOG,*,ERR=1303) "      External data :"
	    WRITE (IDLOG,*,ERR=1303) "      ---------------"
 	    WRITE (IDLOG,861,ERR=1303) FICEXTDATA_AER  
	    WRITE (IDLOG,*,ERR=1303) "      --------------------"
 	    WRITE (IDLOG,*,ERR=1303) "      Mixture definition :"
	    WRITE (IDLOG,*,ERR=1303) "      --------------------"
 	    WRITE (IDLOG,8611,ERR=1303) FICMIXTURE_AER  
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,*,ERR=1303) "      MIE computations :"
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,864,ERR=1303) DIR_MIE 
	    WRITE (IDLOG,865,ERR=1303) FICMIE_LOG 
	    WRITE (IDLOG,*,ERR=1303) "      --------------"
	    WRITE (IDLOG,*,ERR=1303) "      Output files :"
	    WRITE (IDLOG,*,ERR=1303) "      --------------"
 	    WRITE (IDLOG,862,ERR=1303) FICGRANU_LOG 
	    WRITE (IDLOG,863,ERR=1303) FICGRANU
	    
	  ELSE
	 
	    WRITE (IDLOG,*,ERR=1303) " "
           WRITE (IDLOG,*,ERR=1303) 
     &	    "--> SOS_AEROSOLS is not launched"
     	    IF (FICUSER_AER.NE.'NO_USER_AEROSOLS') THEN
               WRITE (IDLOG,*,ERR=1303) "UserFile : ",FICUSER_AER
 	    ENDIF
	  ENDIF ! fin test fichier utilisateur ou AOT=0
    
      ENDIF !Fin Trace   
    

C     Si pas de fichier utilisateur => on lance les traitements
 
      IF (( FICUSER_AER.EQ.'NO_USER_AEROSOLS' )
     &     .OR.( AOT_REF.EQ.0.D+00)) THEN

C --- Lancement de la routine aérosols 

         ITRONC_AER_IN = ITRONC_AER !Si ITRONC_AER=1 mais que le code estime
				 	 !la troncature inutile, seul le paramètre 
					 !ITRONC_AER_IN basculera de 1 à 0, sans
					 !impacter la valeur de ITRONC_AER

         CALL SOS_AEROSOLS(FICANGLES_RES_MIE, 
     &        WA, TA, WAREF_AOT, AOT_REF,	  
     &        ITRONC_AER_IN, IMOD_AER,
     &        RN_WAREF, IN_WAREF, IGRANU, VARGRANU1_MMD_AER,
     &        VARGRANU2_MMD_AER, VARGRANU3_MMD_AER,
     &        IModele_WMO, 
     &        C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &        IModele_SF, RH, 
     &        MODE_PARAM_BILND, 
     &        USER_CV_COARSE, USER_CV_FINE, RTAUct_WAREF,
     &        BMD_CM_MRWAREF, BMD_CM_MIWAREF, 
     &        BMD_CM_MRWAREF, BMD_CM_MIWAREF,
     &        BMD_CM_RMODAL, BMD_CM_VAR,
     &        BMD_FM_MRWAREF, BMD_FM_MIWAREF, 
     &        BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &        BMD_FM_RMODAL, BMD_FM_VAR,	 
     &        FICEXTDATA_AER, FICMIXTURE_AER,
     &        DIR_TMP, DIR_MIE, FICMIE_LOG,
     &        FICGRANU, FICGRANU_LOG, KMAT1,
     &        PIZAER, COEF_TRONCA, IER)
           			  
C --- Gestion cas d'erreur 
				  
         IF (IER.NE.0) THEN ! cas d'erreur 
	  			
             IF (TRACE) THEN
	         WRITE (IDLOG,*,ERR=1303) 
     &	         "   ==> Error in subroutine SOS_AEROSOLS"
 	      ENDIF
	 
	      GOTO 4100   
	 
         ENDIF 	  

      ELSE

C    Avec fichier utilisateur => pas de traitements

         !Fichier matrice de phase des aérosols
         LFIC = INDEX(FICUSER_AER,' ') - 1
         FICGRANU = FICUSER_AER(1:LFIC)
      	  
      ENDIF
      
      
C --- Cas où WAREF != WA_SIMU   (rem: pas le cas si utilisation d'un fichier utilisateur)
C --- Lancement de la routine aérosol pour WA_SIMU

      IF ((AOT_REF.NE.0.D+00).AND.(WA_SIMU.NE.WAREF_AOT)) THEN

C ---    Sauvegarde du coefficient KMAT calculé précédemment
	  KMAT1_WAREF=KMAT1
	  
C ---    Epaisseur optique de simu fixee arbitrairement à 0.1
C ---    pour permettre les calculs aérosols
	  TA=0.1	        

C ---    Gestion de la trace (affichage paramètres de calcul)
      
         IF (TRACE) THEN
	 
           WRITE (IDLOG,*,ERR=1303) " "
	    WRITE (IDLOG,*,ERR=1303) 
     &	     "=> Second Calculation of aerosol parameters"
	    WRITE (IDLOG,*,ERR=1303) "Case WA = WA_SIMU"	         
           WRITE (IDLOG,*,ERR=1303) 
     &	     "--> SOS_AEROSOLS is called with the parameters :"
	    WRITE (IDLOG,807,ERR=1303) FICANGLES_RES_MIE
	    WRITE (IDLOG,841,ERR=1303) WA_SIMU
	    WRITE (IDLOG,842,ERR=1303) TA
	    WRITE (IDLOG,843,ERR=1303) WAREF_AOT
	    WRITE (IDLOG,844,ERR=1303) AOT_REF
	    WRITE (IDLOG,845,ERR=1303) ITRONC_AER
	    WRITE (IDLOG,846,ERR=1303) IMOD_AER
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,*,ERR=1303) "      mono-modal model :"
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
 	    WRITE (IDLOG,847,ERR=1303) RN_WA, IN_WA
	    WRITE (IDLOG,849,ERR=1303) IGRANU
	    WRITE (IDLOG,850,ERR=1303) VARGRANU1_MMD_AER, 
     &	     			          VARGRANU2_MMD_AER, 
     &	     			          VARGRANU3_MMD_AER
	    WRITE (IDLOG,*,ERR=1303) "      -----------"
	    WRITE (IDLOG,*,ERR=1303) "      WMO model :"
	    WRITE (IDLOG,*,ERR=1303) "      -----------"
           WRITE (IDLOG,851,ERR=1303) IModele_WMO
	    WRITE (IDLOG,852,ERR=1303) C_WMO_DL,C_WMO_WS,C_WMO_OC,
     &	    			      C_WMO_SO
	    WRITE (IDLOG,*,ERR=1303) "      ----------------------"
	    WRITE (IDLOG,*,ERR=1303) "      Shettle & Fenn model :"
	    WRITE (IDLOG,*,ERR=1303) "      ----------------------"
	    WRITE (IDLOG,853,ERR=1303) IModele_SF
	    WRITE (IDLOG,854,ERR=1303) RH
	    WRITE (IDLOG,*,ERR=1303) "      -----------------------"
	    WRITE (IDLOG,*,ERR=1303) "      Bi-modal distribution :"
	    WRITE (IDLOG,*,ERR=1303) "      -----------------------"
	    WRITE (IDLOG,855,ERR=1303) MODE_PARAM_BILND
	    WRITE (IDLOG,856,ERR=1303) USER_CV_COARSE, USER_CV_FINE
	    WRITE (IDLOG,857,ERR=1303) RTAUct_WAREF
	    WRITE (IDLOG,*,ERR=1303) "      Coarse mode :"
	    WRITE (IDLOG,858,ERR=1303) BMD_CM_MRWA, BMD_CM_MIWA
	    WRITE (IDLOG,859,ERR=1303) BMD_CM_MRWAREF, BMD_CM_MIWAREF
	    WRITE (IDLOG,860,ERR=1303) BMD_CM_RMODAL, BMD_CM_VAR
	    WRITE (IDLOG,*,ERR=1303) "      Fine mode :"
	    WRITE (IDLOG,858,ERR=1303) BMD_FM_MRWA, BMD_FM_MIWA
	    WRITE (IDLOG,859,ERR=1303) BMD_FM_MRWAREF, BMD_FM_MIWAREF
	    WRITE (IDLOG,860,ERR=1303) BMD_FM_RMODAL, BMD_FM_VAR
	    WRITE (IDLOG,*,ERR=1303) "      ---------------"
 	    WRITE (IDLOG,*,ERR=1303) "      External data :"
	    WRITE (IDLOG,*,ERR=1303) "      ---------------"
 	    WRITE (IDLOG,861,ERR=1303) FICEXTDATA_AER  
	    WRITE (IDLOG,*,ERR=1303) "      --------------------"
 	    WRITE (IDLOG,*,ERR=1303) "      Mixture definition :"
	    WRITE (IDLOG,*,ERR=1303) "      --------------------"
 	    WRITE (IDLOG,8611,ERR=1303) FICMIXTURE_AER  
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,*,ERR=1303) "      MIE computations :"
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,864,ERR=1303) DIR_MIE 
	    WRITE (IDLOG,865,ERR=1303) FICMIE_LOG 
	    WRITE (IDLOG,*,ERR=1303) "      --------------"
	    WRITE (IDLOG,*,ERR=1303) "      Output files :"
	    WRITE (IDLOG,*,ERR=1303) "      --------------"
 	    WRITE (IDLOG,862,ERR=1303) FICGRANU_LOG 
	    WRITE (IDLOG,863,ERR=1303) FICGRANU
	    
         ENDIF !Fin Trace       
            
C ---     Lancement de la routine aérosols 

         ITRONC_AER_IN = ITRONC_AER   !Si ITRONC_AER=1 mais que le code estime
					   !la troncature inutile, seul le paramètre 
					   !ITRONC_AER_IN basculera de 1 à 0, sans
					   !impacter la valeur de ITRONC_AER

         CALL SOS_AEROSOLS(FICANGLES_RES_MIE, 
     &        WA_SIMU, TA, WAREF_AOT, AOT_REF,	  
     &        ITRONC_AER_IN, IMOD_AER,
     &        RN_WA, IN_WA, IGRANU, VARGRANU1_MMD_AER,
     &        VARGRANU2_MMD_AER, VARGRANU3_MMD_AER,
     &        IModele_WMO, 
     &        C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &        IModele_SF, RH, 
     &        MODE_PARAM_BILND, 
     &        USER_CV_COARSE, USER_CV_FINE, RTAUct_WAREF,
     &        BMD_CM_MRWA, BMD_CM_MIWA, 
     &        BMD_CM_MRWAREF, BMD_CM_MIWAREF,
     &        BMD_CM_RMODAL, BMD_CM_VAR,
     &        BMD_FM_MRWA, BMD_FM_MIWA, 
     &        BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &        BMD_FM_RMODAL, BMD_FM_VAR,	 
     &        FICEXTDATA_AER, FICMIXTURE_AER,
     &        DIR_TMP, DIR_MIE, FICMIE_LOG,
     &        FICGRANU, FICGRANU_LOG, KMAT1,
     &        PIZAER, COEF_TRONCA, IER)

C ---    Gestion cas d'erreur 
				  
         IF (IER.NE.0) THEN ! cas d'erreur 
	  			
             IF (TRACE) THEN
	          WRITE (IDLOG,*,ERR=1303) 
     &	          "   ==> Error in subroutine SOS_AEROSOLS"
 	      ENDIF
	 
	      GOTO 4100   
	 
         ENDIF 	 

C ---    Calcul de l'épaisseur optique à WA_SIMU
	  TA = (KMAT1 / KMAT1_WAREF) * AOT_REF
         IF (TA.LT.0.0) GOTO 2516

C ---    Gestion de la trace (paramètres de calcul de AOT à WA_SIMU)

	  IF (TRACE) THEN
	      WRITE (IDLOG,*,ERR=1303) " "
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "-----------------------------------------"
	      WRITE (IDLOG,*,ERR=1303)
     &	      "AEROSOLS RADIATIVE PROPERTIES COMPUTATION"
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "-----------------------------------------"	
	      WRITE (IDLOG,*,ERR=1303) " "
	      WRITE (IDLOG,*,ERR=1303) "=> AOT estimation at WA : "
	      WRITE (IDLOG,*,ERR=1303) "    WAREF   : ", WAREF_AOT
	      WRITE (IDLOG,*,ERR=1303) "    AOT_REF : ", AOT_REF
	      WRITE (IDLOG,*,ERR=1303) "    WA      : ", WA_SIMU         
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "    Sect eff ext (WAREF) : ", KMAT1_WAREF
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "    Sect eff ext (WA)    : ", KMAT1
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "    ===> AOT estimated at WA : ", TA
 	  ENDIF !Fin Trace 

                 			  	    
      ENDIF ! fin cas WAREF != WA_SIMU 
      	      				  
C =========================================================
C --- Préparation des propriétés de surface
C =========================================================     
      				

C --- Détermination répertoire de sortie
      IF (ISURF.EQ.1) THEN
      	   ! Cas glitter
          SURF_OUTDIR(1:LDIRGLIT)=DIR_GLIT(1:LDIRGLIT)
      ENDIF 
      IF (ISURF.EQ.3) THEN
      	   ! Cas Roujean
          SURF_OUTDIR(1:LDIRROUJ)=DIR_ROUJ(1:LDIRROUJ)
      ENDIF       
      IF (ISURF.EQ.4) THEN
      	   ! Cas Rondeaux Herman
          SURF_OUTDIR(1:LDIRRH)=DIR_RH(1:LDIRRH)
      ENDIF       
      IF (ISURF.EQ.5) THEN
      	   ! Cas Breon
          SURF_OUTDIR(1:LDIRBREON)=DIR_BREON(1:LDIRBREON)
      ENDIF       
      IF (ISURF.EQ.6) THEN
      	   ! Cas Nadal
          SURF_OUTDIR(1:LDIRNADAL)=DIR_NADAL(1:LDIRNADAL)
      ENDIF       
      IF (ISURF.EQ.7) THEN
      	   ! Cas Maignan
          SURF_OUTDIR(1:LDIRMAIGNAN)=DIR_MAIGNAN(1:LDIRMAIGNAN)
      ENDIF             

C --- Gestion de la trace
      
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "SURFACE CALCULATION"
	  WRITE (IDLOG,*,ERR=1303) "------------------"
	  
	  !Cas lambertien
	  IF (ISURF.EQ.0) THEN	  
	  
	     WRITE (IDLOG,*,ERR=1303) "Lambertian case "
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "--> SOS_SURFACE is not launched"     
     
         ELSE
	  
	     !Cas Fesnel
	     IF (ISURF.EQ.2) THEN
	     
	        WRITE (IDLOG,*,ERR=1303) "Lambertian case + Fresnel"
	        WRITE (IDLOG,*,ERR=1303) 
     &	        "--> SOS_SURFACE is not launched"
     
     	     ELSE  

		! Cas fichier utilisateur
                IF (FICSURF.NE."DEFAULT") THEN
		
                  WRITE (IDLOG,*,ERR=1303) "User Surface File Case"
                  WRITE (IDLOG,*,ERR=1303) "FICSURF ",FICSURF
	           WRITE (IDLOG,*,ERR=1303) 
     &                  "--> SOS_SURFACE is not launched"
                   
                ELSE            ! Autres cas, lancement de la routine
                   
          	    WRITE (IDLOG,*,ERR=1303) 
     &	  	    "--> SOS_SURFACE is called with the parameters :"
	  	    WRITE (IDLOG,*,ERR=1303) "ISURF ",ISURF
	  	    WRITE (IDLOG,*,ERR=1303) "SURF_IND ",SURF_IND
	  	    WRITE (IDLOG,*,ERR=1303) "WIND ",WIND
	  	    WRITE (IDLOG,*,ERR=1303) "K0 ",K0_ROUJEAN
	  	    WRITE (IDLOG,*,ERR=1303) "K1 ",K1_ROUJEAN 
	  	    WRITE (IDLOG,*,ERR=1303) "K2 ",K2_ROUJEAN 
	  	    WRITE (IDLOG,*,ERR=1303) "ALPHA_NADAL ",ALPHA_NADAL
	  	    WRITE (IDLOG,*,ERR=1303) "BETA_NADAL ",BETA_NADAL
	  	    WRITE (IDLOG,*,ERR=1303) "COEF_C_MAIGNAN ",
     &	  	    COEF_C_MAIGNAN
	  	    WRITE (IDLOG,*,ERR=1303) "FICANGLES_RES_LUM ",
     &	  	    FICANGLES_RES_LUM
	  	    WRITE (IDLOG,*,ERR=1303) "FICSURF_LOG ",FICSURF_LOG
	  	    WRITE (IDLOG,*,ERR=1303) "SURF_OUTDIR ",SURF_OUTDIR
	  	    WRITE (IDLOG,*,ERR=1303) "DIR_ROUJ ",DIR_ROUJ
	     
	     	  ENDIF          ! Fin cas utilisateur
                 
            ENDIF             ! Fin test Fresnel
	     
         ENDIF                !  Fin test lambertien 
	  
      ENDIF                   ! Fin trace
      
C --- Lancement de l'outil de simulation des propriétés de surface
C --- Cet outil n'est pas exécuté si ISURF=0 ou 2 
C --- Ou si un fichier de surface FICSURF est défini par l'utilisateur

      IF ((ISURF.NE.0).AND.(ISURF.NE.2)
     &	    .AND.(FICSURF.EQ."DEFAULT")) THEN
 
         WRITE(6,*) 'SURFACE BRDF/BPDF MATRICES COMPUTATION'   

         IF (TRACE) THEN
            WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "-----------------------------------------------"
	     WRITE (IDLOG,*,ERR=1303)
     &	     "SURFACE BRDF/BPDF MATRICES COMPUTATION"
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "-----------------------------------------------"
            WRITE (IDLOG,*,ERR=1303) " "
         ENDIF
        
         CALL SOS_SURFACE(ISURF, FICANGLES_RES_LUM, SURF_IND, WIND,  
     &        K0_ROUJEAN, K1_ROUJEAN , K2_ROUJEAN,  
     &        ALPHA_NADAL, BETA_NADAL, COEF_C_MAIGNAN,  
     &        FICSURF_LOG,                
     &        SURF_OUTDIR, DIR_ROUJ,  DIR_TMP, IER)
     
      ENDIF

      
      
C --- Gestion du cas d'erreur
      
      IF (IER.NE.0) THEN ! cas d'erreur
      		
         IF (TRACE) THEN
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "   ==> Error in subroutine SOS_SURFACE"
 	  ENDIF
	  GOTO 4200   
	  
      ENDIF
      
C =============================================================
C --- Préparation du profil atmosphérique + calcul luminances
C =============================================================     

      WRITE(6,*) 'PROFILE AND RADIANCE COMPUTATION'   
C --- Gestion de la trace du main
      
      IF (TRACE) THEN
      
C	  Pour le profil atmosphérique
          
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "PROFILE CALCULATION"
	  WRITE (IDLOG,*,ERR=1303) "-------------------"	  
         WRITE (IDLOG,*,ERR=1303) 
     &	  "--> SOS_PROFILE is called with the parameters :"
	  WRITE (IDLOG,*,ERR=1303) "WA_SIMU ",WA_SIMU
	  WRITE (IDLOG,*,ERR=1303) "IPROFIL ",IPROFIL
	  WRITE (IDLOG,*,ERR=1303) "TR ",TR
	  WRITE (IDLOG,*,ERR=1303) "HR ",HR
	  WRITE (IDLOG,*,ERR=1303) "TA ",TA
	  WRITE (IDLOG,*,ERR=1303) "HA ",HA
	  WRITE (IDLOG,*,ERR=1303) "PSURF ",PSURF
	  WRITE (IDLOG,*,ERR=1303) "ZMIN ",ZMIN
	  WRITE (IDLOG,*,ERR=1303) "ZMAX ",ZMAX
	  WRITE (IDLOG,*,ERR=1303) "H2O ",H2O
	  WRITE (IDLOG,*,ERR=1303) "O3 ",O3
	  WRITE (IDLOG,*,ERR=1303) "CO2 ",CO2
	  WRITE (IDLOG,*,ERR=1303) "CH4 ",CH4
	  WRITE (IDLOG,*,ERR=1303) "ABSPROFIL ",ABSPROFIL
	  WRITE (IDLOG,*,ERR=1303) "FICABSPROFIL ",FICABSPROFIL
	  WRITE (IDLOG,*,ERR=1303) "FICPROFIL_LOG ",FICPROFIL_LOG

C	  Pour le calcul des luminances
      
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "--------------------"
	  WRITE (IDLOG,*,ERR=1303)
     &	  "RADIANCE COMPUTATION"
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "--------------------"    
     
         IF (IGMAX.NE.CTE_DEFAULT_IGMAX)	THEN
	     WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303)  
     &	        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            WRITE (IDLOG,*,ERR=1303)  
     &	        "WARNING : The maximal order of interactions has been"
	     WRITE (IDLOG,*,ERR=1303)  
     &	        "          changed by the user  ==> IGMAX =",IGMAX
	     WRITE (IDLOG,*,ERR=1303)  
     &	        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	     WRITE (IDLOG,*,ERR=1303) " "
         ENDIF 

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "--> SOS_OS is called with the parameters :"
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_RES_LUM ",
     &	  FICANGLES_RES_LUM
	  WRITE (IDLOG,*,ERR=1303) "FICGRANU ",FICGRANU
	  WRITE (IDLOG,*,ERR=1303) "ZOUT ",ZOUT
	  WRITE (IDLOG,*,ERR=1303) "IGMAX ",IGMAX
	  WRITE (IDLOG,*,ERR=1303) "IPOLAR ",IPOLAR
	  WRITE (IDLOG,*,ERR=1303) "ITRPHI ",ITRPHI
	  WRITE (IDLOG,*,ERR=1303) "PHIOS ",PHIOS
	  WRITE (IDLOG,*,ERR=1303) "PAS_PHI ",PAS_PHI	  
	  WRITE (IDLOG,*,ERR=1303) "RON ",RON
	  WRITE (IDLOG,*,ERR=1303) "RHO ",RHO
	  WRITE (IDLOG,*,ERR=1303) "ISURF ",ISURF
	  WRITE (IDLOG,*,ERR=1303) "SURF_IND ",SURF_IND
	  WRITE (IDLOG,*,ERR=1303) "WIND ",WIND
	  WRITE (IDLOG,*,ERR=1303) "K0 ",K0_ROUJEAN
	  WRITE (IDLOG,*,ERR=1303) "K1 ",K1_ROUJEAN 
	  WRITE (IDLOG,*,ERR=1303) "K2 ",K2_ROUJEAN 
	  WRITE (IDLOG,*,ERR=1303) "ALPHA_NADAL ",ALPHA_NADAL
	  WRITE (IDLOG,*,ERR=1303) "BETA_NADAL ",BETA_NADAL
	  WRITE (IDLOG,*,ERR=1303) "SURF_OUTDIR ",SURF_OUTDIR
	  WRITE (IDLOG,*,ERR=1303) "FICSURF ",FICSURF
	  WRITE (IDLOG,*,ERR=1303) "FICSOS_RES_BIN ",FICSOS_RES_BIN
	  WRITE (IDLOG,*,ERR=1303) "FICTRANS ",FICTRANS
	  WRITE (IDLOG,*,ERR=1303) "FICFLUX ",FICFLUX
	  WRITE (IDLOG,*,ERR=1303) "FICSOS_LOG ",FICSOS_LOG	  
	     
      ENDIF   ! Fin trace

C ---- Ouverture des fichiers trace des sous fonctions de profils et d'OS

      IF (TRACEPROF) THEN
          OPEN (IDLOGPROF,FILE=FICPROFIL_LOG,ERR=3104)
      ENDIF
       
      IF (TRACEOS) THEN
          OPEN (IDLOGOS,FILE=FICSOS_LOG,ERR=3106)
      ENDIF
      
C --- Calcul de l'épaisseur optique moléculaire

C     Si l'épaisseure optique n'est pas définie, on la calcule
C     en fonction de la pression atmosphérique selon la formule
C     CNES (Perbos, 1982)

      IF (TR.EQ.CTE_NOT_DEFINED_VALUE_DBLE) THEN
      
         TR = (PSURF/CTE_HT_STD_PSURF) * 1.d-4 *       
     &     (84.35/WA_SIMU**4 - 1.225/WA_SIMU**5 + 1.4/WA_SIMU**6)
     	  WRITE(6,*)       
     &   '==> Rayleigh optical thickness from CNES formulation ', TR

         IF (TRACE) THEN            
            WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303)
     &      "Molecular Optical Thickness calculated by"
	     WRITE (IDLOG,*,ERR=1303)
     &      "the CNES formulation (Perbos, 1982), TR= ",TR
            WRITE (IDLOG,*,ERR=1303) 
     &      "for the wavelength (microns)=",WA_SIMU
            WRITE (IDLOG,*,ERR=1303) 
     &      "and surface pressure (mbar)=",PSURF	      
            WRITE (IDLOG,*,ERR=1303) " "
	  ENDIF 

      ENDIF

       
C --- Lancement de la routine de préparation du profil pour l'absorption
C --- Executée seulement si l'absorption est prise en compte IPROFIL=1 et ABSPROFIL!=7
      
      IF ((IPROFIL.EQ.1).AND.(ABSPROFIL.NE.7)) THEN
     	 
         CALL SOS_PREPA_ABSPROFILE(WA_SIMU, NUSTEP, PSURF, H2O, O3, 
     &		           CO2, CH4, ABSPROFIL, FICABSPROFIL, TRACEPROF,
     &		           IDLOGPROF, NU, LAMB1, IABS,
     &		           ABS_USERPROFIL, ALTABS, RO,
     &			    NEXP, KDIS_AI, KDIS_KI, KDIS_KI_H2O,
     &                   TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                   TAB_CONC_H2O,NB_CONC_H2O, IER)
  
 
C ---    Gestion du cas d'erreur
      
         IF (IER.NE.0) THEN ! cas d'erreur
      		
            IF (TRACE) THEN
	        WRITE (IDLOG,*,ERR=1303) 
     &	        "   ==> Error in subroutine SOS_PREPA_ABSPROFILE"
 	     ENDIF
	     GOTO 4300   
	  
         ENDIF
 
C ---    Calcul de la somme des coefficients AIK de pondération des luminances par exponentielle
         sumAIK=0.D+00
         DO IK1=1,NEXP(1,LAMB1)
         DO IK2=1,NEXP(2,LAMB1)
         DO IK3=1,NEXP(3,LAMB1)
         DO IK4=1,NEXP(4,LAMB1)
         DO IK5=1,NEXP(5,LAMB1)
         DO IK6=1,NEXP(6,LAMB1)
         DO IK7=1,NEXP(7,LAMB1)
         DO IK8=1,NEXP(8,LAMB1)
           !Calcul du produit des fonctions de répartition (poids affecté aux luminances)
            AIK=KDIS_AI(IK1,1,LAMB1)*KDIS_AI(IK2,2,LAMB1)*
     &          KDIS_AI(IK3,3,LAMB1)*KDIS_AI(IK4,4,LAMB1)*
     &          KDIS_AI(IK5,5,LAMB1)*KDIS_AI(IK6,6,LAMB1)*
     &	         KDIS_AI(IK7,7,LAMB1)*KDIS_AI(IK8,8,LAMB1)

            sumAIK=sumAIK+AIK
         ENDDO  !IK1
         ENDDO  !IK2
         ENDDO  !IK3
         ENDDO  !IK4
         ENDDO  !IK5
         ENDDO  !IK6
         ENDDO  !IK7
         ENDDO  !IK8

         IF (TRACE) THEN            
            WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303)
     &      "Summation of AIK coefficients:", sumAIK
            WRITE (IDLOG,*,ERR=1303) " "
	  ENDIF 

         !Controle de l'amplitude d'écart à 1.
         IF (ABS(sumAIK-1.D+00).GE.1.D-06) GOTO 3500

     
      ENDIF ! fin test appel SOS_PREPA_ABSPROFILE


C --- Lancement de la routine de préparation pour le calcul des OS
C --- Cette routine lit les fichiers aérosols et angles
       
      CALL SOS_PREPA_OS(FICANGLES_RES_LUM, FICGRANU, 
     &	     ISURF, SURF_IND,WIND, 
     &	     K0_ROUJEAN, K1_ROUJEAN , K2_ROUJEAN ,  
     &	     ALPHA_NADAL, BETA_NADAL, COEF_C_MAIGNAN,
     &	     SURF_OUTDIR, FICSURF, TRACEOS, IDLOGOS,
     &	     LUM_NBMU, N0, PIZ, PIZTR, COEFTRONCA, RMU, GA, OS_NB, 
     &	     IMAT_SURF, IGLI, IFRESNEL, IROUJEAN, IRONDEAUX, IBREON, 
     &	     INADAL, IMAIGNAN,
     &	     AER_ALPHA, AER_BETA, AER_GAMMA, AER_ZETA, 
     &	     IND_ANGOUT_FIN, IER)   
	
C --- Gestion du cas d'erreur
      
      IF (IER.NE.0) THEN ! cas d'erreur
      		
         IF (TRACE) THEN
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "   ==> Error in subroutine SOS_PREPA_OS"
         ENDIF
	  GOTO 4400   
	  
      ENDIF      
   
    
C* ----------------------------------------------------------------------------------------- 
C* Calcul par approche CKD de base: 1 calcul par terme de la décomposition en exponentielles 
C* -----------------------------------------------------------------------------------------    
C --- Le profil et les luminances sont calculés pour tous les coefficients CKD de chaque espèce.
C --- Il y a 8 absorbants => la luminance est calculée pour chaque produit d'exponentielles sur les 8 absorbants
C --- Dans le cas où l'absorption n'est pas prise en compte, le nombre d'exponentielles pour chaque composant vaut 1

      IF (IMODE_CKD_CALCUL.EQ.1) THEN  ! CAS CALCUL LUM POUR EP OPT UNIQUE PAR CKD

         WRITE(6,*)  
     &	  "==> Radiance computation for each CKD term" //
     &   " of the exponential expansion"
         DO IK1=1,NEXP(1,LAMB1)
         DO IK2=1,NEXP(2,LAMB1)
         DO IK3=1,NEXP(3,LAMB1)
         DO IK4=1,NEXP(4,LAMB1)
         DO IK5=1,NEXP(5,LAMB1)
         DO IK6=1,NEXP(6,LAMB1)
         DO IK7=1,NEXP(7,LAMB1)
         DO IK8=1,NEXP(8,LAMB1)


C ---       Affichage des infos à l'écran
            write(6,*) '     IK1 (H2O) ',IK1,'/',NEXP(1,LAMB1)
            write(6,*) '     IK2 (CO2) ',IK2,'/',NEXP(2,LAMB1)
            write(6,*) '     IK3 (O3)  ',IK3,'/',NEXP(3,LAMB1)
            write(6,*) '     IK4 (N2O) ',IK4,'/',NEXP(4,LAMB1)
            write(6,*) '     IK5 (CO)  ',IK5,'/',NEXP(5,LAMB1)
            write(6,*) '     IK6 (CH4) ',IK6,'/',NEXP(6,LAMB1)
            write(6,*) '     IK7 (O2)  ',IK7,'/',NEXP(7,LAMB1)
            write(6,*) '     IK8 (NO2) ',IK8,'/',NEXP(8,LAMB1)
            write(6,*) '     --------------------------------'
	   
C ---       Calcul du produit des fonctions de répartition (poids affecté aux luminances)
            AIK=KDIS_AI(IK1,1,LAMB1)*KDIS_AI(IK2,2,LAMB1)*
     &          KDIS_AI(IK3,3,LAMB1)*KDIS_AI(IK4,4,LAMB1)*
     &          KDIS_AI(IK5,5,LAMB1)*KDIS_AI(IK6,6,LAMB1)*
     &	         KDIS_AI(IK7,7,LAMB1)*KDIS_AI(IK8,8,LAMB1)

C ---       Renormalisation (gestion des cas pour lesquels la somme est un peu différente de 1.)
	     AIK=AIK/sumAIK

            !write(6,*) 'AIK ',AIK

	   
C ---       Calcul du profil d'épaisseurs optiques d'absorption en fonction de l'altitude

            CALL SOS_ABSPROFILE(ABSPROFIL, NU, LAMB1, IABS, 
     &		               ABS_USERPROFIL, ALTABS, RO, NEXP,
     &		               KDIS_KI, KDIS_KI_H2O, 
     &		               IK1,IK2,IK3,IK4,IK5,IK6,IK7,IK8,
     &                       TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                       TAB_CONC_H2O, NB_CONC_H2O,
     &			        TAUABS, TRACEPROF, IDLOGPROF,IER) 
                  


C ---       Gestion du cas d'erreur      
            IF (IER.NE.0) THEN ! cas d'erreur
      		
               IF (TRACE) THEN
                   WRITE (IDLOG,*,ERR=1303) 
     &	            "   ==> Error in subroutine SOS_ABSPROFILE"
 	        ENDIF
	        GOTO 4500 
	  
            ENDIF

       				
            ! Lancement de l'outil de calcul du profil atmo en épaisseurs optiques

            CALL SOS_PROFILE(IPROFIL, TR, HR, TA, HA, ZMIN,  
     &		               ZMAX, ABSPROFIL, ALTABS, TAUABS, 
     &			        TRACEPROF, IDLOGPROF, FICPROFIL_TMP,
     &			        NT, IER)


            ! Gestion du cas d'erreur
            IF (IER.NE.0) THEN ! cas d'erreur
      		
                IF (TRACE) THEN
	             WRITE (IDLOG,*,ERR=1303) 
     &	             "   ==> Error in subroutine SOS_PROFILE"
                ENDIF
	         GOTO 4600   
	  
            ENDIF

 
            !Lancement de l'outil de calcul des luminances

            CALL SOS(FICOS_TMP, FICTRANS, FICPROFIL_TMP,
     &	             NT, ZOUT, IGMAX, IPOLAR, RON,         
     &	             SURF_IND, RHO, IMAT_SURF, IFRESNEL, FICSURF,
     &		     N0, PIZ, PIZTR, COEFTRONCA, 
     &		     RMU, GA, TETAS, OS_NB, LUM_NBMU,
     &		     AER_ALPHA, AER_BETA, AER_GAMMA, AER_ZETA,
     &		     TTOT_TRONC_TMP, TTOT_VRAI_TMP, TAUOUT_TMP, 
     &		     TDIFMUS_TMP, TDIFMUG_TMP, EMOINS_TMP, EPLUS_TMP,
     &               TRACEOS, IDLOGOS, IER)      
      
            ! Gestion du cas d'erreur
            IF (IER.NE.0) THEN ! cas d'erreur
      		
               IF (TRACE) THEN
	            WRITE (IDLOG,*,ERR=1303) 
     &	            "   ==> Error in subroutine SOS"
               ENDIF
	        GOTO 4700   
	  
            ENDIF

            ! Lecture des résultats pour l'octuple ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8
            ! Pondération et aggrégation des éléments du vecteur de Stokes
            ! Pondération et aggrégation des transmissions
            ! Pondération et aggrégation des épaisseurs optiques
	   
            CALL SOS_AGGREGATE(LUM_NBMU, AIK, FICOS_TMP, 
     &           TTOT_TRONC_TMP, TTOT_VRAI_TMP, 
     &           TAUOUT_TMP, TDIFMUS_TMP, TDIFMUG_TMP,
     &           EMOINS_TMP, EPLUS_TMP,
     &           FICOS_AGG_TMP,
     &           FICSOS_RES_BIN,
     &           TTOT_TRONC, TTOT_VRAI,           
     &           TAUOUT, TDIFMUS, TDIFMUG, 
     &           EMOINS, EPLUS, IER)


            ! Gestion du cas d'erreur      
            IF (IER.NE.0) THEN ! cas d'erreur
      		
               IF (TRACE) THEN
	            WRITE (IDLOG,*,ERR=1303) 
     &	            "   ==> Error in subroutine SOS_AGGREGATE"
               ENDIF
	        GOTO 4800   
	  
            ENDIF

C ---    Fin des boucles sur les coeff CKD de chaque élément absorbant
         ENDDO  !IK1
         ENDDO  !IK2
         ENDDO  !IK3
         ENDDO  !IK4
         ENDDO  !IK5
         ENDDO  !IK6
         ENDDO  !IK7
         ENDDO  !IK8

         ! Destruction du fichier temporaire
	  OPEN(UNIT=3,FILE=FICOS_TMP,FORM='UNFORMATTED') 
         CLOSE(3,STATUS='DELETE')    

      ENDIF ! Fin CAS CALCUL LUMINANCE POUR CHAQUE EP OPT CKD


C* ----------------------------------------------------------------------------------------- 
C* Calcul par approche CKD appliquée au calcul de transmission gazeuse uniquement: 
C* ==> Approche simplifiée, approximative, avec 1 calcul d'épaisseur optique d'absorption unique
C* -----------------------------------------------------------------------------------------    
C --- Le profil et les luminances sont calculés pour l'épaisseur optique d'absorption estimée globalement.

      IF (IMODE_CKD_CALCUL.EQ.2) THEN  ! CAS CALCUL LUM POUR EP OPT ESTIMEE PAR CKD

         WRITE(6,*)  
     &	  "==> Radiance computation for the absorption optical " //
     &   "thickness calculated by the CKD expansion approach"
         
         DO I=1,CTE_ABS_NBLEV
            TRSCKD(I)=0.D+00
         ENDDO

         DO IK1=1,NEXP(1,LAMB1)
         DO IK2=1,NEXP(2,LAMB1)
         DO IK3=1,NEXP(3,LAMB1)
         DO IK4=1,NEXP(4,LAMB1)
         DO IK5=1,NEXP(5,LAMB1)
         DO IK6=1,NEXP(6,LAMB1)
         DO IK7=1,NEXP(7,LAMB1)
         DO IK8=1,NEXP(8,LAMB1)

	   
C ---       Calcul du produit des fonctions de répartition (poids affecté aux luminances)
            AIK=KDIS_AI(IK1,1,LAMB1)*KDIS_AI(IK2,2,LAMB1)*
     &          KDIS_AI(IK3,3,LAMB1)*KDIS_AI(IK4,4,LAMB1)*
     &          KDIS_AI(IK5,5,LAMB1)*KDIS_AI(IK6,6,LAMB1)*
     &	         KDIS_AI(IK7,7,LAMB1)*KDIS_AI(IK8,8,LAMB1)

C ---       Renormalisation (gestion des cas pour lesquels la somme est un peu différente de 1.)
	     AIK=AIK/sumAIK
	   
C ---       Calcul du profil d'épaisseurs optiques d'absorption en fonction de l'altitude
            CALL SOS_ABSPROFILE(ABSPROFIL, NU, LAMB1, IABS, 
     &		               ABS_USERPROFIL, ALTABS, RO, NEXP,
     &		               KDIS_KI, KDIS_KI_H2O,
     &		               IK1,IK2,IK3,IK4,IK5,IK6,IK7,IK8,
     &                       TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                       TAB_CONC_H2O, NB_CONC_H2O,
     &			        TAUABS, TRACEPROF, IDLOGPROF,IER) 

C ---       Gestion du cas d'erreur      
            IF (IER.NE.0) THEN ! cas d'erreur
      		
               IF (TRACE) THEN
                   WRITE (IDLOG,*,ERR=1303) 
     &	            "   ==> Error in subroutine SOS_ABSPROFILE"
 	        ENDIF
	        GOTO 4500
	  
            ENDIF

            !Calcul de la transmission par niveau I, pour masse d'air fictive de 1 selon la méthode CKD
            DO I=1,CTE_ABS_NBLEV
               TRSCKD(I)=TRSCKD(I)+AIK*DEXP(-TAUABS(I))
            ENDDO              
	  
C ---    Fin des boucles sur les coeff CKD de chaque élément absorbant
         ENDDO  !IK1
         ENDDO  !IK2
         ENDDO  !IK3
         ENDDO  !IK4
         ENDDO  !IK5
         ENDDO  !IK6
         ENDDO  !IK7
         ENDDO  !IK8

         ! Estimation du profil d'épaisseur optique d'absorption
         DO I=1,CTE_ABS_NBLEV
             TAUABS(I)=-DLOG(TRSCKD(I))
             if (TAUABS(I).LT.0.D+00) TAUABS(I)=0.D+00	!MODIF TEST : pour gestion des arrondis donnant transmission très légèrement sup à 1
         ENDDO              

          ! Lancement de l'outil de calcul du profil atmo en épaisseurs optiques

         CALL SOS_PROFILE(IPROFIL, TR, HR, TA, HA, ZMIN,  
     &		            ZMAX, ABSPROFIL, ALTABS, TAUABS, 
     &			     TRACEPROF, IDLOGPROF, FICPROFIL_TMP,
     &			     NT, IER)

         ! Gestion du cas d'erreur
         IF (IER.NE.0) THEN ! cas d'erreur
      		
              IF (TRACE) THEN
	           WRITE (IDLOG,*,ERR=1303) 
     &	           "   ==> Error in subroutine SOS_PROFILE"
              ENDIF
	       GOTO 4600   
	  
         ENDIF

 
         !Lancement de l'outil de calcul des luminances

         CALL SOS(FICSOS_RES_BIN, FICTRANS, FICPROFIL_TMP,
     &         NT, ZOUT, IGMAX, IPOLAR, RON,         
     &         SURF_IND, RHO, IMAT_SURF, IFRESNEL, FICSURF,
     &         N0, PIZ, PIZTR, COEFTRONCA, 
     &         RMU, GA, TETAS, OS_NB, LUM_NBMU,
     &         AER_ALPHA, AER_BETA, AER_GAMMA, AER_ZETA,
     &         TTOT_TRONC, TTOT_VRAI, TAUOUT, 
     &         TDIFMUS, TDIFMUG, EMOINS, EPLUS,
     &         TRACEOS, IDLOGOS, IER)      
      
         ! Gestion du cas d'erreur
         IF (IER.NE.0) THEN ! cas d'erreur
      		
             IF (TRACE) THEN
	          WRITE (IDLOG,*,ERR=1303) 
     &	          "   ==> Error in subroutine SOS"
             ENDIF
	      GOTO 4700   
	  
         ENDIF

         ! Destruction du fichier de profil temporaire
         OPEN(UNIT=3,FILE=FICPROFIL_TMP,FORM='UNFORMATTED') 
         CLOSE(3,STATUS='DELETE') 

      ENDIF ! Fin CAS CALCUL LUM POUR EP OPT UNIQUE PAR CKD



C* ----------------------------------------------------------------------------------------- 
C* Recomposition de la série de Fourier pour un azimuth de visée donné
C* ==> Génération des fichiers résultats
C* -----------------------------------------------------------------------------------------    
C --- Création des fichiers de sortie (I,Q,U)
      WRITE(6,*) "FINALYSING COMPUTATIONS AND WRITING OUTPUT FILES"


      IF (TRACEOS) THEN
       WRITE(IDLOGOS,*) 
       WRITE(IDLOGOS,*)"***********************************************"
       WRITE(IDLOGOS,*)"Recombination of the radiation field by azimuth"
       WRITE(IDLOGOS,*)"***********************************************"
       WRITE(IDLOGOS,*) 
       IF (ITRPHI.EQ.1) THEN
	  WRITE(IDLOGOS,*) 
     &       "Recombination for the azimuth (degrees) : ",PHIOS
       ELSE
	  WRITE(IDLOGOS,*)
     &      "Recombination for azimuths from 0 to 360 " //
     &      "with a step (degrees): ",PAS_PHI
       ENDIF   
       WRITE(IDLOGOS,*)
      ENDIF    !TRACE
    

      CALL SOS_TRPHI_OPTION(LUM_NBMU, RMU, GA,
     &     FICSOS_RES_BIN, TTOT_TRONC, TAUOUT, ZOUT,
     &     IGLI, N0,WIND, SURF_IND, IFRESNEL,
     &     IROUJEAN, K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN,
     &     IRONDEAUX, IBREON,
     &     INADAL, ALPHA_NADAL, BETA_NADAL,
     &     IMAIGNAN, COEF_C_MAIGNAN,
     &     ITRPHI ,PHIOS, PAS_PHI, IPOLAR,
     &     PHI_FIN, THETA_FIN, 
     &     SCA_UP_FIN, I_UP_FIN, Q_UP_FIN, U_UP_FIN,
     &     POL_ANG_UP_FIN, POL_RATE_UP_FIN, L_POL_UP_FIN,
     &     SCA_DOWN_FIN, I_DOWN_FIN, Q_DOWN_FIN, U_DOWN_FIN,
     &     POL_ANG_DOWN_FIN, POL_RATE_DOWN_FIN, L_POL_DOWN_FIN,     
     &     IER)
           
C --- Gestion du cas d'erreur
      
      IF (IER.NE.0) THEN ! cas d'erreur
      		
         IF (TRACEOS) THEN
	         WRITE (IDLOGOS,*,ERR=1303) 
     &	         "   ==> Error in subroutine SOS_TRPHI_OPTION"
         ENDIF
         GOTO 4900   
	  
      ENDIF

C* ----------------------------------------------------------------------------------------- 
C* Sauvegarde des transmissions dans un fichier
C* -----------------------------------------------------------------------------------------    
      IF (FICTRANS.NE.'NO_OUTPUT') THEN
      
         OPEN(UNIT=20,FILE=FICTRANS,FORM='FORMATTED',ERR=997)

C     --> Transmission directe pour épaisseur optique vraie (non tronquée)
  
	  TDIR_TRONC = 
     &      DEXP(-TTOT_TRONC/DCOS(INCTE_PI*TETAS/180.D+00))
	  TDIR_VRAI = 
     &      DEXP(-TTOT_VRAI/DCOS(INCTE_PI*TETAS/180.D+00))
	  
	  WRITE(20,1005,err=998) TETAS 
	  WRITE(20,1006,err=998) TDIR_VRAI
	  WRITE(20,*,err=998) ' '

C     --> Transmission diffuse TOA --> Sol	 
	
	  WRITE(20,*,err=998)'Diffuse transmittance : TOA -> surface'
	  WRITE(20,1010,err=998)TETAS,(TDIFMUS + TDIR_TRONC - TDIR_VRAI)
	  WRITE(20,*,err=998) ' '
	  
C     --> Transmission diffuse Sol --> TOA : 
C         Boucle sur les valeurs IMUS : N0 = J	pour calcul TOA --> Sol (reciprocite)  
         WRITE(20,*,err=998)'Diffuse transmittance : surface -> TOA'
	  DO J=1,LUM_NBMU
	  
	     TDIR_TRONC = DEXP(-TTOT_TRONC/RMU(J))
	     TDIR_VRAI  = DEXP(-TTOT_VRAI/RMU(J))
	  	      
	     WRITE(20,2010,err=998) DACOS(RMU(J))*180.D+00/INCTE_PI, 
     &	                           (TDIFMUG(J) + TDIR_TRONC - TDIR_VRAI)
         ENDDO	  

	  CLOSE(20)
  
      ENDIF !Fin test si fichier transmission à créer

C* ----------------------------------------------------------------------------------------- 
C* Sauvegarde des flux TOA et BOA dans un fichier
C* -----------------------------------------------------------------------------------------  

!     On récupère les grandeurs TDIR_VRAI, FLUX_DIFF_DOWN et FLUX_DOWN pour les fournir 
!     en sorties de SOS_PROC (même si on ne les enregistre pas dans un fichier optionnel).
!     --> Utile à l'appel de la routine SOS_PROC via binding Python :

C     Transmission directe pour épaisseur optique tronquée
      TDIR_TRONC = DEXP(-TTOT_TRONC/DCOS(INCTE_PI*TETAS/180.D+00))
C     Transmission directe pour épaisseur optique vraie (non tronquée)
      TDIR_VRAI = DEXP(-TTOT_VRAI/DCOS(INCTE_PI*TETAS/180.D+00))
C     Flux diffus descendant
      FLUX_DIFF_DOWN = EMOINS + TDIR_TRONC - TDIR_VRAI
C     Flux total descendant
      FLUX_DOWN = EMOINS + TDIR_TRONC

      IF (FICFLUX.NE.'NO_OUTPUT') THEN
      
          OPEN(UNIT=30,FILE=FICFLUX,FORM='FORMATTED',ERR=987)
	  
	   WRITE(30,1005,err=988)TETAS  	                       
	   WRITE(30,*,err=988) ' '

C     --> Flux diffus descendant au BOA: TOA --> Sol	 
	   WRITE(30,*,err=988)
     &          'Downward fluxes at BOA (normalized by TOA solar flux)'	  
	   WRITE(30,2016,err=988) TDIR_VRAI
	   WRITE(30,2017,err=988) FLUX_DIFF_DOWN
	   WRITE(30,2018,err=988) FLUX_DOWN
	   WRITE(30,*,err=988) ' '
           
C     --> Flux diffus montant au TOA: Sol --> TOA : 
          WRITE(30,*,err=988)
     &    'Upward diffuse flux at TOA (normalized by TOA solar flux):',
     &    EPLUS

          WRITE(30,*,err=988)
          WRITE(30,*,err=988)
          WRITE(30,*,err=988) 'According to the following profile'
          WRITE(30,*,err=988) 'Z(km)    MOT     AOT     GOT     TOTAL'
          DO I = CTE_ABS_NBLEV, 1 , -1
             ZALT = ABS_USERPROFIL(I,1)
             TR_ZALT = TR*DEXP(-ZALT/HR)
             TA_ZALT = TA*DEXP(-ZALT/HA)
             TG_ZALT = TAUABS(CTE_ABS_NBLEV+1-I)
             WRITE(30,2020,err=988) ZALT, TR_ZALT , TA_ZALT, TG_ZALT,
     &                              (TR_ZALT + TA_ZALT + TG_ZALT)
          ENDDO

          CLOSE(30)
	  
      ENDIF !Fin test si fichier de flux a creer
     
C* Fin nominale 
C-------------------
      IF (TRACEPROF) CLOSE(IDLOGPROF)
      IF (TRACEOS)   CLOSE(IDLOGOS)
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &   "--> Nominal end of SOS_PROC processes"
         WRITE (IDLOG,*,ERR=1303) "JOB_STATUS=OK"
         CLOSE(IDLOG)
      ENDIF
      GOTO 9999
      

C =========================================================
C --- Cas d'erreur et retour du status 1 au shell  
C =========================================================  

C --- Définition Arborescence de sortie
 1000 WRITE(6,*) '  SOS_PROC : ERROR_1000 on parameters'
      WRITE(6,*) '     RESROOT parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS_Main.ResRoot keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Working folder of SOS computations'
      GOTO 9998 

 1001 WRITE(6,*) '  SOS_PROC : ERROR_1001 on parameters'
      WRITE(6,*) '     DIR_MIE parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.DirMie keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Storage directory for MIE files'
      GOTO 9998 

 1002 WRITE(6,*) '  SOS_PROC : ERROR_1002 on parameters'
      WRITE(6,*) '     DIR_SURF parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Dir keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Storage directory for SURFACE files'
      GOTO 9998 
      
 1100 WRITE(6,*) '  SOS_PROC : ERROR_1100 : '
      WRITE(6,*) 
     &'  => Working folder name : maybe too long or not defined'
      WRITE(6,*) '     RESROOT parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS_Main.ResRoot keyword parameter)'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1101 WRITE(6,*) '  SOS_PROC : ERROR_1101 : '
      WRITE(6,*) '  => Name of storage directory for MIE files :'
      WRITE(6,*) '     maybe too long or not defined'
      WRITE(6,*) '        DIR_MIE parameter in SOS_PROC'
      WRITE(6,*) '        (-AER.DirMie keyword parameter)'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1102 WRITE(6,*) '  SOS_PROC : ERROR_1102 : '
      WRITE(6,*) '  => Surface storage directory name : '
      WRITE(6,*) '     maybe too long or not defined'
      WRITE(6,*) '        DIR_SURF parameter in SOS_PROC'
      WRITE(6,*) '        (-SURF.Dir keyword parameter)'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1103 WRITE(6,*) '  SOS_PROC : ERROR_1103 : '
      WRITE(6,*) '  => Glitter storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1104 WRITE(6,*) '  SOS_PROC : ERROR_1104 : '
      WRITE(6,*) '  => Roujean storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1105 WRITE(6,*) '  SOS_PROC : ERROR_1105 : '
      WRITE(6,*) '  => Rondeaux storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1106 WRITE(6,*) '  SOS_PROC : ERROR_1106 : '
      WRITE(6,*) '  => Breon storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1107 WRITE(6,*) '  SOS_PROC : ERROR_1107 : '
      WRITE(6,*) '  => Nadal storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

11071 WRITE(6,*) '  SOS_PROC : ERROR_11071 : '
      WRITE(6,*) '  => Maignan storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998  

 1108 WRITE(6,*) '  SOS_PROC : ERROR_1108 : '
      WRITE(6,*) '  => SOS Working folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1109 WRITE(6,*) '  SOS_PROC : ERROR_1109 : '
      WRITE(6,*) '  => Log folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1110 WRITE(6,*) '  SOS_PROC : ERROR_1110 : '
      WRITE(6,*) '  => TMP folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998   

C --- Création arborescence de sortie

 1201 WRITE(6,*) '  SOS_PROC : ERROR_1201 : '
      WRITE(6,*) 
     &'  => Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1202 WRITE(6,*) '  SOS_PROC : ERROR_1202 : '
      WRITE(6,*) 
     &'  => Mie Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1203 WRITE(6,*) '  SOS_PROC : ERROR_1203 : '
      WRITE(6,*) 
     &'  => Surface Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1204 WRITE(6,*) '  SOS_PROC : ERROR_1204 : '
      WRITE(6,*) 
     &'  => Glitter Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1205 WRITE(6,*) '  SOS_PROC : ERROR_1205 : '
      WRITE(6,*) 
     &'  => Roujean Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1206 WRITE(6,*) '  SOS_PROC : ERROR_1206 : '
      WRITE(6,*) 
     &'  => Rondeaux Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1207 WRITE(6,*) '  SOS_PROC : ERROR_1207 : '
      WRITE(6,*) 
     &'  => Breon Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1208 WRITE(6,*) '  SOS_PROC : ERROR_1208 : '
      WRITE(6,*) 
     &'  => Nadal Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

12081 WRITE(6,*) '  SOS_PROC : ERROR_12081 : '
      WRITE(6,*) 
     &'  => Maignan Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998  

 1209 WRITE(6,*) '  SOS_PROC : ERROR_1209 : '
      WRITE(6,*)  
     &'  => SOS Directory results creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1210 WRITE(6,*) '  SOS_PROC : ERROR_1210 : '
      WRITE(6,*) 
     &'  => Log folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1211 WRITE(6,*) '  SOS_PROC : ERROR_1211 : '
      WRITE(6,*) 
     &'  => TMP folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998   

C --- gestion du log
                                           
 1301 WRITE(6,*) '  SOS_PROC : ERROR_1301 : '
      WRITE(6,*) '  => Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 1302 WRITE(6,*) '  SOS_PROC : ERROR_1302 on Log file opening'
      GOTO 9998

 1303 WRITE(6,*) '  SOS_PROC : ERROR_1303 on Log file writing'
      GOTO 9998
            
C --- Controle paramètres simulation

 2100 WRITE(6,*) '  SOS_PROC : ERROR_2100 on parameters'
      WRITE(6,*) '     WA_SIMU parameter is required in SOS_PROC'
      WRITE(6,*) '     (-SOS_Main.Wa keyword parameter)'
      WRITE(6,*) '  ==>  Wavelength for radiance simulation (microns)'
      GOTO 9998 

 2101 WRITE(6,*) '  SOS_PROC : ERROR_2101 on parameter value'
      WRITE(6,*) '     WA_SIMU in SOS_PROC (-SOS_Main.Wa keyword)'
      WRITE(6,*) '     is lower than the minimal available value'
      WRITE(6,*) '     or higher than the maximal available value'
      WRITE(6,*) '         WA_SIMU = ', WA_SIMU
      WRITE(6,*) '         CTE_WAMIN = ', CTE_WAMIN
      WRITE(6,*) '         CTE_WAMAX = ', CTE_WAMAX

C --- Controle des paramètres pour les calculs d'angles

 2200 WRITE(6,*) '  SOS_PROC : ERROR_2200 on parameters'
      WRITE(6,*) '     TETAS parameter is required in SOS_PROC'
      WRITE(6,*) '     (-ANG.Thetas keyword parameter)'
      WRITE(6,*) '  ==> Solar zenithal angle (degrees)'
      GOTO 9998 

 2201 WRITE(6,*) '  SOS_PROC : ERROR_2201 on parameter value'
      WRITE(6,*) '     TETAS in SOS_PROC (-ANG.Thetas keyword)'
      WRITE(6,*) '  --> Solar Zenith Angle out of available range'
      GOTO 9998  	

C --- Controle des paramètres de calculs d'aérosols

 2301 WRITE(6,*) '  SOS_PROC : ERROR_2301 on parameters'
      WRITE(6,*) '     AOT_REF parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.AOTref keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) 
     &'  ==> Aerosol optical thickness at reference wavelength '
      GOTO 9998   

 2302 WRITE(6,*) '  SOS_PROC : ERROR_2302 on parameters'
      WRITE(6,*) '     WAREF_AOT parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.Waref keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*)  
     &'  ==> Reference wavelength for the aerosol optical thickness'
      GOTO 9998   

 2304 WRITE(6,*) '  SOS_PROC : ERROR_2304 on parameters'
      WRITE(6,*) '     IMOD_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.Model keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Type of aerosol model : '
      WRITE(6,*) '         0 : Mono-modal '
      WRITE(6,*) '         1 : WMO multi-modal '   
      WRITE(6,*) '         2 : Shettle & Fenn bi-modal '    
      WRITE(6,*) '         3 : Log-Normal bi-modal '
      WRITE(6,*) '         4 : Phaze functions from an external source'
      WRITE(6,*) '         5 : Multi-modes model defined by user'
      GOTO 9998   

 2305 WRITE(6,*) '  SOS_PROC : ERROR_2305 on parameter value'
      WRITE(6,*) '     IMOD_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.Model keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range.'
      WRITE(6,*) '  ==> Type of aerosol model : '
      WRITE(6,*) '         0 : Mono-modal '
      WRITE(6,*) '         1 : WMO multi-modal '   
      WRITE(6,*) '         2 : Shettle & Fenn bi-modal '    
      WRITE(6,*) '         3 : Log-Normal bi-modal '
      WRITE(6,*) '         4 : Phaze functions from an external source'
      WRITE(6,*) '         5 : Mixture defined by user file'    
      GOTO 9998 
   
 2309 WRITE(6,*) '  SOS_PROC : ERROR_2309 on parameters'
      WRITE(6,*) '     RN_WA and IN_WA parameters in SOS_PROC'
      WRITE(6,*) '     keyword parameter:'
      WRITE(6,*) '      -AER.MMD.MRwa Value (format : F5.3)'
      WRITE(6,*) '      -AER.MMD.MIwa Value (format : F8.5)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) 
     &'  ==> Aerosol refractive index for mono-modal size distribution'
      WRITE(6,*) '      (real and imaginary parts)'
      WRITE(6,*) '      at the wavelength of radiance calculation.'
      GOTO 9998   

 2310 WRITE(6,*) '  SOS_PROC : ERROR_2310 on parameter value'
      WRITE(6,*) '     IN_WA parameter in SOS_PROC'
      WRITE(6,*) '     (AER.MMD.MIwa keyword parameter)'
      WRITE(6,*) '  ==> Imagery part of Aerosol refractive index '
      WRITE(6,*) '      must be a negative value'
      GOTO 9998   

 2311 WRITE(6,*) '  SOS_PROC : ERROR_2311 on parameters'
      WRITE(6,*) '     IGRANU parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.MMD.SDtype keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Model index for mono-modal size distribution :'
      WRITE(6,*) '         1 : LND '
      WRITE(6,*) '         2 : Junge law '
      GOTO 9998    

 2312 WRITE(6,*) '  SOS_PROC : ERROR_2312 on parameter value'
      WRITE(6,*) '     IGRANU parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.MMD.SDtype keyword parameter - Integer)'
      WRITE(6,*) '  ==> Model index for mono-modal size distribution :'
      WRITE(6,*) '         1 : LND '
      WRITE(6,*) '         2 : Junge law '
      GOTO 9998      

 2314 WRITE(6,*) '  SOS_PROC : ERROR_2314 on parameters'
      WRITE(6,*) '     RN_WAREF and IN_WAREF parameters in SOS_PROC'
      WRITE(6,*) '     keyword parameter:'
      WRITE(6,*) '      -AER.MMD.MRwaref Value (format : F5.3)'
      WRITE(6,*) '      -AER.MMD.MIwaref Value (format : F8.5)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Aerosol refractive index '
      WRITE(6,*) '      for mono-modal size distribution'
      WRITE(6,*) '      (real and imaginary parts)'
      WRITE(6,*) '      at the wavelength of reference.'
      GOTO 9998  

23131 WRITE(6,*) '  SOS_PROC : ERROR_23131 on parameters'
      WRITE(6,*) 
     &'     LND_RADIUS_MMD_AER and LND_LNVAR_MMD_AER in SOS_PROC'
      WRITE(6,*) '     (-AER.MMD.LNDradius keyword parameter)'
      WRITE(6,*) '     (-AER.MMD.LNDvar keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  Case IGRANU = 1 (AER.MMD.SDtype = 1)'
      WRITE(6,*) '  ==> Parameters of LND mono-modal aerosol size '
      WRITE(6,*) '      distribution :'
      WRITE(6,*) '         - modal radius (microns) '
      WRITE(6,*) '         - and standard deviation'
      GOTO 9998  	     

23132 WRITE(6,*) '  SOS_PROC : ERROR_23132 on parameters'
      WRITE(6,*) 
     &'     JD_SLOPE_MMD_AER and JD_RMIN_MMD_AER in SOS_PROC'
      WRITE(6,*) '     -AER.MMD.JD.slope keyword parameter (Float)'
      WRITE(6,*) '     -AER.MMD.JD.rmin  keyword parameter (Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  Case IGRANU = 2 (AER.MMD.SDtype = 2)'
      WRITE(6,*) '  ==> Parameters of Junge law mono-modal'
      WRITE(6,*) '      aerosol size distribution :'
      WRITE(6,*) '         - slope'
      WRITE(6,*) '         - minimal radius (microns)'
      GOTO 9998 

23141 WRITE(6,*) '  SOS_PROC : ERROR_23141 on parameters'
      WRITE(6,*) '     ITRONC_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.Tronca keyword parameter)'
      WRITE(6,*) '  --> Defined value of truncation option (0 or 1)'
      GOTO 9998     

 2315 WRITE(6,*) '  SOS_PROC : ERROR_2315 on parameters'
      WRITE(6,*) '     IModele_WMO parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.WMO.Model keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Type of WMO model :'
      WRITE(6,*) '         1 : Continental WMO model. '
      WRITE(6,*) '         2 : Maritime WMO model. '
      WRITE(6,*) '         3 : Urban WMO model. '
      WRITE(6,*) '         4 : WMO model by user definition. '
      GOTO 9998    
         
 2316 WRITE(6,*) '  SOS_PROC : ERROR_2316 on parameter value'
      WRITE(6,*) '     IModele_WMO parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.WMO.Model keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range'
      WRITE(6,*) '  ==> Type of WMO model :'
      WRITE(6,*) '         1 : Continental WMO model. '
      WRITE(6,*) '         2 : Maritime WMO model. '
      WRITE(6,*) '         3 : Urban WMO model. '
      WRITE(6,*) '         4 : WMO model by user definition. '
      GOTO 9998      		  

 2317 WRITE(6,*) '  SOS_PROC : ERROR_2317 on parameters'
      WRITE(6,*) '     C_WMO_DL, C_WMO_WS,C_WMO_OC '
      WRITE(6,*) '     and C_WMO_SO parameters in SOS_PROC'
      WRITE(6,*) '     (-AER.WMO.DL keyword parameter - float)'
      WRITE(6,*) '     (-AER.WMO.WS keyword parameter - float)'
      WRITE(6,*) '     (-AER.WMO.OC keyword parameter - float)'
      WRITE(6,*) '     (-AER.WMO.SO keyword parameter - float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Volumetric concentration (%) for '
      WRITE(6,*) '         WMO Dust-Like particles'
      WRITE(6,*) '         WMO Water Soluble particles'	    
      WRITE(6,*) '         WMO OCeanic particles'
      WRITE(6,*) '         WMO SOot particles'	     
      GOTO 9998 

 2318 WRITE(6,*) '  SOS_PROC : ERROR_2318 on parameters'
      WRITE(6,*) '     IModele_SF parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.SF.Model keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Type of Shettle & Fenn model :'		  
      WRITE(6,*) '         1 : Tropospheric S&F model. '
      WRITE(6,*) '         2 : Urban S&F model. '
      WRITE(6,*) '         3 : Maritime S&F model. '
      WRITE(6,*) '         4 : Coastal S&F model. '
      GOTO 9998    		 

 2319 WRITE(6,*) '  SOS_PROC : ERROR_2319 on parameters'
      WRITE(6,*) '     RH parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.SF.RH keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> % of relative humidity'
      GOTO 9998 

 2320 WRITE(6,*) '  SOS_PROC : ERROR_2320 on parameter value'
      WRITE(6,*) '     IModele_SF parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.SF.Model keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range'
      WRITE(6,*) '  ==> Type of Shettle & Fenn model :'		  
      WRITE(6,*) '         1 : Tropospheric S&F model. '
      WRITE(6,*) '         2 : Urban S&F model. '
      WRITE(6,*) '         3 : Maritime S&F model. '
      WRITE(6,*) '         4 : Coastal S&F model. '
      GOTO 9998    		 

 2321 WRITE(6,*) '  SOS_PROC : ERROR_2321 on parameters'
      WRITE(6,*) '  --> -AER.SF.RH Value    (format : F5.2): '
      WRITE(6,*) '  ==> Relative humidity: between 0 and 99%'
      GOTO 9998   		    
     
 2322 WRITE(6,*) '  SOS_PROC : ERROR_2322 on parameters'
      WRITE(6,*) '     MODE_PARAM_BILND parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.VCdef keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Type of selected method for bi-modal LND'
      WRITE(6,*) '      volumetric concentration definition :'
      WRITE(6,*) '       1 : Use of user volumetric concentrations.'
      WRITE(6,*) '       2 : Use of user ratio AOT_coarse / AOT_total.'
      GOTO 9998   		

 2323 WRITE(6,*) '  SOS_PROC : ERROR_2323 on parameter value'
      WRITE(6,*) '     MODE_PARAM_BILND parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.VCdef keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range'
      WRITE(6,*) '  ==> Type of selected method for bi-modal LND'
      WRITE(6,*) '      volumetric concentration definition :'
      WRITE(6,*) '       1 : Use of user volumetric concentrations.'
      WRITE(6,*) '       2 : Use of user ratio AOT_coarse / AOT_total.'
      GOTO 9998   		    

 2324 WRITE(6,*) '  SOS_PROC : ERROR_2324 on parameters'
      WRITE(6,*) '     USER_CV_COARSE parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.CoarseVC keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Volumetric concentration for bimodal LND'      
      GOTO 9998   		  

 2325 WRITE(6,*) '  SOS_PROC : ERROR_2325 on parameters'
      WRITE(6,*) '     USER_CV_FINE parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.FineVC keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Volumetric concentration for bimodal LND'      
      GOTO 9998   		  

 2326 WRITE(6,*) '  SOS_PROC : ERROR_2326 on parameters'
      WRITE(6,*) '     RTAUct_WAREF parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.RAOT keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> User volumetric concentration for bimodal LND'  
      WRITE(6,*) '      defined by AOT ratio coarse / total'  
      GOTO 9998  

 2327 WRITE(6,*) '  SOS_PROC : ERROR_2327 on parameters'
      WRITE(6,*) '     BMD_CM_MRWA, BMD_CM_MIWA, BMD_CM_RMODAL '
      WRITE(6,*) '     and BMD_CM_VAR in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.CM.MRwa keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.CM.MIwa keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.CM.SDradius keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.CM.SDvar keyword parameter -Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Parameters of the bimode LND coarse mode '
      WRITE(6,*) '      for wavelength -SOS_Main.Wa (WA_SIMU)'
      GOTO 9998  

 2328 WRITE(6,*) '  SOS_PROC : ERROR_2328 on parameters'
      WRITE(6,*) '     BMD_FM_MRWA, BMD_FM_MIWA, BMD_FM_RMODAL '
      WRITE(6,*) '     and BMD_FM_VAR in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.FM.MRwa keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.FM.MIwa keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.FM.SDradius keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.FM.SDvar keyword parameter -Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Parameters of the bimode LND fine mode '
      WRITE(6,*) '      for wavelength -SOS_Main.Wa (WA_SIMU)'
      GOTO 9998  

 2329 WRITE(6,*) '  SOS_PROC : ERROR_2329 on parameters'
      WRITE(6,*) '     BMD_CM_MRWAREF, BMD_CM_MIWAREF,'
      WRITE(6,*) '     BMD_FM_MRWAREF, BMD_FM_MIWAREF in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.CM.MRwaref keyword parameter - Float)'
      WRITE(6,*) '     (-AER.BMD.CM.MIwaref keyword parameter - Float)'
      WRITE(6,*) '     (-AER.BMD.FM.MRwaref keyword parameter - Float)'
      WRITE(6,*) '     (-AER.BMD.FM.MIwaref keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Refractive index at the reference wavelength' 
      WRITE(6,*) '      for the bimode LND coarse and fine modes'
      GOTO 9998     

 2330 WRITE(6,*) '  SOS_PROC : ERROR_2330 on parameters'
      WRITE(6,*) '     FICEXTDATA_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.ExtData keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> User phase function filename'  
      GOTO 9998  

 2331 WRITE(6,*) '  SOS_PROC : ERROR_2331 on parameter value'
      WRITE(6,*)  
     & '  --> Simulation and reference wavelengths must be equal'
      WRITE(6,*) '    WA_SIMU and WAREF_AOT parameters in SOS_PROC'
      WRITE(6,*) '    (-SOS_Main.Wa and -AER.Waref keyword parameters)'
      WRITE(6,*) '  ==> Case of an external User phase function'  
      GOTO 9998  

 2340 WRITE(6,*) '  SOS_PROC : ERROR_2340 on parameters'
      WRITE(6,*) '     FICMIXTURE_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.DefMixture keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) 
     & '  ==> User definition of the mixture of aerosol models'
      GOTO 9998  

 2350 WRITE(6,*) '  SOS_PROC : ERROR_2350 on parameter value'
      WRITE(6,*)  
     & '  --> Simulation and reference wavelengths must be equal'
      WRITE(6,*) '    WA_SIMU and WAREF_AOT parameters in SOS_PROC'
      WRITE(6,*) '    (-SOS_Main.Wa and -AER.Waref keyword parameters)'
      WRITE(6,*) '  ==> Case of an aerosol User File'  
      GOTO 9998  

 2351 WRITE(6,*) '  SOS_PROC : ERROR_2351 on parameter value'
      WRITE(6,*) '     FIC_GRANU parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.ResFile keyword parameter)'
      WRITE(6,*)  
     &'     is defined for the case of using an aerosol User File'
      WRITE(6,*) '     This is inconsistent.'
      WRITE(6,*)  
     &'  --> Parameter -AER.ResFile (FIC_GRANU) must not be defined !'
      GOTO 9998   


C --- Controle des paramètres de surface

 2401 WRITE(6,*) '  SOS_PROC : ERROR_2401 on parameters'
      WRITE(6,*) '     RHO parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Alb keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Surface albedo'
      GOTO 9998 

 2402 WRITE(6,*) '  SOS_PROC : ERROR_2402 on parameter value'
      WRITE(6,*) '     RHO parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Alb keyword parameter - Float)'
      WRITE(6,*) '  --> must be positive '
      WRITE(6,*) '  ==>  Surface albedo'
      GOTO 9998 
 
 2403 WRITE(6,*) '  SOS_PROC : ERROR_2403 on parameters'
      WRITE(6,*) '     ISURF parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Type keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Type of surface model'
      WRITE(6,*) '         0 : lambertian reflectance '
      WRITE(6,*) '         1 : lambertian reflectance + glitter'
      WRITE(6,*) '         2 : lambertian reflectance + fresnel'
      WRITE(6,*) '         3 : lambertian reflectance + Roujean BRDF'
      WRITE(6,*) '         4 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Rondeaux'
      WRITE(6,*) '         5 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Breon'  
      WRITE(6,*) '         6 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Nadal'     
      WRITE(6,*) '         7 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Maignan'             
      GOTO 9998 

 2404 WRITE(6,*) '  SOS_PROC : ERROR_2404 on parameter value'
      WRITE(6,*) '     ISURF parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Type keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range'
      WRITE(6,*) '  ==>  Type of surface model'
      WRITE(6,*) '         0 : lambertian reflectance '
      WRITE(6,*) '         1 : lambertian reflectance + glitter'
      WRITE(6,*) '         2 : lambertian reflectance + fresnel'
      WRITE(6,*) '         3 : lambertian reflectance + Roujean BRDF'
      WRITE(6,*) '         4 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Rondeaux'
      WRITE(6,*) '         5 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Breon'  
      WRITE(6,*) '         6 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Nadal'     
      WRITE(6,*) '         7 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Maignan'             
      GOTO 9998 

 2405 WRITE(6,*) '  SOS_PROC : ERROR_2405 on parameter'
      WRITE(6,*) '     SURF_IND parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Ind keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '         Required if surface type=1 2 4 5 6 or 7'
      GOTO 9998 

 2406 WRITE(6,*) '  SOS_PROC : ERROR_2406 on parameter value'
      WRITE(6,*) '     WIND parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Glitter.Wind keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Wind speed in m/s (required if type=1)'
      GOTO 9998 

24061 WRITE(6,*) '  SOS_PROC : ERROR_24061 on parameter value'
      WRITE(6,*) '     WIND parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Glitter.Wind keyword parameter - Float)'
      WRITE(6,*) '  ==>  Wind speed must be positive (in m/s)'
      GOTO 9998 

 2407 WRITE(6,*) '  SOS_PROC : ERROR_2407 on parameters'
      WRITE(6,*) '     K0_ROUJEAN, K1_ROUJEAN,K2_ROUJEAN in SOS_PROC'
      WRITE(6,*) '     (-SURF.Roujean.K0 keyword parameter - Float)'
      WRITE(6,*) '     (-SURF.Roujean.K1 keyword parameter - Float)'
      WRITE(6,*) '     (-SURF.Roujean.K2 keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'  
      WRITE(6,*) '  ==> Roujean BRDF parameters (required if type>=3)'
      GOTO 9998 

 2408 WRITE(6,*) '  SOS_PROC : ERROR_2408 on parameters'
      WRITE(6,*) '     ALPHA_NADAL and BETA_NADAL in SOS_PROC'
      WRITE(6,*) '     (-SURF.Nadal.Alpha keyword parameter - Float)'
      WRITE(6,*) '     (-SURF.Nadal.Beta keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Nadal BPDF parameters (required if type=6)'
      GOTO 9998 
 
 2409 WRITE(6,*) '  SOS_PROC : ERROR_2409 on parameter value'
      WRITE(6,*) '     ALPHA_NADAL parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Nadal.Alpha keyword parameter - Float)'
      WRITE(6,*) 
     &'  --> parameter value  must be included between 0 and 10' 
      GOTO 9998 

 2410 WRITE(6,*) '  SOS_PROC : ERROR_2410 on parameter value'
      WRITE(6,*) '     BETA_NADAL parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Nadal.Beta keyword parameter - Float)'
      WRITE(6,*) 
     &'  --> parameter value  must be included between 0 and 100' 
      GOTO 9998 

 2411 WRITE(6,*) '  SOS_PROC : ERROR_2411 on parameters'
      WRITE(6,*) '     COEF_C_MAIGNAN parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Maignan.C  keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) 
     &"  ==>  C*exp(-NDVI) coefficient of the Maignan's BPDF model"
      GOTO 9998 

C --- Controle des paramètres de profil

 2502 WRITE(6,*) '  SOS_PROC : ERROR_2502 on parameter value'
      WRITE(6,*) '     TR parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.MOT  keyword parameter - Float)'
      WRITE(6,*) '  --> - must be positive ' 
      WRITE(6,*) '  ==> Molecular Optical Thickness' 
      GOTO 9998 

 2503 WRITE(6,*) '  SOS_PROC : ERROR_2503 on parameters'
      WRITE(6,*) '     HR parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.HR  keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Molecular height scale' 
      WRITE(6,*) '      Case of non null AP.MOT value (TR in SOS_PROC)'       
      GOTO 9998    

 2504 WRITE(6,*) '  SOS_PROC : ERROR_2504 on parameter value'
      WRITE(6,*) '     HR parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.HR  keyword parameter - Float)'
      WRITE(6,*) '  --> must be positive ' 
      WRITE(6,*) '  ==> Molecular height scale' 
      WRITE(6,*) '      Case of non null AP.MOT value (TR in SOS_PROC)' 
      GOTO 9998    
 
 2505 WRITE(6,*) '  SOS_PROC : ERROR_2505 on parameters'
      WRITE(6,*) '     IPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AerProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Aerosol atmospheric profil type' 
      GOTO 9998    
 
 2506 WRITE(6,*) '  SOS_PROC : ERROR_2506 on parameter value'
      WRITE(6,*) '     IPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AerProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) '  --> Min value: 1'    
      WRITE(6,*) '  --> Max value: CTE_NB_PROFIL (constant)'          
      WRITE(6,*) '  CTE_NB_PROFIL: ',CTE_NB_PROFIL
      WRITE(6,*) '  -AP.AerProfile.Type: ',IPROFIL   
      WRITE(6,*) '  ==> Aerosol atmospheric profil type' 
      GOTO 9998    

 2507 WRITE(6,*) '  SOS_PROC : ERROR_2507 on parameters'
      WRITE(6,*) '     HA parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.AerHS.HA  keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Height scale of the aerosol profile'
      WRITE(6,*) 
     &'  ==> Case of non null -AER.AOTref value (AOT_REF in SOS_PROC)'
      GOTO 9998    

 2508 WRITE(6,*) '  SOS_PROC : ERROR_2508 on parameter value'
      WRITE(6,*) '     HA parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.AerHS.HA  keyword parameter - Float)'
      WRITE(6,*) '  --> must be positive '
      WRITE(6,*) '  ==> Height scale of the aerosol profile'
      WRITE(6,*) 
     &'  ==> Case of non null -AER.AOTref value (AOT_REF in SOS_PROC)'
      GOTO 9998    
 
 2509 WRITE(6,*) '  SOS_PROC : ERROR_2509 on parameters'
      WRITE(6,*) '     ZMIN and ZMAX parameters in SOS_PROC'
      WRITE(6,*) '     (-AP.AerLayer.Zmin keyword parameter - Float)'
      WRITE(6,*) '     (-AP.AerLayer.Zmax keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> altitudes of aerosol layer (in km)'
      WRITE(6,*) 
     &'  ==> Case of AP.AerProfile.Type = 2 (IPROFIL in SOS_PROC)'
      GOTO 9998     

 2510 WRITE(6,*) '  SOS_PROC : ERROR_2510 on parameters'
      WRITE(6,*) '     ABSPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AbsProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Absorption profil type'
      WRITE(6,*) '         0 : User profile '
      WRITE(6,*) '         1 : Tropical profile'
      WRITE(6,*) '         2 : MidLat Summer profile'
      WRITE(6,*) '         3 : Midlat Winter profile'
      WRITE(6,*) '         4 : SAS profile'
      WRITE(6,*) '         5 : SAW profile'
      WRITE(6,*) '         6 : USST62'
      WRITE(6,*) '         7 : No absorption'             
      GOTO 9998   
  
 2511 WRITE(6,*) '  SOS_PROC : ERROR_2511 on parameter value'
      WRITE(6,*) '     ABSPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AbsProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) '  ==> Absorption profil type'
      WRITE(6,*) '         0 : User profile '
      WRITE(6,*) '         1 : Tropical profile'
      WRITE(6,*) '         2 : MidLat Summer profile'
      WRITE(6,*) '         3 : Midlat Winter profile'
      WRITE(6,*) '         4 : SAS profile'
      WRITE(6,*) '         5 : SAW profile'
      WRITE(6,*) '         6 : USST62'
      WRITE(6,*) '         7 : No absorption'             
      GOTO 9998    
    
 2512 WRITE(6,*) '  SOS_PROC : ERROR_2512 on parameters'
      WRITE(6,*) '     FICABSPROFIL parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.AbsProfile.UserFile keyword parameter)'
      WRITE(6,*) '     --> must be defined'    
      WRITE(6,*) '  ==> User absorption profile filename'
      WRITE(6,*) 
     &'  ==> Case of AP.AbsProfile.Type = 0 (ABSPROFIL in SOS_PROC)'
      GOTO 9998     

 2513 WRITE(6,*) '  SOS_PROC : ERROR_2513 on parameters'
      WRITE(6,*) '     ABSPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AbsProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) 
     &'  ==> Absorption profil type must be set to 7 (no absorption)'
      WRITE(6,*) 
     &'  ==> Case of -AP.AerProfile.Type = 2 (IPROFIl in SOS_PROC)'
      WRITE(6,*) '  ==> Aerosol layer between two altitudes'
      GOTO 9998     

 2514 WRITE(6,*) '  SOS_PROC : ERROR_2514 on parameters'
      WRITE(6,*) '     NUSTEP parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.SpectralResol keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'   
      WRITE(6,*) '  ==> Spectral resolution of CKD coefficients (cm-1)'
      GOTO 9998  

25141 WRITE(6,*) '  SOS_PROC : ERROR_25141 on parameter value'
      WRITE(6,*) '     NUSTEP parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.SpectralResol keyword parameter - Float)'
      WRITE(6,*) '  ==> Spectral resolution of CKD coefficients'
      WRITE(6,*) '  --> Available values are 1, 5 or 10 cm-1'
      GOTO 9998

 2515 WRITE(6,*) '  SOS_PROC : ERROR_2515 on parameters'
      WRITE(6,*) '     IMODE_CKD_CALCUL parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.AbsModeCKD  keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'   
      WRITE(6,*) '      1 : Fine calculation of the impact of gaz'
      WRITE(6,*) '          absorption on radiance by the CKD method.'
      WRITE(6,*) '      2 : Simplified calculation of the impact of'
      WRITE(6,*) '           gaz absorption radiance.'
      GOTO 9998  
      
 2516 WRITE(6,*) '  SOS_PROC : ERROR_2516 on AOT extrapolated to the',
     &           ' simulation wavelength'
      WRITE(6,*) '  --> AOT = ',TA, ' < 0'
      WRITE(6,*) 
     &'  ==> Check the value for -AER.AOTref (AOT_REF in SOS_PROC)'
      WRITE(6,*) 
     &'  ==> Check the parameters defining the aerosol model'
      GOTO 9998           

C --- Controle des paramètres SOS
 2604 WRITE(6,*) '  SOS_PROC : ERROR_2604 on parameter value'
      WRITE(6,*) '     IGMAX parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.IGmax keyword parameter - Integer)'
      WRITE(6,*) '  --> must be higher than 1' 
      GOTO 9998 

 2605 WRITE(6,*) '  SOS_PROC : ERROR_2605 on parameters'
      WRITE(6,*) '     ITRPHI parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined' 
      WRITE(6,*) '  ==> Azimut mode  '
      GOTO 9998   

 2606 WRITE(6,*) '  SOS_PROC : ERROR_2606 on parameter value'
      WRITE(6,*) '     ITRPHI parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View keyword parameter - Integer)'
      WRITE(6,*) '  --> must be equal to 1 or 2' 
      WRITE(6,*) '     1: Azimut plan mode' 
      WRITE(6,*) '     2: Polar diagram azimut mode'       
      GOTO 9998 
 
 2607 WRITE(6,*) '  SOS_PROC : ERROR_2607 on parameters'
      WRITE(6,*) '     PHIOS parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View.Phi keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined' 
      WRITE(6,*) 
     &'  ==> Relative azimut required when SOS.View=1 ',
     &' (ITRPHI=1 in SOS_PROC)'  
      GOTO 9998     
 
 2608 WRITE(6,*) '  SOS_PROC : ERROR_2608 on parameters'
      WRITE(6,*) '     PAS_PHI parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View.Dphi keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined' 
      WRITE(6,*) '  ==> Azimut sampling step required when SOS.View=2'
      WRITE(6,*) '      (ITRPHI=2 in SOS_PROC)'    
      GOTO 9998   

 2609 WRITE(6,*) '  SOS_PROC : ERROR_2609 on parameter value'
      WRITE(6,*) '     PAS_PHI parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View.Dphi keyword parameter - Integer)'
      WRITE(6,*) '  --> must be a positive value'
      WRITE(6,*) '  ==> Case SOS.View = 2 (ITRPHI in SOS_PROC)'      
      WRITE(6,*) '      Polar diagram azimut mode'      
      GOTO 9998     
   
 2610 WRITE(6,*) '  SOS_PROC : ERROR_2610 on parameter value'
      WRITE(6,*) '     IPOLAR parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.Ipolar keyword parameter - Integer)'
      WRITE(6,*) '  -->  must be equal to 0 or 1 ' 
      GOTO 9998 

 2611 WRITE(6,*) '  SOS_PROC : ERROR_2611 on parameter value'
      WRITE(6,*) '     ZOUT parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.OutputAlt keyword parameter - Float)'
      WRITE(6,*) '  -->  must be equal to:' 
      WRITE(6,*) 
     &'       -1 : standard output (Up TOA, Down BOA)'
      WRITE(6,*) 
     &'       Zout : output level altitude (in km) (Up & Down at Zout)'
      WRITE(6,*) 
     &'          for 0 < Zout < CTE_TOA_ALT (constant): TOA altitude' 
      WRITE(6,*) 
     &'          with CTE_TOA_ALT =',CTE_TOA_ALT            
      GOTO 9998 
   
C --- Définition des fichiers de trace

 3001 WRITE(6,*) '  SOS_PROC : ERROR_3001 : '
      WRITE(6,*) '  => Angles Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	
   
 3002 WRITE(6,*) '  SOS_PROC : ERROR_3002 : '
      WRITE(6,*) '  => Mie Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3003 WRITE(6,*) '  SOS_PROC : ERROR_3003 : '
      WRITE(6,*) '  => Aerosol Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3004 WRITE(6,*) '  SOS_PROC : ERROR_3004 : '
      WRITE(6,*) '  => Profil Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	
   
 3005 WRITE(6,*) '  SOS_PROC : ERROR_3005 : '
      WRITE(6,*) '  => Surface Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3006 WRITE(6,*) '  SOS_PROC : ERROR_3006 : '
      WRITE(6,*) '  => SOS Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

C --- Définition des fichiers de sortie

 3007 WRITE(6,*) '  SOS_PROC : ERROR_3007 : '
      WRITE(6,*) '  => Radiance Angles filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3008 WRITE(6,*) '  SOS_PROC : ERROR_3008 : '
      WRITE(6,*) '  => Mie Angles filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3009 WRITE(6,*) '  SOS_PROC : ERROR_3009 : '
      WRITE(6,*) '  => Aerosol Phaze function filename: maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	          	

 3011 WRITE(6,*) '  SOS_PROC : ERROR_3011 : '
      WRITE(6,*) '  => Radiance binary filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

30111 WRITE(6,*) '  SOS_PROC : ERROR_30111 : '
      WRITE(6,*) '  => Cte CTE_REMOVE_PREVIOUS_BIN_FILE is set to 0:'
      WRITE(6,*) '    --> Remove this file' 
      WRITE(6,*) '    --> or change the CTE_RESULT directory' 
      WRITE(6,*) ''   
      GOTO 9998
        
 3012 WRITE(6,*) '  SOS_PROC : ERROR_3012 : '
      WRITE(6,*) '  => Transmittance filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998
      
30121 WRITE(6,*) '  SOS_PROC : ERROR_30121 : '
      WRITE(6,*) '  => Flux filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998   

C --- Ouverture des fichiers de trace

 3104 WRITE(6,*) '  SOS_PROC : ERROR_3104 : '
      WRITE(6,*) '  => Error opening Profile log file: ',FICPROFIL_LOG
      GOTO 9998

 3106 WRITE(6,*) '  SOS_PROC : ERROR_3106 : '
      WRITE(6,*) '  => Error opening SOS log file: ',FICSOS_LOG
      GOTO 9998          
    
C --- Lancement des routines
 3500 WRITE(6,*) '  SOS_PROC : ERROR_3500 : '
      WRITE(6,*) '  => The summation of AIK coefs is not 1 : ',sumAIK
      WRITE(6,*) '  => Check the CKD coefficients'
      GOTO 9998          
  
C --- Lancement des routines

 4000 WRITE(6,*) '  SOS_PROC : ERROR_4000 on subroutine SOS_ANGLES'
      GOTO 9998 

 4100 WRITE(6,*) '  SOS_PROC : ERROR_4100 on subroutine SOS_AEROSOLS'
      GOTO 9998 

 4200 WRITE(6,*) '  SOS_PROC : ERROR_4200 on subroutine SOS_SURFACE'
      GOTO 9998 

 4300 WRITE(6,*) '  SOS_PROC : ERROR_4300 on subroutine 
     &            SOS_PREPA_ABSPROFILE'
      GOTO 9998 

 4400 WRITE(6,*) '  SOS_PROC : ERROR_4400 on subroutine SOS_PREPA_OS'
      GOTO 9998 

 4500 WRITE(6,*) '  SOS_PROC : ERROR_4500 on subroutine SOS_ABSPROFILE'
      GOTO 9998 

 4600 WRITE(6,*) '  SOS_PROC : ERROR_4600 on subroutine SOS_PROFIL'
      GOTO 9998 

 4700 WRITE(6,*) '  SOS_PROC : ERROR_4700 on subroutine SOS'
      GOTO 9998 

 4800 WRITE(6,*) '  SOS_PROC : ERROR_4800 on subroutine SOS_AGGREGATE'
      GOTO 9998 

 4900 WRITE(6,*) 
     &  '  SOS_PROC : ERROR_4900 on subroutine SOS_TRPHI_OPTION'
      GOTO 9998 

C -- Ouverture / écriture des fichiers Transmissions et Flux
  997 WRITE(6,*) '  SOS_PROC : ERROR_997 : '
      WRITE(6,*) '  => ERROR on transmission file opening'
      GOTO 9998         

  998 WRITE(6,*) '  SOS_PROC : ERROR_998 : '
      WRITE(6,*) '  => ERROR on transmission file writing '
      GOTO 9998         

  987 WRITE(6,*) '  SOS_PROC : ERROR_987 : '
      WRITE(6,*) '  => ERROR on flux file opening'
      GOTO 9998         

  988 WRITE(6,*) '  SOS_PROC : ERROR_988 : '
      WRITE(6,*) '  => ERROR on flux file writing'
      GOTO 9998         

      CALL FLUSH(6)             !Force l'envoie du buffer à l'écran

       
 9998 IER=1
 
      
 9999 RETURN



C =========================================================
C --- Formats
C =========================================================            
C


C --- Formats pour le fichier Trace : entrées de SOS_AEROSOLS
  807 FORMAT(28h       FICANGLES_RES_LUM :  ,A)
  841 FORMAT(45h       Simulation wavelength (Wa)          : ,F12.5)
  842 FORMAT(45h       AOT                                 : ,F12.5) 
  843 FORMAT(45h       Reference wavelength (Waref)        : ,F12.5) 
  844 FORMAT(45h       Reference AOT                       : ,F12.5) 
  845 FORMAT(45h       Truncature option                   : ,I4) 
  846 FORMAT(45h       Type of aerosols model              : ,I4)
  847 FORMAT(45h       MMD : Refractive index for Wa       : ,
     &      (F12.3,2X,F12.5))
c$$$  848 FORMAT(45h       MMD : Refractive index for Waref    : ,
c$$$     &      (F12.3,2X,F12.5))  
  849 FORMAT(45h       Type of size distribution           : ,I4)
  850 FORMAT(45h       Size distribution parameters        : ,
     &      (F12.5,2X,F12.5,2X,F12.5))    
  851 FORMAT(45h       Type of WMO model                   : ,I4)  
  852 FORMAT(45h       WMO component conc. (DL, WS, OC, SO): ,
     &      4(F12.5,2X))    
  853 FORMAT(45h       Type of Shettle & Fenn model        : ,I4)    
  854 FORMAT(45h       S&F model : relative humidity       : ,F12.5)
  855 FORMAT(45h       BMD : Choice of mixture description : ,I4) 
  856 FORMAT(45h       BMD : user conc. Coarse / Fine      : ,
     &      (F12.5,2X,F12.5))
  857 FORMAT(45h       BMD : AOT coarse / AOTtot for Waref : ,F12.5)   
  858 FORMAT(45h         Refractive index for Wa           : ,
     &      (F12.3,2X,F12.5))      
  859 FORMAT(45h         Refractive index for Waref        : ,
     &      (F12.3,2X,F12.5))      
  860 FORMAT(45h         Size distr parameters             : ,
     &      (F12.5,2X,F12.5))  
  861 FORMAT(45h       User's file defining aerosols IOP   : ,A)  
 8611 FORMAT(45h       User's file for aerosols mixture    : ,A)  
  862 FORMAT(45h       Log file                            : ,A)  
  863 FORMAT(45h       Res file                            : ,A) 
  864 FORMAT(45h       Storage directory for MIE files     : ,A)  
  865 FORMAT(45h       Mie computations log file           : ,A)    

C --- Formats pour le fichier Transmissions et Flux : 
1005  FORMAT(20hSolar Zenith Angle :,1x,f7.3)
1006  FORMAT(36hDirect transmission TOA -> surface :,1x,f8.4)
1010  FORMAT(12h    thetas =,1x,f6.3,3x,12htd(thetas) =,1x,f7.4)  
2010  FORMAT(12h    thetav =,1x,f6.3,3x,12htd(thetav) =,1x,f7.4)
2016  FORMAT(34h   - Downward direct flux at BOA :,1x,f9.5)
2017  FORMAT(34h   - Downward diffuse flux at BOA:,1x,f9.5)
2018  FORMAT(34h   ==> Downward total flux at BOA:,1x,f9.5)
2020  FORMAT(f7.2,2x,f7.4,1x,f7.4,1x,f7.4,1x,f7.4) 

     
      END     !FIN DU PROGRAMME SOS_PROC
