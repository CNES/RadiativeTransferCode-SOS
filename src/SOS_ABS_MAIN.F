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
C* FICHIER: SOS_ABS_MAIN.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Script de lancement du code des OS
C*           - gère les entrées / sorties et options de traitements,
C*           - exécute successivement les différents modules des OS
C*           - les luminances sont calculées pour chaque produit de la
C*             décomposition de la méthode K-Correlated
C*
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmosphérique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, B.Lafrance, X.Lenot
C*
C*
C* MOD:VERSION:1.0: DATE: 06/06/2013
C*    Création de la fonction (remplace ancien script main_SOS.ksh)
C*
C* MOD:VERSION:1.1: DATE: 22/02/2016
C*   - Ajout d'un test sur l'existence du fichier FICSOS_RES_BIN (-SOS.ResBin)
C*     pour éviter en cas de relance d'un traitement de poursuivre la complétude 
C*     de ce fichier.
C*
C*    - Suppression des formats à la lecture des paramètres d'entrée en faveur 
C*      de READ(PROG_ARG,*,ERR=NNN)
C*      ==> Permet une compilation sous gfortran et évite les mauvaises conversions 
C*          si la valeur d'entrée ne respecte pas strictement le format attendu 
C*          (par exemple, saisie d'un entier N au lieu d'un réel N. 
C*           aboutissant à une valeur interprétée fausse). 
C*
C*    - Ajustements mineurs pour le respect strict des 72 colonnes.
C*    - Corrections mineures de syntaxe de messages d'erreur 
C*      (AP.HA corrigé en AP.AerHS.HA)
C*    - Pas de fichier LOG aérosols par défaut.
C*
C*    - Tableau NEXP passe de Double précision à Integer
C*
C*    - Modification d'un controle du TA (AOT) après extrapolation spectrale
C*
C*    - Gestion du cas où la pression n'est pas définie en entrée (affectation 
C*      de la valeur par défaut)
C*      
C*    - Correction du passage des indices de réfraction à la routine SOS_AEROSOLS
C*      pour la longueur d'onde de référence (premier appel de SOS_AEROSOLS).
C*
C*    - Correction de la gestion du passage de l'indice d'option de troncature 
C*      à SOS_AEROSOLS pour éviter de basculer à zéro l'option si la troncature est
C*      annulée par le code.
C*
C*    - Ajustement pour correctement localiser les fichiers FICUP_ANGLES_USER 
C*      et FICDOWN_ANGLES_USER dans l'arborescence des fichiers résultats.
C*
C*    - Modification de la gestion des variances pour les granulométrie LND:
C*      abandon de l'usage d'un SIG = log10(variance) au profit de l'utilisation
C*      d'une formulation classique du modèle LND en fonction de la variance.
C*		N(r) = exp(-0.5 *(log10(r/rm)/sig)**2) / r
C*	devient 
C*		N(r) = exp(-0.5 *(log(r/rm)/sig)**2) / r
C*      ==> Impacte les valeurs des paramètres en entrée de SOS_AEROSOLS et donc
C*          fournis en ligne de commandes à SOS_ABS_MAIN.   
C*
C*
C* MOD:VERSION:1.2: DATE: 10/12/2018
C*   - Correction pour ne plus avoir de calcul systématique du fichier des 
C*     transmission même quand l'option n'est pas sollicitée (génération d'un
C*     fichier NO_OUTPUT non souhaité).
C*
C*   - Suppression de la lecture du paramère AP.ResFile définissant le nom du
C*     fichier de profil généré car le profil est propre à chaque terme de la 
C*     décomposition CKD. Le stocker induit en erreur sur son contenu.
C*     Utilisation d'un fichier de profil temporaire 
C*     (nom par défaut, fichier effacé en fin de traitement).
C*
C*   - Par défaut le fichier FICABSPROFIL est nommé "NO_USER_ABS_PROFILE_FILE"
C*
C*
C* MOD:VERSION:1.3: DATE: 01/08/2019
C*   - Correction pour ne plus avoir de conflit, entre le programme principal 
C*     SOS_ABS_MAIN et la routine SOS_TRPHI_OPTION, concernant la déclaration 
C*     du nom d'un fichier de sortie de luminance pour les angles utilisateurs 
C*     quand l'option n'est pas sollicitée (i.e quand il n'y a pas de définition
C*     de -SOS.ResFileUp.UserAng ou -SOS.ResFileDown.UserAng mais que 
C*     -ANG.Rad.UserAngFile est défini).
C*     ==> L'initialisation du nom des fichiers FICUP_ANGLES_USER 
C*         et FICDOWN_ANGLES_USER passe de "NO_USER_ANGLES" à "NO_OUTPUT".
C*     ==> Ajustement du test pour réaliser la complétude du nom de ces fichiers.
C*
C*
C* MOD:VERSION:2.0: DATE: 20/04/2021
C*   - Ajustement pour usage de nouvelles données CKD
C*     ==> Suppression de la distinction des absorptions par raies fines ou continuums.
C*     ==> Suppression des tables AA et BB 
C*     ==> Ajustement des tables XKI et AL
C*     ==> Ajout des constantes de définition des dimensions de tables pour coefs CKD
C*     ==> Ajout du nouveau paramètre d'entrée NUSTEP donnant le pas spectral souhaité
C*         des fichiers de coefficient CKD à lire.
C*     ==> Abandon des continuums SO2 (absorption au-delà de 3.7 microns)
C*         et N2 (absorption pour longueurs d'onde inférieures à 0.33 microns)
C*     ==> Ajustements pour traiter le continuum NO2 à partir de coefficients CKD
C*     ==> Suppression du test de non simulation si la transmission sur la colonne
C*         atmosphérique est trop faible, de sorte à préserver les diffusions dans
C*         les hautes couches atmosphériques.
C*   - Le mot clé AP.H2o devient AP.H2O
C*   - Introduction de la possibilité d'utiliser le modèle de BPDF de Maignan.
C*   - Désactivation de l'option d'utilisation du modèle de BPDF de Nadal.
C*   - Ajout d'un test pour avertir si un paramètre lu ne fait pas parti des 
C*     paramètres possibles.
C*   - Ajout d'un répertoire de dépôt des fichiers temporaires de traitement (DIR_TMP)
C*     pour permettre plusieurs lancements simultanés depuis le même répertoire
C*     ==> fichier temporaire FICPROFIL_TMP : fichier de profil pour un traitement unitaire de la routine SOS.
C*     ==> fichier temporaire FICOS_TMP : résultat d'un traitement unitaire de la routine SOS.
C*     ==> fichier temporaire FICOS_AGG_TMP : aggrégation des résultats des traitements unitaires. 
C*   - Clarification de l'utilisation des modèles mono-modaux
C*     ==> Suppression de la possibilité de définir le fichier de Mie à utiliser (le code les définit par lui-même).
C*         La paramètre -AER.MMD.Mie.Filename n'est plus proposé.
C*     ==> Suppression des paramètres -SDparam1 et -SDparam2 des modèles mono-modaux,
C*         en faveur de nouveaux paramètres plus explicites:
C*         -AER.MMD.LNDradius :  Modal radius of the LND.
C*         -AER.MMD.LNDvar : Variance of the LND for a mono-modal distribution.
C*         -AER.MMD.JD.slope : Slope of the Junge's law.
C*         -AER.MMD.JD.rmin  : Minimal radius of the Junge's law.
C*         -AER.MMD.JD.rmax  : Maximal radius of the Junge's law (optional)
C*     ==> Ajout d'une valeur par défaut du rayon maximal d'un modèle de Junge
C*         (SOS_DEFAULT_AER_JUNGE_RMAX)
C*     ==> Adaptation de l'appel de la routine SOS_AEROSOLS pour passer l'argument associé à -AER.MMD.JD.rmax
C*     ==> Suppression du paramètre -AER.MMD.Mie.AlphaMax 
C*         et du passage d'argument ALPHAF à la routine SOS_AEROSOLS.
C*         La valeur maximale du paramètre de taille pour les calculs de Mie est définie par la routine SOS_AEROSOLS.
C*    - Ajout d'une option de traitement IMOD=5 pour utiliser la définition 
C*      d'un mélange de modes fournis dans un fichier par l'utilisateur.
C*      ==> Introduction du paramètre: FICMIXTURE_AER qui donne la location du fichier.
C*      ==> Adaptation de l'appel de la routine SOS_AEROSOLS.
C*    - Ajout de paramètres pour définir la concentration en surface du dioxyde de carbone (AP.CO2) 
C*      et du méthane (AP.CH4).
C*
C* MOD:VERSION:2.1: 31/08/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Suppression des arguments TRACEOS et IDLOGOS à l'appel de SOS_AGGREGATE
C*    - Conversion des traces en Anglais.
C*
C* MOD:VERSION:2.2: 24/11/2022
C*    - Refonte du programme SOS_ABS_MAIN pour introduire la routine SOS_PROC permettant un appel du code par Python.
C*
C* MOD:VERSION:2.3: 05/04/2023
C*    - Usage de la constante CTE_TOA_ALT pour définir l'altitude max (TOA)
C*    
C*  MOD:VERSION:2.4: 19/04/2023
C*    - Ajout de la variable IND_ANGOUT_FIN aux paramètres de sortie de SOS_PROC    
C*
C*  MOD:VERSION:2.5: 10/08/2023
C*    - Ajout du coeffficient de troncature de la fonction de phase aérosols en paramètre de sortie 
C*      de la routine SOS_PROC
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Définition des constantes   
C---------------------------------------------------------------------------- 
C Constantes utilisées :
C    CTE_LENDIR : Longueur des noms de répertoires.
C    CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_LENKEYWORD : Taille maximale des Keywords pour le passage des arguments.
C    CTE_HT_STD_PSURF: Pression atmosphèrique standard
C    CTE_NOT_DEFINED_VALUE_INT: valeur par défaut pour les entiers
C    CTE_NOT_DEFINED_VALUE_DBLE: valeur par défaut pour les réels
C    CTE_DEFAULT_FICANGLES_RES_LUM: nom par défaut du fichier angles pour le calcul des luminances
C    CTE_DEFAULT_FICANGLES_RES_MIE: nom par défaut du fichier angles pour le calcul des aérosols
C    CTE_DEFAULT_FICGRANU: nom par défaut du fichier granulométrie 
C    CTE_DEFAULT_RESBIN: nom par défaut du fichier binaire de résultats
C    CTE_DEFAULT_RESUP: nom par défaut du fichier de luminances montantes
C    CTE_DEFAULT_RESDOWN: nom par défaut du fichier de luminances descendantes
C    CTE_OS_NBMU_MAX  : Nombre maximal d'angles de Gauss positifs pour les calculs de luminance
C    CTE_DEFAULT_AER_JUNGE_RMAX : Rayon maximal (en microns) d'une granulomètrie du modèle de Junge
C    CTE_DEFAULT_IGMAX : Nombre maximal d'ordres de diffusion par défaut

C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define INCTE_IDLOG_MAIN 77

C==============================================================================
C PROGRAMME: SOS_ABS_MAIN
C ==========  
C
C Programme principal du code des OS avec absorption 
C
C
C Donnees en entrée du programme
C ------------------------------
C     -->Les paramètres/arguments du programme
C
C        Les arguments sont passés sous la forme d'une liste de couples (-Keyword Value)
C
C        Définition des mots clef, statut obligatoire / option, format de la valeur associée:
C
C 	 ------- Arguments propres au programme principal	
C 	 ------------------------------------------------
C
C         -SOS_Main.ResRoot ==> Obligatoire 
C             RESROOT   (CHARACTER*CTE_LENDIR) : répertoire de travail du traitement 
C						       (chemin complet)
C 	
C         -SOS_Main.Log  ==> Optionnel :  
C            FICMAIN_LOG (CHARACTER*CTE_LENFIC2) : Nom du fichier Trace du programme principal.
C                                          	  (nom fichier + extension) 
C						         ==> Pas de  trace si non défini.
C
C         -SOS_Main.Wa ==> Obligatoire 
C	     PROG_ARG  (CHARACTER*10)  : Longueur d'onde de simulation de la luminance 
C                                       en microns.	
C                                       (format de décodage : libre --> Float)
C
C
C
C 	 ------- Arguments propres à la définition des angles		
C 	 ----------------------------------------------------	
C
C    Keyword  : -ANG.Rad.NbGauss  
C    ---------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_NBMU_LUM is used : see SOS.h )       
C       Value :  Number of Gauss angles to be used for radiance computations 
C       Value format : No format constraint --> read as Integer
C
C    Keyword  : -ANG.Rad.UserAngFile  
C    -------------------------------   
C       Status : Optional      
C       Value :  Filename of the complementary list of user's angles to complete the 
C                ANG.Rad.NbGauss angles (complete path). 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -ANG.Thetas  
C    --------------------------   
C       Status : Required     
C       Solar zenithal angle (degrees) 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -ANG.Rad.ResFile  
C    ---------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_FICANGLES_RES_LUM is used : see SOS.h )      
C       Value :  Filename of list of angles used to BRDF/BPDF and radiance computations. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -ANG.Aer.NbGauss  
C    ---------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_NBMU_MIE is used : see SOS.h )       
C       Value :  Number of Gauss angles to be used for phaze matrix computations 
C       Value format : No format constraint --> read as Integer
C
C    Keyword  : -ANG.Aer.UserAngFile 
C    -------------------------------   
C       Status : Optional      
C       Value :  Filename of the complementary list of user's angles to complete the 
C                ANG.Aer.NbGauss angles (complete path). 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -ANG.Aer.ResFile  
C    ----------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_FICANGLES_RES_MIE is used : see SOS.h )       
C       Value :  Filename of list of angles used to the matrix phase function computations. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -ANG.Log
C    -------------------   
C       Status : Optional
C       Value :  log filename for ANGLES computations. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C
C 	 ------- Arguments propres à la définition du modèle de surface		
C 	 ---------------------------------------------------------------	
C
C    Keyword  : -SURF.Dir
C    -------------------   
C       Status : Conditional (required if a surface BRDF/BPDF is simulated)
C       Value :  Storage directory for SURFACE files producted by SOS_SURFACE computations (complete path). 
C       Value format : Encoding as string (CTE_LENDIR characters max : see SOS.h).
C
C    Keyword  : -SURF.Type  
C    ---------------------   
C       Status : Required      
C       Value :  Kind of surface 
C      		0 : lambertian surface for an albedo RHO.
C			1 : lambertian surface + glitter (roughness water).
C			2 : lambertian surface + Fresnel's reflexion (plane water).
C		       3 : lambertian surface + Roujean's BRDF.
C			4 : lambertian surface + Roujean's BRDF + Rondeaux's BPDF.
C			5 : lambertian surface + Roujean's BRDF + Breon's BPDF .
C			6 : lambertian surface + Roujean's BRDF + Nadal's BPDF.
C       Value format : No format constraint --> read as Integer.
C
C    Keyword  : -SURF.Log  
C    --------------------   
C       Status : Optional    
C       Value :  log filename for SOS_SURFACE computations.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -SURF.Ind  
C    --------------------   
C       Status : Associated to  -SURF.Type 1 2 4 5 6 or 7     
C       Value :  Surface / atmosphere refractive index. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Glitter.Wind  
C    -----------------------------   
C       Status : Associated to  -SURF.Type 1    
C       Value :  Wind velocity (m/s). 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Roujean.K0  
C    ---------------------------   
C       Status : Associated to  -SURF.Type 3 4 5 6 or 7      
C       Value :  Parameter K0 of Roujean's BRDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Roujean.K1  
C    --------------------------- 
C       Status : Associated to  -SURF.Type 3 4 5 6 or 7      
C       Value :  Parameter K1 of Roujean's BRDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Roujean.K2  
C    ---------------------------   
C       Status : Associated to  -SURF.Type 3 4 5 6 or 7      
C       Value :  Parameter K2 of Roujean's BRDF model. 
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -SURF.Nadal.Alpha 
C    ----------------------------   
C       Status : Associated to  -SURF.Type 6     
C       Value :  Parameter Alpha of Nadal's BPDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Nadal.Beta  
C    ---------------------------   
C       Status : Associated to  -SURF.Type 6      
C       Value :  Parameter Beta of Nadal's BPDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Maignan.C   
C    ---------------------------   
C       Status : Associated to  -SURF.Type 7      
C       Value :  C*exp(-NDVI) parameter of the Maignan's surface BPDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Alb 
C    --------------------   
C       Status : Required      
C       Value :  Surface albedo for the wavelength of radiance simulation
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.File  
C    ---------------------   
C       Status : Optional. Associated to  -SURF.Type 1 3 4 5 6 or 7     
C       Value :  User filename for the surface reflexion matrix (complete path).
C                --> Glitter or others BRDF/BPDF. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C 	 ------- Arguments propres à la définition des modèles d'aerosols		
C 	 ----------------------------------------------------------------	
C
C    Keyword  : -AER.Waref  
C    ---------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Reference wavelength for the aerosol optical thickness (microns).
C       Value format : No format constraint --> read as Float
C	  
C    Keyword  : -AER.AOTref  
C    ----------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Aerosol optical thickness for the reference wavelength.
C       Value format : No format constraint --> read as Float
C     		       
C    Keyword  : -AER.Tronca 
C    ----------------------   
C       Status : Default value = 1 (for common case with AER.AOT > 0)     
C       Value :  Option for aerosol phaze function troncation 
C                (1 by default to apply a truncation, 0 to skip the truncation)
C       Value format : No format constraint --> read as Integer.
C
C    Keyword  : -AER.Log
C    -------------------   
C       Status : Optional (for common case with AER.AOT > 0)     
C       Value :  log filename for SOS_AEROSOLS computations (complete path).
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -AER.MieLog  
C    ----------------------   
C       Status : Optional  
C       Value :  log filename for MIE computations (complete path).
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -AER.DirMie
C    ----------------------   
C       Status : Required      
C       Value :  Storage directory for MIE files producted by SOS_AEROSOLS computations (complete path). 
C       Value format : Encoding as string (CTE_LENDIR characters max : see SOS.h).
C
C    Keyword  : -AER.ResFile
C    -----------------------   
C       Status : Optional (if not defined, neither AER.UserFile, 
C                          the default value CTE_DEFAULT_FICGRANU is used : see SOS.h )   
C       Value :  Filename of the result SOS_AEROSOLS computations. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -AER.UserFile
C    -----------------------   
C       Status : Optional (required if AER.ResFile is not defined)      
C       Value :  Filename of the user's aerosol filename (complete path). 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -AER.Model 
C    ----------------------   
C       Status : Required (for common case with AER.AOT > 0)     
C       Value :  Type of size distribution (mono or multi-mode)
C      		        0 : Mono-modal models.
C			        1 : WMO models.
C			        2 : Shettle & Fenn models.
C		               3 : LND bi-modal models.
C                            4 : Use of external phase functions
C       Value format : No format constraint --> read as Integer.
C
C					
C    Keyword  : -AER.MMD.MRwa  
C    ------------------------ 
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
C       Value :  Real part of the refractive index at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float
C	  
C    Keyword  : -AER.MMD.MIwa  
C    ------------------------ 
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
C       Value :  Imaginary part of the refractive index for the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -AER.MMD.MRwaref  
C    ------------------------ 
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
C       Value :  Real part of the refractive index at the reference wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float
C	  
C    Keyword  : -AER.MMD.MIwaref  
C    ------------------------ 
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
C       Value :  Imaginary part of the refractive index for the reference wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float
C	  
C    Keyword  : -AER.MMD.SDtype 
C    --------------------------   
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)  
C       Value :  Type of mono-modal size distribution 
C			1 : LND.
C			2 : Junge's law.
C       Value format : No format constraint --> read as Integer.			
C	    
C    Keyword  : -AER.MMD.LNDradius  
C    ----------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 1 (for common case with AER.AOT > 0)     
C       Value :  Modal radius of the LND (microns).
C			-AER.MMD.SDtype 2 : minimal radius of the Junge's law (microns).
C       Value format : No format constraint --> read as Float	  
C	    
C    Keyword  : -AER.MMD.LNDvar  
C    -------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 1 (for common case with AER.AOT > 0)     
C       Value :  Variance of the LND (microns) for a mono-modal distribution.
C       Value format : No format constraint --> read as Float	  
C	  
C    Keyword  : -AER.MMD.JD.slope 
C    ---------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 2 (for common case with AER.AOT > 0)     
C       Value :  Slope of the Junge's law.
C       Value format : No format constraint --> read as Float	  
C	  
C    Keyword  : -AER.MMD.JD.rmin 
C    ---------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 2 (for common case with AER.AOT > 0)     
C       Value :  Minimal radius of the Junge's law (microns).
C       Value format : No format constraint --> read as Float	  
C	  
C    Keyword  : -AER.MMD.JD.rmax 
C    ---------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 2 (for common case with AER.AOT > 0)     
C       Value :  Maximal radius of the Junge's law (microns). 
C                If not defined : the default value CTE_DEFAULT_AER_JUNGE_RMAX is used
C       Value format : No format constraint --> read as Float	  
C  
C    Keyword  : -AER.WMO.Model  
C    ------------------------- 
C       Status : Associated to -AER.Model 1 (for common case with AER.AOT > 0)     
C       Value :  Type of WMO model.
C      		        1 : Continental WMO model.
C			        2 : Maritime WMO model.
C			        3 : Urban WMO model.
C		               4 : WMO model by user définition.
C       Value format : No format constraint --> read as Integer.	  
C
C    Keyword  : -AER.WMO.DL  
C    ---------------------- 
C       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
C       Value :  Volumetric concentration (%) for WMO Dust-Like particles.
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -AER.WMO.WS  
C    ---------------------- 
C       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
C       Value :  Volumetric concentration (%) for WMO Water Soluble particles.
C       Value format : Encoding as real F9.5 (10 characters max).
C
C    Keyword  : -AER.WMO.OC  
C    ---------------------- 
C       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
C       Value :  Volumetric concentration (%) for WMO OCeanic particles.
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AER.WMO.SO  
C    ---------------------- 
C       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
C       Value :  Volumetric concentration (%) for WMO SOot particles.
C       Value format : No format constraint --> read as Float.	  
C	  
C    Keyword  : -AER.SF.Model  
C    ------------------------ 
C       Status : Associated to -AER.Model 2 (for common case with AER.AOT > 0)     
C       Value :  Type of Shettle & Fenn model.
C      		1 : Tropospheric S&F model.
C			2 : Urban S&F model.
C			3 : Maritime S&F model.
C		       4 : Coastal S&F model.
C       Value format : No format constraint --> read as Integer.	  
C
C    Keyword  : -AER.SF.HR  
C    --------------------- 
C       Status : Associated to -AER.Model 2 (for common case with AER.AOT > 0)     
C       Value :  Relative humidity (%) for Shettle & Fenn model.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.VCdef 
C    ---------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)  
C       Value :  Type of selected method for bi-modal LND volumetric concentration définition 
C			1 : Use of user volumetric concentrations.
C			2 : Use of user ratio AOT_coarse / AOT_total.
C       Value format : No format constraint --> read as Integer.			
C
C    Keyword  : -AER.BMD.CoarseVC  
C    ---------------------------- 
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 1 (for common case with AER.AOT > 0)     
C       Value :  User volumetric concentration for the "LND coarse mode".
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AER.BMD.FineVC  
C    -------------------------- 
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 1 (for common case with AER.AOT > 0)     
C       Value :  User volumetric concentration for the "LND fine mode".
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AER.BMD.RAOT  
C    ------------------------ 
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)     
C       Value :  User value of the ration AOT_coarse / AOT_total for the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  	  
C    Keyword  : -AER.BMD.CM.MRwa  
C    ---------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Real part of the refractive index for the "LND coarse mode" 
C                at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.MIwa  
C    ---------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Imaginary part of the refractive index for the "LND coarse mode" 
C                at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.MRwaref  
C    ------------------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Real part of the refractive index for the "LND coarse mode" 
C                at the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.MIwaref  
C    -----------------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Imaginary part of the refractive index for the "LND coarse mode" 
C                at the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.SDradius  
C    -------------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Modal radius of the "LND coarse mode" (microns).
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.SDvar  
C    ----------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Log10 of the "LND coarse mode" variance.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.MRwa  
C    ---------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Real part of the refractive index for the "LND fine mode" 
C                at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.MIwa  
C    ---------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Imaginary part of the refractive index for the "LND fine mode" 
C                at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.MRwaref  
C    ------------------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Real part of the refractive index for the "LND fine mode" 
C                at the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.MIwaref  
C    -----------------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Imaginary part of the refractive index for the "LND fine mode" 
C                at the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.SDradius  
C    -------------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Modal radius of the "LND fine mode" (microns).
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.SDvar  
C    ----------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Log10 of the "LND fine mode" variance.
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AER.ExtData
C    -----------------------   
C       Status : Associated to -AER.Model 4 (for common case with AER.AOT > 0)       
C       Value :  Filename (complete path) of user's external phaze functions 
C                and radiative parameters (extinction and scattering coefficients)
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C       Special requirement : the reference aerosol wavelength (-AER.Waref) and 
C                             the radiance simulation wavelength (-SOS.Wa) have to be equal.	 
C
C    Keyword  : -AER.DefMixture
C    --------------------------   
C       Status : Associated to -AER.Model 5 (for common case with AER.AOT > 0)       
C       Value :  Filename (complete path) of user's file defining the mixture of aerosol models.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C    
C
C 	 ------- Arguments propres à la définition des profils		
C 	 ------------------------------------------------------------------	
C
C    Keyword  : -AP.Log
C    ------------------   
C       Status : Optional
C       Value :  log filename for SOS_PROFIL computations.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -AP.MOT  
C    ------------------   
C       Status : Optional      
C       Value :  Molecular optical thickness for the wavelength of radiance simulation
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.HR  
C    ------------------   
C       Status : Associated to -AP.MOT >= 0.    
C       Value :  Molecular optical thickness for the wavelength of radiance simulation
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.AerProfile.Type  
C    -------------------   
C       Status : Required      
C       Value :  Type of aerosol profile
C      		     1 : Aerosol profile described by an exponential decrease using a heigth scale
C		     2 : Homogeneous aerosol layer between two altitudes.
C       Value format : No format constraint --> read as Integer.
C
C    Keyword  : -AP.AerHS.HA  
C    -----------------------   
C       Status : Associated to -AP.Type 1 and -AP.AOT >= 0.   
C       Value :  Aerosol heigth scale (km).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.AerLayer.Zmin  
C    ----------------------------   
C       Status : Associated to -AP.Type 2      
C       Value :  Minimal altitude of the aerosol layer (km).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.AerLayer.Zmax  
C    ----------------------------   
C       Status : Associated to -AP.Type 2      
C       Value :  Maximal altitude of the aerosol layer (km).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.Psurf  
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)      
C       Value :  Surface Pressure (mbar).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.H2O  
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)	 
C       Value :  Water vapor content (g/cm2).
C       Value format : No format constraint --> read as Float.
C  
C    Keyword  : -AP.O3  
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)	 
C       Value :  Ozone content (Dobson).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.CO2 
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)	 
C       Value :  Carbon dioxyde concentration at surface level (ppmv).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.CH4  
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)	 
C       Value :  Methane concentration at surface level (ppmv).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.AbsProfile.Type  
C    ----------------------------   
C       Status : Required	 
C       Value :  Type of asborption atmospheric profile
C			 0:User
C			 1:TROPICAL
C			 2:MLS
C			 3:MLW
C			 4:SAS
C			 5:SAW
C			 6:USST62
C			 7:No absorption
C       Value format : No format constraint --> read as Integer.
C 
C
C    Keyword  : -AP.SpectralResol  
C    ----------------------------   
C       Status : Conditional (required in case of simulation of absorptions, i.e. if AP.AbsProfile.Type != 7)	 
C       Value :  Spectral resolution in cm-1 of CKD datafiles
C       Value format : No format constraint --> read as Float.
C 

C    Keyword  : -AP.AbsProfile.UserFile  
C    ----------------------------   
C       Status : Required in association to  AP.AbsProfile.Type = 0	 
C       Value :  Filename of user's absorption profile file (complete path)
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h)	 
C
C
C 	 ------- Arguments propres aux calculs d'Ordres Successifs 
C 	 ------- et à la production des fichiers résultats
C 	 ---------------------------------------------------------	
C
C    Keyword  : -SOS.ResBin  
C    ----------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_RESBIN is used : see SOS.h)        
C       Value :  Filename of the binary file resulting from SOS computations . 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -SOS.ResFileUp  
C    -------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_RESUP is used : see SOS.h)     
C       Value :  Filename of the ascii file resulting from SOS computations  
C                --> upward radiance.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -SOS.ResFileDown  
C    ---------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_RESDOWN is used : see SOS.h)      
C       Value :  Filename of the ascii file resulting from SOS computations  
C                --> downward radiance.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -SOS.ResFileUp.UserAng 
C    ---------------------------------  
C       Status : Optional  (requires the use of -SOS.UserAng.File)       
C       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
C                giving radiance values only for user's angles 
C                --> upward radiance.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted     
C
C    Keyword  : -SOS.ResFileDown.UserAng 
C    -----------------------------------  
C       Status : Optional  (requires the use of -SOS.UserAng.File)       
C       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
C                giving radiance values only for user's angles 
C                --> downward radiance.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted     
C
C    Keyword  : -SOS.Log
C    -------------------   
C       Status : Optional
C       Value :  log filename for SOS computations.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -SOS.AbsModeCKD  
C    --------------------------   
C       Status : Required      
C       Value :  Number defining the CKD mode calculation
C                1 : Computation of the radiance from CKD expansion
C                2 : Computation of the absorption optical thickness from CKD expansion 
C       Value format : No format constraint --> read as Integer
C
C    Keyword  : -SOS.Trans
C    ---------------------  
C       Status : Optional 
C       Value :  Filename of the file given information about the simulated 
C                atmospheric transmissions (will be located in the repertory Working_dir/SOS).
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -SOS.Flux
C    ---------------------  
C       Status : Optional 
C       Value :  Filename of the file given information about the fluxes at BOA and TOA
C                (will be located in the repertory Working_dir/SOS).
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -SOS.OutputAlt  
C    ---------------------------   
C       Status : Optionnal   
C       Value :  Altitude (in km) for radiance results.
C                  If not defined : -1 is the default value attributed 
C                  -1 : standard results 
C			--> Upward radiance at TOA
C                       --> Downward radiance at surface level
C                  z  : results for the altitude z of the atmospheric profile
C			(0 <= z <= CTE_TOA_ALT)
C			--> Upward and downward radiance are given for the
C			    altitude z of the profile.
C       Value format : No format constraint --> read as Float.
C			
C    Keyword  : -SOS.Ipolar 
C    ----------------------   
C       Status : Optionnal   
C       Value :  Profile level for radiance results.
C                  If not defined : 1 is the default value attributed
C                  0 : radiance calculation without taking account of polarization. 
C                  1 : radiance calculation with polarization processes.
C       Value format : No format constraint --> read as Integer.
C
C    Keyword  : -SOS.IGmax  
C    ---------------------   
C       Status : Optionnal      
C       Value :  Maximal order of interaction (scattering & surface reflexion). 
C       Value format : No format constraint --> read as Integer.   	 	 	
C
C    Keyword  : -SOS.View  
C    --------------------   
C       Status : Required      
C       Value :  Index for output type
C                  1 : radiance for viewing directions in a constante azimutal plane. 
C                  2 : radiance for viewing directions covering half the space, with a constant step on azimut.  
C       Value format : No format constraint --> read as Integer.   	 	 	
C		 	 	 
C    Keyword  : -SOS.View.Phi  
C    ------------------------   
C       Status : Associated to -SOS.View 1      
C       Value :  Relative azimut angle (degrees).  
C       Value format : No format constraint --> read as Float.  	 	 	
C
C    Keyword  : -SOS.View.Dphi  
C    -------------------------   
C       Status : Associated to -SOS.View 2      
C       Value :  Constant step on azimut (degrees). Integer value.  
C       Value format : No format constraint --> read as Integer.   
C  					
C       
C Résultats fournis par le programme
C ----------------------------------
C   Sorties de SOS_PROC :
C      - Champs de luminance normalisée par défaut émergeant au TOA et descendant en surface
C        ou montant et descendant pour le niveau du profil sélectionné.
C      - Données de polarisation :
C        Paramètres Q et U du vecteur de Stokes, angle de polarisation, intensité polarisée, taux de polarisation.
C      - Transmissions directes et diffuses, descendantes et montantes.
C      - Flux direct et diffus descendant au niveau de la surface, flux diffus montant au TOA.
C      - Liste d'indice des angles utilisateurs pour filtrage des données de luminance.
C      - Coefficient de troncature de la fonction de phase des aérosols.
C
C   -->Contenu du fichier résultat pour le champ montant
C       Fichier Ascii, format d'un enregistrement :
C          2(2X,F7.2),2X,3(E13.6,2X),2(F7.2,2X),E13.6
C       ou 3(2X,F7.2),1X,3(E13.6,2X),2(1X,F7.2),E13.6
C
C       2 cas : ITRPHI = 1, 
C               on obtient le plan de visées pour PHIos .
C               Entête de 31 lignes.
C               Nombre d'enregistrements : 
C                   un pour chaque valeur de -NBMU a NBMU (angles utiles)
C               Chaque enregistrement contient 8 valeurs :
C                   angle de visée(en degrés), angle de diffusion (en degrés),
C                   paramètres de Stokes I,Q,U (en sr-1), 
C                   angle de polarisation (en degrés), taux de polarisation (en %),
C                   intensité de polarisation (en sr-1).
C
C               ITRPHI = 2, 
C               on obtient un diagramme polaire en fonction de l'angle 
C               de visée (teta) et de l'angle azimutal (phi) 
C               Entête de 26 lignes.
C               Nombre d'enregistrements :
C                   un pour chaque valeur de 1 a NBMU (angles utiles)
C               Chaque enregistrement contient 9 valeurs :
C                   coordonnees phi, teta (en degrés), angle de diffusion (en degrés),
C                   paramètres de Stokes I,Q,U (en sr-1), 
C                   angle de polarisation (en degrés), taux de polarisation (en %),
C                   intensité de polarisation (en sr-1).
C
C   -->Contenu du fichier résultat pour le champ descendant
C      Fichier Ascii, format d'un enregistrement (F7.2,2X,3(F13.6,2X)) 
C      ou (2(F7.2,2X),3(F13.6,2X))
C      Meme description que pour le fichier résultat du champ montant
C


C
C   Fichiers de trace (optionnels)
C
C
C Variable d'environnement
C ------------------------
C   SOS_ABS_ROOT : Root of the SOS-ABS code directory
C
C
C Common utilisé:
C --------------
C  Aucun
C
C
C Cas d'erreur : 
C ------------        
C     -Decodage des paramètres incorrect
C     -Erreur dans un sous-programme
C     -Erreur à l'ouverture du fichier Trace
C     -Paramètres non valables
C
C     Affichage d'un message à l'écran, arrêt du programme et retour du status 1 au shell
C   
C
C  Sous programmes utilisés:
C --------------------------
C  Ce programme fait appel à la routine:
C     - SOS_PROC
C     - SOS_OUTPUT_HEADER
C     - SOS_OUTPUT_HEADER_POLAR_DIAG
C
C============================================================================


      PROGRAM SOS_ABS_MAIN
      
      IMPLICIT NONE
      

C Constantes 
C -----------------            
      DOUBLE PRECISION CONVDEGRAD         !Coefficient de conversion des degrés en radians.
      
      INTEGER*2 IDLOG			! Numéro identifiant du fichier Trace	    
      PARAMETER(IDLOG=INCTE_IDLOG_MAIN)
 



     	            
C* Définition des variables                         
C*-----------------------------------------------------------------
    
C-- Variables liées aux paramètres d'entrée 
      INTEGER*2 NBARGS		 ! Nombre d'arguments en entrée 
      INTEGER*4 IARGC     	        ! Type de la fonction IARGC
      INTEGER*4 NPARAM 		 ! Numéro du paramètre en entrée 
                                    ! (type INTEGER*4 obligatoire pour le f77)
	
      CHARACTER*10 PROG_ARG,DARG,IARG   ! Arguments du programme
					
      CHARACTER*CTE_LENKEYWORD KEYWORD  ! Mot clef du passage d'arguments.
      INTEGER*2 KEYWORD_VALID	 ! Indicateur de validité du mot clé lu
                                    
      
C-- Variables liées à l'arborescence de travail

      INTEGER*4 LDIR     	       ! longueur de la chaine de caractères RESROOT
      INTEGER*4 LDIRSOS     	! longueur de la chaine de caractères DIR_SOS
      INTEGER*4 LFIC    	       ! longueur de la chaine de caractères d'un nom de fichier

      CHARACTER*CTE_LENDIR RESROOT      ! Répertoire de sortie 
      CHARACTER*CTE_LENDIR DIR_SURF     ! Répertoire de sortie des calculs de surface
      CHARACTER*CTE_LENDIR DIR_MIE      ! Répertoire de sortie des calculs de Mie
      CHARACTER*CTE_LENDIR DIR_SOS      ! Répertoire de sortie des ordres successifs
      CHARACTER*CTE_LENDIR DIR_LOG      ! Répertoire de sortie des fichiers traces

C-- Variables globales

      DOUBLE PRECISION WA_SIMU	! Longueur d'onde de simulation des luminances (en microns)
      
      CHARACTER*CTE_LENFIC2 FICMAIN_LOG ! Nom du fichier trace principal 

      INTEGER*4 I,J,K,JJ,IP ! Indices de tableau

      INTEGER*4 IER        	! code d'erreur =0 si pas d'erreur, =-1 sinon

      LOGICAL TRACE     	! = vrai si écriture dans le fichier trace

                  
C---- Variables liées à la définition des angles  
      		
      DOUBLE PRECISION TETAS	! Angle zénithal solaire en degrés (0 < tetas < 90).
      				       ! (angle complémentaire : mus= -cos(tetas))
 
      INTEGER*4 NBMU_GAUSS_LUM     ! Nombre d'angles de Gauss utilisés pour les calculs de luminance      
      INTEGER*4 NBMU_GAUSS_MIE     ! Nombre d'angles de Gauss utilisés pour les calculs de Mie      

      !Calculs de luminance :
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_LUM  ! chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*CTE_LENFIC2 FICANGLES_RES_LUM   ! chemin et nom du fichier des angles utilisés 

      !Calculs de fonctions de phase :
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_MIE  ! chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*CTE_LENFIC2 FICANGLES_RES_MIE   ! chemin et nom du fichier des angles utilisés 
  
      CHARACTER*CTE_LENFIC2 FICANGLOG           ! nom du fichier trace


C---- Variables liées au modèle d'aérosols
 
      DOUBLE PRECISION AOT_REF	! Epaisseur optique des aérosols à la longueur d'onde de référence (non tronquée).
      DOUBLE PRECISION WAREF_AOT	! Longueur d'onde de référence
    
   	!Cas mono-mode aérosols :         		      		
      DOUBLE PRECISION RN_WA		! partie réelle de l'indice de réfraction à la longueur d'onde de simulation
      DOUBLE PRECISION IN_WA		! partie imaginaire de l'indice de réfraction à la longueur d'onde de simulation
      DOUBLE PRECISION RN_WAREF	! partie réelle de l'indice de réfraction à la longueur d'onde de référence
      DOUBLE PRECISION IN_WAREF	! partie imaginaire de l'indice de réfraction à la longueur d'onde de référence    
      DOUBLE PRECISION LND_RADIUS_MMD_AER ! Rayon modal de la LND.	
      DOUBLE PRECISION LND_LNVAR_MMD_AER  ! Variance de la LND. 	
      DOUBLE PRECISION JD_SLOPE_MMD_AER   ! Pente de la loi de Junge. 
      DOUBLE PRECISION JD_RMIN_MMD_AER    ! Rayon minimal de la loi de Junge.
      DOUBLE PRECISION JD_RMAX_MMD_AER    ! Rayon maximal de la loi de Junge.
    	!Proportion volumique des constituants WMO
      DOUBLE PRECISION C_WMO_DL	    	! Valeur utilisateur pour particule DL du modèle WMO
      DOUBLE PRECISION C_WMO_WS	    	! Valeur utilisateur pour particule WS du modèle WMO
      DOUBLE PRECISION C_WMO_OC	    	! Valeur utilisateur pour particule OC du modèle WMO
      DOUBLE PRECISION C_WMO_SO	    	! Valeur utilisateur pour particule SO du modèle WMO
      DOUBLE PRECISION RH                 ! % d'humidite relative des modèles de Shettle & Fenn

	!Paramètres pour un modèle bimodal de LND
      DOUBLE PRECISION USER_CV_COARSE	!Valeur utilisateur de la concentration volumique de la "LND coars
      DOUBLE PRECISION USER_CV_FINE	!Valeur utilisateur de la concentration volumique de la "LND fine"
      DOUBLE PRECISION RTAUct_WAREF       ! Rapport des épaisseurs optique AOT_coarse / AOT_total 
      					       ! pour la longueur d'onde de référence.    
      DOUBLE PRECISION BMD_CM_MRWA	       ! Valeur lue pour RN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MIWA	       ! Valeur lue pour IN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MRWAREF	! Valeur lue pour RN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MIWAREF	! Valeur lue pour RN_LND_WAREF(2) : fine mode 
      DOUBLE PRECISION BMD_CM_RMODAL	! Valeur lue pour RMODAL_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_VAR	       ! Valeur lue pour VAR_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MRWA	       ! Valeur lue pour RN_LND(2) : fine mode    					 
      DOUBLE PRECISION BMD_FM_MIWA	       ! Valeur lue pour IN_LND(2) : fine mode .
      DOUBLE PRECISION BMD_FM_MRWAREF	! Valeur lue pour IN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MIWAREF	! Valeur lue pour IN_LND_WAREF(2) : fine mode      
      DOUBLE PRECISION BMD_FM_RMODAL	! Valeur lue pour RMODAL_LND(2) : fine mode        
      DOUBLE PRECISION BMD_FM_VAR	       ! Valeur lue pour VAR_LND(2) : fine mode 
 			     
      DOUBLE PRECISION COEF_TRONCA	! Coefficient de troncature de la fonction de phase des aérosols	
                                  
      INTEGER*4 ITRONC_AER   	! Option de troncature de la fonction de phase des aérosols 
					! (1 pour la réaliser). 
      INTEGER*2 IMOD_AER   		! Granulometrie mono ou multi-modale.
      INTEGER*4 IGRANU	   	! Numéro du modèle de granulométrie 
      INTEGER*2 IModele_WMO 	! Type du modèle WMO.
      INTEGER*2 IModele_SF  	! Type du modèle Shettle & Fenn.
      INTEGER*2 MODE_PARAM_BILND   ! Indice du choix de description du mélange de 2 LND.
                      
      CHARACTER*CTE_LENFIC2 FICMIE_LOG      ! nom du fichier TRACE de MIE  
      CHARACTER*CTE_LENFIC2 FICGRANU_LOG    ! nom du fichier TRACE de GRANU
      CHARACTER*CTE_LENFIC2 FICGRANU        ! nom du fichier résultat GRANU  
      CHARACTER*CTE_LENFIC2 FICEXTDATA_AER  ! chemin et nom du fichier utilisateur donnant la
      						  ! fonction de phase et les coefficients de diffusion
						  ! et d'extinction.
      CHARACTER*CTE_LENFIC2 FICUSER_AER     ! chemin et nom du fichier aérosols dénommé par l'utilisateur
      CHARACTER*CTE_LENFIC2 FICMIXTURE_AER  ! chemin et nom du fichier de définition d'un mélange de modèles d'aérosols par l'utilisateur  
				
C --- Variables liées à la définition du profil atmosphérique 
  
      DOUBLE PRECISION TR     ! Epaisseur optique rayleigh.
      DOUBLE PRECISION HR     ! Echelle de hauteur du profil moléculaire (km)
      DOUBLE PRECISION HA     ! Echelle de hauteur du profil d'aérosols (km)
      DOUBLE PRECISION ZMIN   ! Altitude minimale de la couche molécules + aérosols.
      DOUBLE PRECISION ZMAX   ! Altitude maximale de la couche molécules + aérosols.
      DOUBLE PRECISION PSURF  ! Pression atmosphérique de surface (mbar).
      DOUBLE PRECISION H2O    ! Contenu intégré en vapeur d'eau (g/cm2).
      DOUBLE PRECISION O3     ! Contenu intégré en ozone (Dobson).
      DOUBLE PRECISION CO2    ! Concentration du dioxyde de carbone en surface (ppmv).
      DOUBLE PRECISION CH4    ! Concentration du méthane en surface (ppmv).

      DOUBLE PRECISION TDIR_VRAI   !Transmission directe pour la direction d'incidence solaire
      				   	!de l'atmosphere vraie

      INTEGER*2 IPROFIL	  ! Type du profil aérosol
      INTEGER*2 ABSPROFIL     ! Type du profil en absorption
        
      CHARACTER*CTE_LENFIC2 FICABSPROFIL    ! Fichier profil en absorption fourni par utilisateur
      CHARACTER*CTE_LENFIC2 FICPROFIL_LOG   ! Chemin et nom du fichier Trace des profils par ordre CKD
      
 
C --- Variables liées à la définition du choix des données de coefficients CKD
        
      DOUBLE PRECISION  NUSTEP  ! Pas du nombre d'onde (cm-1) dans le fichier de coeff CKD
                
C --- Variables liées à la définition du modèle de surface 
  
      INTEGER*2 ISURF	       ! Type de surface   

      DOUBLE PRECISION RHO	       ! Albédo de surface.        
      DOUBLE PRECISION WIND	       ! Vitesse du vent (m/s).    
      DOUBLE PRECISION SURF_IND	! Indice de réfraction du dioptre.	   
      DOUBLE PRECISION K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN 	! Paramètres du modèle de Roujean     
      DOUBLE PRECISION ALPHA_NADAL	! Paramètre alpha du modèle de Nadal.
      DOUBLE PRECISION BETA_NADAL	! Paramètre beta du modèle de Nadal.   
      DOUBLE PRECISION COEF_C_MAIGNAN  ! Paramètre C*exp(-NDVI) du modèle de Maignan.   
       
      CHARACTER*CTE_LENFIC2 FICSURF_LOG   ! nom du fichier trace SURFACE.
      CHARACTER*CTE_LENFIC2 FICSURF       ! nom d'un fichier contenant une matrice de réflexion existante
      

C --- Variables liées aux calculs SOS

      DOUBLE PRECISION ZOUT	! Type de sortie :
				!      -1 =>  Sortie standard: Luminance émergeante au TOA
				!                              + descendante au niveau du sol.
				!       z =>  Luminance montante et descendante
				!             pour l'altitude z (0 <= z <= CTE_TOA_ALT).
			
      INTEGER*4 IGMAX	! Ordre maximal de diffusion.
      
      INTEGER*4 IPOLAR	! Option de simulation sans polarisation.
      			       ! 0 : calculs sans polarisation.
			       ! 1 : calculs normaux avec polarisation.

      INTEGER*4 IMODE_CKD_CALCUL ! Option de mode d'application de la décomposition CKD
                                 ! 1 : Calcul le plus fin avec un calcul de luminance par ordre de la décomposition
                                 !     puis recomposition selon la pondération CKD.
                                 ! 2 : Calcul pour une épaisseur optique d'absorption estimée par l'approche CKD
                                 !     --> Un unique calcul de luminance

      INTEGER*4 ITRPHI	! Indice de type de sortie.
			       !   1 =>  Plan de visée pour une valeur d'azimut.
			       !   2 =>  Diagramme polaire en fonction de la direction
			       !         d'azimut et de visée.
			
      INTEGER*4 PAS_PHI	! Pas sur l'azimut en degrés (pour ITRPHI=2)
      			
      DOUBLE PRECISION PHIOS ! Azimut relatif en degrés (pour ITRPHI=1).
       
      CHARACTER*CTE_LENFIC2 FICSOS_LOG ! chemin et nom du fichier Trace
      CHARACTER*CTE_LENFIC2 FICUP	    ! Nom du fichier résultat pour le champ montant.					
      CHARACTER*CTE_LENFIC2 FICDOWN    ! Nom du fichier résultat pour le champ descendant.
      CHARACTER*CTE_LENFIC2 FICUP_ANGLES_USER ! Nom du fichier résultat pour le champ
      				                  ! montant (pour les angles utilisateur).					
      CHARACTER*CTE_LENFIC2 FICDOWN_ANGLES_USER	! Nom du fichier résultat pour le champ
      				                     ! descendant (pour les angles utilisateur).			
      CHARACTER*CTE_LENFIC2 FICTRANS    ! Nom du fichier optionnel de sortie des transmissions.		
      CHARACTER*CTE_LENFIC2 FICFLUX     ! Nom du fichier optionnel de sortie des flux up TOA et down BOA.

      CHARACTER*CTE_LENFIC2 FICSOS_RES_BIN  ! Nom du fichier résultat des OS (pour un traitement unitaire)

   
C --- Variables internes
      
C -   Pour les OS

      INTEGER*4 LUM_NBMU     ! Nombre d'angles effectifs des simulations

      INTEGER*4 IND_ANGOUT_FIN(0:CTE_OS_NBMU_MAX) 
                        ! Indices des tableaux d'angles pour lesquels il doit y avoir 
      			   ! une sortie dans le fichier résultat des angles utilisateur 
			   ! (1 : pour une sortie requise)		 
                        				 

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
         
      DOUBLE PRECISION FLUX_DIFF_DOWN	!Flux diffus descendant au niveau de la surface 
      DOUBLE PRECISION FLUX_DOWN		!Flux total descendant au niveau de la surface

      DOUBLE PRECISION EPLUS    ! Flux diffus émergeant au TOA, normalisé par le flux solaire incident au TOA.
    
C       
      INTEGER*4 IPHI        ! Indice de boucle sur l'azimut relatif.
   
      DOUBLE PRECISION XI   ! Valeur du paramètre de Stokes I pour un angle
      DOUBLE PRECISION XQ   ! Valeur du paramètre de Stokes Q pour un angle
      DOUBLE PRECISION XU   ! Valeur du paramètre de Stokes U pour un angle  

      DOUBLE PRECISION XAN  ! Angle de la polarisation (degrés)	
      DOUBLE PRECISION LPOL ! Intensité polarisée normalisée à l'éclairement solaire au TOA 	       
      DOUBLE PRECISION TPOL ! Taux de polarisation en %
      DOUBLE PRECISION PHI1, PHI2 ! Azimut relatif (en radians).
   
   
C =========================================
C --- Initialisations 
C =========================================

C --- Initialisation de l'indicateur d'erreur de traitement
C----------------------------------------------------------
      IER=0
	 
C --- Initialisation des paramètres optionnels avec leur valeur par défaut 
C-----------------------------------------------------------------------------

      FICMAIN_LOG="NO_LOG_FILE"

C --- Initialisation des valeurs par défaut des paramètres pour les angles

C     Nombre d'angles de Gauss pour les calculs de luminance : non défini     
      NBMU_GAUSS_LUM=CTE_NOT_DEFINED_VALUE_INT
      
C     Nombre d'angles de Gauss pour les calculs de Mie : non défini      
      NBMU_GAUSS_MIE=CTE_NOT_DEFINED_VALUE_INT
      
C     Absence de fichier d'angles "utilisateur" pour les calculs de luminance  : non défini     
      FICANGLES_USER_LUM="NO_USER_ANGLES"
      
C     Absence de fichier d'angles "utilisateur" pour les calculs de Mie : non défini      
      FICANGLES_USER_MIE="NO_USER_ANGLES"
 
C     Pas de fichier trace de la définition des angles
      FICANGLOG = "NO_LOG_FILE"
       
C     Fichiers résultats
      FICANGLES_RES_LUM = CTE_DEFAULT_FICANGLES_RES_LUM
      FICANGLES_RES_MIE = CTE_DEFAULT_FICANGLES_RES_MIE

C --- Initialisation des paramètres obligatoires
      WA_SIMU=CTE_NOT_DEFINED_VALUE_DBLE   
      TETAS=CTE_NOT_DEFINED_VALUE_DBLE   
      AOT_REF=CTE_NOT_DEFINED_VALUE_DBLE

      RESROOT="UNDEFINED_REPERTORY"
      DIR_SURF="UNDEFINED_REPERTORY"
      DIR_MIE="UNDEFINED_REPERTORY"
      DIR_SOS="UNDEFINED_REPERTORY"
      DIR_LOG="UNDEFINED_REPERTORY"

      IMOD_AER=CTE_NOT_DEFINED_VALUE_INT
      WAREF_AOT=CTE_NOT_DEFINED_VALUE_DBLE
      RHO=CTE_NOT_DEFINED_VALUE_DBLE
      ISURF=CTE_NOT_DEFINED_VALUE_INT
      IPROFIL=CTE_NOT_DEFINED_VALUE_INT



C --- Initialisation des valeurs par défaut des paramètres pour les aérosols

c     Paramètres pour le model mono-modal
      RN_WA=CTE_NOT_DEFINED_VALUE_DBLE
      IN_WA=CTE_NOT_DEFINED_VALUE_DBLE
      RN_WAREF=CTE_NOT_DEFINED_VALUE_DBLE
      IN_WAREF=CTE_NOT_DEFINED_VALUE_DBLE 
      IGRANU=CTE_NOT_DEFINED_VALUE_INT
      LND_RADIUS_MMD_AER=CTE_NOT_DEFINED_VALUE_DBLE	
      LND_LNVAR_MMD_AER=CTE_NOT_DEFINED_VALUE_DBLE
      JD_SLOPE_MMD_AER=CTE_NOT_DEFINED_VALUE_DBLE
      JD_RMIN_MMD_AER=CTE_NOT_DEFINED_VALUE_DBLE
      JD_RMAX_MMD_AER=CTE_DEFAULT_AER_JUNGE_RMAX		!Rmax par défaut  

C     Paramètres model WMO
      C_WMO_DL=CTE_NOT_DEFINED_VALUE_DBLE	
      C_WMO_WS=CTE_NOT_DEFINED_VALUE_DBLE	
      C_WMO_OC=CTE_NOT_DEFINED_VALUE_DBLE     
      C_WMO_SO=CTE_NOT_DEFINED_VALUE_DBLE	
      IModele_WMO=CTE_NOT_DEFINED_VALUE_INT

C     Paramètres S&F
      IModele_SF=CTE_NOT_DEFINED_VALUE_INT          
      RH=CTE_NOT_DEFINED_VALUE_DBLE

C     Paramètres bimode LND
      MODE_PARAM_BILND = CTE_NOT_DEFINED_VALUE_INT 
      USER_CV_COARSE = CTE_NOT_DEFINED_VALUE_DBLE 
      USER_CV_FINE 	= CTE_NOT_DEFINED_VALUE_DBLE  
      RTAUct_WAREF 	= CTE_NOT_DEFINED_VALUE_DBLE  
      BMD_CM_MRWA 	= CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_CM_MIWA 	= CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_CM_MRWAREF = CTE_NOT_DEFINED_VALUE_DBLE  
      BMD_CM_MIWAREF = CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_CM_RMODAL 	= CTE_NOT_DEFINED_VALUE_DBLE  
      BMD_CM_VAR 	= CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_FM_MRWA	= CTE_NOT_DEFINED_VALUE_DBLE  
      BMD_FM_MIWA	= CTE_NOT_DEFINED_VALUE_DBLE
      BMD_FM_MRWAREF = CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_FM_MIWAREF = CTE_NOT_DEFINED_VALUE_DBLE
      BMD_FM_RMODAL	= CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_FM_VAR 	= CTE_NOT_DEFINED_VALUE_DBLE 	 

C     Paramètres données utilisateur
      FICEXTDATA_AER = "NO_USER_AEROSOLS_PHAZE_FCT" 

C     Paramètres mélange de modèles par utilisateur
      FICMIXTURE_AER = "NO_USER_AEROSOLS_MIXTURE" 
      
C     Nom du fichier d'aerosols utilisateur
      FICUSER_AER="NO_USER_AEROSOLS"	 
              
C     Nom par défaut du fichier des propriétés radiatives des aérosols
      FICGRANU = CTE_DEFAULT_FICGRANU
     
C     Pas de fichier trace des calculs de Mie pour les aérosols
      FICMIE_LOG = "NO_LOG_FILE"

C     Pas de fichier trace de la définition des matrices de phase aérosols
      FICGRANU_LOG = "NO_LOG_FILE" 

C     Troncature de la pointe avant de la fonction de phase aérosols
      ITRONC_AER = 1
 
C --- Initialisation des valeurs par défaut des paramètres de surface
      
C     Pas de fichier trace de la définition des matrices de phase aérosols
      FICSURF_LOG = "NO_LOG_FILE"

C     Indice de réfraction
      SURF_IND=CTE_NOT_DEFINED_VALUE_DBLE
      
C     Paramètres glitter
      WIND=CTE_NOT_DEFINED_VALUE_DBLE
      
C     Paramètres Roujean
      K0_ROUJEAN=CTE_NOT_DEFINED_VALUE_DBLE
      K1_ROUJEAN=CTE_NOT_DEFINED_VALUE_DBLE
      K2_ROUJEAN=CTE_NOT_DEFINED_VALUE_DBLE 
      
C     Paramètres Nadal
      ALPHA_NADAL=CTE_NOT_DEFINED_VALUE_DBLE
      BETA_NADAL=CTE_NOT_DEFINED_VALUE_DBLE   

C     Paramètres Maignan
      COEF_C_MAIGNAN=CTE_NOT_DEFINED_VALUE_DBLE     

C     Fichier résultat de BRDF/BPDF
      FICSURF="DEFAULT"
      		
C --- Initialisation des valeurs par défaut des paramètres de profil
     
C     Epaisseur optique moléculaire
      TR=CTE_NOT_DEFINED_VALUE_DBLE
           
C     Echelles de hauteur
      HA=CTE_NOT_DEFINED_VALUE_DBLE
      HR=CTE_NOT_DEFINED_VALUE_DBLE
      
C     Pression (valeur par défaut)
      PSURF=CTE_HT_STD_PSURF
 
C     Altitudes min et max du profil aérosol
      ZMIN=CTE_NOT_DEFINED_VALUE_DBLE
      ZMAX=CTE_NOT_DEFINED_VALUE_DBLE     
 
C     Profil d'absorption
      ABSPROFIL=CTE_NOT_DEFINED_VALUE_INT

C     Contenus en gaz absorbants
      H2O=CTE_NOT_DEFINED_VALUE_DBLE
      O3=CTE_NOT_DEFINED_VALUE_DBLE
      CO2=CTE_NOT_DEFINED_VALUE_DBLE
      CH4=CTE_NOT_DEFINED_VALUE_DBLE

C     Fichier de profil de gaz absorbants
      FICABSPROFIL = "NO_USER_ABS_PROFILE_FILE"             

C     Données CKD
      NUSTEP=CTE_NOT_DEFINED_VALUE_DBLE
      IMODE_CKD_CALCUL=CTE_NOT_DEFINED_VALUE_INT
  
C     Pas de fichier trace 
      FICPROFIL_LOG = "NO_LOG_FILE"

 
            
C --- Initialisation des valeurs par défaut des paramètres SOS
      
C     Pas de fichier trace 
      FICSOS_LOG = "NO_LOG_FILE"

C     Option de coupure de la polarisation : optionnelle (calcul avec polarisation par défaut)
      IPOLAR=1   

C     Nombre maximal de diffusion
      IGMAX=CTE_DEFAULT_IGMAX

      
C     Niveau de sortie des simulations de luminance : optionnel  
      ZOUT=-1.0

C     Fichiers de sortie sur angles utilisateurs
      FICUP_ANGLES_USER="NO_OUTPUT"
      FICDOWN_ANGLES_USER="NO_OUTPUT"            

C     Fichier de sortie des transmissions
      FICTRANS="NO_OUTPUT"

C     Fichier de sortie des flux
      FICFLUX="NO_OUTPUT"

C     Fichiers de résultats
      FICSOS_RES_BIN=CTE_DEFAULT_RESBIN
      FICUP=CTE_DEFAULT_RESUP
      FICDOWN=CTE_DEFAULT_RESDOWN
   
C     Définition des azimut
      ITRPHI=CTE_NOT_DEFINED_VALUE_INT  
      PHIOS=CTE_NOT_DEFINED_VALUE_DBLE
      PAS_PHI=CTE_NOT_DEFINED_VALUE_INT    
  
             		       
C =========================================================
C --- Lecture des paramètres/arguments d'entrée et decodage
C =========================================================           				
      
C --- Nombre d'arguments en entrée
      NBARGS = IARGC()
      
C --- Boucle de lecture des couples d'arguments       
      DO NPARAM=1,NBARGS, 2
         CALL GETARG(NPARAM,KEYWORD)
         KEYWORD_VALID=0
		
C ------- Controle élémentaire de validité du keyword	  
         IF (KEYWORD(1:1).NE.'-') GOTO 100

C ------- Arguments propres au programme principal

	 IF (KEYWORD.EQ."-SOS_Main.ResRoot") THEN
	     CALL GETARG (NPARAM+1,RESROOT)		!Nom du répertoire de travail.
            KEYWORD_VALID=1
	 ENDIF	

	 IF (KEYWORD.EQ."-SOS_Main.Log") THEN
	     CALL GETARG (NPARAM+1,FICMAIN_LOG)		!Nom du fichier Trace global.
            KEYWORD_VALID=1
	 ENDIF		 
		
	 IF (KEYWORD.EQ."-SOS_Main.Wa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)     !Longueur d'onde de simulation de la luminance
	     READ(PROG_ARG,*,ERR=101) WA_SIMU    
            KEYWORD_VALID=1
	 ENDIF	
	 
C ------- Arguments propres à la définition des angles	
	
	 IF (KEYWORD.EQ."-ANG.Rad.NbGauss") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Nb d'angles de Gauss pour les simulations de luminance
	     READ(PROG_ARG,*,ERR=201) NBMU_GAUSS_LUM
            KEYWORD_VALID=1
	 ENDIF	 
	 	 
	 IF (KEYWORD.EQ."-ANG.Rad.UserAngFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_USER_LUM)	!Nom du fichier des angles "utilisateurs" (chemin complet)
            KEYWORD_VALID=1			              !pour les calculs de luminance
	 ENDIF	  	 
	 		 
	 IF (KEYWORD.EQ."-ANG.Thetas") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Angle zenithal solaire (degrés)
	     READ(PROG_ARG,*,ERR=202) TETAS   
            KEYWORD_VALID=1
	 ENDIF	

	 IF (KEYWORD.EQ."-ANG.Rad.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_RES_LUM)	!Nom du fichier résultat (chemin complet)
            KEYWORD_VALID=1	     				!pour les calculs de luminance
	 ENDIF	  	
				
	 IF (KEYWORD.EQ."-ANG.Aer.NbGauss") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Nb d'angles de Gauss pour les calculs de Mie
	     READ(PROG_ARG,*,ERR=203) NBMU_GAUSS_MIE
            KEYWORD_VALID=1
	 ENDIF	 
	 	 
	 IF (KEYWORD.EQ."-ANG.Aer.UserAngFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_USER_MIE)	!Nom du fichier des angles "utilisateurs" (chemin complet)
            KEYWORD_VALID=1		        	       !pour les calculs de Mie
	 ENDIF	  	 
	 	 	 
	 IF (KEYWORD.EQ."-ANG.Aer.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_RES_MIE)	!Nom du fichier résultat (chemin complet)
            KEYWORD_VALID=1	     				!pour les calculs de Mie
	 ENDIF	  	
	 	 
	 IF (KEYWORD.EQ."-ANG.Log") THEN
	     CALL GETARG (NPARAM+1,FICANGLOG)		!Nom du fichier Trace.
            KEYWORD_VALID=1
	 ENDIF		

C ------- Arguments propres à la définition des modèles d'aérosols

	 IF (KEYWORD.EQ."-AER.Waref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)     	!Longueur d'onde de référence 
	     READ(PROG_ARG,*,ERR=301) WAREF_AOT 	!pour l'ep. opt aérosol.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.AOTref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !AOT totale pour la longueur d'onde WAREF 
	     READ(PROG_ARG,*,ERR=302) AOT_REF
            KEYWORD_VALID=1
	 ENDIF 	 

	 IF (KEYWORD.EQ."-AER.Tronca") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)       	!Indice de troncature
	     READ(PROG_ARG,*,ERR=303) ITRONC_AER
            KEYWORD_VALID=1
	 ENDIF	 
	 
	 IF (KEYWORD.EQ."-AER.Log") THEN
	     CALL GETARG (NPARAM+1,FICGRANU_LOG) !Fichier Trace du programme Aérosols
            KEYWORD_VALID=1
	 ENDIF	 	  	 

	 IF (KEYWORD.EQ."-AER.MieLog") THEN
	     CALL GETARG (NPARAM+1,FICMIE_LOG)	!Fichier Trace des calculs de Mie   
            KEYWORD_VALID=1
	 ENDIF
         
         IF (KEYWORD.EQ."-AER.DirMie") THEN
	     CALL GETARG (NPARAM+1,DIR_MIE)	!Répertoire de stockage des fichiers de Mie
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICGRANU)	!Fichier résultat
            KEYWORD_VALID=1
	 ENDIF	 	  	
	 	 	  
	 IF (KEYWORD.EQ."-AER.Model") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Type de modèle d'aérosols
	     READ(PROG_ARG,*,ERR=304) IMOD_AER
            KEYWORD_VALID=1
	 ENDIF	 				   	 	  	  

	 IF (KEYWORD.EQ."-AER.MMD.MRwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!partie réelle de l'indice de réfraction pour modèle mono-mode
	     READ(PROG_ARG,*,ERR=306) RN_WA
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.MMD.MIwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!partie imaginaire de l'indice de réfraction pour modèle mono-mode
	     READ(PROG_ARG,*,ERR=307) IN_WA
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.MMD.MRwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Partie réelle de l'indice de réfraction pour modèle 
	     READ(PROG_ARG,*,ERR=308) RN_WAREF		!aérosol mono-mode à la longueur d'onde de référence 
            KEYWORD_VALID=1                      	
	 ENDIF	 					

	 IF (KEYWORD.EQ."-AER.MMD.MIwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Partie réelle de l'indice de réfraction pour modèle
	     READ(PROG_ARG,*,ERR=309) IN_WAREF	 	!daérosol mono-mode à la longueur d'onde de référence 
            KEYWORD_VALID=1                      	
	 ENDIF	 					
	 
	 IF (KEYWORD.EQ."-AER.MMD.SDtype") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!indice de type de granulométrie pour modèle mono-mode
	     READ(PROG_ARG,*,ERR=310) IGRANU
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.MMD.LNDradius") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Rayon modal (microns) de granulométrie LND pour aérosol mono-mode
	     READ(PROG_ARG,*,ERR=311) LND_RADIUS_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.LNDvar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Variance de granulométrie LND pour aérosol mono-mode
	     READ(PROG_ARG,*,ERR=312) LND_LNVAR_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.JD.slope") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Pente de loi de Junge pour modèle aérosol mono-mode
	     READ(PROG_ARG,*,ERR=3121) JD_SLOPE_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.JD.rmin") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Rayon minimal de la loi de Junge (microns) pour modèle aérosol mono-mode
	     READ(PROG_ARG,*,ERR=3122) JD_RMIN_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.JD.rmax") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Rayon maximal de la loi de Junge (microns) pour modèle aérosol mono-mode
	     READ(PROG_ARG,*,ERR=3123) JD_RMAX_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF							
	 
	 IF (KEYWORD.EQ."-AER.WMO.Model") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Indice de choix du modèle WMO
	     READ(PROG_ARG,*,ERR=313) IModele_WMO 
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.DL") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Pourcentage volumique des particules DL du modèle WMO
	     READ(PROG_ARG,*,ERR=314) C_WMO_DL
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.WMO.WS") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)    	!Pourcentage volumique des particules WS du modèle WMO
	     READ(PROG_ARG,*,ERR=314) C_WMO_WS
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.OC") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   	!Pourcentage volumique des particules OC
	     READ(PROG_ARG,*,ERR=314) C_WMO_OC
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.SO") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   	!Pourcentage volumique des particules SO
	     READ(PROG_ARG,*,ERR=314) C_WMO_SO
            KEYWORD_VALID=1
	 ENDIF	
	 
	 IF (KEYWORD.EQ."-AER.SF.Model") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Indice de choix du modèle S&F
	     READ(PROG_ARG,*,ERR=315) IModele_SF  
            KEYWORD_VALID=1
	 ENDIF 	
	 
	 IF (KEYWORD.EQ."-AER.SF.RH") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Humidite relative pour modèle S&F
	     READ(PROG_ARG,*,ERR=316) RH 
            KEYWORD_VALID=1
	 ENDIF 	 
	 	 
	 IF (KEYWORD.EQ."-AER.BMD.VCdef") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Indice de type de description du mélange bimodal LND
	     READ(PROG_ARG,*,ERR=317) MODE_PARAM_BILND 
            KEYWORD_VALID=1
	 ENDIF 	 	 

	 IF (KEYWORD.EQ."-AER.BMD.CoarseVC") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Concentration volumique du mode "coarse"
	     READ(PROG_ARG,*,ERR=318) USER_CV_COARSE
            KEYWORD_VALID=1
	 ENDIF 	 

	 IF (KEYWORD.EQ."-AER.BMD.FineVC") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Concentration volumique du mode "fine"
	     READ(PROG_ARG,*,ERR=319) USER_CV_FINE 
            KEYWORD_VALID=1
	 ENDIF 	 	 

	 IF (KEYWORD.EQ."-AER.BMD.RAOT") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	 !AOT_coarse / AOT_tot pour la longueur d'onde WAREF
	     READ(PROG_ARG,*,ERR=320) RTAUct_WAREF 
            KEYWORD_VALID=1
	 ENDIF 	 
	 

	 IF (KEYWORD.EQ."-AER.BMD.CM.MRwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=321) BMD_CM_MRWA
            KEYWORD_VALID=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.CM.MIwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=322) BMD_CM_MIWA
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.CM.MRwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=323) BMD_CM_MRWAREF
            KEYWORD_VALID=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.CM.MIwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=324) BMD_CM_MIWAREF
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.CM.SDradius") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=325) BMD_CM_RMODAL
            KEYWORD_VALID=1
	 ENDIF 

	 IF (KEYWORD.EQ."-AER.BMD.CM.SDvar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=326) BMD_CM_VAR
            KEYWORD_VALID=1
	 ENDIF 
	 	 
	 IF (KEYWORD.EQ."-AER.BMD.FM.MRwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=327) BMD_FM_MRWA
            KEYWORD_VALID=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.FM.MIwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=328) BMD_FM_MIWA	
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.FM.MRwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=329) BMD_FM_MRWAREF
            KEYWORD_VALID=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.FM.MIwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=330) BMD_FM_MIWAREF
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.FM.SDradius") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=331) BMD_FM_RMODAL  
            KEYWORD_VALID=1
	 ENDIF 

	 IF (KEYWORD.EQ."-AER.BMD.FM.SDvar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=332) BMD_FM_VAR	 
            KEYWORD_VALID=1
	 ENDIF 	
	 
	 IF (KEYWORD.EQ."-AER.ExtData") THEN
	     CALL GETARG (NPARAM+1,FICEXTDATA_AER) 
            KEYWORD_VALID=1
	 ENDIF 	 	 	 	 	 	 		  	 	 	 	 	 	 	 		  
	 
	 IF (KEYWORD.EQ."-AER.DefMixture") THEN
	     CALL GETARG (NPARAM+1,FICMIXTURE_AER) 
            KEYWORD_VALID=1
	 ENDIF 	
	 	 	 
	 IF (KEYWORD.EQ."-AER.UserFile") THEN
	     CALL GETARG (NPARAM+1,FICUSER_AER) 
            KEYWORD_VALID=1
	 ENDIF
	 
C ------- Arguments propres à la définition du profil atmosphérique 
	 IF (KEYWORD.EQ."-AP.Log") THEN
	     CALL GETARG (NPARAM+1,FICPROFIL_LOG) !Fichier Trace du programme
            KEYWORD_VALID=1
	 ENDIF	 	  	 
	 	 
	 IF (KEYWORD.EQ."-AP.MOT") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Epaisseur optique moléculaire
	     READ(PROG_ARG,*,ERR=401) TR
            KEYWORD_VALID=1
	 ENDIF	
	 	 
	 IF (KEYWORD.EQ."-AP.HR") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Echelle de hauteur moléculaire (km)
	     READ(PROG_ARG,*,ERR=402) HR
            KEYWORD_VALID=1
	 ENDIF	

	 IF (KEYWORD.EQ."-AP.AerHS.HA") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Echelle de hauteur des aérosols (km)
	     READ(PROG_ARG,*,ERR=403) HA
            KEYWORD_VALID=1
	 ENDIF	

	 IF (KEYWORD.EQ."-AP.AerProfile.Type") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Type de profil des aérosols simule
	     READ(PROG_ARG,*,ERR=404) IPROFIL 
            KEYWORD_VALID=1
	 ENDIF	
	 	 	 	 
	 IF (KEYWORD.EQ."-AP.AerLayer.Zmin") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Altitude minimale de la couche d'aérosols (km)
	     READ(PROG_ARG,*,ERR=405) ZMIN
            KEYWORD_VALID=1
	 ENDIF		 
	 	 
	 IF (KEYWORD.EQ."-AP.AerLayer.Zmax") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Altitude maximale de la couche d'aérosols (km)
	     READ(PROG_ARG,*,ERR=406) ZMAX
            KEYWORD_VALID=1
	 ENDIF	
	 
	 IF (KEYWORD.EQ."-AP.Psurf") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Pression à la surface
	     READ(PROG_ARG,*,ERR=407) PSURF
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.H2O") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Contenu en vapeur d'eau
	     READ(PROG_ARG,*,ERR=408) H2O
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.O3") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Contenu en ozone
	     READ(PROG_ARG,*,ERR=409) O3
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.CO2") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Concentration en dioxyde de carbone en surface
	     READ(PROG_ARG,*,ERR=4091) CO2
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.CH4") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Concentration en méthane en surface
	     READ(PROG_ARG,*,ERR=4092) CH4
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.AbsProfile.Type") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Type de profil pour l'absorption
	     READ(PROG_ARG,*,ERR=410) ABSPROFIL 
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.AbsProfile.UserFile") THEN
	     CALL GETARG (NPARAM+1,FICABSPROFIL)	!Fichier profil d'absorption Utilisateur
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.SpectralResol") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Résolution spectrale des fichiers CKD à utiliser (cm-1)
	     READ(PROG_ARG,*,ERR=411) NUSTEP 
            KEYWORD_VALID=1
	 ENDIF
     
     
C ------- Arguments propres au model de surface

	 IF (KEYWORD.EQ."-SURF.Type") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
	     READ(PROG_ARG,*,ERR=501) ISURF
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Dir") THEN
	     CALL GETARG (NPARAM+1,DIR_SURF)
            KEYWORD_VALID=1
	 ENDIF	  	

	 IF (KEYWORD.EQ."-SURF.Log") THEN
	     CALL GETARG (NPARAM+1,FICSURF_LOG)
            KEYWORD_VALID=1
	 ENDIF	  	 

	 IF (KEYWORD.EQ."-SURF.Ind") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
	     READ(PROG_ARG,*,ERR=502) SURF_IND
            KEYWORD_VALID=1
	 ENDIF
	  
	 IF (KEYWORD.EQ."-SURF.Glitter.Wind") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=503) WIND    
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Roujean.K0") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=504) K0_ROUJEAN    
            KEYWORD_VALID=1
	 ENDIF
	  
	 IF (KEYWORD.EQ."-SURF.Roujean.K1") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=505) K1_ROUJEAN    
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Roujean.K2") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
             READ(PROG_ARG,*,ERR=506) K2_ROUJEAN    
            KEYWORD_VALID=1
	 ENDIF
	  	 
	 IF (KEYWORD.EQ."-SURF.Nadal.Alpha") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=507) ALPHA_NADAL    
            KEYWORD_VALID=1
	 ENDIF	 

        IF (KEYWORD.EQ."-SURF.Nadal.Beta") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=508) BETA_NADAL	
            KEYWORD_VALID=1
        ENDIF	  

        IF (KEYWORD.EQ."-SURF.Maignan.C") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=5081) COEF_C_MAIGNAN	
            KEYWORD_VALID=1
        ENDIF	

	 IF (KEYWORD.EQ."-SURF.Alb") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Albédo de surface
	     READ(PROG_ARG,*,ERR=509) RHO
            KEYWORD_VALID=1
	 ENDIF 	 	 

	 IF (KEYWORD.EQ."-SURF.File") THEN
	     CALL GETARG (NPARAM+1,FICSURF)	!Nom du fichier de Matrice de réflexion
            KEYWORD_VALID=1
	 ENDIF
	 
C ------- Arguments propres aux calculs d'Ordres Successifs et à la production
C         des fichiers résultats

	 IF (KEYWORD.EQ."-SOS.Log") THEN
	     CALL GETARG (NPARAM+1,FICSOS_LOG)	!Nom du fichier Trace.
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-SOS.ResBin") THEN
	     CALL GETARG (NPARAM+1,FICSOS_RES_BIN) !Nom du fichier OS binaire.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileUp") THEN
	     CALL GETARG (NPARAM+1,FICUP)	!Nom du fichier de sortie pour le champ UP.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileDown") THEN
	     CALL GETARG (NPARAM+1,FICDOWN) !Nom du fichier de sortie pour le champ DOWN.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileUp.UserAng") THEN
            CALL GETARG (NPARAM+1,FICUP_ANGLES_USER) !Nom du fichier de sortie pour le champ UP 
            KEYWORD_VALID=1     !aux angles "utilisateur".	     						     
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileDown.UserAng") THEN
	     CALL GETARG (NPARAM+1,FICDOWN_ANGLES_USER)	!Nom du fichier de sortie pour le champ DOWN. 
            KEYWORD_VALID=1                             !aux angles "utilisateur".	     		
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.Trans") THEN
	     CALL GETARG (NPARAM+1,FICTRANS)  !Nom du fichier de sortie des valeurs de transmission.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.Flux") THEN
	     CALL GETARG (NPARAM+1,FICFLUX)  !Nom du fichier de sortie des valeurs de flux
            KEYWORD_VALID=1
	 ENDIF
          	 	 		 
         IF (KEYWORD.EQ."-SOS.OutputAlt") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Altitude de sortie
            READ(PROG_ARG,*,ERR=601) ZOUT  
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.IGmax") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Ordre maximal d'interaction
            READ(PROG_ARG,*,ERR=602) IGMAX
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.Ipolar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Option de coupure de la polarisation
            READ(PROG_ARG,*,ERR=603) IPOLAR
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.View") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Option de visée
            READ(PROG_ARG,*,ERR=604) ITRPHI
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.View.Phi") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Angle azimutal relatif 
            READ(PROG_ARG,*,ERR=605) PHIOS 
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.View.Dphi") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Pas en angles azimutaux
            READ(PROG_ARG,*,ERR=606) PAS_PHI     
            KEYWORD_VALID=1
         ENDIF	 		 

         IF (KEYWORD.EQ."-SOS.AbsModeCKD") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Option de coupure de la polarisation
            READ(PROG_ARG,*,ERR=608) IMODE_CKD_CALCUL
            KEYWORD_VALID=1
         ENDIF
	 

         IF (KEYWORD_VALID.EQ.0) THEN
	     WRITE(6,*) 
	     WRITE(6,*)	     
     &      "**************************************************"
	     WRITE(6,*) "WARNING:"
	     WRITE(6,*) "The defined parameter:",KEYWORD
	     WRITE(6,*) "is not part of the SOS-ABS parameters."
	     WRITE(6,*) 
     &      "==> Its definition does not act on the simulation."
	     WRITE(6,*) 
     &      "==> Be sure not to confuse the keyword parameters!"
	     WRITE(6,*)
     &      "**************************************************"
	     WRITE(6,*) 
         ENDIF
	 		
      ENDDO	!Fin de boucle de lecture des arguments


C =========================================================
C --- Définition de l'arborescence de travail
C =========================================================           				

C ---- Controle de la définition de l'espace de travail (obligatoire)
      IF (RESROOT.EQ."UNDEFINED_REPERTORY") GOTO 1000
      
      LDIR = INDEX(RESROOT,' ') - 1
      IF((LDIR.LE.0).OR.(LDIR.GT.CTE_LENDIR)) GOTO 1100
		
C ---- Définition de l'arborescence de résultats SOS
      LDIRSOS= LDIR + 4  	
      IF(LDIRSOS.GT.CTE_LENDIR) GOTO 1108
      DIR_SOS(1:LDIRSOS) = RESROOT(1:LDIR)//'/'//'SOS'	


C =========================================================
C --- Définition du fichier trace principal
C =========================================================        

       !Fichier Trace principal
       IF (FICMAIN_LOG.EQ."NO_LOG_FILE") THEN
	   TRACE=.FALSE.
       ELSE
          TRACE=.TRUE.	
       ENDIF ! Fin Trace 
 
 
      
C -------------------------------------------------------------------  
C --- Controle des paramètres propres aux calculs d'Ordres Successifs 
C     et à la production des fichiers résultats
C -------------------------------------------------------------------  
C  
C ---- Test d'existence et de validité des paramètres obligatoires
      
C     Option de visée : obligatoire  
      IF (ITRPHI.EQ.CTE_NOT_DEFINED_VALUE_INT) THEN
          GOTO 2605
      ELSE
	  
          !Si choix ITRPHI=1 : il faut définir l'azimut relatif
          IF (ITRPHI.EQ.1) THEN
	      IF (PHIOS.EQ.CTE_NOT_DEFINED_VALUE_DBLE) THEN
	         GOTO 2607
	      ELSE
	        !Pour le cas d'une définition du plan azimutal d'observation 
	        !on n'a pas besoin de définir PAS_PHI car il n'est pas utilisé par le code SOS, 
	        !mais l'argument est attendu en entrée SOS_TRPHI_OPTION
	        ! ==> on définit une valeur quelleconque (entière)
	        PAS_PHI=-9999
	     ENDIF
	  ENDIF
	  
	  !Si choix ITRPHI=2 : il faut définir le pas en azimut
	  IF (ITRPHI.EQ.2) THEN
	     IF (PAS_PHI.EQ.CTE_NOT_DEFINED_VALUE_INT) THEN
	        GOTO 2608
	     ELSE
	        !Pour le cas d'une définition du champ d'observation sur pas en azimut 
	        !on n'a pas besoin de définir PHIOS car il n'est pas utilisé par le code SOS, 
	        !mais l'argument est attendu en entrée SOS_TRPHI_OPTION
	        ! ==> on définit une valeur quelleconque
	        PHIOS=-9999.0
	     ENDIF
	  ENDIF
      ENDIF ! ITRPHI
 

C =========================================================
C --- Gestion de l'arborescence des fichiers de sortie
C =========================================================     
                              
C ---- Fichiers de sortie
         	     	
      !Fichier de luminances montantes
      LFIC = INDEX(FICUP,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3013
      FICUP = DIR_SOS(1:LDIRSOS)//'/'//FICUP(1:LFIC)

      IF(FICUP_ANGLES_USER.NE."NO_OUTPUT") THEN
	  LFIC = INDEX(FICUP_ANGLES_USER,' ') - 1
         IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3013
         FICUP_ANGLES_USER 
     &    = DIR_SOS(1:LDIRSOS)//'/'//FICUP_ANGLES_USER(1:LFIC)
      ENDIF

      !Fichier de luminances descendantes
      LFIC = INDEX(FICDOWN,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3014
      FICDOWN = DIR_SOS(1:LDIRSOS)//'/'//FICDOWN(1:LFIC)

      IF(FICDOWN_ANGLES_USER.NE."NO_OUTPUT") THEN
	  LFIC = INDEX(FICDOWN_ANGLES_USER,' ') - 1
         IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3013
         FICDOWN_ANGLES_USER 
     &    = DIR_SOS(1:LDIRSOS)//'/'//FICDOWN_ANGLES_USER(1:LFIC)
      ENDIF


C     Lancement de la subroutine SOS_PROC
      CALL SOS_PROC(RESROOT, FICMAIN_LOG, WA_SIMU, NBMU_GAUSS_LUM,
     &     FICANGLES_USER_LUM, TETAS, FICANGLES_RES_LUM,
     &     NBMU_GAUSS_MIE, FICANGLES_USER_MIE, FICANGLES_RES_MIE,
     &     FICANGLOG, WAREF_AOT, AOT_REF, ITRONC_AER, FICGRANU_LOG,
     &     FICMIE_LOG, DIR_MIE, FICGRANU, IMOD_AER, 
     &     RN_WA, IN_WA, RN_WAREF, IN_WAREF, 
     &     IGRANU, LND_RADIUS_MMD_AER, LND_LNVAR_MMD_AER,
     &     JD_SLOPE_MMD_AER, JD_RMIN_MMD_AER, JD_RMAX_MMD_AER,
     &     IModele_WMO, C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &     IModele_SF, RH, MODE_PARAM_BILND, USER_CV_COARSE,
     &     USER_CV_FINE, RTAUct_WAREF, BMD_CM_MRWA, BMD_CM_MIWA,
     &     BMD_CM_MRWAREF, BMD_CM_MIWAREF, BMD_CM_RMODAL, BMD_CM_VAR,
     &     BMD_FM_MRWA, BMD_FM_MIWA, BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &     BMD_FM_RMODAL, BMD_FM_VAR, FICEXTDATA_AER, FICMIXTURE_AER,
     &     FICUSER_AER, FICPROFIL_LOG, TR, HR, HA, IPROFIL, ZMIN, ZMAX,
     &     PSURF, H2O, O3, CO2, CH4, ABSPROFIL, FICABSPROFIL, NUSTEP,
     &     ISURF, DIR_SURF, FICSURF_LOG, SURF_IND, WIND, K0_ROUJEAN,
     &     K1_ROUJEAN, K2_ROUJEAN, ALPHA_NADAL, BETA_NADAL,
     &     COEF_C_MAIGNAN, RHO, FICSURF, FICSOS_LOG, FICSOS_RES_BIN,
     &     FICTRANS, FICFLUX, ZOUT, IGMAX, IPOLAR, ITRPHI, PHIOS,
     &     PAS_PHI, IMODE_CKD_CALCUL,
     &     IER, TRACE,
     &     LUM_NBMU,IND_ANGOUT_FIN, PHI_FIN, THETA_FIN, 
     &     SCA_UP_FIN, I_UP_FIN, Q_UP_FIN, U_UP_FIN,
     &     POL_ANG_UP_FIN, POL_RATE_UP_FIN, L_POL_UP_FIN,
     &     SCA_DOWN_FIN, I_DOWN_FIN, Q_DOWN_FIN, U_DOWN_FIN,
     &     POL_ANG_DOWN_FIN, POL_RATE_DOWN_FIN, L_POL_DOWN_FIN,
     &     TDIR_VRAI, FLUX_DIFF_DOWN, FLUX_DOWN, EPLUS,
     &     COEF_TRONCA)


C --- Gestion du cas d'erreur
      
      IF (IER.NE.0) THEN        ! cas d'erreur
         WRITE (6,*)
     &        "   ==> Error in subroutine SOS_PROC"
         GOTO 9998
      ENDIF

C --- Fichier Trace      
      IF (TRACE) THEN
         OPEN (IDLOG,FILE=FICMAIN_LOG,STATUS='OLD',
     &         ACCESS='APPEND',ERR=1302)
         WRITE (IDLOG,*,ERR=1303) " "

    	  WRITE (IDLOG,*,ERR=1303) "FICUP_ANGLES_USER ",
     &	  FICUP_ANGLES_USER
    	  WRITE (IDLOG,*,ERR=1303) "FICDOWN_ANGLES_USER ",
     &	  FICDOWN_ANGLES_USER
 	  WRITE (IDLOG,*,ERR=1303) "FICUP ",FICUP
    	  WRITE (IDLOG,*,ERR=1303) "FICDOWN ",FICDOWN
      ENDIF

C* Ouverture des fichiers 
C------------------------------------------------------
C     Fichiers résultats pour le champ montant et le champ descendant
C     Si le fichier existe déjà, il est écrasé par le nouveau fichier
      OPEN(UNIT=14,FILE=FICUP,FORM='FORMATTED',ERR=991)
      
      OPEN(UNIT=15,FILE=FICDOWN,FORM='FORMATTED',ERR=991)

C     Fichiers résultats  pour le champ montant et le champ descendant
C     Si le fichier existe déjà, il est écrasé par le nouveau fichier  
C     --> Cas des fichiers optionnels contenant les résultats pour les
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

         PHI1 = PHIOS
         PHI2 = PHIOS+180.D+00

         !Ecriture des informations des entêtes de fichier
         !------------------------------------------------
         IF (ZOUT.EQ.-1) THEN   

            CALL SOS_OUTPUT_HEADER(14, 1, PHI1, PHI2, CTE_TOA_ALT) 
            CALL SOS_OUTPUT_HEADER(15, 2, PHI1, PHI2, 0.D+00)

            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER(24, 1, PHI1, PHI2, CTE_TOA_ALT) 
            ENDIF

            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER(25, 2, PHI1, PHI2, 0.D+00) 
            ENDIF
         
         ELSE   
            CALL SOS_OUTPUT_HEADER(14, 1, PHI1, PHI2, ZOUT)
            CALL SOS_OUTPUT_HEADER(15, 2, PHI1, PHI2, ZOUT)

            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER(24, 1, PHI1, PHI2, ZOUT) 
            ENDIF
            
            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER(25, 2, PHI1, PHI2, ZOUT) 
            ENDIF

         ENDIF

         !Demi plan PHIos + 180 degrés : Angles de visée négatifs par convention
         !----------------------------------------------------------------------

         DO J=LUM_NBMU,1,-1
	 
            JJ=J-1
	    	    
            !Enregistrement pour le champ montant
            !------------------------------------
            XI = I_UP_FIN(0, JJ)
            XQ = Q_UP_FIN(0, JJ)
            XU = U_UP_FIN(0, JJ)
            XAN = POL_ANG_UP_FIN(0, JJ)
            TPOL = POL_RATE_UP_FIN(0, JJ)
            LPOL = L_POL_UP_FIN(0, JJ)
      
            !Cas nominal : champ montant pour tous les angles sont à sortir
            WRITE(14,55,err=992) 
     &         -THETA_FIN(JJ),SCA_UP_FIN(0, JJ),XI,XQ,XU,XAN,TPOL,LPOL

            !Cas optionnel : champ montant pour les angles utilisateurs uniquement
            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
                  WRITE(24,55,err=992) 
     &                -THETA_FIN(JJ),SCA_UP_FIN(0,JJ),XI,XQ,XU,XAN,
     &                 TPOL,LPOL
               ENDIF
            ENDIF  

           
            !Enregistrement pour le champ descendant
            !---------------------------------------
            XI = I_DOWN_FIN(0, JJ)
            XQ = Q_DOWN_FIN(0, JJ)
            XU = U_DOWN_FIN(0, JJ)
            XAN = POL_ANG_DOWN_FIN(0, JJ)
            TPOL = POL_RATE_DOWN_FIN(0, JJ)
            LPOL = L_POL_DOWN_FIN(0, JJ)

            !Cas nominal : champ descendant pour tous les angles 
            WRITE(15,55,err=992)  
     &           -THETA_FIN(JJ),SCA_DOWN_FIN(0, JJ),XI,XQ,XU,XAN,
     &           TPOL,LPOL   	     

            !Cas optionnel : champ descendant pour les angles utilisateurs uniquement
            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
                  WRITE(25,55,err=992) 
     &                 -THETA_FIN(J-1),SCA_DOWN_FIN(0, JJ),XI,XQ,XU,
     &                 XAN,TPOL,LPOL
               ENDIF
            ENDIF 	    
         ENDDO

         DO J=1,LUM_NBMU
	 
            JJ = J -1

            !Enregistrement pour le champ montant
            !------------------------------------
            XI = I_UP_FIN(1, JJ)
            XQ = Q_UP_FIN(1, JJ)
            XU = U_UP_FIN(1, JJ)
            XAN = POL_ANG_UP_FIN(1, JJ)
            TPOL = POL_RATE_UP_FIN(1, JJ)
            LPOL = L_POL_UP_FIN(1, JJ)

            !Cas nominal : champ montant pour tous les angles 
            WRITE(14,55,err=992) 
     &            THETA_FIN(JJ),SCA_UP_FIN(1, JJ),XI,XQ,XU,XAN,
     &            TPOL,LPOL
     
            !Cas optionnel : champ montant pour les angles utilisateurs uniquement
            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
                IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
	            WRITE(24,55,err=992) 
     &                  THETA_FIN(JJ),SCA_UP_FIN(1, JJ),XI,XQ,XU,
     &                  XAN,TPOL,LPOL
                ENDIF
	     ENDIF  

          
            !Enregistrement pour le champ descendant
            !---------------------------------------
            XI = I_DOWN_FIN(1, JJ)
            XQ = Q_DOWN_FIN(1, JJ)
            XU = U_DOWN_FIN(1, JJ)
            XAN = POL_ANG_DOWN_FIN(1, JJ)
            TPOL = POL_RATE_DOWN_FIN(1, JJ)
            LPOL = L_POL_DOWN_FIN(1, JJ)
	    
            !Cas nominal : champ descendant pour tous les angles sont à sortir
	     WRITE(15,55,err=992)  
     &            THETA_FIN(JJ),SCA_DOWN_FIN(1, JJ),XI,XQ,XU,XAN,
     &            TPOL,LPOL   	     

            !Cas optionnel : champ descendant pour les angles utilisateurs uniquement
            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
                IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
	            WRITE(25,55,err=992) 
     &                  THETA_FIN(JJ),SCA_DOWN_FIN(1, JJ),XI,XQ,XU,XAN,
     &                  TPOL,LPOL
                ENDIF
	     ENDIF 	    

         ENDDO	    
	 
	    
      ENDIF  !ITRPHI=1
      


C* Sortie d'un diagramme polaire
C---------------------------------
      IF(ITRPHI.EQ.2) THEN        
	 
         !Ecriture des informations des entêtes de fichier
         !------------------------------------------------
         IF (ZOUT.EQ.-1) THEN   

            CALL SOS_OUTPUT_HEADER_POLAR_DIAG(14, 1, CTE_TOA_ALT) 
            CALL SOS_OUTPUT_HEADER_POLAR_DIAG(15, 2, 0.D+00)

            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER_POLAR_DIAG(24, 1, CTE_TOA_ALT) 
            ENDIF

            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER_POLAR_DIAG(25, 2, 0.D+00) 
            ENDIF
         
         ELSE   
            CALL SOS_OUTPUT_HEADER_POLAR_DIAG(14, 1, ZOUT)
            CALL SOS_OUTPUT_HEADER_POLAR_DIAG(15, 2, ZOUT)

            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER_POLAR_DIAG(24, 1,  ZOUT) 
            ENDIF

            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER_POLAR_DIAG(25, 2, ZOUT) 
            ENDIF

         ENDIF

         IP=0
         DO IPHI=0,360,PAS_PHI

            DO J=1,LUM_NBMU

               JJ=J-1
	       
               !Enregistrement pour le champ montant
               !------------------------------------
               XI = I_UP_FIN(IP, JJ)
               XQ = Q_UP_FIN(IP, JJ)
               XU = U_UP_FIN(IP, JJ)
               XAN = POL_ANG_UP_FIN(IP, JJ)
               TPOL = POL_RATE_UP_FIN(IP, JJ)
               LPOL = L_POL_UP_FIN(IP, JJ)

               !Cas nominal : champ montant pour tous les angles   
               WRITE(14,56,err=992) 
     &              PHI_FIN(IP),THETA_FIN(JJ),SCA_UP_FIN(IPHI, JJ),
     &              XI,XQ,XU,XAN,TPOL, LPOL
     
               !Cas optionnel : champ montant pour les angles utilisateurs uniquement
               IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
                  IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
                     WRITE(24,56,err=992) 
     &                    PHI_FIN(IP),THETA_FIN(JJ),
     &                    SCA_UP_FIN(IP, JJ),XI,XQ,XU, XAN,TPOL,LPOL
                     
                  ENDIF
               ENDIF  

               !Cas nominal : champ descendant pour tous les angles  
               XI = I_DOWN_FIN(IP, JJ)
               XQ = Q_DOWN_FIN(IP, JJ)
               XU = U_DOWN_FIN(IP, JJ)
               XAN = POL_ANG_DOWN_FIN(IP, JJ)
               TPOL = POL_RATE_DOWN_FIN(IP, JJ)
               LPOL = L_POL_DOWN_FIN(IP, JJ) 
               WRITE(15,56,err=992) 
     &              PHI_FIN(IP),THETA_FIN(JJ),SCA_DOWN_FIN(IP, JJ),
     &              XI,XQ,XU,XAN,TPOL,LPOL

	 
               !Cas optionnel : champ descendant pour les angles utilisateurs uniquement
                IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
                  IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
                      WRITE(25,56,err=992) 
     &                    PHI_FIN(IP),THETA_FIN(JJ),
     &                    SCA_DOWN_FIN(IP, JJ),XI,XQ,XU,XAN,TPOL,LPOL
                    ENDIF
	        ENDIF  

      
            ENDDO !MU
	     
	     IP=IP+1
	     
	  ENDDO !PHI

      ENDIF   !ITRPHI=2


      
C* Fermeture fichier
C-------------------
      CLOSE(14)
      CLOSE(15)
      CLOSE(21)

      IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') CLOSE(24)
      IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') CLOSE(25)	 
      
C* Fin nominale 
C-------------------
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &   "--> Nominal end of SOS_ABS_MAIN processes"
         WRITE (IDLOG,*,ERR=1303) "JOB_STATUS=OK"
         CLOSE(IDLOG)
      ENDIF
      
      GOTO 9999
      

C =========================================================
C --- Cas d'erreur et retour du status 1 au shell  
C =========================================================  

C --- Format des inputs

 100  WRITE(6,*) '  SOS_MAIN : ERROR_100 on argument loading'
      WRITE(6,*) '  --> Check and correct the keyword : ', KEYWORD
      GOTO 9998  

C --- Paramètres généraux

 101  WRITE(6,*) '  SOS_MAIN : ERROR_101 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS_Main.Wa value : Wavelength for '
      WRITE(6,*) '       radiance simulation (microns) '
      GOTO 9998

C --- Paramètres pour les calculs d'angles

 201  WRITE(6,*) '  SOS_MAIN : ERROR_201 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -ANG.Rad.NbGauss value : Number of Gauss'
      WRITE(6,*) '       angles for radiance computations (I4).' 
      GOTO 9998    

 202  WRITE(6,*) '  SOS_MAIN : ERROR_202 on argument loading'
      WRITE(6,*) '  --> Check required formats for -ANG.Thetas value :'
      WRITE(6,*) '      Solar zenithal angle (degrees).'
      GOTO 9998  

 203  WRITE(6,*) '  SOS_MAIN : ERROR_203 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -ANG.Aer.NbGauss value : Number of Gauss'
      WRITE(6,*) '       angles for Mie computations (I4).' 
      
C --- Paramètres aérosols

 301  WRITE(6,*) '  SOS_MAIN : ERROR_301 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.Waref : Reference wavelength for the '
      WRITE(6,*) '       aerosol optical thickness .' 
      GOTO 9998       

 302  WRITE(6,*) '  SOS_MAIN : ERROR_302 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.AOTref : Aerosol optical thickness at '
      WRITE(6,*) '       reference wavelength (F9.5).' 
      GOTO 9998   
     
 303  WRITE(6,*) '  SOS_MAIN : ERROR_303 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.Tronca : Option of aerosols phase '
      WRITE(6,*) '       function truncature (I4).' 
      GOTO 9998   

 304  WRITE(6,*) '  SOS_MAIN : ERROR_304 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.Model : Type of aerosols model (I2).' 
      GOTO 9998                 

 306  WRITE(6,*) '  SOS_MAIN : ERROR_306 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.MRwa : Mono-modal aerosol model -'
      WRITE(6,*) '       Real part of the aerosol refractive index '
      WRITE(6,*) '       for the simulation wavelength (F5).' 
      GOTO 9998            

 307  WRITE(6,*) '  SOS_MAIN : ERROR_307 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.MIwa : Mono-modal aerosol model -'
      WRITE(6,*) '      Imaginary part of the aerosol refractive index'
      WRITE(6,*) '      for the simulation wavelength (F8).' 
      GOTO 9998             

 308  WRITE(6,*) '  SOS_MAIN : ERROR_308 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.MRwaref : Mono-modal aerosol model -'
      WRITE(6,*) '       Real part of the aerosol refractive index '
      WRITE(6,*) '       for the aerosol reference wavelength (F5).' 
      GOTO 9998            

 309  WRITE(6,*) '  SOS_MAIN : ERROR_309 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.MIwaref : Mono-modal aerosol model -'
      WRITE(6,*) '      Imaginary part of the aerosol refractive index'
      WRITE(6,*) '      for the aerosol reference wavelength (F8).' 
      GOTO 9998             

 310  WRITE(6,*) '  SOS_MAIN : ERROR_310 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.SDtype (I4) : Mono-modal aerosol model'
      WRITE(6,*) '       Model index : 1 LND '
      WRITE(6,*) '                     2 Junge Law '
      GOTO 9998            

 311  WRITE(6,*) '  SOS_MAIN : ERROR_311 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.LNDradius : Mono-modal aerosol model -'
      WRITE(6,*) '       Modal radius (microns) of the LND (Float)'
      GOTO 9998           
            
 312  WRITE(6,*) '  SOS_MAIN : ERROR_312 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.LNDvar : Mono-modal aerosol model -'
      WRITE(6,*) '       Variance of the LND (Float)'
      GOTO 9998           

 3121 WRITE(6,*) '  SOS_MAIN : ERROR_3121 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.JD.slope : Aerosol Junge model - '
      WRITE(6,*) '       Slope (Float)'
      GOTO 9998  

 3122 WRITE(6,*) '  SOS_MAIN : ERROR_3122 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.JD.rmin : Aerosol Junge model - '
      WRITE(6,*) '       Minimal radius in microns (Float)'
      GOTO 9998  

 3123 WRITE(6,*) '  SOS_MAIN : ERROR_3123 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.JD.rmax : Aerosol Junge model - '
      WRITE(6,*) '       Maximal radius in microns (Float)'
      GOTO 9998  
	              
 313  WRITE(6,*) '  SOS_MAIN : ERROR_313 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.WMO.Model : Type of WMO Model (I2)'
      GOTO 9998       	

 314  WRITE(6,*) '  SOS_MAIN : ERROR_314 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      Volumetric concentration of WMO user model :'
      WRITE(6,*) '      -AER.WMO.DL :  Dust-Like component (F9.5)'
      WRITE(6,*) '      -AER.WMO.WS :  Water Soluble component (F9.5)'
      WRITE(6,*) '      -AER.WMO.OC :  Oceanic component (F9.5)'
      WRITE(6,*) '      -AER.WMO.SO :  Soot component (F9.5)'
      GOTO 9998      	 

 315  WRITE(6,*) '  SOS_MAIN : ERROR_315 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.SF.Model : Type of Shettle & Fenn Model'
      WRITE(6,*) '       (I2)'
      GOTO 9998       	

 316  WRITE(6,*) '  SOS_MAIN : ERROR_316 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.SF.RH : % of humidity relative (F5.2)'
      GOTO 9998       	

 317  WRITE(6,*) '  SOS_MAIN : ERROR_317 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.VCdef : Type of selected method for'
      WRITE(6,*) '       bi-modal LND volumetric concentration'
      WRITE(6,*) '       definition (I2)'
      GOTO 9998       	
	
 318  WRITE(6,*) '  SOS_MAIN : ERROR_318 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CoarseVC : User volumetric '
      WRITE(6,*) '       concentration for "LND coarse mode" (F8.5)'
      GOTO 9998       	 
	
 319  WRITE(6,*) '  SOS_MAIN : ERROR_319 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FineVC : User volumetric '
      WRITE(6,*) '       concentration for "LND fine mode" (F8.5)'
      GOTO 9998       

 320  WRITE(6,*) '  SOS_MAIN : ERROR_320 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.RAOT : User value of the ration '
      WRITE(6,*) '       AOT_coarse / AOT_total '
      WRITE(6,*) '       for the aerosol reference wavelength (F9.5)'
      GOTO 9998       
      
 321  WRITE(6,*) '  SOS_MAIN : ERROR_321 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.MRwa : '
      WRITE(6,*) '       Real part of the refractive '
      WRITE(6,*) '       index for the "LND coarse mode" '
      WRITE(6,*) '       at the wavelength of aerosol properties '
      WRITE(6,*) '       calculation (F5.3)'
      GOTO 9998       	 		

 322  WRITE(6,*) '  SOS_MAIN : ERROR_322 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.MIwa : '
      WRITE(6,*) '       Imaginary part of the refractive '
      WRITE(6,*) '       index for the "LND coarse mode" '
      WRITE(6,*) '       at the wavelength of aerosol properties '
      WRITE(6,*) '       calculation (F8.5)'
      GOTO 9998       
	 
 323  WRITE(6,*) '  SOS_MAIN : ERROR_323 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.MRwaref : '
      WRITE(6,*) '       Real part of the refractive '
      WRITE(6,*) '       index for the "LND coarse mode" '
      WRITE(6,*) '       at the aerosol reference wavelength'
      WRITE(6,*) '       calculation (F5.3)'
      GOTO 9998       	 		

 324  WRITE(6,*) '  SOS_MAIN : ERROR_324 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.MIwaref : '
      WRITE(6,*) '       Imaginary part of the refractive '
      WRITE(6,*) '       index for the "LND coarse mode" '
      WRITE(6,*) '       at the aerosol reference wavelength'
      WRITE(6,*) '       calculation (F8.5)'
      GOTO 9998  

 325  WRITE(6,*) '  SOS_MAIN : ERROR_325 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.SDradius : '
      WRITE(6,*) '       Modal radius of the "LND coarse mode"'
      WRITE(6,*) '       in microns (F9.5)'
      GOTO 9998       	 		
           
 326  WRITE(6,*) '  SOS_MAIN : ERROR_326 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.SDvar : '
      WRITE(6,*) '       "LND coarse mode" variance (F9.5)'
      GOTO 9998       	 			

 327  WRITE(6,*) '  SOS_MAIN : ERROR_327 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.MRwa : '
      WRITE(6,*) '       Real part of the refractive '
      WRITE(6,*) '       index for the "LND fine mode" '
      WRITE(6,*) '       at the wavelength of aerosol properties'
      WRITE(6,*) '       calculation (F5.3)'
      GOTO 9998       	 		

 328  WRITE(6,*) '  SOS_MAIN : ERROR_328 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.MIwa : '
      WRITE(6,*) '       Imaginary part of the refractive '
      WRITE(6,*) '       index for the "LND fine mode" '
      WRITE(6,*) '       at the wavelength of aerosol properties'
      WRITE(6,*) '       calculation (F8.5)'
      GOTO 9998       
	 
 329  WRITE(6,*) '  SOS_MAIN : ERROR_329 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.MRwaref : '
      WRITE(6,*) '       Real part of the refractive '
      WRITE(6,*) '       index for the "LND fine mode" '
      WRITE(6,*) '       at the aerosol reference wavelength'
      WRITE(6,*) '       calculation (F5.3)'
      GOTO 9998       	 		

 330  WRITE(6,*) '  SOS_MAIN : ERROR_330 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.MIwaref : '
      WRITE(6,*) '       Imaginary part of the refractive '
      WRITE(6,*) '       index for the "LND fine mode" '
      WRITE(6,*) '       at the aerosol reference wavelength'
      WRITE(6,*) '       calculation (F8.5)'
      GOTO 9998  

 331  WRITE(6,*) '  SOS_MAIN : ERROR_331 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.SDradius : '
      WRITE(6,*) '       Modal radius of the "LND fine mode"'
      WRITE(6,*) '       in microns (F9.5)'
      GOTO 9998       	 		
           
 332  WRITE(6,*) '  SOS_MAIN : ERROR_332 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.SDvar : '
      WRITE(6,*) '       "LND fine mode" variance (F9.5)'
      GOTO 9998       	 				

C --- Paramètres de profil

 401  WRITE(6,*) '  SOS_MAIN : ERROR_401 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.MOT : Molecular atmos optical thickness '
      WRITE(6,*) '    (for the radiance simulation wavelength) (F9.5).'
      GOTO 9998   

 402  WRITE(6,*) '  SOS_MAIN : ERROR_402 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.HR : Height scale of the molecular '
      WRITE(6,*) '        atmospheric profile (F9).' 
      GOTO 9998   

 403  WRITE(6,*) '  SOS_MAIN : ERROR_403 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AerHS.HA : Height scale of the aerosol '
      WRITE(6,*) '       atmospheric profile (F9).' 
      GOTO 9998   

 404  WRITE(6,*) '  SOS_MAIN : ERROR_404 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AerProfile.Type : type of '
      WRITE(6,*) '        atmospheric profile (I2).' 
      GOTO 9998   

 405  WRITE(6,*) '  SOS_MAIN : ERROR_405 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AerLayer.Zmin : Minimal altitude (km)'
      WRITE(6,*) '        of aerosol layer (F9.5).' 
      GOTO 9998   

 406  WRITE(6,*) '  SOS_MAIN : ERROR_406 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AerLayer.Zmax : Maximal altitude (km)'
      WRITE(6,*) '       of aerosol layer (F9.5).' 
      GOTO 9998   

 407  WRITE(6,*) '  SOS_MAIN : ERROR_407 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.Psurf : Atmospheric pressure at sea '
      WRITE(6,*) '       level (mbar) (F9.5).' 
      GOTO 9998   

 408  WRITE(6,*) '  SOS_MAIN : ERROR_408 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.H2O : Water vapor content (g/cm2) (F9) '
      GOTO 9998   

 409  WRITE(6,*) '  SOS_MAIN : ERROR_409 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.O3 : Ozone content (Dobson) (F9) '
      GOTO 9998   

 4091 WRITE(6,*) '  SOS_MAIN : ERROR_4091 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.CO2 : Carbon dioxyde concentration'
      WRITE(6,*) '                at surface level (ppmv) (F9)'
      GOTO 9998   

 4092 WRITE(6,*) '  SOS_MAIN : ERROR_4092 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.CH4 : Methane concentration'
      WRITE(6,*) '                at surface level (ppmv) (F9)'
      GOTO 9998   

 410  WRITE(6,*) '  SOS_MAIN : ERROR_410 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AbsProfile.Type : ',
     &		 'Absorption profile type (I2)'
      GOTO 9998   

 411  WRITE(6,*) '  SOS_MAIN : ERROR_411 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.SpectralResol : ',
     &		 'Spectral resolution in cm-1 of CKD datafiles (Float)'
      GOTO 9998   

C --- Paramètres de surface

 501  WRITE(6,*) '  SOS_MAIN : ERROR_501 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Type : Type of surface model (I2) '
      GOTO 9998   

 502  WRITE(6,*) '  SOS_MAIN : ERROR_502 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Ind : Refraction indice (F5.3) '
      GOTO 9998   

 503  WRITE(6,*) '  SOS_MAIN : ERROR_503 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Glitter.Wind : Wind speed (m/s) (F4) '
      GOTO 9998   

 504  WRITE(6,*) '  SOS_MAIN : ERROR_504 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Roujean.K0 : Roujean K0 parameter (F7.3)'
      GOTO 9998   

 505  WRITE(6,*) '  SOS_MAIN : ERROR_505 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Roujean.K1 : Roujean K1 parameter (F7.3)'
      GOTO 9998   

 506  WRITE(6,*) '  SOS_MAIN : ERROR_506 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Roujean.K2 : Roujean K2 parameter (F7.3)'
      GOTO 9998   

 507  WRITE(6,*) '  SOS_MAIN : ERROR_507 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Nadal.Alpha : Nadal Alpha parameter (F6)'
      GOTO 9998   
 
 508  WRITE(6,*) '  SOS_MAIN : ERROR_508 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Nadal.Beta : Nadal Beta parameter (F6) '
      GOTO 9998   
 
 5081 WRITE(6,*) '  SOS_MAIN : ERROR_5081 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Maignan.C : C*exp(-NDVI) parameter (F6) '
      WRITE(6,*) "      of the Maignan's surface BPDF model "
      GOTO 9998 

 509  WRITE(6,*) '  SOS_MAIN : ERROR_509 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Alb : Surface Albedo'
      GOTO 9998   

C --- Paramètres SOS

 601  WRITE(6,*) '  SOS_MAIN : ERROR_601 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.OutputAlt : Output Altitude in km (F9) '
      GOTO 9998   

 602  WRITE(6,*) '  SOS_MAIN : ERROR_602 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.IGmax : maximal order of interaction (I4)'
      GOTO 9998	 

 603  WRITE(6,*) '  SOS_MAIN : ERROR_603 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '   -SOS.Ipolar : Polarization activation option (I2)'
      GOTO 9998	 

 604  WRITE(6,*) '  SOS_MAIN : ERROR_604 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.View : Output Type (I4)'
      GOTO 9998	 

 605  WRITE(6,*) '  SOS_MAIN : ERROR_605 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.View.Phi : Relative azimut angle (F9)'
      GOTO 9998	 

 606  WRITE(6,*) '  SOS_MAIN : ERROR_606 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.View.Dphi : Azimut sampling step (I4)'
      GOTO 9998	 

 608  WRITE(6,*) '  SOS_MAIN : ERROR_608 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.AbsModeCKD : Selected ',
     & 		 ' mode for application of CKD method factor (I4)'
      GOTO 9998	 

C --- Existence des paramètres d'arborescence

 1000 WRITE(6,*) '  SOS_MAIN : ERROR_1000 on parameters'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS_Main.ResRoot Value ' 
      WRITE(6,*) '  ==>  Working folder of SOS computations'
      GOTO 9998 

C --- Définition Arborescence de sortie

 1100 WRITE(6,*) '  SOS_MAIN : ERROR_1100 : '
      WRITE(6,*) '  => Working folder name : maybe too long'
      WRITE(6,*) '     -SOS_Main.ResRoot parameter'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	             	             	                 	              	

 1108 WRITE(6,*) '  SOS_MAIN : ERROR_1108 : '
      WRITE(6,*) '  => SOS Working folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	               	  

 1109 WRITE(6,*) '  SOS_MAIN : ERROR_1109 : '
      WRITE(6,*) '  => Main log folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998     

C --- gestion du log
                                           
 1301 WRITE(6,*) '  SOS_MAIN : ERROR_1301 : '
      WRITE(6,*) '  => Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 1302 WRITE(6,*) '  SOS_MAIN : ERROR_1302 on Log file opening'
      GOTO 9998

 1303 WRITE(6,*) '  SOS_MAIN : ERROR_1303 on Log file writing'
      GOTO 9998
            

C --- Controle des paramètres de profil

 2515 WRITE(6,*) '  SOS_MAIN : ERROR_2515 on parameters'
      WRITE(6,*) '  --> -SOS.AbsModeCKD value must be defined'
      WRITE(6,*) '      1 : Fine calculation of the impact of gaz'
      WRITE(6,*) '          absorption on radiance by the CKD method.'
      WRITE(6,*) '      2 : Simplified calculation of the impact of'
      WRITE(6,*) '           gaz absorption radiance.'
      GOTO 9998  
      

C --- Controle des paramètres SOS

 2605 WRITE(6,*) '  SOS_MAIN : ERROR_2605 on parameters'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.View  (I4) '
      WRITE(6,*) '  ==> Azimut mode  '
      GOTO 9998     

 2607 WRITE(6,*) '  SOS_MAIN : ERROR_2607 on parameters'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.View.Phi  (F9) '
      WRITE(6,*) '  ==> Relative azimut (required when SOS.View=1) '      
      GOTO 9998     
 
 2608 WRITE(6,*) '  SOS_MAIN : ERROR_2608 on parameters'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.View.Dphi  (Integer) '
      WRITE(6,*) '  ==> Azimut sampling step'
      WRITE(6,*) '      (required when SOS.View=2)'    
      GOTO 9998     

C --- Définition des fichiers de sortie

 3013 WRITE(6,*) '  SOS_MAIN : ERROR_3013 : '
      WRITE(6,*) '  => Upward radiance filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3014 WRITE(6,*) '  SOS_MAIN : ERROR_3014 : '
      WRITE(6,*) '  => Downward radiance filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	
     
      
C* Cas d'erreur : retour du status 1 au programme appelant
C------------------------------------------------------
  991 WRITE(6,*) "SOS_TRPHI : Error while opening a file"
      GOTO 9998

  992 WRITE(6,*) 
     &  "SOS_TRPHI : Error while reading or writing on a file"
      GOTO 9998

                    
 9998 IF (TRACE) THEN
         !OPEN (IDLOG,FILE=FICMAIN_LOG,STATUS='OLD',ERR=1302)
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) "--> Processing error"
	  WRITE (IDLOG,*,ERR=1303) "JOB_STATUS=ERROR"
	  CLOSE(IDLOG)
      ENDIF 

      CALL FLUSH(6)        !Force l'envoie du buffer à l'écran

      CALL EXIT(1)    
      
 9999 STOP



C =========================================================
C --- Formats
C =========================================================            
C
     
C --- Formats pour les fichiers SOS_Up et SOS_Down
55    FORMAT(2(2X,F7.2),2X,3(E13.6,2X),2(F7.2,2X),E13.6)
56    FORMAT(3(2X,F7.2),1X,3(E13.6,2X),2(1X,F7.2),E13.6)

     
      END   	!FIN DU PROGRAMME SOS_MAIN
