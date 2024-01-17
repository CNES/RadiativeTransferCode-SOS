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
C* FICHIER: SOS_AEROSOLS.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module gère le calcul des proprietes radiatives de 
C*       différentes granulométrieF d'aérosols.
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et développement : CS GROUP France, B.Lafrance, X.Lenot
C*
C* DATE: 30/04/2002
C*
C* MOD:VERSION:1.1: Mise a niveau du code des OS du LOA
C* MOD:VERSION:1.2: Mise à jour des commentaires.
C*
C* MOD:VERSION:2.0: 24/08/2008 
C*       - Introduction d'une modelisation bi-LND.
C*       - Conversion des sorties KMAT1 et KMAT2 de la routine SOS_GRANU
C*         en sections efficaces d'extinction et de diffusion (micron^2)
C*         au lieu de coefficients d'extinction et de diffusion (/mic)
C*         qui n'étaient pas normalisés à une particule.
C*       - Mise à jour des commentaires.
C*
C*
C* MOD:VERSION:2.1: 08/09/2008 
C*       - Correction de dépassements de zone d'indentation pour compilation f77.
C*       - Suppression de la sortie de boucle pour fin de fichier de MIE atteinte.
C*         Car si c'est le cas, il s'agit d'une erreur : le fichier de MIE est 
C*         incomplet (probable manque d'espace disque à sa génération).
C*         Ce cas est alors correctement géré en erreur.
C*       - Correction pour le calcul des proportions de composants des modèles bimodaux
C*         Dubovik : cas du calcul à partir des rapports d'épaisseur optique des deux modes
C*         (fin et grossier) et de l'épaisseur optique totale pour une longueur
C*         d'onde de référence ==> utilisation de l'indice de réfraction des 
C*         particules pour la longueur d'onde de référence (uniquement).
C*
C* MOD:VERSION:3.0: 19/11/2009 
C*       - Modification du passage des arguments : 
C*           - passage par couples (mot clef, valeur). 
C*           - adaptation de la gestion des cas d'erreurs. 
C*       - Conversion des messages d'erreurs en anglais. 
C*       - Introduction du paramètre CTE_MIE_ALPHAMIN qui précise la valeur minimale
C*         du paramètre de taille ALPHAO des calculs de Mie.
C*
C* MOD:VERSION:4.0: 10/01/2010 
C*   - Suppression de la constante CTE_OS_FICANGLE (de SOS.h) 
C*     au profit d'un passage par argument du nom du fichier des angles (produit par SOS_ANGLES).
C*
C*   - Modification de la taille des tableaux de luminance :
C*     ==> Passage de CTE_MIE_NBMU à CTE_MIE_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'à MIE_NBMU égal au nombre effectif d'angles utilises.
C*                                   
C*   - Modification de la taille des tableaux des fonctions de phase :
C*     ==> Passage de CTE_MIE_NBMU à CTE_MIE_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'à MIE_NBMU égal au nombre effectif d'angles utilises.
C*
C*   - Modification de la taille des tableaux de decomposition en polynomes de Legendre :
C*     ==> Passage de CTE_OS_NB à CTE_OS_NB_MAX.
C*     ==> Boucle sur les ordre jusqu'à OS_NB (lu dans le fichier des angles)
C*
C*   - CTE_MIE : Ajout des arguments d'entrée : NBMU_MIE (nb angles utiles) et XMU (tableau des angles)
C*
C*   - SOS_NOMFIC_MIE : Ajout des arguments d'entrée : NBMU_GAUSS (nb angles de Gauss) 
C*                                                     et FICANGLES_USER (fichier d'angles utilisateurs)
C*
C*   - SOS_GRANU : 
C*         * Ajout de l'argument d'entrée : NBMU_MIE (nb angles utiles)
C*         * Modification de la ligne de lecture  des tableaux 
C*           de fonctions de phase par valeur de paramètre de taille Alpha
C*           ==> limitation aux bornes utiles -MIE_NBMU:MIE_NBMU
C*
C*   - SOS_DECOMPO_LEGENDRE : 
C*         * Ajout des arguments d'entrée : NBMU_MIE (nb angles utiles) et OS_NB
C*         * Boucles limitées à NBMU_MIE au lieu CTE_MIE_NBMU_MAX
C*  
C*   - Suppression de l'utilisation de la variable d'environnement SOS_RACINE_MU
C*     ==> devient SOS_ABS_FIC car utilisee seulement pour les fichiers WMO et S&F.
C*
C*   - Ajout d'une option de traitement IMOD=4 pour ingérer une matrice de phase externe
C*     --> Fonctions de phase externes et coefficients de diffusion et d'extinction
C*     ==> Introduction des routines : 
C*            	- SOS_INTERPO_SPLINT
C*     		- SOS_SPLINE
C*    		- SOS_SPLINT
C*     ==> Ajout de la constante (SOS.h) :  CTE_MAXNB_ANG_EXT
C*
C* MOD:VERSION:5.0: 06/06/2013
C*     
C*   -  Les paramètres d'entrée sont fournis sous la forme de liste de paramètres
C*      et non plus sous la forme de couples avec mot cle associe => tous les
C*	paramètres deviennent donc obligatoires
C*
C*   -  Modification de la valeur de Pi, declaree en constante dans le fichier
C*      (INCTE_PI remplace CTE_PI)
C* 
C* MOD:VERSION:6.0: 02/03/2016
C*
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*      (incluant au décodage des constante) : requis pour compilation gfortran 
C*
C*    - Changement de la manière de lire les fichiers formatés 
C*      (Angles, fichier de données externes pour les aérosols)
C*      pour être compatible avec le compilateur gfortran.
C*
C*    - Modification du format du fichier GRANU résultat pour permettre sa lecture
C*      par le code SOS_PREPA_OS.F compilé avec gfortran et pour clarifier son contenu
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
C* MOD:VERSION:6.1: 12/06/2020
C* 
C*    - Ajout du paramètre DIR_TMP qui définit le répertoire temporaire des 
C*      calculs de Mie.
C*
C*    - Définition 'un fichier de Mie temporaire FICMIE_TMP qui est généré 
C*      sous DIR_TMP avant d'être déplacé sous le répertoire de stockage
C*      des fichiers de Mie (DIRS) avec le nom formalisé du fichier (FICMIE).
C*
C*    - Ajout de l'argument VARGRANU3_MMD qui donne le rayon max d'une
C*      granulomètrie suivant la loi de Junge.
C*      On renomme VARGRANU1 en VARGRANU1_MMD et VARGRANU2 en VARGRANU2_MMD
C*      ==> Adaptation de l'appel de CTE_GRANU
C*
C*    - Le paramètre de taille maximal des calculs de Mie (ALPHAF)  
C*      pour un modèle mono-modal n'est plus fourni mais est calculé par la routine.
C*
C*    - Ajout d'une option de traitement IMOD=5 pour utiliser la définition 
C*      d'un mélange de modes fournis dans un fichier par l'utilisateur.
C*      ==> Introduction de paramètre: FICMIXTURE_AER qui donne la location du fichier.
C*      ==> Implémentation du bloc de traitement pour IMOD=5.
C*
C*
C* MOD:VERSION:6.2: 31/08/2022
C*    - Suppression du passage par la variable d'environnement SOS_ABS_FIC.
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Conversion des traces en Anglais.
C*    - Introduire des constantes CTE_AER_MU1_TRONCA et CTE_AER_MU2_TRONCA 
C*      pour définir les angles permettant d'estimer la troncature de la pointe
C*      de diffusion avant de la fonction de phase des aérosols.
C*    - The information about the units for extinction and scattering coefficients,
C*      and for the number of particles, are clarified in the log file.
C*
C* MOD:VERSION:6.3: 08/08/2023
C*    - Ajout d'un test sur la somme des taux d'épaisseur optique définissant un mélange d'aérosols
C*       ==> Normalisation appliquée si l'écart à 1 est inférieur au seuil CTE_GAP_TOLER_SUM_RATES
C*    - Clarification du nom de la variable donnant le coefficient de troncature (COEF_TRONCA)
C*      et passage de cette variable en paramètre de sortie de la routine SOS_AEROSOLS.
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    CTE_LENDIR : Longueur des noms de répertoires.
C    CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_LENKEYWORD : Taille maximale des Keywords pour le passage des arguments.
C    CTE_LENCOM: longueur des chaines de commande système
C    INCTE_PI : Valeur de PI.
C    CTE_MIE_NBMU_MAX : Nombre maximal d'angles positifs.
C    CTE_MIE_ALPHAMIN : Valeur minimale par defaut du paramètre de taille 
C                       Alpha des calculs de Mie.
C    CTE_ALPHAMAX_WMO_DL : Valeur limite du paramètre de taille pour le calcul
C                          des fichiers du composant Dust-Like de la WMO.
C    CTE_ALPHAMAX_WMO_WS : Idem pour le composant Water Soluble.
C    CTE_ALPHAMAX_WMO_OC : Idem pour le composant Oceanic.
C    CTE_ALPHAMAX_WMO_SO : Idem pour le composant Soot.
C    CTE_ALPHAMAX_SF_SR : Idem pour le composant Small Rural de Shettle & Fenn.
C    CTE_ALPHAMAX_SF_SU : Idem pour le composant Small Urban de Shettle & Fenn.
C    CTE_COEF_NRMAX : Constante pour la determination du paramètre de taille max.
C    CTE_OS_NB_MAX :  Valeur maximale pour l'ordre limite du développement en 
C                     polynomes de Legendre.
C    CTE_PH_SEUIL_TRONCA : Valeur seuil pour la troncature.
C    CTE_AER_DATAWMO : Nom du fichier contenant les paramètres WMO (rm,sig,...).
C    CTE_AER_DATASF  : Nom du fichier contenant rayon modal et log variance S&F.
C    CTE_AER_XX_SF   : Nom des fichiers contenant les indices de réfraction S&F
C		       pour chaque constituant XX = SR, LR, SU, LU, OM.
C    CTE_MAXNB_ANG_EXT : Taille maximale des tableaux de fonctions de phase externesC
C    CTE_AER_MU1_TRONCA : cos de l'angle Theta1 pour l'estimation du segment de
C                         troncature de la fonction de phase des aérosols.
C    CTE_AER_MU2_TRONCA : cos de l'angle Theta2 pour l'estimation du segment de
C                         troncature de la fonction de phase des aérosols.
C
C
C Constantes spécifiques :
C    LENLINE  : Longueur maximale d'une ligne de données lue dans un fichier
C    CTE_NOT_DEFINED_VALUE_DBLE: valeur par défaut pour les reels
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define LENLINE 1000

C==============================================================================
C PROGRAMME: SOS_AEROSOLS
C ==========  
C  Modèles simulés :
C  ----------------
C      Ce programme calcule les proprietes radiatives d'une granulométrie 
C      d'aérosols mono-modale, multi-modale (de la WMO) ou bi-modale (S&F).
C
C      Les constituants de base des modèles de la WMO sont:
C      dust-like, oceanic, water-soluble et soot. 
C      Quatre modèles WMO types peuvent etre simules: continental, maritime, 
C      urbain et une composition definie par l'utilisateur.
C
C      Les constituants de base des modèles de Shettle & Fenn sont:
C      small rural, large rural, small urban, large urban et oceanic. 
C      Cinq modèles S&F types peuvent etre simules: tropospherique, urbain,  
C      maritime et cotier.
C
C      Si imodalite = 0 => prise en compte d'une seule granulométrie d'aérosols
C			   (Log Normale ou loi de Junge).
C
C      Si imodalite = 1 => cas de granulométries multi-modales, utilisation  
C                          des quatre modèles standards WMO.
C
C      Si imodalite = 2 => cas de granulométries bi-modales, utilisation  
C                          des cinq modèles standards S&F.
C
C      Si imodalite = 3 => cas de granulométries LND bi-modales.
C
C      Si imodalite = 4 => cas de l'utilisation d'une matrice de phase externe 
C                         (fonctions de phase en intensite et polarisation).
C
C      Si imodalite = 5 => cas de l'utilisation d'un fichier utilisateur 
C                          définissant un mélange de modèles d'aérosols.
C
C
C  Calculs de Mie et application :
C  -----------------------------
C      Les calculs de MIE donnent les proprietes radiatives d'un composant. 
C      Ils sont realises par la procédure "SOS_MIE".
C      Les fichiers de MIE sont stockes sous un répertoire specifique. 
C
C      Le pas sur le paramètre de taille alpha est identique est :
C         si        alpha  <= 0.1  ===> pas = 0.0001
C         si 0.1  < alpha  <= 1.0  ===> pas = 0.001
C         si 1.0  < alpha  <= 10   ===> pas = 0.01
C         si 10   < alpha  <= 30   ===> pas = 0.05
C         si 30   < alpha  <= 100  ===> pas = 0.1
C         si 100  < alpha          ===> pas = 1.0
C
C      A partir des calculs de Mie, on determine la fonction de phase 
C      P(mu) et les fonctions de phase polarisee Q(mu) et U(mu), pour une 
C      longueur d'onde donnée et une distribution en taille des particules 
C      donnée pour chaque modèle.
C
C
C      Une troncature est realisee sur la fonction de phase P(mu) si l'option
C      itronc = 1 est selectionnee.
C
C      Les fonctions de phase sont ensuite developpees dans la base des
C      polynomes de Legendre. On obtient :
C         - Pour P(mu) => les coefficients beta(k) de 0 a OS_NB,
C         - Pour Q(mu) => les coefficients gamma(k) de 0 a OS_NB,
C         - Pour U(mu) => les coefficients delta(k) de 0 a OS_NB.
C
C      On calcule les coefficients alp(k) et zeta(k) en fonction des beta(k)
C      et delta(k). Ils serviront pour l'expression de la matrice de phase
C      decomposee en series de Fourier (dans le code des OS).
C
C 
C  Remarque :
C  --------    
C      Si l'épaisseur optique des aérosols est nulle, il y a creation 
C      d'un fichier Aérosols avec des valeurs mises a zero.
C
C
C Donnees en entrée du programme
C ------------------------------
C     -->Les paramètres/arguments du programme
C
C             FICANGLES (CHARACTER) : (E) nom complet du fichier des angles utilises 
C                                     pour les calculs de luminance et BRDF
C                 		     (répertoire + nom fichier + extension)
C
C      	TA   (DOUBLE)  : (E) Epaisseur optique des aérosols (non tronquée).
C                                (format de decodage en réel F9.5)
C      		       
C	       WA   (DOUBLE)  : (E) Longueur d'onde de simulation en microns.	
C                                (format de decodage en réel F9.5)
C
C             WAREF (DOUBLE) : (E) longueur d'onde de référence WAREF.
C                                (format de decodage en réel F9.5)
C	 
C      	TAWAREF (DOUBLE) : (E) AOT totale pour la longueur d'onde WAREF.
C                                  (format de decodage en réel F9.5)
C     		       
C	       ITRONC  (I4)   : (E) Option de troncature  (1 pour la realiser).
C                                (format de decodage en entier I4)
C
C	       IMOD  (I4)  : (E) Granulométrie mono-modal ou multi-modale.
C                               (format de decodage en entier I2)
C      		        	0 : Modèles mono-mode
C				       1 : Modèles WMO
C				       2 : Modèles Shettle & Fenn
C		        	       3 : Modèles bimodaux de LND
C                       	       4 : Utilisation de fonctions de phase externes
C                                  5 : Utilisation d'un fichier définissant un mélange de modes
C
C	       RN   (DOUBLE)  : (E) ==> Associe a IMOD = 0
C				   Partie réelle de l'indice de réfraction pour modèle mono-modal
C                                  a la longueur d'onde de simulation des luminances (WA)
C                                 (format de decodage en réel F5.3)	
C
C             IN   (DOUBLE)  : (E) ==> Associe a IMOD = 0
C				   Partie imaginaire de l'indice de réfraction pour modèle mono-modal
C                                  a la longueur d'onde de simulation des luminances (WA)
C                                (format de decodage en réel F8.5)				
C	  
C             IGRANU  (I4)   : (E) ==> Associe a IMOD = 0
C			      	   Indice de type de granulométrie pour modèle mono-modal
C                          		1 : LND (paramètres : rayon modal en microns et variance). 
C                          		2 : Loi de Junge (paramètres : rayon minimal  
C                              		en microns et puissance appliquee au rayon).					
C				   (format de decodage en entier I4)
C
C             VARGRANU1_MMD  (DOUBLE) (E) :  1er paramètre de granulométrie pour modèle mono-modal
C             VARGRANU2_MMD  (DOUBLE) (E) :  2nd paramètre de granulométrie pour modèle mono-modal
C             VARGRANU3_MMD  (DOUBLE) (E) :  3eme paramètre de granulométrie pour modèle mono-modal	
C	 		   If LND       : modal radius and standard deviation
C	  		   If Junge law : min radius, slope and max radius
C
C             IMODELE_WMO   (I4)  : (E) ==> Associe a IMOD = 1
C				    Indice de choix du modèle WMO
C				   (format de decodage en entier I2)
C      		        		1 : Modèle WMO Continental
C					2 : Modèle WMO Maritime
C					3 : Modèle WMO Urban
C					4 : Modèle WMO "user definition"
C
C             C_WMO_DL (DOUBLE)  : (E) ==> Associe a IMOD = 1
C				       Pourcentage volumique C des particules DL (WMO).
C				      (format de decodage en réel F9.5)
C
C             C_WMO_WS (DOUBLE)  : (E) ==> Associe a IMOD = 1
C				       Pourcentage volumique C des particules WS (WMO).
C				      (format de decodage en réel F9.5)
C
C             C_WMO_OC (DOUBLE)  : (E) ==> Associe a IMOD = 1
C				       Pourcentage volumique C des particules OC (WMO).
C				      (format de decodage en réel F9.5)
C
C             C_WMO_SO (DOUBLE)  : (E) ==> Associe a IMOD = 1
C				       Pourcentage volumique C des particules SO (WMO).
C				       (format de decodage en réel F9.5)
C
C             IMODELE_SF (DOUBLE) : (E) ==> Associe a IMOD = 2
C					Indice de choix du modèle Shettle & Fenn
C				 	(format de decodage en entier I2)
C      		        		1 : Modèle SF Tropospheric
C					       2 : Modèle SF Urban
C					       3 : Modèle SF Maritime
C					       4 : Modèle SF Coastal
C
C             RH  (DOUBLE) : (E) ==> Associe a IMOD = 2
C				 Humidite relative pour modèle Shettle & Fenn
C				(format de decodage en réel F5.2)
C
C             MODE_PARAM_BILND (I4) : (E) ==> Associe a IMOD = 3
C					  Indice de type de description du mélange bimodal LND.                 
C                               	 (format de decodage en entier I2)
C      		        		 1 : Utilisation des concentrations volumiques de l'utilisateur
C					        2 : Utilisation du rapport AOT_coarse / AOT_totale
C
C             USER_CV_COARSE (DOUBLE) : (E) ==> Associe a IMOD = 3 et MODE_PARAM_BILND = 1
C      	     					Concentration volumique de l'utilisateur pour le mode "LND coarse". 
C                                		       (format de decodage en réel F8.5)
C	 
C             USER_CV_FINE (DOUBLE) : (E) ==> Associe a IMOD = 3 et MODE_PARAM_BILND = 1
C      	     			      	  Concentration volumique de l'utilisateur pour le mode "LND fine". 
C                                	         (format de decodage en réel F8.5)
C	 
C             RTAUct_WAREF (DOUBLE) : (E) ==> Associe a IMOD = 3 et MODE_PARAM_BILND = 2
C      	     			     	  rapport AOT_coarse / AOT_tot pour la longueur d'onde WAREF.
C                                    	  (format de decodage en réel F9.5)
C
C             BMD_CM_MRwa  (DOUBLE) : (E) ==> Associe a IMOD = 3 
C      	     				  Partie réelle de l'indice de réfraction
C                                   	  des particules du mode "LND coarse" pour WA de simulation
C                               	         (format de decodage en réel F5.3)
C
C             BMD_CM_MIwa  (DOUBLE) : (E) ==> Associe a IMOD = 3 
C      	     				  Partie imaginaire de l'indice de réfraction
C                                   	  des particules du mode "LND coarse" pour WA de simulation
C                               	         (format de decodage en réel F8.5)
C
C             BMD_CM_MRwaref  (DOUBLE) : (E) ==> Associe a IMOD = 3  
C      	     				     Partie réelle de l'indice de réfraction
C                                   	     des particules du mode "LND coarse" pour WAREF
C                               	            (format de decodage en réel F5.3)
C
C             BMD_CM_MIwaref (DOUBLE) : (E) ==> Associe a IMOD = 3 
C      	     				    Partie imaginaire de l'indice de réfraction
C                                   	    des particules du mode "LND coarse" pour WAREF
C                               	           (format de decodage en réel F8.5)
C
C             BMD_CM_RMODAL (DOUBLE) : (E)  ==> Associe a IMOD =  3 
C      	     				    Rayon modal de la "LND coarse"
C                               	           (format de decodage en réel F9.5)
C
C             BMD_CM_VAR (DOUBLE) : (E) ==> Associe a IMOD = 3 
C      	     				    Variance de la "LND coarse"
C                               	    (format de decodage en réel F9.5)
C
C             BMD_FM_MRwa  (DOUBLE) : (E) ==> Associe a IMOD = 3 
C      	     				  Partie réelle de l'indice de réfraction
C                                   	  des particules du mode "LND fine" pour WA de simulation
C                               	  (format de decodage en réel F5.3)
C
C             BMD_FM_MIwa  (DOUBLE) : (E) ==> Associe a IMOD = 3 
C      	     				  Partie imaginaire de l'indice de réfraction
C                                   	  des particules du mode "LND fine" pour WA de simulation
C                               	         (format de decodage en réel F8.5)
C
C             BMD_FM_MRwaref  (DOUBLE) : (E) ==> Associe a IMOD = 3  
C      	     				     Partie réelle de l'indice de réfraction
C                                   	     des particules du mode "LND fine" pour WAREF
C                               	            (format de decodage en réel F5.3)
C
C             BMD_FM_MIwaref (DOUBLE) : (E) ==> Associe a IMOD = 3 
C      	     				    Partie imaginaire de l'indice de réfraction
C                                   	    des particules du mode "LND fine" pour WAREF
C                               	           (format de decodage en réel F8.5)
C
C             BMD_FM_RMODAL (DOUBLE) : (E)  ==> Associe a IMOD =  3 
C      	     				    Rayon modal de la "LND fine"
C                               	           (format de decodage en réel F9.5)
C
C             BMD_FM_VAR (DOUBLE) : (E) ==> Associe a IMOD = 3 
C      	     				    Variance de la "LND fine"
C                               	           (format de decodage en réel F9.5)
C
C             FICEXTDATA (CHARACTER) : (E) ==> Associe a IMOD = 4 
C                     			   Nom du fichier de données externes de fonctions de phase
C                     			   et paramètres radiatifs (coefficients extinction et diffusion) 
C                    			   (répertoire + nom fichier + extension)
C
C             FICMIXTURE_AER (CHARACTER) : (E) ==> Associe a IMOD = 5 
C                     			   Nom du fichier de définition d'un mélange de modes aérosols 
C                    			   (répertoire + nom fichier + extension)
C
C             DIR_TMP  (CHARACTER)  : (E) répertoire temporaire de génération des fichiers
C
C             DIRS   (CHARACTER)  : (E) répertoire de stockage des fichiers MIE
C
C             FICLOGMIE (CHARACTER) : (E) Fichier trace MIE  
C				              (0 ou chemin complet au fichier Trace)
C                 		              (répertoire + nom fichier + extension)	
C
C             FICGRANU (CHARACTER) : (E) nom complet du fichier résultat GRANU
C                 		             (répertoire + nom fichier + extension)
C 
C             FICLOGGRANU (CHARACTER) : (E) Fichier trace du programme Aérosols 
C					        (0 ou chemin complet au fichier Trace)
C                 			        (répertoire + nom fichier + extension)
C
C      	KMAT1 (DOUBLE) : (S) Section efficace d'extinction (micron^2). Non tronqué.
C
C      	PIZ (DOUBLE) : (S) Albédo de simple diffusion.C
C
C      	COEF_TRONCA (DOUBLE) : Coefficient de troncature de la fonction de phase (S) 
C
C      	IER (I4) : (S) code d'erreur =0 si pas d'erreur, =-1 sinon
C 
C                          
C
C    -->Fichiers
C       Le fichier MIE
C       Le fichier Trace Granu
C       Le fichier Trace Mie
C       Le fichier des angles et paramètres internes pour les limites des développements.
C       Le fichier des fonctions de phase externes (entrée optionnelle)
C       Le fichier de définition d'un mélange de modes aérosols (entrée optionnelle)
C
C Resultats fournis par le programme
C ----------------------------------
C   Fichier résultat GRANU contenant :
C	- coefficient d'extinction et de diffusion, Beta(0)/3
C	- coefficient de troncature
C	- albédo de simple diffusion tronquée
C	- les paramètres de decomposition de la fonction de phase en 
C	  polynomes de Legendre.
C
C Variable d'environnement
C ------------------------
C   SOS_ABS_ROOT : qui définit l'arborescence d'accès au répertoire principal du code SOS_ABS
C   Cette variable d'environnement permet de localiser les fichiers sous $SOS_ABS_ROOT/fic 
C
C
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier de MIE
C       Fichier binaire non formatte, cree par SOS_MIE
C       1 enregistrement d'entete avec
C                Indice de réfraction (partie réelle)
C                Indice de réfraction (partie imaginaire)
C                Paramètre de taille limite
C                Nombre d'angles de Gauss positifs
C
C      1 enregistrement pour chaque valeur du paramètre de taille avec:
C        - Paramètre de taille (réel),
C        - Coefficient d'efficacite en extinction (réel),
C        - Coefficient d'efficacite en diffusion (réel),
C        - Facteur d'anisotropie (double précision),
C        - Tableau de la fonction de phase Mie(mu) (réel),
C        - Tableau de la fonction de phase polarisee Qmie(mu) (réel),
C        - Tableau de la fonction de phase polarisee Umie(mu) (réel).
C
C    --> Les tableaux ont une taille utile limitée à MIE_NBMU : 
C        ils sont crees et lus avec cette limitation.
C
C
C    -->Contenu du fichier des angles:
C       Fichier ascii formatte, cree par SOS_ANGLES.
C       Il contient l'information sur les angles à utiliser 
C       et sur les ordres des développements adaptes : 
C		NB_TOTAL_ANGLES :  nb d'angles total à utiliser
C               NB_GAUSS_ANGLES :  nb d'angles de Gauss
C               ANGLES_USERFILE :  Nom du fichier des angles utilisateur
C               INTERNAL_OS_NB :   Ordre des développements des fonctions de phase 
C				   en Polynomes de Legendre à utiliser.
C               INDEX COS_ANGLE WEIGHT : lignes donnant le numéro d'angle, son cosinus et son poids
C                                        au format I4,1X,2D21.14
C 
C
C    -->Contenu du fichier des fonctions de phase externes:  
C       Fichier ascii formaté, cree par l'utilisateur.
C       Le format de ce fichier s'inspire de celui de l'outil DLS d'Oleg Dubovik (LOA)
C       fournissant les propriétés de diffusion primaire de particules
C       sphéroides homogènes orientées de facon aléatoire.
C       Les premières lignes ont dû être adaptées.
C 
C          1ere ligne : EXTINCTION_COEF : Valeur 
C              Valeurs : Section efficace d'extinction (micron²)
C              Format  : Commentaire + séparateur « : » + valeur (réel non formaté)
C                
C          2eme ligne : SCATTERING_COEF : Valeur 
C              Valeurs : Section efficace de diffusion (micron²)
C              Format  : Commentaire + séparateur « : » + valeur (réel non formaté)
C
C	  3eme ligne : NB_LINES : nbAng
C              Valeur : Nombre d'angles décrivant les fonctions de phase.
C              Format : Commentaire + séparateur « : » + valeur (entier non formaté)
C              
C          4eme ligne : Ligne de commentaire 
C          
C          Lignes 5 à 5+nbAng-1
C               Valeurs : ANGLE  F11  -F12/F11  F22/F1   F33/F11
C                           ANGLE 	: angle en degrés.
C                           F11 	: fonction de phase en intensité
C 			   -F12/F11 	: rapport des fonctions -F12 et F11
C 			   F22/F11 	: rapport des fonctions F22 et F11
C 			   F33/F11 	: rapport des fonctions F33 et F11 
C               ----------------------   ===========================================
C               NB CONVENTION DE SIGNE : F12 est négative pour la diffusion Rayleigh.
C                                        = convention des OS
C               ----------------------   ===========================================
C 		 Format  : Valeurs réelles, non formatées
C 
C          Pas de lecture des lignes suivantes.	
C
C 
C
C    -->Contenu du fichier de définition d'un mélange de modèles d'aérosols:  
C       Fichier ascii formaté, créé par l'utilisateur.
C       
C          1ere ligne : NUMBER OF AEROSOLS MODES: Valeur 
C              Valeurs : Nombre de blocs de définition de modes
C              Format  : Commentaire + séparateur « : » + valeur (entier non formaté)
C               
C          Par bloc de données d'un mode: 
C                 Bloc ligne 1 : SIZE DISTRIBUTION  : Valeur 
C                 Valeurs : Type de granulomètrie (LND ou JUNGE)
C                 Format  : Commentaire + séparateur « : » + valeur (chaine de caractères)
C
C                 Cas d'une granulométrie LND:
C	              Bloc ligne 2 : MODAL RADIUS  : Valeur
C                    Valeur : Rayon modal de la LND  (en microns)
C                    Format : Commentaire + séparateur « : » + valeur (réel non formaté)
C              
C	              Bloc ligne 3 : STANDARD DEVIATION  : Valeur
C                     Valeur : Variance de la LND (en microns)
C                     Format : Commentaire + séparateur « : » + valeur (réel non formaté)
C  
C                 Cas d'une granulométrie JUNGE:
C	               Bloc ligne 2 : SLOPE : Valeur
C                     Valeur : Pente de la distribution de Junge 
C                     Format : Commentaire + séparateur « : » + valeur (réel non formaté)
C              
C	               Bloc ligne 3 : MIN RADIUS : Valeur
C                     Valeur : Rayon minimal de la distribution de Junge (en microns)
C                     Format : Commentaire + séparateur « : » + valeur (réel non formaté)
C              
C	               Bloc ligne 4 : MAX RADIUS : Valeur
C                     Valeur : Rayon maximal de la distribution de Junge (en microns)
C                     Format : Commentaire + séparateur « : » + valeur (réel non formaté)                         
C    
C                 Bloc ligne + 1 : REFRACTIVE INDEX at WA_SIMU - Real part  : Valeur 
C                 Valeurs : Partie réelle de l'indice de réfraction à la longueur d'onde de simulation de la luminance
C                 Format  : Commentaire + séparateur « : » + valeur (chaine de caractères)
C    
C                 Bloc ligne + 2 : REFRACTIVE INDEX at WA_SIMU - Imaginary part  : Valeur 
C                 Valeurs : Partie imaginaire (<0) de l'indice de réfraction à la longueur d'onde de simulation de la luminance
C                 Format  : Commentaire + séparateur « : » + valeur (chaine de caractères)
C    
C                 Bloc ligne + 3 : REFRACTIVE INDEX at WA_SIMU - Real part  : Valeur 
C                 Valeurs : Partie réelle de l'indice de réfraction à la longueur d'onde de référence de l'AOT
C                 Format  : Commentaire + séparateur « : » + valeur (chaine de caractères)
C    
C                 Bloc ligne + 4 : REFRACTIVE INDEX at WA_SIMU - Imaginary part  : Valeur 
C                 Valeurs : Partie imaginaire (<0) de l'indice de réfraction à la longueur d'onde de référence de l'AOT
C                 Format  : Commentaire + séparateur « : » + valeur (chaine de caractères)
C    
C                 Bloc ligne + 5 : PROPORTION TO THE TOTAL AOT at WA_REF  : Valeur 
C                 Valeurs : Proportion du mode (entre 0 et 1) à la contribution à l'AOT à la longueur d'onde de référence
C                 Format  : Commentaire + séparateur « : » + valeur (chaine de caractères)
C
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
C     - SOS_NOM_FICMIE
C     - SOS_MIE
C     - SOS_INIT_PARAMWMO
C     - SOS_INIT_PARAMSF
C     - SOS_INTERPOL
C     - SOS_GRANU
C     - SOS_DECOMPO_LEGENDRE
C     - SOS_INTERPO_SPLINT
C     - SOS_SPLINE
C     - SOS_SPLINT
C
C============================================================================
      SUBROUTINE SOS_AEROSOLS(FICANGLES, 
     &	    	        WA, TA, WAREF, TA_WAREF,	  
     &	  	        ITRONC, IMOD,
     &	  	        RN, IN, IGRANU, 
     &	  	        VARGRANU1_MMD, VARGRANU2_MMD, VARGRANU3_MMD,
     &                IMODELE_WMO, 
     &	  	        C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &	  	        IMODELE_SF, RH, 
     &	  	        MODE_PARAM_BILND, 
     &	  	        USER_CV_COARSE, USER_CV_FINE, RTAUct_WAREF,
     &	  	        BMD_CM_MRWA, BMD_CM_MIWA, 
     &	  	        BMD_CM_MRWAREF, BMD_CM_MIWAREF,
     &	  	        BMD_CM_RMODAL, BMD_CM_VAR,
     &	  	        BMD_FM_MRWA, BMD_FM_MIWA, 
     &	  	        BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &	  	        BMD_FM_RMODAL, BMD_FM_VAR,	 
     &	  	        FICEXTDATA, FICMIXTURE_AER,
     &			 DIR_TMP, DIRS, FICLOGMIE,
     &	  	        FICGRANU, FICLOGGRANU, KMAT1,
     &			 PIZ, COEF_TRONCA, IER)
 
      
      IMPLICIT NONE

C Constantes AEROSOLS
C -------------------

      INTEGER CTE_MODALITE	!Valeur max du paramètre modalite 
      PARAMETER (CTE_MODALITE=4)

      INTEGER CTE_NB_MODELE	!Nombre de modèles WMO ou Shettle & Fenn possibles.
      PARAMETER (CTE_NB_MODELE=4)

      INTEGER CTE_NB_WMO	    !Nombre de types d'aérosols de la WMO utilisés:
      PARAMETER (CTE_NB_WMO=4)  !DL,WS,OC,SO.
      
      INTEGER CTE_NB_SF	    !Nombre de types d'aérosols du modèle Shettle & Fenn 
      PARAMETER (CTE_NB_SF=5)   !utilises:SR,LR,SU,LU,OM.
          
      INTEGER CTE_NB_AERO	    !Initialisation avec le plus grand nombre de types 
      PARAMETER (CTE_NB_AERO=5) !d'aérosols (cas Shettle & Fenn).	
     

C* Definition des variables                         
C*-----------------------------------------------------------------


      DOUBLE PRECISION XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)  
                                                  !Cosinus des angles de Gauss (et utilisateur).
      DOUBLE PRECISION XHR(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)  
                                                  !Poids des angles de Gauss.				

      DOUBLE PRECISION TA	! Epaisseur optique des aérosols (non tronquée).

      DOUBLE PRECISION WA	! Longueur d'onde en microns selectionnee.	
            						  						  
      DOUBLE PRECISION C(CTE_NB_WMO)	 !Proportion volumique des constituants WMO
      DOUBLE PRECISION C_WMO_DL	    	 ! --> Valeur utilisateur pour particule DL
      DOUBLE PRECISION C_WMO_WS	    	 ! --> Valeur utilisateur pour particule WS
      DOUBLE PRECISION C_WMO_OC	    	 ! --> Valeur utilisateur pour particule OC
      DOUBLE PRECISION C_WMO_SO	    	 ! --> Valeur utilisateur pour particule SO
      DOUBLE PRECISION V(CTE_NB_WMO)	 !Concentration volumique des constituants WMO
      DOUBLE PRECISION N(CTE_NB_WMO)	 ! C/V
      DOUBLE PRECISION NTOT		        ! Somme des coefficients de ponderation 
                                           ! des composants aérosols.
      DOUBLE PRECISION ALPHAF_PARTICULE(4) ! Valeur Alpha Max pour chaque
      					        ! composant de la WMO (DL, WS, OC, SO)
					   
      DOUBLE PRECISION Ni(CTE_NB_SF)	 !densite de particules (Shettle & Fenn)					   
      DOUBLE PRECISION RH                  ! % d'humidite relative des modèles de Shettle & Fenn



 	!   Modèle mono-mode
	!------------------- 
 	! Indice de réfraction pour modèle mono-mode     
      DOUBLE PRECISION RN	! partie réelle de l'indice de réfraction
      DOUBLE PRECISION IN	! partie imaginaire de l'indice de réfraction
 	! Paramètres de granulométrie pour modèle mono-mode   
      DOUBLE PRECISION VARGRANU1_MMD    ! Paramètre 1 de la granulométrie.
      DOUBLE PRECISION VARGRANU2_MMD    ! Paramètre 2 de la granulométrie.
      DOUBLE PRECISION VARGRANU3_MMD    ! Paramètre 3 de la granulométrie.



      DOUBLE PRECISION ALPHAO   ! valeur de alpha minimale des calculs de Mie		
      DOUBLE PRECISION ALPHAF   ! valeur de alpha maximale des calculs de Mie			

      DOUBLE PRECISION SOMME_NR	! nombre de particules : integration sur les rayons de NR * PR.  
      
      DOUBLE PRECISION MR(CTE_NB_AERO)   ! partie réelle de l'indice de réfraction
      					      ! pour chaque composant d'aérosols
      DOUBLE PRECISION MI(CTE_NB_AERO)   ! partie imaginaire de l'indice de réfraction.
      					      ! pour chaque composant d'aérosols.
 



      DOUBLE PRECISION VARGRANU1I(CTE_NB_AERO)   ! Paramètre 1 de la granulométrie
      						       ! pour chaque modèle de constituant.
      DOUBLE PRECISION VARGRANU2I(CTE_NB_AERO)   ! Paramètre 2 de la granulométrie
      						       ! pour chaque modèle de constituant.      
	      				

	!Paramètres pour un modèle bimodal de LND
      INTEGER*2 MODE_PARAM_BILND   !Indice du choix de description du mélange de 2 LND.
      DOUBLE PRECISION CVI(2)	!Concentration volumique de la "LND coarse" (1) et 
      				      	!de la "LND fine" (2)
      DOUBLE PRECISION USER_CV_COARSE	!Valeur utilisateur de la concentration volumique de la "LND coarse"
      DOUBLE PRECISION USER_CV_FINE	!Valeur utilisateur de la concentration volumique de la "LND fine"	
      DOUBLE PRECISION RTAUct_WAREF	!Rapport des épaisseurs optique AOT_coarse / AOT_total 
      DOUBLE PRECISION WAREF		!Longueur d'onde du rapport AOT_coarse / AOT_total
      DOUBLE PRECISION TA_WAREF		!Epaisseur optique aérosol totale 
      					       !pour la longueur d'onde de référence.
					
      DOUBLE PRECISION RN_LND(2)	!Partie réelle de l'indice de réfraction de la 
      					!"LND coarse" (1) et de la "LND fine" (2) pour la longeur d'onde WA.					
      DOUBLE PRECISION BMD_CM_MRWA	! Valeur lue pour RN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MRWA	! Valeur lue pour RN_LND(2) : fine mode    
       
      DOUBLE PRECISION IN_LND(2) 	!Partie imaginaire de l'indice de réfraction de la 
      					!"LND coarse" (1) et de la "LND fine" (2) pour la longeur d'onde WA.					
      DOUBLE PRECISION BMD_CM_MIWA	! Valeur lue pour IN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MIWA	! Valeur lue pour IN_LND(2) : fine mode 
      
      DOUBLE PRECISION RN_LND_WAREF(2)	!Partie réelle de l'indice de réfraction de la 
      					       !"LND coarse" (1) et de la "LND fine" (2) pour la longueur d'onde WAREF.
      DOUBLE PRECISION BMD_CM_MRWAREF	! Valeur lue pour RN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MIWAREF	! Valeur lue pour RN_LND_WAREF(2) : fine mode 
      
      DOUBLE PRECISION IN_LND_WAREF(2) 	!Partie imaginaire de l'indice de réfraction de la 
      					       !"LND coarse" (1) et de la "LND fine" (2) pour la longueur d'onde WAREF.
      DOUBLE PRECISION BMD_FM_MRWAREF	! Valeur lue pour IN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MIWAREF	! Valeur lue pour IN_LND_WAREF(2) : fine mode      
      
      DOUBLE PRECISION RMODAL_LND(2)  	!Rayon modal des "LND coarse" (1) et "fine" (2).
      DOUBLE PRECISION BMD_CM_RMODAL	! Valeur lue pour RMODAL_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_RMODAL	! Valeur lue pour RMODAL_LND(2) : fine mode  
      
      DOUBLE PRECISION VAR_LND(2)	!Variance des "LND coarse" (1) et "fine" (2).
      DOUBLE PRECISION BMD_CM_VAR	! Valeur lue pour VAR_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_VAR	! Valeur lue pour VAR_LND(2) : fine mode 
          
      
      DOUBLE PRECISION KMAT1c      ! Section efficace d'extinction (micron^2) pour le mode 
      					! "coarse" a la longueur d'onde WAREF. Non tronqué.
      DOUBLE PRECISION KMAT1f	! Section efficace d'extinction (micron^2) pour le mode 
      					! "fine" a la longueur d'onde WAREF. Non tronqué.
   

  
      DOUBLE PRECISION KMAT1I	! Section efficace d'extinction (micron^2)
      				       ! pour un constituant I. Non tronqué.
      DOUBLE PRECISION KMAT2I 	! Section efficace de diffusion (micron^2)
      				       ! pour un constituant I.  Non tronqué.       
      DOUBLE PRECISION TAUEXTI     ! Epaisseurs optiques d'extinction d'un mode I des aérosols 
      DOUBLE PRECISION TAUSCAI     ! Epaisseurs optiques de diffusion d'un mode I des aérosols

      DOUBLE PRECISION KMAT1	! Section efficace d'extinction (micron^2). Non tronqué.
      DOUBLE PRECISION KMAT2 	! Section efficace de diffusion (micron^2). Non tronqué.      
      DOUBLE PRECISION COEF_TRONCA ! Coefficient de troncature.
      DOUBLE PRECISION PIZ	       ! Albédo de simple diffusion.
      DOUBLE PRECISION PIZTR	! Albédo de simple diffusion apres troncature.
      
      
      DOUBLE PRECISION P11I(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Fonction de phase calculée à partir des fonctions de 
		         !la théorie de MIE pour le constituant I.
      DOUBLE PRECISION P12I(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P12 de la matrice de phase pour le constituant I.
      DOUBLE PRECISION P33I(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P33 de la matrice de phase pour le constituant I.
      
      DOUBLE PRECISION P11(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Fonction de phase du modèle simule.
      DOUBLE PRECISION P12(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P12 de la matrice de phase.
      DOUBLE PRECISION P22(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P22 de la matrice de phase.
      DOUBLE PRECISION P33(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P33 de la matrice de phase.
		       
      DOUBLE PRECISION TTT(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Copie de la fonction de phase (P11) avant troncature.
     
      DOUBLE PRECISION Z1	!Valeur de BETA11(0) avant normalisation.
      DOUBLE PRECISION TETA !Angle de diffusion.
      
                  
      DOUBLE PRECISION ALP(0:CTE_OS_NB_MAX)	! Paramètres alpha(k).
      DOUBLE PRECISION BETA11(0:CTE_OS_NB_MAX)	! Paramètres beta11(k).
      DOUBLE PRECISION BETA22(0:CTE_OS_NB_MAX)	! Paramètres beta22(k).
      DOUBLE PRECISION GAMMA12(0:CTE_OS_NB_MAX)	! Paramètres gamma12(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Paramètres zeta(k).
      DOUBLE PRECISION DELTA33(0:CTE_OS_NB_MAX)	! Paramètres delta33(k).
     
	    
      DOUBLE PRECISION X, RMAX   ! Paramètres pour le calcul de AlphaMax 
      				     ! selon la taille des particules.
               
      DOUBLE PRECISION VAL	! Variable pour la lecture de valeurs inutilisees
      
      DOUBLE PRECISION ANG_EXT(CTE_MAXNB_ANG_EXT)     ! Table des angles des fonctions de phase externes
      DOUBLE PRECISION MU_EXT(CTE_MAXNB_ANG_EXT)	     ! Table des cosinus des angles des fonctions de phase externes   
      DOUBLE PRECISION F11_EXT(CTE_MAXNB_ANG_EXT)     ! Fonction de phase P11 de données externes
      DOUBLE PRECISION F12_EXT(CTE_MAXNB_ANG_EXT)     ! Fonction de phase P11 de données externes
      DOUBLE PRECISION F22_EXT(CTE_MAXNB_ANG_EXT)     ! Fonction de phase P11 de données externes
      DOUBLE PRECISION F33_EXT(CTE_MAXNB_ANG_EXT)     ! Fonction de phase P11 de données externes
      DOUBLE PRECISION F12sF11_EXT(CTE_MAXNB_ANG_EXT) ! Fonction de phase P12/P11 de données externes
      DOUBLE PRECISION F22sF11_EXT(CTE_MAXNB_ANG_EXT) ! Fonction de phase P22/P11 de données externes
      DOUBLE PRECISION F33sF11_EXT(CTE_MAXNB_ANG_EXT) ! Fonction de phase P33/P11 de données externes   
      
      INTEGER*2 IMOD	       ! Indice du cas de modélisation des aérosols.
      INTEGER*2 IMODELE_WMO ! Type du modèle WMO.
      INTEGER*2 IMODELE_SF  ! Type du modèle Shettle & Fenn.

      INTEGER*4 NBMU_GAUSS  ! Nombre d'angles de Gauss utilises 
      INTEGER*4 MIE_NBMU   ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB   ! Ordre des développements des fonctions de phase 
			   ! en Polynomes de Legendre.      
      INTEGER*4 IGRANU	   ! Numéro du modèle de granulométrie
      INTEGER*4 ITRONC	   ! Option de troncature  (1 pour la realiser). 
      INTEGER*4 LDIRTMP        ! longueur de la chaine de caractères DIRTMP 
      INTEGER*4 LDIRS	   ! longueur de la chaine de caractères DIRS
      INTEGER*4 LFICS	   ! longueur de la chaine de caractères FICS  
      INTEGER*4 LFICTMP	   ! longueur de la chaine de caractères FICMIE_TMP  
      INTEGER*4 LFIC           ! longueur de chaine de caractères
      INTEGER*4 IFIN           ! nb de caractères significatifs d'une chaine           
      INTEGER*4 LCOM           ! longueur de la chaine de caractères COMMAND     
      INTEGER*4 I	      ! Indice de tableau 
      INTEGER*4 J	      ! Indice de tableau 
      INTEGER*4 K          ! Indice de tableau        
      INTEGER*4 IER        ! code d'erreur =0 si pas d'erreur, =-1 sinon
      INTEGER*4 IARGC      ! Type de la fonction IARGC
      INTEGER*4 NPARAM     ! Numéro du paramètre en entrée (type INTEGER*4 obligatoire pour le f77)
      INTEGER*4 NB_FICEXT_ANGLES ! Nombre d'angles du fichier de fonctions de phase externe
  
      INTEGER*4 NB_MODES_MELANGE   ! Nombre de modes aérosols pour un mélange défini par l'utilisateur
      INTEGER*4 NMOD			! Numéro d'un mode aérosol
    
      CHARACTER*CTE_LENDIR DIR_TMP        ! chemin d'acces aux fichiers MIE en cours de génération
      CHARACTER*CTE_LENDIR DIRS	       ! chemin d'acces au fichier MIE
      CHARACTER*CTE_LENFIC1 FICS		! nom du fichier MIE 
      CHARACTER*CTE_LENFIC2 FICMIE	       ! chemin et nom du fichier MIE
      CHARACTER*CTE_LENFIC2 FICMIE_TMP    ! chemin et nom du fichier MIE en cours de génération
      CHARACTER*CTE_LENCOM  COMMAND       ! Chaine de commande pour appel système.
      CHARACTER*CTE_LENFIC2 FICLOGMIE           ! nom du fichier TRACE de MIE  
      CHARACTER*CTE_LENFIC2 FICLOGGRANU         ! nom du fichier TRACE de GRANU
      CHARACTER*CTE_LENFIC2 FICGRANU	      ! nom du fichier résultat GRANU  
      CHARACTER*CTE_LENFIC2 FICANGLES           ! chemin et nom du fichier des angles effectifs    
      CHARACTER*CTE_LENFIC2 FICANGLES_USER      ! chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*CTE_LENFIC2 FICEXTDATA          ! chemin et nom du fichier de données externes de fonctions
      						      ! de phase et paramètres radiatifs (coefficients extinction et
						      ! diffusion)
 
      CHARACTER*CTE_LENFIC2 FICMIXTURE_AER      ! chemin et nom du fichier de définition d'un mélange de modèles 
                                                ! d'aérosols par l'utilisateur  
     
      CHARACTER*CTE_LENDIR DIRFIC	             ! chemin d'acces au fichier des paramètres des modèles WMO et S&F
      CHARACTER*CTE_LENFIC2 FICDATAWMO          ! chemin et nom du fichier Data_WMO
      
      CHARACTER*15 PARTICULE(CTE_NB_AERO)       ! nom des différentes particules de la WMO
      					             ! ou des modèles S&F	

      CHARACTER*LENLINE LIGNE_TEXTE             ! Ligne de caractères pour décodage de paramètres      
      CHARACTER*1000 CVAL	                    ! Chaine de caractères pour la lecture de valeurs inutilisees

      CHARACTER*10 CMODEL       ! Chaine de caractères précisant le nom d'un mode aérosols (LND ou JUNGE)    
 
      INTEGER*4 MEL_MODEL_TYPE(CTE_MAX_NB_MODE_MIXTURE)           ! Tableau listant les types de modèles d'un mélange (coef IGRANU)   
      DOUBLE PRECISION MEL_MODEL_PARAM(CTE_MAX_NB_MODE_MIXTURE,3) ! Tableau listant les paramètres des modèles d'un mélange    
      DOUBLE PRECISION MEL_MODEL_INDEX_WAREF(CTE_MAX_NB_MODE_MIXTURE,2)  ! Tableau listant les indices de réfraction des modèles d'un mélange à WAREF
      DOUBLE PRECISION MEL_MODEL_INDEX_WASIMU(CTE_MAX_NB_MODE_MIXTURE,2) ! Tableau listant les indices de réfraction des modèles d'un mélange à WASIMU
      DOUBLE PRECISION MEL_MODEL_AOT(CTE_MAX_NB_MODE_MIXTURE)     ! Tableau listant les AOT par modèle d'un mélange
      DOUBLE PRECISION AOT_TAUX					   ! Proportion à l'AOT totale du modèle d'un mélange
                                                                  ! à la longueur d'onde de référence.
      DOUBLE PRECISION SOM_TAUX_AOT				   ! Somme des taux des différents modes
      DOUBLE PRECISION COEF_ALPHA(CTE_MAX_NB_MODE_MIXTURE)        ! Tableau listant les coefficients de pondération alpha par modèle 
      DOUBLE PRECISION SOM_COEF_ALPHA                             ! Somme des coefficients avant normalisation
      DOUBLE PRECISION RN_WAREF_LU, IN_WAREF_LU			   ! Indice de réfraction (partie réelle et partie imaginaire)
                                                                  ! lu dans le fichier utilisateur de définition d'un mélange
                                                                  ! à la longueur d'onde de référence.
      DOUBLE PRECISION RN_WASIMU_LU, IN_WASIMU_LU		   ! Idem à la longueur d'onde de simulation des luminances.

      LOGICAL TRACE     ! = vrai,si écriture dans le fichier trace
      LOGICAL EX	   ! Verifie l'existence d'un fichier


C*-------------
C* Traitements                        
C*-------------    
 
 
      
C* Initialisation                     
C*-----------------------------------------------------------------
      IER=0

C --- Dénomination du fichier de MIE temporaire 
C     pour sa création dans l'espace de travail avant déplacement 
C     dans l'espace d'archivage.
      LDIRTMP = INDEX(DIR_TMP,'/TMP') +3
      FICMIE_TMP=DIR_TMP(1:LDIRTMP)//'/'//"MIE_FILE_BUILDING_TMP"   
      LFICTMP = INDEX(FICMIE_TMP,' ') - 1
	
C --- Affectation de la valeur minimale du paramètre de taille 
C     (passage de la constante en double précision)
      ALPHAO = CTE_MIE_ALPHAMIN      
      
      IF ((ALPHAO.LT.1D-4).OR.(ALPHAO.GE.10D-00)) GOTO 1010
     	
C* Ouverture du fichier Trace 
C------------------------------------------------------
      IF (FICLOGGRANU.EQ.'NO_LOG_FILE') THEN
	 TRACE=.FALSE.
      ELSE
        TRACE=.TRUE.
        OPEN (88,FILE=FICLOGGRANU,ERR=911)
        
        WRITE (88,*,ERR=921) ""
	 WRITE (88,*,ERR=921) ""
	 WRITE (88,*,ERR=921) "****************************"	 
	 WRITE (88,*,ERR=921) "SOS_AEROSOLS :  "
	 WRITE (88,*,ERR=921) "****************************"
	 WRITE (88,*,ERR=921) ""
	 WRITE (88,*,ERR=921) "Reference wavelength :", WAREF
	 WRITE (88,*,ERR=921) "Wavelength for radiance simulation:", WA
	 WRITE (88,*,ERR=921) ""
	 
      ENDIF
      
C* Lecture des tables d'angles et de poids utiles aux simulations
C* -----------------------------------------------------------------     
      OPEN(UNIT=10,FILE=FICANGLES,STATUS='OLD',ERR=930)

      READ(10,'(a)',ERR=931) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) MIE_NBMU
      
      READ(10,'(a)',ERR=931)  LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) NBMU_GAUSS
      
      READ(10,'(a)',ERR=931)  LIGNE_TEXTE
      !On ne conserve que la partie après le séparateur ':'
      IFIN=INDEX(LIGNE_TEXTE,':')	
      LIGNE_TEXTE=LIGNE_TEXTE(IFIN+1:LENLINE)	
      !On recherche la fin de ligne utile ' ' pour réduire la taille de la chaine
      IFIN=INDEX(LIGNE_TEXTE,' ')
      FICANGLES_USER=LIGNE_TEXTE(1:IFIN)
      
      READ(10,'(a)',ERR=931)  LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NB
      
      READ(10,*,ERR=931) CVAL
      DO J=1,MIE_NBMU
	  READ(10,*,ERR=931) CVAL, XMU(J), XHR(J)
      ENDDO 
      
      CLOSE(10)
      
      DO J=1,MIE_NBMU
         XMU(-J) = -XMU(J)
	 XHR(-J) =  XHR(J)
      ENDDO 
                
        
      IF (TRACE) THEN
        WRITE (88,*,ERR=921) "USED ANGLES"
	 WRITE (88,*,ERR=921) "****************"
	 WRITE (88,*,ERR=921) " "
	 WRITE (88,*,ERR=921) "Total number of angles : ", MIE_NBMU
	 WRITE (88,*,ERR=921) "Number of Gauss angles : ", NBMU_GAUSS
	 WRITE (88,*,ERR=921) " "
        IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	    WRITE(88,*,ERR=921)  
     &	       "A file of user's angles is used :",
     &         FICANGLES_USER
     	    WRITE(88,*,ERR=921) "Number of user's angles :",
     &         (MIE_NBMU-NBMU_GAUSS)
	 ELSE
	    WRITE(88,*,ERR=921)  "No file of user's angles"
	 ENDIF
        WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) "Max order NB : ", OS_NB
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,err=921) "List of used angles :"
	 WRITE (88,*,err=921) "  Index, mu, angle (degrees), weigth"
	 DO 1000 J=-MIE_NBMU,MIE_NBMU
	    IF (J.EQ.0) GOTO 1000
	    IF (J.LT.0) THEN 
	        WRITE(88,222,ERR=921) J, XMU(J),
     &	                     -DACOS(-XMU(J))*180./INCTE_PI, XHR(J)
           ELSE
	        WRITE(88,222,ERR=921) J, XMU(J),
     &	                      DACOS(XMU(J))*180./INCTE_PI, XHR(J)
	    ENDIF	    
1000    CONTINUE
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) " "  	
      ENDIF !Fichier Trace	 


C* Récupération du nom du répertoire contenant les fichiers de paramètres WMO 
C* et Shettle & Fenn
C----------------------------------------------------------------------------
      CALL GETENV('SOS_ABS_ROOT',DIRFIC)
      IFIN=INDEX(DIRFIC,' ')
      IFIN=IFIN-1
      IF (IFIN.LE.0) GOTO 925
      
      IF ((IFIN+4).GT.CTE_LENDIR) GOTO 926
      
      DIRFIC(IFIN+1:IFIN+4)='/fic'
      IFIN=IFIN+4
      IF (IFIN.LE.0) IFIN=CTE_LENDIR
            
C* Initialisation des proportions de constituants des modèles WMO et Shettle & Fenn
C----------------------------------------------------------------------------------      
      DO I=1,CTE_NB_WMO
	  C(I)=0.
	  N(I)=0.
      ENDDO
      
       DO I=1,CTE_NB_SF
	  Ni(I)=0.
      ENDDO

      NTOT=0.
            
C* Initialisation des tableaux de paramètres de la matrice de phase
C--------------------------------------------------------------------
      DO 101 K=0,OS_NB
         BETA11(K)=0.D+00
	  BETA22(K)=0.D+00
         GAMMA12(K)=0.D+00
         DELTA33(K)=0.D+00
         ALP(K)=0.D+00
         ZETA(K)=0.D+00
  101 CONTINUE

C* Initialisation des tableaux de fonction de phase
C--------------------------------------------------
      DO 102 J=-MIE_NBMU,MIE_NBMU
         P11(J)=0.
         P12(J)=0.
         P22(J)=0.
	  P33(J)=0.
         TTT(J)=0.
  102 CONTINUE

C* Initialisation des sections efficaces d'extinction et de diffusion
C--------------------------------------------------------------------
      KMAT1=0.
      KMAT2=0.
                        
C* Calcul des paramètres de la matrice de phase des aérosols
C* et de l'albédo de simple diffusion.
C------------------------------------------------------------
      IF (TA.EQ.0.D+00) THEN	!Pas de calcul sur la granulométrie
      
      	PIZ = 0.D+00
	GOTO 1111	!Enregistrement du fichier résultat GRANU
      
      ENDIF

      !Poursuite des traitements pour TA > 0
      !=====================================
      
      	
C* Longueur du nom du répertoire des fichiers de Mie
C---------------------------------------------------
      LDIRS = INDEX(DIRS,' ')
      LDIRS = LDIRS - 1
      IF(LDIRS.LE.0) LDIRS = 100
      

C ===================
C* Calcul mono-modal
C ===================

      IF(IMOD.EQ.0) THEN
	 	 	
      
C* Valeurs de RN et IN forcées au format F5.3 et F8.5
C  -----------------------------------------
	 RN = DNINT(RN*1000.D+00) / 1000.D+00
	 IN = -DNINT(-IN*100000.D+00) / 100000.D+00 	 


        IF (IGRANU.EQ.1) THEN	! Cas LND
	 
	     !!!!!IGRANU_MMD_CALL = 1	! Indice LND pour routine OSOAA_GRANU
	     
	     X = VARGRANU2_MMD   
	     RMAX = VARGRANU1_MMD * DEXP(X*X) *
     &	            DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))
	     
	 ELSE	! Cas loi Junge pour aérosols
	 
	     !!!!!IGRANU_MMD_CALL = 3	! Indice Junge pour aérosols pour OSOAA_GRANU

	     RMAX = VARGRANU3_MMD
	     
	 ENDIF ! Fin cas IGRANU_MMD = 2


	 ALPHAF 
     &	 =REAL(100+100*DINT(2.*INCTE_PI*RMAX / (100.*CTE_WAMIN)))

	 IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009 


        IF (TRACE) THEN
           WRITE(88,*,ERR=921) "MONO-MODAL MODEL"
	    WRITE(88,*,ERR=921) "****************"
	    WRITE(88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) "Wavelength in microns : ",WA
	    WRITE(88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) 
     & 	 "Refractive index (real part) : ",RN
	    WRITE(88,*,ERR=921) 
     & 	 "Refractive index (imaginary part) : ",IN
	    IF(IGRANU.EQ.1) THEN
               WRITE(88,*,ERR=921) "Size distribution : Log-normal"
	        WRITE(88,*,ERR=921) "  Modal radius (microns) : " , 
     &	        VARGRANU1_MMD 
	        WRITE(88,*,ERR=921) "  Variance (sig) : ", VARGRANU2_MMD
	        WRITE(88,*,ERR=921) "  sig / ln10  : " , 
     &         VARGRANU2_MMD  / DLOG(10.D+00)
               WRITE(88,*,ERR=921) 
     &	        "  Calculated Rmax limit (mic) : ", RMAX
	    ELSE
	       WRITE(88,*,ERR=921) "Size distribution : Junge's law"
              WRITE(88,*,ERR=921) "  Slope : " , VARGRANU2_MMD 
              WRITE(88,*,ERR=921) "  Minimal radius (microns) : ", 
     &	       VARGRANU1_MMD
              WRITE(88,*,ERR=921) "  Maximal radius (microns) : ", RMAX
	    ENDIF   

           WRITE(88,*,ERR=921) " "
           WRITE(88,*,ERR=921) "  AlphaF limit : ", ALPHAF
	    WRITE(88,*,ERR=921) " "	
	    WRITE(88,*,ERR=921) " "  		
         ENDIF !Fichier Trace	 
 

     	 
C* Constitution du nom du fichier MIE
C  ----------------------------------

	 CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &                      RN,IN,ALPHAO,ALPHAF,FICS) 	     

	    
	 FICMIE = DIRS(1:LDIRS)//'/'//FICS	 
	 LFICS=INDEX(FICMIE,' ')-1

	 IF (TRACE) THEN 
	    WRITE (88,*,ERR=921) "MIE file : ", 
     &	    FICMIE(1:LFICS)	
	 ENDIF	
	 
	 !Affichage ecran
	 WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)

	 
C* Si le fichier MIE existe déjà, le programme n'est pas lancé
C  -----------------------------------------------------------	 
        INQUIRE(FILE=FICMIE,EXIST=EX)
	 
        IF (EX) THEN  
            WRITE(6,*)'  ==> The file has already been calculated'
	  
	 ELSE   
      
C* Calcul du fichier de MIE
C  ------------------------ 
            CALL SOS_MIE(MIE_NBMU,XMU,XHR,RN,IN,ALPHAO,ALPHAF,
     &	                 FICMIE_TMP,FICLOGMIE,IER)
            IF (IER.NE.0) GOTO 993   ! cas d'erreur

            COMMAND='mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
            LCOM = LCOM + (4+LFICTMP)          
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICMIE,EXIST=EX)
            IF (EX) THEN !Si le fichier de MIE existe on détruit le fichier temporaire (généré par autre traitement)
                OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                CLOSE(3,STATUS='DELETE') 
            ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de MIE
               CALL SYSTEM(COMMAND)  
            ENDIF

        ENDIF
	 
C* Calcul des sections efficaces d'extinction et de diffusion 
C  et des fonctions de phase I, Q, U
C  ----------------------------------------------------------          
        CALL SOS_GRANU(FICMIE,IGRANU,VARGRANU1_MMD,
     &                 VARGRANU2_MMD,VARGRANU3_MMD,
     &                 WA,MIE_NBMU,XMU,TRACE,
     &                 KMAT1,KMAT2,SOMME_NR,
     &                 P11,P12,P33,IER)
        IF (IER.NE.0) GOTO 995   ! cas d'erreur 

C* Cas de particules sphériques : P22(J) = P11(J)
C  -----------------------------------------------------
        DO 1021 J=-MIE_NBMU,MIE_NBMU
           P22(J)=P11(J)
 1021   CONTINUE	 
  
C* Calcul des paramètres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				      P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				      ALP,BETA11,BETA22,
     &                            GAMMA12,DELTA33,ZETA,IER)  
        IF (IER.NE.0) GOTO 996   ! cas d'erreur

      ENDIF	!Fin test sur Modalite IMOD=0 (Mono-modal)



C ====================================================
C* Calcul multi-modal avec les constituants de la WMO
C ====================================================
	    
      IF(IMOD.EQ.1) THEN

	 
C* Initialisation
C  --------------
        ! Nom des différentes particules de la WMO
	 PARTICULE(1)='DUST LIKE'
        PARTICULE(2)='WATER SOLUBLE'
        PARTICULE(3)='OCEANIC'
        PARTICULE(4)='SOOT'

        ALPHAF_PARTICULE(1) = CTE_ALPHAMAX_WMO_DL
	 ALPHAF_PARTICULE(2) = CTE_ALPHAMAX_WMO_WS
	 ALPHAF_PARTICULE(3) = CTE_ALPHAMAX_WMO_OC
	 ALPHAF_PARTICULE(4) = CTE_ALPHAMAX_WMO_SO
	 
	 IGRANU =   1	! distribution Log Normal par defaut
	    
C* Initialisation des paramètres de la WMO
C  ---------------------------------------
	 FICDATAWMO=DIRFIC(1:IFIN)//'/'//CTE_AER_DATAWMO

	 CALL SOS_INIT_PARAMWMO(FICDATAWMO,WA,VARGRANU1I,
     &                          VARGRANU2I,MR,MI,V,IER)
        IF (IER.NE.0) GOTO 994   ! cas d'erreur 

C* Definition de la proportion volumique de chaque constituant
C  pour le modèle selectionné
C  -----------------------------------------------------------	 
        IF(IMODELE_WMO.EQ.1) THEN	!Modèle Continental
	    C(1) = 0.70
	    C(2) = 0.29
	    C(4) = 0.01      
	    
        ELSE IF(IMODELE_WMO.EQ.2) THEN 	!Modèle Maritime
	    C(2) = 0.05
	    C(3) = 0.95
	 
        ELSE IF(IMODELE_WMO.EQ.3) THEN	!Modèle Urbain
	    C(1) = 0.17
	    C(2) = 0.61
	    C(4) = 0.22
	       
        ELSE IF(IMODELE_WMO.EQ.4) THEN	!Modèle défini par l'utilisateur
	 
C*          Controle de la récupération des proportions volumiques 
C*          definies par l'utilisateur
C           ------------------------------------------------------
	    
	    C(1) = C_WMO_DL    !Pourcentage volumique des particules DL
	    C(2) = C_WMO_WS    !Pourcentage volumique des particules WS
	    C(3) = C_WMO_OC    !Pourcentage volumique des particules OC
	    C(4) = C_WMO_SO    !Pourcentage volumique des particules SO	       	    

        ENDIF	!Fin définition des modèles
	    
        DO 103 I=1,CTE_NB_WMO	    
            N(I) = C(I)/V(I)
	     NTOT = NTOT + N(I)
  103   CONTINUE


        IF (TRACE) THEN
          WRITE (88,*,ERR=921) "WMO MODEL"
	   WRITE (88,*,ERR=921) "*********"
	   WRITE (88,*,ERR=921) " "
	   IF(IMODELE_WMO.EQ.1) 
     &	      WRITE(88,*,ERR=921) "CONTINENTAL Model"
	   IF(IMODELE_WMO.EQ.2) WRITE(88,*,ERR=921) "MARITIME Model"
	   IF(IMODELE_WMO.EQ.3) WRITE(88,*,ERR=921) "URBAN Model"
	   WRITE (88,*,ERR=921) " "
	   WRITE (88,*,ERR=921) "  Coefficients alphaK of the mixture:"
	   DO I=1,CTE_NB_WMO
	      WRITE (88,*,ERR=921) "    ",PARTICULE(I)," : ", N(I)/NTOT
          ENDDO
	   WRITE (88,*,ERR=921) " "
	   IF(IMODELE_WMO.EQ.4) THEN
	     WRITE (88,*,ERR=921) "USER Model "
	     DO I=1,4
	        WRITE (88,*,ERR=921) "    ",PARTICULE(I)," : ",
     &		                     100*C(I)," % volume"
            ENDDO
	   ENDIF 
	   WRITE (88,*,err=921) "Wavelength in microns : ",WA
	   WRITE (88,*,ERR=921) " " 
	   WRITE (88,*,ERR=921) " "  		
        ENDIF !Fichier Trace
      
      

C* Pour chaque type d'aérosols de la WMO
C  -------------------------------------	       
        DO 105 I=1,CTE_NB_WMO
	    
	    IF(C(I).EQ.0.) GOTO 105
		
C* Constitution du nom du fichier MIE
C  ----------------------------------
	    ALPHAF = ALPHAF_PARTICULE(I)

           CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                        MR(I),MI(I),ALPHAO,ALPHAF,FICS)       
	    FICMIE = DIRS(1:LDIRS)//'/'//FICS
           LFICS=INDEX(FICMIE,' ')-1
	    
	    IF (TRACE) THEN 
	       WRITE (88,*,ERR=921) "Fichier MIE : ", 
     &	       FICMIE(1:LFICS)	
	    ENDIF	
	 
	    !Affichage ecran
	    WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
C* Si le fichier MIE existe déjà, le programme n'est pas lance
C  -----------------------------------------------------------		              
           INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	    IF (EX) THEN  
               WRITE(6,*)'  ==> The file has already been calculated'
	  
	    ELSE   
      
C* Calcul du fichier de MIE
C  ------------------------ 
              CALL SOS_MIE(MIE_NBMU,XMU,XHR,MR(I),MI(I),ALPHAO,ALPHAF,
     &                      FICMIE_TMP,FICLOGMIE,IER)
    	       IF (IER.NE.0) GOTO 993   ! cas d'erreur
               
              COMMAND=
     &        'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
              LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
              LCOM = LCOM + (4+LFICTMP)
              IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020
              
              INQUIRE(FILE=FICMIE,EXIST=EX)
              IF (EX) THEN !Si le fichier de MIE existe on détruit le fichier temporaire (généré par autre traitement)
                 OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                 CLOSE(3,STATUS='DELETE') 
              ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de MIE
                 CALL SYSTEM(COMMAND)  
              ENDIF
   
           ENDIF
	 
C* Calcul des sections efficaces d'extinction et de diffusion 
C  et des fonctions de phase I, Q, U du constituant I
C  ---------------------------------------------------------- 
	    IF (TRACE) THEN
	       WRITE(88,*,err=921) PARTICULE(I)
	       WRITE(88,*,err=921)'============='
	       WRITE(88,*,err=921) ' '
	    ENDIF
	       
           CALL SOS_GRANU(FICMIE,IGRANU,VARGRANU1I(I),
     &                    VARGRANU2I(I),CTE_NOT_DEFINED_VALUE_DBLE,
     &                    WA,MIE_NBMU,XMU,TRACE,
     &			     KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
           IF (IER.NE.0) GOTO 995   ! cas d'erreur  
	    

C* Calcul des sections efficaces d'extinction et de diffusion du mélange
C  ---------------------------------------------------------------------	     
	    KMAT1 = KMAT1 + (N(I)/NTOT)*KMAT1I
	    KMAT2 = KMAT2 + (N(I)/NTOT)*KMAT2I
	    
	    DO 104 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) + (N(I)/NTOT)*P11I(J)*KMAT2I
	       P12(J) = P12(J) + (N(I)/NTOT)*P12I(J)*KMAT2I
	       P33(J) = P33(J) + (N(I)/NTOT)*P33I(J)*KMAT2I
  104	    CONTINUE  !Fin boucle sur les angles de gauss
    
  105	 CONTINUE	 !Fin boucle sur les constituants de la WMO


C* Normalisation des fonctions de phase du mélange
C  -----------------------------------------------	    
	 DO 106 J=-MIE_NBMU,MIE_NBMU   
	    P11(J) = P11(J) / KMAT2
	    P12(J) = P12(J) / KMAT2
	    P33(J) = P33(J) / KMAT2
  106	 CONTINUE

C* Cas de particules sphériques : P22(J) = P11(J)
C  -----------------------------------------------------
        DO 1022 J=-MIE_NBMU,MIE_NBMU
           P22(J)=P11(J)
 1022   CONTINUE	 
 
C* Calcul des paramètres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------  
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				      P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				      ALP,BETA11,BETA22,
     &                            GAMMA12,DELTA33,ZETA,IER)  
        IF (IER.NE.0) GOTO 996   ! cas d'erreur


      ENDIF	!Fin test sur Modalite IMOD=1 (WMO) 
  
       
C =====================================================================
C* Calcul bi-modal avec les constituants du modèle Shettle & Fenn (S&F)
C =====================================================================
      
      IF(IMOD.EQ.2) THEN
      	 
C* Initialisation
C  --------------
         ! Nom des différentes particules du modèle S&F
	  PARTICULE(1)='SMALL_RURAL'
         PARTICULE(2)='LARGE RURAL'
         PARTICULE(3)='SMALL_URBAN'
         PARTICULE(4)='LARGE_URBAN'
	  PARTICULE(5)='OCEANIC_MODEL'

	  IGRANU =   1	! distribution Log Normal par defaut
	    
C* Initialisation des paramètres du modèle S&F
C  -------------------------------------------
	  CALL SOS_INIT_PARAMSF(DIRFIC,WA,RH,VARGRANU1I,
     &                         VARGRANU2I,MR,MI,IER)
         IF (IER.NE.0) GOTO 997   ! cas d'erreur 

C* Definition de la densite de chaque constituant pour le 
C  modèle selectionne
C  -----------------------------------------------------------	 
         IF(IMODELE_SF.EQ.1) THEN		!Modèle Tropospherique
	    Ni(1) = 1.0     
	    
         ELSE IF(IMODELE_SF.EQ.2) THEN 	!Modèle Urbain
	    Ni(3) = 0.999875
	    Ni(4) = 0.000125
	 
         ELSE IF(IMODELE_SF.EQ.3) THEN	!Modèle Maritime
	    Ni(1) = 0.99
	    Ni(5) = 0.01
	       
         ELSE IF(IMODELE_SF.EQ.4) THEN	!Modèle cotier
	    Ni(1) = 0.995
	    Ni(5) = 0.005
	       
         ENDIF	!Fin definition des modèles S&F


      IF (TRACE) THEN
        WRITE (88,*,ERR=921) "SHETTLE & FENN MODEL"
	 WRITE (88,*,ERR=921) "********************"
	 WRITE (88,*,ERR=921) " "
	 IF(IMODELE_SF.EQ.1) WRITE(88,*,ERR=921) "TROPOSPHERIC Model"
	 IF(IMODELE_SF.EQ.2) WRITE(88,*,ERR=921) "URBAN Model"
	 IF(IMODELE_SF.EQ.3) WRITE(88,*,ERR=921) "MARITIME Model"
	 IF(IMODELE_SF.EQ.4) WRITE(88,*,ERR=921) "COASTAL Model"
	 WRITE (88,*,ERR=921) " "
	 WRITE (88,*,ERR=921) "  Coefficients alphaK of the mixture:"
	 DO I=1,CTE_NB_SF
	    WRITE (88,*,ERR=921) "    ",PARTICULE(I)," : ", Ni(I)
        ENDDO
        WRITE (88,*,ERR=921) " "
	 WRITE (88,*,err=921) "Percentage of relative humidity : ",RH
	 WRITE (88,*,err=921) "Wavelength in microns : ",WA 
	 WRITE (88,*,ERR=921) " "   
	 WRITE (88,*,ERR=921) " "		
      ENDIF !Fichier Trace

      
C* Pour chaque type d'aérosols du modèle S&F
C  ----------------------------------------     
	 DO 109 I=1,CTE_NB_SF
	    
	    IF(Ni(I).EQ.0.) GOTO 109
		
C* Constitution du nom du fichier MIE
C  ----------------------------------	   
           IF (I.EQ.1) ALPHAF = CTE_ALPHAMAX_SF_SR
	    IF (I.EQ.3) ALPHAF = CTE_ALPHAMAX_SF_SU
	    IF ((I.NE.1).AND.(I.NE.3)) THEN
	        X = VARGRANU2I(I)  
	        RMAX = VARGRANU1I(I) * DEXP(X*X) *
     &	               DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))
     	
	        ALPHAF = REAL( 100 +
     &		100*DINT(2.*INCTE_PI*RMAX / (100.*WA)) )
	    ENDIF

	    IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009 

           CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                        MR(I),MI(I),ALPHAO,ALPHAF,FICS)       
	    FICMIE = DIRS(1:LDIRS)//'/'//FICS
           LFICS=INDEX(FICMIE,' ')-1
	    
	    IF (TRACE) THEN 
	       WRITE (88,*,ERR=921) "MIE file : ", 
     &	       FICMIE(1:LFICS)	
	    ENDIF	
	 
	    !Affichage ecran
	    WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
C* Si le fichier MIE existe déjà, le programme n'est pas lancé
C  -----------------------------------------------------------           	              
           INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	    IF (EX) THEN  
               WRITE(6,*)'  ==> The file has already been calculated'
	  
	    ELSE   
      
C* Calcul du fichier de MIE
C  ------------------------ 
              CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &	                    MR(I),MI(I),ALPHAO,ALPHAF,
     &                     FICMIE_TMP,FICLOGMIE,IER)
     	       IF (IER.NE.0) GOTO 993   ! cas d'erreur
       
              COMMAND=
     &        'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
              LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
              LCOM = LCOM + (4+LFICTMP)
              IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

              INQUIRE(FILE=FICMIE,EXIST=EX)
              IF (EX) THEN !Si le fichier de MIE existe on détruit le fichier temporaire (généré par autre traitement)
                 OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                 CLOSE(3,STATUS='DELETE') 
              ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de MIE
                 CALL SYSTEM(COMMAND)  
              ENDIF

            ENDIF
	 
C* Calcul des sections efficaces d'extinction et de diffusion 
C  et des fonctions de phase I, Q, U du constituant I
C  ---------------------------------------------------------- 
	    IF (TRACE) THEN
	       WRITE(88,*,err=921) PARTICULE(I)
	       WRITE(88,*,err=921)'============='
	       WRITE(88,*,err=921) ' '
	       IF ((I.NE.1).AND.(I.NE.3)) THEN
	         WRITE(88,*,err=921) 'RMAX (microns) :',RMAX
		  WRITE(88,*,err=921) 'Wavelength (microns) : ',WA
		  WRITE(88,*,err=921) 'Max size parameter : ',ALPHAF
	       ENDIF
	    ENDIF
	       
           CALL SOS_GRANU(FICMIE,IGRANU,VARGRANU1I(I),
     &                    VARGRANU2I(I),CTE_NOT_DEFINED_VALUE_DBLE,
     &                    WA,MIE_NBMU,XMU,TRACE,
     &			     KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
           IF (IER.NE.0) GOTO 995   ! cas d'erreur  
	     
C* Calcul des sections efficaces d'extinction et de diffusion du mélange
C  ---------------------------------------------------------------------
	    KMAT1 = KMAT1 + Ni(I)*KMAT1I
	    KMAT2 = KMAT2 + Ni(I)*KMAT2I
	    
	    DO 108 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) + Ni(I)*P11I(J)*KMAT2I
	       P12(J) = P12(J) + Ni(I)*P12I(J)*KMAT2I
	       P33(J) = P33(J) + Ni(I)*P33I(J)*KMAT2I
  108	    CONTINUE  !Fin boucle sur les angles de gauss
    
  109	 CONTINUE	 !Fin boucle sur les constituants du modèle Shettle et Fenn
	    

C* Normalisation des fonctions de phase du mélange
C  -----------------------------------------------
	 DO 110 J=-MIE_NBMU,MIE_NBMU   
	    P11(J) = P11(J) / KMAT2
	    P12(J) = P12(J) / KMAT2
	    P33(J) = P33(J) / KMAT2
  110	 CONTINUE

C* Cas de particules sphériques : P22(J) = P11(J)
C  -----------------------------------------------------
        DO 1023 J=-MIE_NBMU,MIE_NBMU
           P22(J)=P11(J)
 1023   CONTINUE	 
 
C* Calcul des paramètres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------  
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				      P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				      ALP,BETA11,BETA22,
     &                            GAMMA12,DELTA33,ZETA,IER)  
        IF (IER.NE.0) GOTO 996   ! cas d'erreur	 

      	 
      ENDIF	!Fin test sur Modalite IMOD=2 (Shettle&Fenn) 



C =====================================================================
C* Calcul des proprietes radiatives d'une granulométrie LND bimodale
C =====================================================================
      
      IF(IMOD.EQ.3) THEN
	
      
C*       Récupération de l'information de la description du mélange
C        ----------------------------------------------------------      
C        Utilisation des concentrations volumiques de l'utilisateur
C         --> Controle que les paramètres utilisateurs sont renseignés pour l'option
         IF (MODE_PARAM_BILND.EQ.1) THEN 
	         
	     CVI(1) = USER_CV_COARSE !mode "coarse"		 
	     CVI(2) = USER_CV_FINE   !mode "fine"
	  
         ENDIF
       	    
C*       Récupération de l'information d'indice et granulométrie des 2 LND
C        -----------------------------------------------------------------         	 
      
         RN_LND(1) = BMD_CM_MRWA  !mode "coarse"
         RN_LND(2) = BMD_FM_MRWA  !mode "fine"
      
         IN_LND(1) = BMD_CM_MIWA  !mode "coarse"
         IN_LND(2) = BMD_FM_MIWA  !mode "fine"
      
         RMODAL_LND(1) = BMD_CM_RMODAL  !mode "coarse"
         RMODAL_LND(2) = BMD_FM_RMODAL  !mode "fine"
      
         VAR_LND(1) = BMD_CM_VAR  !mode "coarse"
         VAR_LND(2) = BMD_FM_VAR  !mode "fine"
      
      
         IF (MODE_PARAM_BILND.EQ.2) THEN 
	  
	    RN_LND_WAREF(1) = BMD_CM_MRWAREF !mode "coarse"
	    RN_LND_WAREF(2) = BMD_FM_MRWAREF !mode "fine"
	  
	    IN_LND_WAREF(1) = BMD_CM_MIWAREF !mode "coarse"
	    IN_LND_WAREF(2) = BMD_FM_MIWAREF !mode "fine"
	    
         ENDIF
 
         

         DO 111 I=1,2
	 
C* Controle le signe de l'argument IN
C  -----------------------------------------
	    IF (IN_LND(I).GT.0.D+00) GOTO 10100	 

C* Valeurs de RN et IN forcées au format F5.3 et F8.5
C  -----------------------------------------
	    RN_LND(I) = DNINT(RN_LND(I)*1000.D+00) / 1000.D+00
	    IN_LND(I) = -DNINT(-IN_LND(I)*100000.D+00) / 100000.D+00 
	    
111      CONTINUE	 
	 

         IF (TRACE) THEN
           WRITE (88,*,ERR=921) "LOG-NORMAL BIMODAL MODEL"
	    WRITE (88,*,ERR=921) "************************"
	    WRITE (88,*,ERR=921) " "
	    WRITE (88,*,ERR=921) "Wavelength in microns : ", WA
	    WRITE (88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) "Size distribution : log-normal coarse"
	    WRITE(88,*,ERR=921) "                    -----------------"
	    WRITE(88,*,ERR=921) 
     &     "   Refractive index (real part) : ",RN_LND(1)
	    WRITE(88,*,ERR=921) 
     &     "   Refractive index (imaginary part) : ",IN_LND(1)
	    WRITE(88,*,ERR=921) 
	    WRITE(88,*,ERR=921) "   Modal radius (microns) : " , 
     &	    RMODAL_LND(1) 
	    WRITE(88,*,ERR=921) "   Standard deviation : " , VAR_LND(1)
	    WRITE(88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) "Size distribution : log-normal fine"
	    WRITE(88,*,ERR=921) "                    ---------------"
	    WRITE(88,*,ERR=921) 
     &     "   Refractive index (real part) : ",RN_LND(2)
	    WRITE(88,*,ERR=921) 
     &     "   Refractive index (imaginary part) : ",IN_LND(2)
	    WRITE(88,*,ERR=921) 
	    WRITE(88,*,ERR=921) "   Modal radius (microns) : " , 
     &	    RMODAL_LND(2) 
	    WRITE(88,*,ERR=921) "   Standard deviation : " , VAR_LND(2) 
	    WRITE(88,*,ERR=921) " "
	 
	    IF (MODE_PARAM_BILND.EQ.1) THEN 
	       WRITE(88,*,ERR=921)  
     &         "Proportion of components defined by", 
     &	       " volume concentrations :"		    
	       WRITE(88,*,ERR=921) 
     & 	"   Volume concentration for LND coarse mode : ",CVI(1)
	       WRITE(88,*,ERR=921) 
     & 	"   Volume concentration for LND fine mode   : ",CVI(2)
	       WRITE(88,*,ERR=921) " "
	    ENDIF 
	  
	    IF (MODE_PARAM_BILND.EQ.2) THEN 
	       WRITE(88,*,ERR=921) 
     &         "Proportion of components defined by the", 
     &	       "by the ratio AOT_coarse / AOT_totale :"
     	       WRITE(88,*,ERR=921) 
     &	       "   Reference wavelength (microns) : ",WAREF
     	       WRITE(88,*,ERR=921) 
     &	       "   --> Total aerosols extinction optical thickness : ",
     &         TA_WAREF	    
	       WRITE(88,*,ERR=921) 
     &         "   --> ratio AOT coarse / totale : ",RTAUct_WAREF
              WRITE(88,*,ERR=921) " "
	       WRITE(88,*,ERR=921) 
     &	       "   --> Refractive index of the coarse mode "
	       WRITE(88,*,ERR=921) 
     & 	       "          - Real part : ",RN_LND_WAREF(1)
	       WRITE(88,*,ERR=921) 
     & 	       "          - Imaginary part : ",IN_LND_WAREF(1)
	       WRITE(88,*,ERR=921) 
     &	       "   --> Refractive index of the fine mode "
	       WRITE(88,*,ERR=921) 
     & 	       "          - Real part : ",RN_LND_WAREF(2)
	       WRITE(88,*,ERR=921) 
     & 	       "          - Imaginary part : ",IN_LND_WAREF(2)
	       WRITE(88,*,ERR=921) " "
	    ENDIF  	        	    
	 
         ENDIF !Fichier Trace	
      	    
	 
 
C* Initialisation
C  --------------
 	  IGRANU =   1	! distribution Log Normal par defaut


C* Calcul des coefficients CVI pour une expression du
C  mélange en AOT_coarse / AOT_totale
C  --------------------------------------------------
         IF (MODE_PARAM_BILND.EQ.2) THEN
	    !Utilisation du rapport AOT_coarse / AOT_totale
 
	    IF (TRACE) THEN
	       WRITE (88,*,ERR=921) " "
	       WRITE (88,*,ERR=921) 
     &	       "Calculation of the mixture coefficients"
	       WRITE (88,*,ERR=921) 
     &	       "======================================="
               WRITE (88,*,ERR=921) "  "
	    ENDIF
	    
C* Pour chacun des deux modes (1 : coarse, 2 : fine)
C  -------------------------------------------------  
           DO 114 I=1,2  
		
C* Constitution du nom du fichier MIE (pour WAREF)
C  -----------------------------------------------	               
	       X = VAR_LND(I)  
	       RMAX = RMODAL_LND(I) * DEXP(X*X) *
     &	              DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))

              ALPHAF = REAL( 100 +
     &	       	 100*DINT(2.*INCTE_PI*RMAX / (100.*WAREF)) )  

	       IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009 	      
        
              CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                           RN_LND_WAREF(I),IN_LND_WAREF(I),
     &	                           ALPHAO,ALPHAF,FICS)       
	       FICMIE = DIRS(1:LDIRS)//'/'//FICS
              LFICS=INDEX(FICMIE,' ')-1
	       
	       IF (TRACE) THEN 
	          WRITE (88,*,ERR=921) "MIE file: ", 
     &	          FICMIE(1:LFICS)	
	       ENDIF	
	 
	       !Affichage ecran
	       WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
C* Si le fichier MIE existe déjà, le programme n'est pas lancé
C  -----------------------------------------------------------           	              
              INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	       IF (EX) THEN  
                  WRITE(6,*)'  ==> The file has already been calculated'
	  
	       ELSE   
      
C* Calcul du fichier de MIE (pour WAREF)
C  ------------------------------------ 
                 CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &		                RN_LND_WAREF(I),IN_LND_WAREF(I),
     &                        ALPHAO,ALPHAF,
     &                        FICMIE_TMP,FICLOGMIE,IER)
     	          IF (IER.NE.0) GOTO 993   ! cas d'erreur

                 COMMAND=
     &           'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
                 LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
                 LCOM = LCOM + (4+LFICTMP)
                 IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

                 INQUIRE(FILE=FICMIE,EXIST=EX)
                 IF (EX) THEN !Si le fichier de MIE existe on détruit le fichier temporaire (généré par autre traitement)
                    OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                    CLOSE(3,STATUS='DELETE') 
                 ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de MIE
                    CALL SYSTEM(COMMAND)  
                 ENDIF
         
              ENDIF
	 
C* Calcul des sections efficaces d'extinction et de diffusion 
C  et des fonctions de phase I, Q, U du constituant I pour WAREF
C  ------------------------------------------------------------- 
            
	       IF (TRACE) THEN
	          IF (I.EQ.1) THEN
	             WRITE(88,*,err=921) "   LND coarse mode"
		      WRITE(88,*,err=921) '   ==============='
	          ELSE !I=2 : fine mode
	             WRITE(88,*,err=921) "   LND fine mode" 
		      WRITE(88,*,err=921) '   =============' 
	          ENDIF	  
	          WRITE(88,*,err=921) ' '
	          WRITE(88,*,err=921) '      RMAX (microns) :',RMAX
	          WRITE(88,*,err=921) 
     &	          '      Wavelength (microns) : ',WAREF
	          WRITE(88,*,err=921) 
     &	          '      Max size parameter : ',ALPHAF	  
	          WRITE(88,*,err=921) ' '     
	       ENDIF
	       
              CALL SOS_GRANU(FICMIE,IGRANU,RMODAL_LND(I),
     &                    VAR_LND(I),CTE_NOT_DEFINED_VALUE_DBLE,
     &                    WAREF,MIE_NBMU,XMU,TRACE,
     &			     KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
              IF (IER.NE.0) GOTO 995   ! cas d'erreur  

    	     
C* Stockage des coefficients d'extinction de chaque mode
C  ---------------------------------------------------------------
	       IF (I.EQ.1)  KMAT1c = KMAT1I	!mode coarse
	       IF (I.EQ.2)  KMAT1f = KMAT1I	!mode fine
      	    	        
  114	    CONTINUE	 !Fin boucle sur les mode coarse et fine de la LND bimodale


C* Calcul des coefficients CVI de chaque mode
C  ---------------------------------------------------------------
           CVI(1)=   (RTAUct_WAREF * TA_WAREF) / KMAT1c
	    CVI(2)=   ( (1.-RTAUct_WAREF)*TA_WAREF) / KMAT1f

	    IF (TRACE) THEN
	       WRITE(88,*,err=921) 
     &	       '      Estimate : CV coarse mode (non-normalized):',
     &        CVI(1)
	       WRITE(88,*,err=921) 
     &	       '      Estimate : CV fine mode (non-normalized):',CVI(2)
	       WRITE(88,*,err=921) ' '  
	       WRITE(88,*,err=921) ' '   
	    ENDIF
			    
	  ENDIF !Fin test si utilisation du rapport AOT_coarse / AOT_tot

		
C* Normalisation à 1 de la somme des coefficients CVI de chaque mode
C  -----------------------------------------------------------------
         NTOT = CVI(1) + CVI(2)
	  CVI(1) = CVI(1) / NTOT
         CVI(2) = CVI(2) / NTOT
	
         IF (TRACE) THEN
	      WRITE(88,*,err=921)
     &	      '      CV coarse mode (normalized):', CVI(1)
	      WRITE(88,*,err=921)
     &	      '      CV fine mode (normalized):', CVI(2)	  
	      WRITE(88,*,err=921) ' '  
	      WRITE(88,*,err=921) ' '   
	     WRITE (88,*,ERR=921) 
     &	     "Calculation of the radiative properties of the mixture"
	     WRITE (88,*,ERR=921) 
     &	     "======================================================"
            WRITE (88,*,ERR=921) "  "
	  ENDIF
	    
C* Pour chacun des deux modes (1 : coarse, 2 : fine)
C  ------------------------------------------------- 
	        
	  DO 112 I=1,2
	    	    
	    IF(CVI(I).EQ.0.) GOTO 112
		
C* Constitution du nom du fichier MIE (pour la longueur d'onde de simulation)
C  ----------------------------------	               
	    X = VAR_LND(I)   
	    RMAX = RMODAL_LND(I) * DEXP(X*X) *
     &	           DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))
     	
	    ALPHAF = REAL( 100 +
     &	      100*DINT(2.*INCTE_PI*RMAX / (100.*WA)) )

	    IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009 

           CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     & 	                 RN_LND(I),IN_LND(I),ALPHAO,ALPHAF,FICS)     
	    FICMIE = DIRS(1:LDIRS)//'/'//FICS
           LFICS=INDEX(FICMIE,' ')-1
	    
	    IF (TRACE) THEN 
	        WRITE (88,*,ERR=921) "MIE file: ", 
     &	        FICMIE(1:LFICS)	
	    ENDIF	
	 
	    !Affichage ecran
	    WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
C* Si le fichier MIE existe déjà, le programme n'est pas lancé
C  -----------------------------------------------------------           	              
           INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	    IF (EX) THEN  
               WRITE(6,*)'  ==> The file has already been calculated'
	  
	    ELSE   
      
C* Calcul du fichier de MIE  (pour la longueur d'onde de simulation)
C  ------------------------ 
              CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &                     RN_LND(I),IN_LND(I),ALPHAO,ALPHAF,
     &                     FICMIE_TMP,FICLOGMIE,IER)
     	       IF (IER.NE.0) GOTO 993   ! cas d'erreur
 
              COMMAND=
     &        'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
              LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !On positionne la recherche après le dernier espace vide du début de chaine
              LCOM = LCOM + (4+LFICTMP)
              IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

              INQUIRE(FILE=FICMIE,EXIST=EX)
              IF (EX) THEN !Si le fichier de MIE existe on détruit le fichier temporaire (généré par autre traitement)
                 OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                 CLOSE(3,STATUS='DELETE') 
              ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de MIE
                 CALL SYSTEM(COMMAND)  
              ENDIF
        
           ENDIF
	 
C* Calcul des sections efficaces d'extinction et de diffusion 
C  et des fonctions de phase I, Q, U du constituant I
C   (pour la longueur d'onde de simulation)
C  ---------------------------------------------------------- 
	    IF (TRACE) THEN
	       IF (I.EQ.1) THEN
	          WRITE(88,*,err=921) "   1) LND coarse mode"
		   WRITE(88,*,err=921) '   =================='
	       ELSE !I=2 : fine mode
	          WRITE(88,*,err=921) "   2) LND fine mode" 
		   WRITE(88,*,err=921) '   ================' 
	       ENDIF	  
	       WRITE(88,*,err=921) ' '
	       WRITE(88,*,err=921) '      RMAX (microns) :',RMAX
	       WRITE(88,*,err=921) 
     &	       '      Wavelength (microns) : ',WA
	       WRITE(88,*,err=921) 
     &	       '      Max size parameter : ',ALPHAF	  
	       WRITE(88,*,err=921) ' '     
	    ENDIF
	       
           CALL SOS_GRANU(FICMIE,IGRANU,RMODAL_LND(I),
     &                    VAR_LND(I),CTE_NOT_DEFINED_VALUE_DBLE,
     &                    WA,MIE_NBMU,XMU,TRACE,
     &			     KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
           IF (IER.NE.0) GOTO 995   ! cas d'erreur  

		 
C* Calcul des coefficients d'extinction et de diffusion du mélange
C  ---------------------------------------------------------------
	    KMAT1 = KMAT1 + CVI(I)*KMAT1I
	    KMAT2 = KMAT2 + CVI(I)*KMAT2I
	    
	    DO 113 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) + CVI(I)*P11I(J)*KMAT2I
	       P12(J) = P12(J) + CVI(I)*P12I(J)*KMAT2I
	       P33(J) = P33(J) + CVI(I)*P33I(J)*KMAT2I	       
  113	    CONTINUE  !Fin boucle sur les angles de gauss
    
  112	  CONTINUE	 !Fin boucle sur les mode coarse et fine de la LND bimodale
	    
	 		
		
C* Normalisation des fonctions de phase du mélange
C  -----------------------------------------------
	DO 120 J=-MIE_NBMU,MIE_NBMU   
	    P11(J) = P11(J) / KMAT2
	    P12(J) = P12(J) / KMAT2
	    P33(J) = P33(J) / KMAT2
  120	CONTINUE

C* Cas de particules sphériques : P22(J) = P11(J)
C  -----------------------------------------------------
       DO 1024 J=-MIE_NBMU,MIE_NBMU
          P22(J)=P11(J)
 1024  CONTINUE	 
 
C* Calcul des paramètres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------  
	CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				     P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				     ALP,BETA11,BETA22,
     &                           GAMMA12,DELTA33,ZETA,IER)     
       IF (IER.NE.0) GOTO 996   ! cas d'erreur	 

          	 
      ENDIF	!Fin test sur Modalite IMOD=3 (LND bimodale) 




C ==============================================================================
C* Données externes :
C* --> Fonctions de phase externes et coefficients de diffusion et d'extinction
C*
C*     NB : Développé sur la base des sorties de l'outil Oleg Dubovik (LOA)
C*          qui fournit les propriétés de diffusion primaire de particules
C*          spheroides homogenes orientees de facon aleatoire.
C*          Un paramètre précisant le nombre d'angles decrivant les fonctions
C*          à toutefois du être ajouté.
C*
C*          Bien faire attention au fait qu'on lit un rapport -F12 / F11.
C*          Importance du signe : F12 < 0 pour la diffusion rayleigh.

C ==============================================================================

      IF(IMOD.EQ.4) THEN
	 
    	
C* Lecture du fichier de données externes
C  -----------------------------------------------------------------
         OPEN(UNIT=20,FILE=FICEXTDATA,STATUS='OLD',ERR=940)

	  READ(20,'(a)',ERR=941) LIGNE_TEXTE
         IFIN=INDEX(LIGNE_TEXTE,':')
         READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) KMAT1
	 
         READ(20,'(a)',ERR=941) LIGNE_TEXTE
         IFIN=INDEX(LIGNE_TEXTE,':')
         READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) KMAT2
	 
	  READ(20,'(a)',ERR=941) LIGNE_TEXTE
         IFIN=INDEX(LIGNE_TEXTE,':')
         READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) NB_FICEXT_ANGLES



	  !Controle du nombre d'angles par rapport a la taille des tableaux
	  IF (NB_FICEXT_ANGLES.GT.CTE_MAXNB_ANG_EXT) GOTO 950
	
	  READ(20,*,ERR=941) 

	  DO J=1,NB_FICEXT_ANGLES
            READ(20,*,ERR=941,END=942) 
     &	                 ANG_EXT(J), F11_EXT(J), F12sF11_EXT(J),
     &	                 F22sF11_EXT(J), F33sF11_EXT(J)
         ENDDO 
	
	  CLOSE(20)
	 
	  IF (TRACE) THEN
            WRITE (88,*,ERR=921) "EXTERNAL DATA (USER'S DATA)"
	     WRITE (88,*,ERR=921) "***************************"
	     WRITE (88,*,ERR=921) " "
	     WRITE (88,*,ERR=921)  "Filename : ",FICEXTDATA
	     WRITE (88,*,ERR=921) " "  		 
            WRITE (88,*,ERR=921) "Extinction coef : ", KMAT1
	     WRITE (88,*,ERR=921) "Scattering coef  : ", KMAT2
	     WRITE (88,*,ERR=921) "Nb of angles : ", NB_FICEXT_ANGLES
	     WRITE (88,*,ERR=921) 
     &	     "  ANGLE     F11         -F12/F11      F22/F11      F33/F11"
	     DO J=1,NB_FICEXT_ANGLES
               WRITE (88,402,ERR=921) ANG_EXT(J), F11_EXT(J), 
     &	       			   F12sF11_EXT(J), F22sF11_EXT(J), 
     &	                               F33sF11_EXT(J)
            ENDDO 
	     WRITE (88,*,ERR=921) " "
         ENDIF !Fichier Trace
	 

C* Passage aux fonctions F12, F22 et F33 pour les mu des angles 
C* decrivant les fonctions externes
C  -----------------------------------------------------------------	 
	  DO J=1,NB_FICEXT_ANGLES
	     MU_EXT(J)  = DCOS(ANG_EXT(J)*INCTE_PI/180.D+00)
            F12_EXT(J) = -F12sF11_EXT(J)*F11_EXT(J)
	     F22_EXT(J) = F22sF11_EXT(J)*F11_EXT(J)
	     F33_EXT(J) = F33sF11_EXT(J)*F11_EXT(J)
         ENDDO 
	 
C* Passage aux fonctions F11, F12, F22 et F33 pour les mu de Gauss 
C  -----------------------------------------------------------------	

C*      P11 : L'appel de la routine se fait en passant XMU et P11
C*            qui poitent sur le premier element utile
C*            (car le domaine couvert est plus grand que celui utile
C*             et qu'il y a un debut de table a indice negatif)
C*      ---	
	  CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F11_EXT,
     &	                          2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                           P11(-MIE_NBMU), IER)	
	  IF (IER.NE.0) GOTO 945   ! cas d'erreur
	
C*      P12 :
C*      ---	
	  CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F12_EXT,
     &	                          2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                           P12(-MIE_NBMU), IER)	
	  IF (IER.NE.0) GOTO 945   ! cas d'erreur
	
C*      P22 :
C*      ---	
	  CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F22_EXT,
     &	                          2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                           P22(-MIE_NBMU), IER)	
	  IF (IER.NE.0) GOTO 945   ! cas d'erreur
		
C*      P33 :
C*      ---	
	  CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F33_EXT,
     &	                          2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                           P33(-MIE_NBMU), IER)		
	  IF (IER.NE.0) GOTO 945   ! cas d'erreur 
	
	!C* Ajouter module de recomposition pour fichier TRACE
	!C* Entetes de fichiers
	!AJUSTER DEFINITION DE LA LIMITE DE TABLE F11, etc. en dur ==> dans SOS.h
	!Attention signe : -F12sF11
	 
	  IF (TRACE) THEN
	    WRITE (88,*,ERR=921) " Functions F11, F12, F22 and F33 for"
	    WRITE (88,*,ERR=921) " the original angles of the datafile"
	    WRITE (88,*,ERR=921) " ==================================="
	    WRITE (88,*,ERR=921) 
     &	    "  ANGLE     F11          F12          F22          F33"
	    DO J=1,NB_FICEXT_ANGLES
              WRITE (88,402,ERR=921) ANG_EXT(J), F11_EXT(J), 
     &	       			  F12_EXT(J), F22_EXT(J), F33_EXT(J)
           ENDDO 
	    WRITE (88,*,ERR=921) " " 
	    WRITE (88,*,ERR=921) " Functions  F11, F12, F22 and F33"
	    WRITE (88,*,ERR=921) " for the simulation angles "
	    WRITE (88,*,ERR=921) " ================================" 
	    DO J=-MIE_NBMU, MIE_NBMU
	       TETA=DACOS(XMU(J))*180./INCTE_PI
              WRITE (88,402,ERR=921) TETA,P11(J),P12(J),P22(J),P33(J)
           ENDDO 	  		
         ENDIF !Fichier Trace	 
	 
	 	
	
	
	 
C* Calcul des paramètres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------
	  CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				       P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				       ALP,BETA11,BETA22,
     &                             GAMMA12,DELTA33,ZETA,IER)  
         IF (IER.NE.0) GOTO 996   ! cas d'erreur

      ENDIF	!Fin test sur Modalite IMOD=4 (données externes)
      
      
      


C ==============================================================================
C* Définition d'un mélange de modèles par l'utilisateur 
C ==============================================================================

      IF(IMOD.EQ.5) THEN

         SOM_TAUX_AOT = 0.	 

    	  LFIC=INDEX(FICMIXTURE_AER,' ')-1

C* Lecture du fichier de définition du mélange de modes aérosols
C  -----------------------------------------------------------------
         OPEN(UNIT=20,FILE=FICMIXTURE_AER,STATUS='OLD',ERR=960)
	 
	  READ(20,'(a)',ERR=961) LIGNE_TEXTE
         IFIN=INDEX(LIGNE_TEXTE,':')
         READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) NB_MODES_MELANGE
	
         DO 122 NMOD=1,NB_MODES_MELANGE 

            !Lecture du type de modèle et controle par rapport aux modèles disponibles
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
            READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) CMODEL
            IFIN=INDEX(CMODEL,' ')-1
            CMODEL=CMODEL(1:IFIN)
 
            IF (CMODEL.EQ."LND") THEN
               MEL_MODEL_TYPE(NMOD)=1
            ELSE IF (CMODEL.EQ."JUNGE") THEN 
               MEL_MODEL_TYPE(NMOD)=2
            ELSE  !Ni LND, ni Junge
               GOTO 962
            ENDIF

            !Lecture des paramètre du modèle
            IF (CMODEL.EQ."LND") THEN

               READ(20,'(a)',ERR=961) LIGNE_TEXTE
               IFIN=INDEX(LIGNE_TEXTE,':')
               READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) 
     &         MEL_MODEL_PARAM(NMOD,1)			!LND: MODAL RADIUS

               READ(20,'(a)',ERR=961) LIGNE_TEXTE
               IFIN=INDEX(LIGNE_TEXTE,':')
               READ(LIGNE_TEXTE(IFIN+1:LENLINE),*)
     &         MEL_MODEL_PARAM(NMOD,2)			!LND: STANDARD DEVIATION

            ELSE 
               IF (CMODEL.EQ."JUNGE") THEN

                  READ(20,'(a)',ERR=961) LIGNE_TEXTE
                  IFIN=INDEX(LIGNE_TEXTE,':')
                  READ(LIGNE_TEXTE(IFIN+1:LENLINE),*)
     &            MEL_MODEL_PARAM(NMOD,2)			!JUNGE: SLOPE

                  READ(20,'(a)',ERR=961) LIGNE_TEXTE
                  IFIN=INDEX(LIGNE_TEXTE,':')
                  READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) 
     &            MEL_MODEL_PARAM(NMOD,1)		       !JUNGE: MIN RADIUS

                  READ(20,'(a)',ERR=961) LIGNE_TEXTE
                  IFIN=INDEX(LIGNE_TEXTE,':')
                  READ(LIGNE_TEXTE(IFIN+1:LENLINE),*)
     &            MEL_MODEL_PARAM(NMOD,3)			!JUNGE: MAX RADIUS

               ENDIF
            ENDIF

            !Lecture des indices de réfraction 
            !à la longueur d'onde des simulations de luminance
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) RN_WASIMU_LU
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) IN_WASIMU_LU

            !Lecture des indices de réfraction 
            !à la longueur d'onde de référence des AOT
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) RN_WAREF_LU
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) IN_WAREF_LU

            MEL_MODEL_INDEX_WAREF(NMOD,1) = RN_WAREF_LU
            MEL_MODEL_INDEX_WAREF(NMOD,2) = IN_WAREF_LU
            MEL_MODEL_INDEX_WASIMU(NMOD,1) = RN_WASIMU_LU
            MEL_MODEL_INDEX_WASIMU(NMOD,2) = IN_WASIMU_LU
 
            !Lecture de l'AOT à WAREF
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) AOT_TAUX

            MEL_MODEL_AOT(NMOD) = TA_WAREF*AOT_TAUX

            SOM_TAUX_AOT = SOM_TAUX_AOT + AOT_TAUX
 122     CONTINUE !Fin de boucle sur les modèles du fichier utilisateur
	
	  CLOSE(20)

	 
	  IF (TRACE) THEN
            WRITE (88,*,ERR=921) "AEROSOLS MIXTURE"
	     WRITE (88,*,ERR=921) "****************"
	     WRITE (88,*,ERR=921) " "
	     WRITE (88,*,ERR=921)  "Filename : ",FICMIXTURE_AER(1:LFIC)
	     WRITE (88,*,ERR=921) " "  		 
            WRITE (88,*,ERR=921) "Number of modes :", NB_MODES_MELANGE
            DO NMOD=1,NB_MODES_MELANGE
	        WRITE (88,*,ERR=921) ""
	        WRITE (88,*,ERR=921) " - Mode: ",NMOD
               IF (MEL_MODEL_TYPE(NMOD).EQ.1) THEN
	           WRITE (88,*,ERR=921) 
     &            "   - Type of size distribution : LND"
	           WRITE(88,*,ERR=921) 
     &            "        Modal radius (microns) : " , 
     &	           MEL_MODEL_PARAM(NMOD,1)
	           WRITE(88,*,ERR=921)  
     &            "        Variance (sig) : ", 
     &	           MEL_MODEL_PARAM(NMOD,2)
               ELSE
	           WRITE (88,*,ERR=921) 
     &            "   - Type of size distribution : JUNGE"
	           WRITE(88,*,ERR=921)  
     &            "        Slope : " , MEL_MODEL_PARAM(NMOD,2)
	           WRITE(88,*,ERR=921)  
     &            "        Minimal radius (microns) : ",
     &	           MEL_MODEL_PARAM(NMOD,1)
	           WRITE(88,*,ERR=921)  
     &            "        Maximal radius (microns) : ", 
     &	           MEL_MODEL_PARAM(NMOD,3)
               ENDIF
	        WRITE (88,*,ERR=921) 
     &         "   - Wavelength WA for aerosols simu (microns) : ", WA
	        WRITE (88,*,ERR=921) 
     &         "   - Refractive index (mr+j*mi) at WA: ",
     &         MEL_MODEL_INDEX_WASIMU(NMOD,1), 
     &         MEL_MODEL_INDEX_WASIMU(NMOD,2)
	        WRITE (88,*,ERR=921) 
     &         "   - Reference wavelength for AOT (mic.) ==> WAREF: ",
     &         WAREF
	        WRITE (88,*,ERR=921) 
     &         "   - Refractive index (mr+j*mi) at WAREF: ",
     &         MEL_MODEL_INDEX_WAREF(NMOD,1), 
     &         MEL_MODEL_INDEX_WAREF(NMOD,2)

	        WRITE (88,*,ERR=921) 
     &         "   - AOT ratio at WAREF: ",MEL_MODEL_AOT(NMOD)/TA_WAREF
	        WRITE (88,*,ERR=921) 
     &         "   - AOT at WAREF: ",MEL_MODEL_AOT(NMOD)
            ENDDO
	     WRITE (88,*,ERR=921) " "
            WRITE (88,*,ERR=921)  
     &      "Sommation of AOT ratios:",SOM_TAUX_AOT
	     WRITE (88,*,ERR=921) " "
         ENDIF !Fichier Trace


         !Test que la somme des taux de proportion d'AOT sur les différents modes vaut 1	 
         IF (DABS(SOM_TAUX_AOT-1.D+00).GT.CTE_GAP_TOLER_SUM_RATES) THEN
            IF (TRACE) THEN
               WRITE (88,*,ERR=921)   
     &         "The sum of AOT rates is not equal to 1"
               WRITE (88,*,ERR=921)  
     &         "Difference to 1 exceeding the threshold ", 
     &         "CTE_GAP_TOLER_SUM_RATES:", CTE_GAP_TOLER_SUM_RATES
	        WRITE (88,*,ERR=921) "Processing stopped"
            ENDIF 
            GOTO 963	!Ecart à 1 trop important ==> sortie en erreur
         ELSE
            IF (SOM_TAUX_AOT.NE.1.D+00) THEN
                 DO NMOD=1,NB_MODES_MELANGE 
                    MEL_MODEL_AOT(NMOD) =  
     &              MEL_MODEL_AOT(NMOD) / SOM_TAUX_AOT
                 ENDDO
                 IF (TRACE) THEN
                    WRITE (88,*,ERR=921)  
     &              "==> Normalisation of AOT ratios ",
     &              "so that their sum is strictly equal to 1"
                    DO NMOD=1,NB_MODES_MELANGE	      
	                WRITE (88,*,ERR=921) "    - Mode: ",NMOD
	                WRITE (88,*,ERR=921)                
     &                 "      Normalized AOT ratio at WAREF: ",                
     &                 MEL_MODEL_AOT(NMOD)/TA_WAREF
	                WRITE (88,*,ERR=921) 
     &                 "      ==> AOT at WAREF: ",MEL_MODEL_AOT(NMOD)
                    ENDDO
	             WRITE (88,*,ERR=921) " "
                ENDIF
            ENDIF  
         ENDIF

         


C* Boucles de calcul sur les modes:
C  --------------------------------  
           !----------------
           !CALCULS A WAREF: Estimation des coefficients de pondération Alpha
           !----------------
           SOM_COEF_ALPHA = 0.
		
           DO 124 NMOD=1,NB_MODES_MELANGE
		
              IF (MEL_MODEL_TYPE(NMOD).EQ.1) THEN   !Cas granu LND           
	          X = MEL_MODEL_PARAM(NMOD,2)
	          RMAX = MEL_MODEL_PARAM(NMOD,1) * DEXP(X*X) *
     &	                 DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))

	       ELSE					  !Cas granu JUNGE
	          RMAX = MEL_MODEL_PARAM(NMOD,3)
	       ENDIF

              ALPHAF = REAL( 100 +
     &	       	100*DINT(2.*INCTE_PI*RMAX / (100.*CTE_WAMIN)) )

	       IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009   
 
              IF (TRACE) THEN
	          WRITE(88,*,err=921) ' '
	          WRITE(88,*,err=921) '      RMAX (microns) :',RMAX
	          WRITE(88,*,err=921) 
     &	          '      Maximal size parameter : ',ALPHAF	  
	          WRITE(88,*,err=921) ' '     
              ENDIF

              !Constitution du nom du fichier MIE (pour indices de réfraction du Mode
              !                                    à la longueur d'onde de référence
              !                                    et limite ALPHAF estimée à WAMIN)    
              CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                           MEL_MODEL_INDEX_WAREF(NMOD,1),
     &	                           MEL_MODEL_INDEX_WAREF(NMOD,2),
     &	                           ALPHAO,ALPHAF,FICS)       
	       FICMIE = DIRS(1:LDIRS)//'/'//FICS
              LFICS=INDEX(FICMIE,' ')-1
	       
	       IF (TRACE) THEN 
	           WRITE (88,*,ERR=921) "Calculations for WA = ", WAREF
	           WRITE (88,*,ERR=921) "MIE file: ", FICMIE(1:LFICS)
	       ENDIF	
	 
	       !Affichage écran
	       WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
              !Si le fichier MIE existe déjà, le programme n'est pas lancé     
              !-----------------------------------------------------------     
              INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	       IF (EX) THEN  
                  WRITE(6,*)
     &	           '  ==> The file has already been calculated'
	  
	       ELSE   
      
                 !Calcul du fichier de MIE (pour limite ALPHAF)
                 !--------------------------------------------- 
                 CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &	                       MEL_MODEL_INDEX_WAREF(NMOD,1),
     &	                       MEL_MODEL_INDEX_WAREF(NMOD,2),
     &                        ALPHAO,ALPHAF,
     &                        FICMIE_TMP,FICLOGMIE,IER)
     	          IF (IER.NE.0) GOTO 993   ! cas d'erreur

                 COMMAND=
     &           'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
                 LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   
                 LCOM = LCOM + (4+LFICTMP)
                 IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

                 INQUIRE(FILE=FICMIE,EXIST=EX)
                 IF (EX) THEN !Si le fichier de MIE existe on détruit le fichier temporaire (généré par autre traitement)
                    OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                    CLOSE(3,STATUS='DELETE') 
                 ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de MIE
                    CALL SYSTEM(COMMAND)  
                 ENDIF
         
              ENDIF
	 
              ! Calcul des sections efficaces d'extinction et de diffusion 
              ! et des fonctions de phase I, Q, U du constituant I pour WAREF
              ! ------------------------------------------------------------- 
 	       
              !Appel de SOS_GRANU pour le modèle IGRANU=MEL_MODEL_TYPE(NMOD)
              !et les paramètres VARGRANU1 = MEL_MODEL_PARAM(NMOD,1)
              !et VARGRANU2 = MEL_MODEL_PARAM(NMOD,2)
              CALL SOS_GRANU(FICMIE,MEL_MODEL_TYPE(NMOD),
     &                       MEL_MODEL_PARAM(NMOD,1),
     &                       MEL_MODEL_PARAM(NMOD,2),
     &                       MEL_MODEL_PARAM(NMOD,3),
     &                       WAREF,MIE_NBMU,XMU,TRACE,
     &			        KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)

              IF (IER.NE.0) GOTO 995   ! cas d'erreur  

    	     
              ! Calcul des coefficients Alpha avant normalisation
              ! Le calcul doit être fait à WAREF car l'AOT de référence est à cette longueur d'onde
              ! ---------------------------------------------------------------------	
              TAUEXTI = MEL_MODEL_AOT(NMOD)
              COEF_ALPHA(NMOD) = TAUEXTI/KMAT1I
              SOM_COEF_ALPHA = SOM_COEF_ALPHA + COEF_ALPHA(NMOD)
   
  124	    CONTINUE	 !Fin boucle sur les modes

           ! Finalisation du calcul des coefficients de pondération des modes
           ! ----------------------------------------------------------------
           DO NMOD=1,NB_MODES_MELANGE
              COEF_ALPHA(NMOD) = COEF_ALPHA(NMOD) / SOM_COEF_ALPHA
           ENDDO

            
           IF (TRACE) THEN
	       WRITE(88,*,err=921) ' '
	       WRITE(88,*,err=921) '      Coef alpha :'
              DO NMOD=1,NB_MODES_MELANGE
	          WRITE(88,*,err=921) 
     &	          '      Mode, coef alpha : ',NMOD, COEF_ALPHA(NMOD)
	       ENDDO  
	       WRITE(88,*,err=921) ' '     
           ENDIF



           !----------------
           !CALCULS A WASIMU: Application du mélange des modes
           !----------------
           KMAT1 = 0.
           KMAT2 = 0.

           DO 125 NMOD=1,NB_MODES_MELANGE
		
              IF (MEL_MODEL_TYPE(NMOD).EQ.1) THEN   !Cas granu LND           
	          X = MEL_MODEL_PARAM(NMOD,2)
	          RMAX = MEL_MODEL_PARAM(NMOD,1) * DEXP(X*X) *
     &	                 DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))

	       ELSE					  !Cas granu JUNGE
	          RMAX = MEL_MODEL_PARAM(NMOD,3)
	       ENDIF

              ALPHAF = REAL( 100 +
     &	       	100*DINT(2.*INCTE_PI*RMAX / (100.*CTE_WAMIN)) )

	       IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009   
 
              IF (TRACE) THEN
	          WRITE(88,*,err=921) ' '
	          WRITE(88,*,err=921) '      RMAX (microns) :',RMAX
	          WRITE(88,*,err=921) 
     &	          '      Maximal size parameter : ',ALPHAF	  
	          WRITE(88,*,err=921) ' '     
              ENDIF

              IF (WA.EQ.WAREF) THEN	!Si WA == WAREF cela signigie que la routine est appelée pour WAREF
					!==> On doit faire les calculs de Mie pour les indices à WAREF
                 RN = MEL_MODEL_INDEX_WAREF(NMOD,1) 
                 IN = MEL_MODEL_INDEX_WAREF(NMOD,2) 
              ELSE !WASIMU		!Si WA != WAREF cela signigie que la routine est appelée pour WA
					!==> On doit faire les calculs de Mie pour les indices à WA
                 RN = MEL_MODEL_INDEX_WASIMU(NMOD,1) 
                 IN = MEL_MODEL_INDEX_WASIMU(NMOD,2) 
              ENDIF

              !Constitution du nom du fichier MIE (pour indices de réfraction du Mode
              !                                    à la longueur d'onde de simulation
              !                                    et limite ALPHAF estimée à WAMIN)    
              CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                           RN,IN,ALPHAO,ALPHAF,FICS)       
	       FICMIE = DIRS(1:LDIRS)//'/'//FICS
              LFICS=INDEX(FICMIE,' ')-1
	       
	       IF (TRACE) THEN 
	           WRITE (88,*,ERR=921) "Calculations for WA = ", WA
	           WRITE (88,*,ERR=921) "MIE file: ", FICMIE(1:LFICS)
	       ENDIF	
	 
	       !Affichage écran
	       WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
              !Si le fichier MIE existe déjà, le programme n'est pas lancé     
              !-----------------------------------------------------------     
              INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	       IF (EX) THEN  
                  WRITE(6,*)
     &	           '  ==> The file has already been calculated'
	  
	       ELSE   
      
                 !Calcul du fichier de MIE (pour limite ALPHAF)
                 !--------------------------------------------- 
                 CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &	                       RN,IN,ALPHAO,ALPHAF,
     &                        FICMIE_TMP,FICLOGMIE,IER)
     	          IF (IER.NE.0) GOTO 993   ! cas d'erreur

                 COMMAND=
     &           'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
                 LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   
                 LCOM = LCOM + (4+LFICTMP)
                 IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

                 INQUIRE(FILE=FICMIE,EXIST=EX)
                 IF (EX) THEN !Si le fichier de MIE existe on détruit le fichier temporaire (généré par autre traitement)
                    OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                    CLOSE(3,STATUS='DELETE') 
                 ELSE !Sinon on le positionne dans l'espace de stockage des fichiers de MIE
                    CALL SYSTEM(COMMAND)  
                 ENDIF
         
              ENDIF
	 
              ! Calcul des sections efficaces d'extinction et de diffusion 
              ! et des fonctions de phase I, Q, U du constituant I pour WAREF
              ! ------------------------------------------------------------- 
            	       
              !Appel de SOS_GRANU pour le modèle IGRANU=MEL_MODEL_TYPE(NMOD)
              !et les paramètres VARGRANU1 = MEL_MODEL_PARAM(NMOD,1)
              !et VARGRANU2 = MEL_MODEL_PARAM(NMOD,2)
              CALL SOS_GRANU(FICMIE,MEL_MODEL_TYPE(NMOD),
     &                       MEL_MODEL_PARAM(NMOD,1),
     &                       MEL_MODEL_PARAM(NMOD,2),
     &                       MEL_MODEL_PARAM(NMOD,3),
     &                       WA,MIE_NBMU,XMU,TRACE,
     &			        KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)

              IF (IER.NE.0) GOTO 995   ! cas d'erreur  

    	     
              ! Calcul des termes pour les sections efficaces d'extinction et de diffusion du mélange
              ! Calcul des fonctions de la matrice de phase du mélange
              ! ---------------------------------------------------------------------	
              KMAT1 = KMAT1 + COEF_ALPHA(NMOD)*KMAT1I
              KMAT2 = KMAT2 + COEF_ALPHA(NMOD)*KMAT2I
	    
	       DO 134 J=-MIE_NBMU,MIE_NBMU	     
	          P11(J) = P11(J) + P11I(J)*COEF_ALPHA(NMOD)*KMAT2I
	          P12(J) = P12(J) + P12I(J)*COEF_ALPHA(NMOD)*KMAT2I
	          P33(J) = P33(J) + P33I(J)*COEF_ALPHA(NMOD)*KMAT2I
  134	       CONTINUE  !Fin boucle sur les angles de gauss
             
  125	    CONTINUE	 !Fin boucle sur les modes

           ! Finalisation du calcul des fonctions de la matrice de phase du mélange
	    ! ----------------------------------------------------------------------	     
	    DO 135 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) / KMAT2
	       P12(J) = P12(J) / KMAT2
	       P33(J) = P33(J) / KMAT2
  135	    CONTINUE  !Fin boucle sur les angles de gauss

C* Cas de particules sphériques : P22(J) = P11(J)
C  -----------------------------------------------------
       DO 1034 J=-MIE_NBMU,MIE_NBMU
          P22(J)=P11(J)
 1034  CONTINUE
	 
C* Calcul des paramètres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------
	  CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				       P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				       ALP,BETA11,BETA22,
     &                             GAMMA12,DELTA33,ZETA,IER)  
         IF (IER.NE.0) GOTO 996   ! cas d'erreur

      ENDIF	!Fin test sur Modalite IMOD=5 (mélange utilisateur)
      
      
            
      
      
C ======================================
C* Calculs communs pour tous les modèles
C ======================================
	 
C* Calcul de l'albédo de simple diffusion
C-------------------------------------------
      PIZ=KMAT2/KMAT1

C* Troncature de l'albédo de simple diffusion des aérosols
C---------------------------------------------------------
      PIZTR = PIZ*(1.-COEF_TRONCA/2.)/(1.-PIZ*COEF_TRONCA/2.)      
    
C* Ecriture du fichier trace
C-------------------------------------------
      IF (TRACE) THEN
        WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) '=============================='
        WRITE(88,*,err=921) '* ROUTINE AEROSOLS : RESULTS *'
	 WRITE(88,*,err=921) '=============================='
        WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) ' '
	
	 IF(ITRONC.EQ.1) THEN
	      WRITE(88,*,err=921) '  Calculation with truncation'
	 ELSE
	      WRITE(88,*,err=921) '  Calculation without truncation'	   
	 ENDIF
	 
        WRITE(88,*,err=921) ' '
            	 
	 IF((COEF_TRONCA.NE.0.D+00).OR.(IMOD.NE.0)) THEN
	    WRITE(88,*,err=921) '  Truncation coefficient :',COEF_TRONCA
	    WRITE(88,*,err=921) 
     &	    '  Extinction cross-section (mic^2) - no trunca. applied:'
     &          ,KMAT1
	    WRITE(88,*,err=921) 
     &	    '  Scattering cross-section (mic^2) - no trunca. applied: '
     &          ,KMAT2
	    IF (COEF_TRONCA.NE.0.D+00) THEN
	       WRITE(88,*,err=921) 
     &		 '  Single scattering albedo '//
     &         ' - adjusted to the truncation: ',PIZTR 
	       WRITE(88,*,err=921) 
     &		 '  Single scattering albedo - no trunca. applied:',
     &             PIZTR / (1.-COEF_TRONCA/2.+PIZTR*COEF_TRONCA/2.)      
	       WRITE(88,*,err=921)
     &		 '  Beta11(0) non-normalized '//
     &         ' - adjusted to the truncation: ',Z1  
              WRITE(88,*,err=921) 
     &	        '  Asymmetry factor '//
     &         ' - adjusted to the truncation: ',BETA11(1)/3.
              WRITE(88,*,err=921) 
     &	        '  Asymmetry factor - no trunca. applied: ',
     &          COEF_TRONCA/2.+(1.-COEF_TRONCA/2.)*BETA11(1)/3.     
	    ELSE
	       WRITE(88,*,err=921) 
     & 	 ' Single scattering albedo - no trunca. applied: ',PIZTR
	       WRITE(88,*,err=921) 
     & 	'  Beta11(0) non-normalized - no trunca. applied: ',Z1
              WRITE(88,*,err=921) 
     &	       '  Asymmetry factor - no trunca. applied: ',BETA11(1)/3.
     	    ENDIF !Fin test si troncature   
           WRITE(88,*,err=921) ' '
	    WRITE(88,*,err=921) '  1st column : Scattering angle '
	    WRITE(88,*,err=921) 
     &	       '  2nd column : Truncated phase function P11' 	 
	    WRITE(88,*,err=921)
     &        '  3rd column : Original phase function P11 '
	    WRITE(88,*,err=921) 
     &        '  4th column : Truncated polarized phase function '
           WRITE(88,*,err=921) '		     P12(mu) truncated'
	    WRITE(88,*,err=921)
     &        '  5th column : Original polarized phase function '
           WRITE(88,*,err=921) '		     P12(mu) '
	    WRITE(88,*,err=921)
     &        '  6th column : Truncated polarized phase function '
           WRITE(88,*,err=921) '                 P33(mu) truncated'
	    WRITE(88,*,err=921) 
     &        '  7th column : Original polarized phase function '
           WRITE(88,*,err=921) '                 P33(mu) '
	    WRITE(88,*,err=921) ' '	    
           DO 44 J=-MIE_NBMU,MIE_NBMU
	       TETA=DACOS(XMU(J))*180./INCTE_PI
	       WRITE(88,10,err=921) TETA,P11(J),TTT(J),
     &                              P12(J)*P11(J)/TTT(J),P12(J),
     &                              P33(J)*P11(J)/TTT(J),P33(J)
   44      CONTINUE
        END IF
    
        WRITE(88,*,err=921) ' '  
	 WRITE(88,35,err=921)'I','ALP(I)','BETA11(I)','GAMMA12(I)',
     &                       'DELTA33(I)','ZETA(I)'  
        DO 45 K=0,OS_NB
	    WRITE(88,25,err=921) K,ALP(K),BETA11(K),GAMMA12(K),
     &	                         DELTA33(K),ZETA(K)
   45   CONTINUE
   	 WRITE(88,*,err=921) ' '
	 	   	   	   	   	   
      END IF

C* Enregistrement des paramètres de granulométrie dans le fichier Resultat
C  -----------------------------------------------------------------------
 1111 OPEN(1,FILE=FICGRANU,ERR=913)

      ! Section efficace d'extinction
      WRITE(1,40,ERR=923) KMAT1

      ! Section efficace de diffusion
      WRITE(1,41,ERR=923) KMAT2
     
      ! Facteur d'asymétrie (non tronqué)
      WRITE(1,42,ERR=923) 
     &	 ( COEF_TRONCA/2.+(1.-COEF_TRONCA/2.)*BETA11(1)/3. )

      ! Coefficient de troncature      
      WRITE(1,46,ERR=923) COEF_TRONCA
      
      ! Albédo de simple diffusion (éventuellement tronqué)      
      WRITE(1,47,ERR=923) PIZTR

      ! Info sur le contenu des colonnes suivantes 
      WRITE(1,39,ERR=923) 
      WRITE(1,48,ERR=923) OS_NB
      WRITE(1,49,ERR=923)

      ! Paramètres de la matrice de phase (éventuellement tronqués)       
      DO K=0,OS_NB
         WRITE(1,50,ERR=923) ALP(K),BETA11(K),GAMMA12(K),ZETA(K)
      ENDDO

 

C* Fermeture fichier
C-------------------
      CLOSE(99)	! Fermeture du fichier trace MIE
      CLOSE(88)      ! Fermeture du fichier trace GRANU
      CLOSE(1)	! Fermeture du fichier résultat GRANU
           
C* Fin nominale 
C-------------------
      GOTO 9999
            
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
 925  WRITE(6,*) '  SOS_AEROSOLS : ERROR_925: '
      WRITE(6,*) '  => Error while getting SOS_ABS_ROOT variable'
      GOTO 9998

 926  WRITE(6,*) '  SOS_AEROSOLS : ERROR_926: '
      WRITE(6,*) '  => Length of $SOS_ABS_ROOT/fic is too long' 
      WRITE(6,*) '  => CTE_LENDIR has to be adjusted.' 
      GOTO 9998

  911 WRITE(6,*) '  ERROR on logfile opening for GRANU'
      GOTO 9998

  930 WRITE(6,*) '  ERROR on angles file opening for SOS_AEROSOLS'
      GOTO 9998

  931 WRITE(6,*) '  ERROR on angles file reading for SOS_AEROSOLS'
      GOTO 9998

  940 WRITE(6,*) '  ERROR on external aerosols IOP file opening'
      WRITE(6,*) '  for SOS_AEROSOLS'
      GOTO 9998

  950 WRITE(6,*) '  ERROR : The number of angles for the external phase'
      WRITE(6,*) '  functions description is higher than the defined '
      WRITE(6,*) '  size of internal tables'
      WRITE(6,*) '  --> Check the value of CTE_MAXNB_ANG_EXT in SOS.h'
      WRITE(6,*) '      compared to the value associated to NbAng '
      WRITE(6,*) '      in the external phase functions file'
      GOTO 9998

  941 WRITE(6,*) '  ERROR on external phaze functions file reading'
      WRITE(6,*) '  for SOS_AEROSOLS'
      GOTO 9998

  942 WRITE(6,*) '  ERROR on external phaze functions file reading'
      WRITE(6,*) '  for SOS_AEROSOLS : End of file encountered'
      WRITE(6,*) '  --> Check if the value associated to NbAng is '
      WRITE(6,*) '      correct, in the external phase functions file'
      GOTO 9998

  945 WRITE(6,*) '  ERROR while interpolating external phaze functions'
      WRITE(6,*) '  for SOS_AEROSOLS'
      WRITE(6,*) '  Check that the abscissas of the phaze functions'
      WRITE(6,*) '  are correct : do not use twice same values'
      GOTO 9998

  960 WRITE(6,*) '  ERROR while opening the user file for the aerosols'
      WRITE(6,*) '  mixture definition (SOS_AEROSOLS) :'
      WRITE(6,*) FICMIXTURE_AER(1:LFIC)
      GOTO 9998

  961 WRITE(6,*) '  ERROR while reading the user file for the aerosols'
      WRITE(6,*) '  mixture definition (SOS_AEROSOLS) :'
      WRITE(6,*) FICMIXTURE_AER(1:LFIC)
      GOTO 9998

  962 WRITE(6,*) '  Non expected model name in the file for the'
      WRITE(6,*) '  mixture definition (SOS_AEROSOLS) :'
      WRITE(6,*) '  --> Wrong name : ', CMODEL
      WRITE(6,*) '  in file : ', FICMIXTURE_AER(1:LFIC)
      WRITE(6,*) '  Available models are : LND, JUNGE'
      GOTO 9998

  963 WRITE(6,*) "  The summation of AOT rates is not equal to 1"
      WRITE(6,*) "  in the user's file:", FICMIXTURE_AER(1:LFIC)
      WRITE(6,*) "  The summation is:", SOM_TAUX_AOT
      WRITE(6,*) "  --> Check the AOT rate by mode"
      WRITE(6,*) "  --> Or adjust the value of the threshold"
      WRITE(6,*) "      CTE_GAP_TOLER_SUM_RATES in SOS.h"
      GOTO 9998

  913 WRITE(6,*) '  ERROR on GRANU result file opening'
      GOTO 9998

  921 WRITE(6,*) '  ERROR on logfile writing for GRANU'
      GOTO 9998

  922 WRITE(6,*) '  ERROR on Gauss angles file reading'
      GOTO 9998

  923 WRITE(6,*) '  ERROR on GRANU result file writing'
      GOTO 9998

  993 WRITE(6,*) '  ERROR on subroutine SOS_MIE'
      GOTO 9998

  994 WRITE(6,*) '  ERROR on subroutine SOS_INIT_PARAMWMO'
      GOTO 9998

  995 WRITE(6,*) '  ERROR on subroutine SOS_GRANU'
      GOTO 9998

  996 WRITE(6,*) '  ERROR on subroutine SOS_DECOMPO_LEGENDRE'
      GOTO 9998

  997 WRITE(6,*) '  ERROR on subroutine SOS_INIT_PARAMSF'
      GOTO 9998

 1006 WRITE(6,*) '  ERROR on IMOD value : out of available range.'
      GOTO 9998

 1007 WRITE(6,*) '  ERROR on TA value : wrong negative value.'
      GOTO 9998

 1008 WRITE(6,*) '  ERROR on WA value : '
      WRITE(6,*) '  --> wavelenth has to be defined in ]0.2, 4.0[ mic.'
      GOTO 9998

 1009 WRITE(6,*) '  ERROR : CTE_MIE_ALPHAMIN > ALPHAF'
      WRITE(6,*) '  or ALPHAF > 1D5 (max value) '
      WRITE(6,*) '  ALPHAF = ', ALPHAF
      WRITE(6,*) 
     &'  --> check the definition of the size distribution model'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
      GOTO 9998

 1010 WRITE(6,*) '  ERROR : CTE_MIE_ALPHAMIN F6.4'
      WRITE(6,*) '          ALPHAF F9.2'
      WRITE(6,*) '  --> check SOS.h and argument value & compile again'
      GOTO 9998  

 1020 WRITE(6,*) '  ERROR : mv command'
      WRITE(6,*) '  => The command to move the temporary MIE file', 
     &           '  to the final depository is maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998
  
10100 WRITE(6,*) '  ERROR : --> imaginary part of refractive'
      WRITE(6,*) '              indexes have to be negative'

 9998 IER=-1
      GOTO 9999          
      
       
C* Format
C---------   
   10 FORMAT(1X,F6.2,6(2X,E12.4))
   25 FORMAT(1X,I3,5(2X,E13.5))
   35 FORMAT(3X,A1,5X,A6,9X,A9,6X,A10,5X,A10,5X,A7)
  
   
   39 FORMAT(33h---------------------------------)
   40 FORMAT(38hEXTINCTION CROSS SECTION (mic^2)     :,E13.5)
   41 FORMAT(38hSCATTERING CROSS SECTION (mic^2)     :,E13.5) 
   42 FORMAT(38hASYMMETRY FACTOR (no truncation)     :,E13.5)   
   46 FORMAT(38hTRUNCATION COEFFICIENT               :,F9.5)  
   47 FORMAT(38hSINGLE SCATTERING ALBEDO (truncation):,F9.5) 
   48 FORMAT(36hPHASE MATRIX COEFFICIENTS FOR K=0 TO,I4)
   49 FORMAT(55hALPHA(K)        BETA11(K)       GAMMA12(K)      ZETA(K))  
   50 FORMAT(E15.8,3(1X,E15.8))  

  402 FORMAT(F7.2,4E14.5)
     
  222 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X)
  
 9999 RETURN
      END   	!FIN DU PROGRAMME SOS_AEROSOLS

 
 
 
C==============================================================================
C PROCEDURE: SOS_NOM_FICMIE
C ==========
C      Cette procédure genere le nom du fichier de Mie en fonction de :
C      	- l'indice de réfraction : partie réelle RN
C	- l'indice de réfraction : partie imaginaire IN
C	- le paramètre de taille minimal ALPHAO
C	- le paramètre de taille maximal ALPHAF
C       - le nombre d'angles de Gauss utilises NBMU_GAUSS
C       - le nom du fichier d'angles utilisateur FICANGLES_USER
C          (en supprimant son arborescence) 
C
C      On impose d'avoir :
C         RN au format F5.3 : precision a 1.E-O3 pres.
C         IN au format F8.5 : precision a 1.E-O5 pres pour valeur negative.
C         ALPHAO au format F6.4 : 1.E-04 <= ALPHA0 < 10
C         ALPHAF au format F8.2 : 1.E-02 <= ALPHA0 < 1.E+05
C         NBMU_GAUSS au format I4 : 1 < NBMU_GAUSS < 9999
C
C      Le nom du fichier de Mie sera sous la forme:
C
C      MIE-RN-IN-ALPHAO-ALPHAF-MUnb-NomFicUser avec 
C      - RN au format F5.3 (ex : 1.400)
C      - IN (oppose a la valeur vraie negative) au format F7.5 (ex : 0.00800)
C      - ALPHAO au format F6.4 (ex : 0.0001)
C      - ALPHAF au format F8.2 (ex : 00500.00)
C      - MUnb le nombre d'angles de Gauss (ex : MU40)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit RN en une chaine de caractères CRN(1:4) 
C      Pour cela on passe du Format F5.3 * 1000 => A4
C           Le premier caractère CRN(1) est le chiffre unite de RN,
C           Le second caractère CRN(2) est la premiere decimale de RN,
C           etc.
C
C
C Description des paramètres entrée / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	RN (double)	: (E) partie réelle de l'indice de réfraction
C	IN (double)	: (E) partie imaginaire de l'indice de réfraction
C       ALPHAO (double) : (E) paramètre de taille minimal
C       ALPHAF (double) : (E) paramètre de taille maximal
C	FICS (CHARACTER*CTE_LENFIC1) : (S) nom du fichier de Mie (sans arborescence)	
C
C Constante utilise:
C ------------------
C    CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Aucun
C==============================================================================
      SUBROUTINE SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &                          RN,IN,ALPHAO,ALPHAF,FICS) 

      IMPLICIT NONE
         
      DOUBLE PRECISION RN	! partie réelle de l'indice de réfraction.
      DOUBLE PRECISION IN	! partie imaginaire de l'indice de réfraction.
      DOUBLE PRECISION ALPHAO   ! valeur de alpha à partir de laquelle 
      				! on commence le calcul de Mie.				
      DOUBLE PRECISION ALPHAF   ! valeur de alpha pour laquelle on arrete le calcul.

      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! chemin et nom du fichier des angles "utilisateurs"				
      CHARACTER*CTE_LENFIC1 FICS	    ! nom du fichier MIE 	
          
      CHARACTER*4 CRN		! Chaines de caractères
      CHARACTER*6 CIN 
      CHARACTER*5 CALPHAO   
      CHARACTER*7 CALPHAF
      CHARACTER*4 CNBMU
      
      INTEGER*2 I
      INTEGER*4 IDEB, IDEB_CH, IFIN, IFIN_CH  ! position de debut et fin d'une chaine de caractères  
      INTEGER*4 NBMU_GAUSS  ! Nombre d'angles de Gauss utilises 
      			      	 
	    
C --------- Conversion en chaines de caractères des paramètres d'entrée	   
      WRITE(CRN,'(I4)') int(RN*1000)  		
      WRITE(CIN,'(I6)') int(-IN*100000)		!Valeur positive
      
      WRITE(CALPHAO,'(I5)') int(ALPHAO*10000)
      WRITE(CALPHAF,'(I7)') int(ALPHAF*100)
      
      WRITE(CNBMU,'(I4)') NBMU_GAUSS

C --------- Les blancs sont remplaces par des 0	      
      DO I=1,4
         IF( CRN(I:I).EQ.' ') CRN(I:I) = '0'
      ENDDO
      
      DO I=1,6
         IF( CIN(I:I).EQ.' ') CIN(I:I) = '0'
      ENDDO
      
      DO I=1,5
         IF( CALPHAO(I:I).EQ.' ') CALPHAO(I:I) = '0'
      ENDDO
      
      DO I=1,7
         IF( CALPHAF(I:I).EQ.' ') CALPHAF(I:I) = '0'
      ENDDO
      	 
	 
      FICS = 'MIE'
      
C -------- Construction de la chaine de caractères pour l'indice RN	    
      FICS(4:4) = CRN(1:1)     ! Chiffre des unités
      FICS(5:5) = '.'	       ! Virgule
      FICS(6:8) = CRN(2:4)     ! Trois premieres decimales   
      FICS(9:9) = '-'
      
C -------- Construction de la chaine de caractères pour l'indice IN	    
      FICS(10:10) = CIN(1:1)     ! Chiffre des unités
      FICS(11:11) = '.'	         ! Virgule
      FICS(12:16) = CIN(2:6)     ! Cinq premieres decimales   	    
      FICS(17:17) = '-'	
      
C -------- Construction de la chaine de caractères pour AlphaO	    
      FICS(18:18) = CALPHAO(1:1)     ! Chiffre des unités
      FICS(19:19) = '.'	             ! Virgule
      FICS(20:23) = CALPHAO(2:5)     ! Quatre premieres decimales   	    
      FICS(24:24) = '-'	

C -------- Construction de la chaine de caractères pour AlphaF	    
      FICS(25:29) = CALPHAF(1:5)     ! Chiffres des unités
      FICS(30:30) = '.'	             ! Virgule
      FICS(31:32) = CALPHAF(6:7)     ! Deux premieres decimales   	    
      FICS(33:33) = '-'	            

C -------- Construction de la chaine de caractères pour MUgauss
C --------Structuration pour avoir NBMU_GAUSS a gauche (sans "blancs")	            
      
      FICS(34:35) = 'MU'
      
      IF (NBMU_GAUSS.LT.10) THEN
          FICS(36:36) = CNBMU(4:4)	!MUn
      ELSE 
          IF (NBMU_GAUSS.LT.100) THEN
              FICS(36:37) = CNBMU(3:4)  !MUnn
	  ELSE 
	     IF (NBMU_GAUSS.LT.1000) THEN
                 FICS(36:38) = CNBMU(2:4) !MUnnn
	     ELSE !NBMU_GAUSS < 10000
	     	 FICS(36:39) = CNBMU(1:4) !MUnnnn
	     ENDIF !NBMU_GAUSS < 1000
	  ENDIF !NBMU_GAUSS < 100
      ENDIF !NBMU_GAUSS < 10
  
       
C -------- Construction de la chaine de caractères pour NomFicUser
C -------- Structuration pour avoir NBMU_GAUSS a gauche (sans "blancs")	            


      IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN

C       Determination de la position du nom de fichier sans son arborescence          
	  IFIN_CH=INDEX(FICANGLES_USER,' ')
         IFIN_CH=IFIN_CH-1
	  IF (IFIN_CH.LE.0) IFIN_CH=CTE_LENDIR
	  
	  IDEB_CH=0
	  	  
10	  IDEB=INDEX(FICANGLES_USER(IDEB_CH:IFIN_CH),'/')
	  IDEB_CH=IDEB_CH+IDEB
         IF (IDEB.NE.0) GOTO 10	  

C       Ajout du nom de fichier d'angles "utilisateur" au nom du fichier de MIE    
         IFIN=INDEX(FICS,' ')  
	  FICS(IFIN:IFIN) = '-'
	  
	  IFIN=IFIN+1
	  FICS(IFIN:IFIN+IFIN_CH-IDEB_CH+1) 
     &	    = FICANGLES_USER(IDEB_CH:IFIN_CH)
	  
      ENDIF 
     	
	
      END !SOS_NOM_FICMIE
 
 
 
 
C==============================================================================
C PROCEDURE: SOS_INIT_PARAMWMO
C ==========
C      Cette procédure initialise l'ensemble des paramètres de chaque
C      constituant de la WMO pour la longueur d'onde de simulation :
C      	- le rayon modal
C	- la variance
C	- concentration volumique
C	- partie réelle de l'indice de réfraction
C	- partie imaginaire
C
C      Interpolation lineaire des indices pour la longueur d'onde d'entrée.
C
C      Force les indices résultats au format F5.3 pour la partie réelle
C      et F8.5 pour la partie imaginaire (valeur negative).
C
C Description des paramètres entrée / sortie
C ------------------------------------------
C	FICDATAWMO (CHARACTER*CTE_LENFIC2)	: (E) chemin et nom du fichier Data_WMO
C	WA (double)	: (E) longueur d'onde selectionnee
C       VARGRANU1I(CTE_NB_WMO) (double) : (S) rayon modal de chaque constituant
C       VARGRANU2I(CTE_NB_WMO) (double) : (S) log decimal de la variance de 
C					  chaque constituant
C	MR(CTE_NB_WMO) (double)	: (S) partie réelle de l'indice de réfraction
C       MI(CTE_NB_WMO) (double) : (S) partie imaginaire de l'indice de réfraction
C	V(CTE_NB_WMO) (double)	: (S) concentration volumique micron3/cm3
C	IER (I4)	: (S) code d'erreur = 0 si pas d'erreur, sinon =1	
C
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier des données de la WMO
C      Fichier Ascii, format d'un enregistrement :  et 
C      1ere ligne : 4(1X,F9.5)
C                   valeur du rayon modal (mic) 
C		    pour les constituants DL,WS,OC,SO
C      2eme ligne : 4(1X,F9.5)
C                   log decimal de la variance 
C		    pour les constituants DL,WS,OC,SO  
C      3eme ligne : 2(1X,F9.5,1X,E12.5)
C                   concentration volumique (mic**3/cm**3) 
C		    pour les constituants DL,WS,OC,SO        
C      4eme ligne a la fin : 9(1X,F9.5)
C                            longueur d'onde, partie réelle et imaginaire
C			     de l'indice de réfraction de chaque constituant.   
C
C
C Constantes utilisees :
C --------------------
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------    
C     -Ouverture d'un fichier
C       Erreur si le fichier des données de la WMO n'existe pas
C     -Lecture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel a la routine: SOS_INTERPOL
C
C
C  MOD:VERSION 2.0: 02/03/2016
C      Conversion de la variance du fichier de données WMO (en log décimal)
C      en variance conforme à la définition classique des LND: VAR_WMO * LOG(10)
C
C==============================================================================
      SUBROUTINE SOS_INIT_PARAMWMO(FICDATAWMO,WA,VARGRANU1I,
     &				   VARGRANU2I,MR,MI,V,IER)

      IMPLICIT NONE
      
C Constantes AEROSOLS
C -------------------

      INTEGER CTE_NB_WMO	!Nombre de types d'aérosols de la WMO utilises: DL,WS,OC,SO
      PARAMETER (CTE_NB_WMO=4)

      INTEGER CTE_NB_AERO	!Dimension des tableaux résultats 
      PARAMETER (CTE_NB_AERO=5)	!(4 composants WMO, 5 composants Shettle & Fenn)	
      
            
C* Definition des variables       
C*-----------------------------------------

      CHARACTER*CTE_LENFIC2 FICDATAWMO	! chemin et nom du fichier Data_WMO
      
      DOUBLE PRECISION WA	! Longueur d'onde en microns selectionnee.	
      DOUBLE PRECISION WA1	! Longueur d'onde en microns lue dans le
      DOUBLE PRECISION WA2	! fichier Data_WMO.
      
      DOUBLE PRECISION VARGRANU1I(CTE_NB_AERO)   ! Paramètre 1 de la granulométrie
      						       ! pour chaque modèle de constituant.
      DOUBLE PRECISION VARGRANU2I(CTE_NB_AERO)   ! Paramètre 2 de la granulométrie
      						       ! pour chaque modèle de constituant.
      
      DOUBLE PRECISION MR(CTE_NB_AERO)   ! partie réelle de l'indice de réfraction
      					      ! pour chaque modèle d'aérosols.
      DOUBLE PRECISION MI(CTE_NB_AERO)   ! partie imaginaire de l'indice de réfraction.
      					      ! pour chaque modèle d'aérosols.
      DOUBLE PRECISION MR1(CTE_NB_WMO)  !
      DOUBLE PRECISION MI1(CTE_NB_WMO)  ! utilisé pour le stockage de 
      DOUBLE PRECISION MR2(CTE_NB_WMO)  ! valeurs intermediaires
      DOUBLE PRECISION MI2(CTE_NB_WMO)  !
      
      DOUBLE PRECISION V(CTE_NB_WMO)	!Concentration volumique des constituants
      
      INTEGER*4 I	   !Indice.       
      INTEGER*4 IER        ! code d'erreur =0 si pas d'erreur, =-1 sinon
      
C* Initialisation
C  --------------
       IER=0	 
       
       DO I=1,CTE_NB_WMO
          VARGRANU1I(I)=0.
          VARGRANU2I(I)=0.
	   MR(I)=0.
    	   MI(I)=0.
	   V(I)=0.
       ENDDO

C* Ouverture du fichier contenant les données de la WMO
C  ---------------------------------------------------	 
       OPEN(8,FILE=FICDATAWMO,STATUS='OLD',ERR=991)

C* Lecture du rayon modal
C  ---------------------- 	 
       READ(8,333,err=992) VARGRANU1I(1),VARGRANU1I(2),
     &                     VARGRANU1I(3),VARGRANU1I(4)
     	  
C* Lecture du log decimal de la variance
C  ------------------------------------- 	 
       READ(8,333,err=992) VARGRANU2I(1),VARGRANU2I(2),
     &                     VARGRANU2I(3),VARGRANU2I(4)      

C* Conversion en variance pour LND compatible avec la définition classique
C  ------------------------------------- 	 
       DO I=1,4
          VARGRANU2I(I) = VARGRANU2I(I) * DLOG(10.D+00)
       ENDDO

C* Lecture de la concentration volumique
C  ------------------------------------- 	 
       READ(8,444,err=992) V(1),V(2),V(3),V(4) 

C* Determination des indices de réfraction par interpolation lineaire
C  ------------------------------------------------------------------

C* Initialisation avec les premieres valeurs du fichier
C  ---------------------------------------------------- 
       READ(8,555,END=999,err=992)WA1,MR1(1),MI1(1),MR1(2),MI1(2),
     &                            MR1(3),MI1(3),MR1(4),MI1(4)

C* Pour chaque ligne du fichier, on recherche la longueur d'onde
C  juste superieure a WA selectionnee
C  --------------------------------------------------------------
  882  READ(8,555,END=999,err=992)WA2,MR2(1),MI2(1),MR2(2),MI2(2),
     &                            MR2(3),MI2(3),MR2(4),MI2(4)
     
       ! WA1<=WA<WA2, interpolation des indices	     
       IF((WA1.LE.WA).AND.(WA.LE.WA2)) THEN
       
	  DO I=1,CTE_NB_WMO	  
	     CALL SOS_INTERPOL(MR1(I),MR2(I),WA1,WA2,WA,MR(I))
	     CALL SOS_INTERPOL(MI1(I),MI2(I),WA1,WA2,WA,MI(I))
	     
C******** Conversion aux formats F5.3 et F8.3	   
	     MR(I) =  DNINT( MR(I)*1000.D+00) / 1000.D+00
	     MI(I) = -DNINT(-MI(I)*100000.D+00) / 100000.D+00	    
	  ENDDO
	  GOTO 881
	       
       ELSE
       	 
          WA1=WA2	! Reinitialisation: ligne1 recoit ligne2
	  DO I=1,CTE_NB_WMO
	     MR1(I)=MR2(I)
	     MI1(I)=MI2(I)
     	  ENDDO
	  GOTO 882	! passage a la ligne suivante
	       
       ENDIF

  999  CONTINUE
  	    
  881  CONTINUE
       CLOSE(8)

C* Fin nominale 
C-------------------
       GOTO 9999
              
C* Cas d'erreur et retour du status -1 au programme appelant
C----------------------------------------------------------
  991  WRITE(6,*) 'Error while opening the WMO datafile'
       WRITE(6,*) 'FICDATAWMO :',FICDATAWMO
       IER=-1
       GOTO 9999
  992  WRITE(6,*) 'Error while reading the WMO datafile'
       IER=-1
       GOTO 9999

C* Format
C---------   
  333  FORMAT(4(1X,F9.5))
  444  FORMAT(1X,F9.5,1X,E12.5,1X,F9.5,1X,E12.5)
  555  FORMAT(9(1X,F9.5))

       
 9999  CONTINUE
 
       RETURN
       END	!FIN DE LA PROCEDURE SOS_INIT_PARAMWMO


C==============================================================================
C PROCEDURE: SOS_INIT_PARAMSF
C ==========
C      Cette procédure initialise l'ensemble des paramètres de chaque
C      constituant des modèles de Shettle et Fenn pour la longueur d'onde 
C      et l'humidite relative selectionnees :
C      	- le rayon modal
C	- la variance
C	- partie réelle de l'indice de réfraction
C	- partie imaginaire
C
C Description des paramètres
C --------------------------
C	DIRFIC (CHARACTER*CTE_LENDIR)  : (E) chemin d'acces au fichier de données
C	WA (double)	      : (E) longueur d'onde selectionne
C	RH (double)	      : (E) humidite relative
C       VARGRANU1I(CTE_NB_AERO) (double) : (S) rayon modal de chaque constituant
C       VARGRANU2I(CTE_NB_AERO) (double) : (S) log decimal de la variance de 
C					  chaque constituant
C	MR(CTE_NB_AERO) (double) : (S) partie réelle de l'indice de réfraction
C       MI(CTE_NB_AERO) (double) : (S) partie imaginaire de l'indice de réfraction
C	IER (I4) : (S) code d'erreur = 0 si pas d'erreur, sinon =1	
C
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier des paramètres de granulométrie des composants
C      de Shettle et Fenn :
C      Fichier Ascii, format d'un enregistrement:
C      1ere ligne : 5(1X,F9.5)
C		    log decimal de la variance 
C		    pour les constituants SR,LR,SU,LU et OM    
C      2eme ligne a la fin : 1X,F5.2,5(1X,F9.5)
C			     % d'humidite relative, rayon modal
C			     de chaque constituant.   
C
C   -->Contenu des fichiers IRefrac_XX (avec XX = SR, LR ...)
C      Fichier Ascii, format d'un enregistrement 17(1X,F9.5)
C      pour chaque ligne: longeur d'onde, partie réelle et partie imaginaire 
C			  de l'indice de réfraction pour 8 humidites 	
C			  relatives : 0%, 50%, 70%, 80%, 90%, 95%, 98%, 99%.
C
C
C Constantes utilisees :
C --------------------
C    CTE_LENDIR : Longueur des noms de répertoires.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_AER_DATASF  : Nom du fichier contenant rayon modal et log variance S&F.
C    CTE_AER_XX_SF   : Nom des fichiers contenant les indices de réfraction S&F
C		       pour chaque constituant XX = SR, LR, SU, LU, OM.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------    
C     -Ouverture d'un fichier
C       Erreur si l'un des fichiers de données n'existe pas
C     -Lecture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel a la routine: SOS_INTERPOL
C
C
C  MOD:VERSION 2.0: 02/03/2016
C      Conversion de la variance du fichier de données S&F (en log décimal)
C      en variance conforme à la définition classique des LND: VAR_WMO * LOG(10)
C
C==============================================================================
      SUBROUTINE SOS_INIT_PARAMSF(DIRFIC,WA,RH,VARGRANU1I,
     &                         	  VARGRANU2I,MR,MI,IER)

      IMPLICIT NONE
      
C Constantes AEROSOLS
C -------------------

      INTEGER CTE_NB_SF	        !Nombre des composants des modèles Shettle & Fenn
      PARAMETER (CTE_NB_SF=5)

      INTEGER CTE_NB_RH		!Nombre de valeurs d'humidite relative
      PARAMETER (CTE_NB_RH=8)
 
      INTEGER CTE_NB_AERO	!Dimension des tableaux résultats 
      PARAMETER (CTE_NB_AERO=5)	!(4 composants WMO, 5 composants Shettle & Fenn)	
      
            
C* Definition des variables       
C*-----------------------------------------

      CHARACTER*CTE_LENDIR DIRFIC	! Chemin d'acces aux fichiers S&F
      CHARACTER*CTE_LENFIC2 FICDATASF	! Nom du fichier Data_SF
      CHARACTER*CTE_LENFIC2 FIC_SR_SF 	! Fichier IRefrac_SR
      CHARACTER*CTE_LENFIC2 FIC_LR_SF 	! Fichier IRefrac_LR
      CHARACTER*CTE_LENFIC2 FIC_SU_SF 	! Fichier IRefrac_SU
      CHARACTER*CTE_LENFIC2 FIC_LU_SF 	! Fichier IRefrac_LU
      CHARACTER*CTE_LENFIC2 FIC_OM_SF 	! Fichier IRefrac_OM
      
      DOUBLE PRECISION WA	! Longueur d'onde en microns selectionnee.	
      DOUBLE PRECISION WA1	! Longueur d'onde en microns lue dans les fichiers
      DOUBLE PRECISION WA2	
      
      DOUBLE PRECISION RH	! Humidite relative (%)
      DOUBLE PRECISION RH1	! Humidite relative lue dans les fichiers
      DOUBLE PRECISION RH2	!
      
      DOUBLE PRECISION VARGRANU1I(CTE_NB_AERO)  ! Paramètre 1 de la granulométrie
      						      ! pour chaque modèle de constituant.
      DOUBLE PRECISION VARGRANU2I(CTE_NB_AERO)  ! Paramètre 2 de la granulométrie
      						      ! pour chaque modèle de constituant.
      
      DOUBLE PRECISION MR(CTE_NB_AERO)  ! partie réelle de l'indice de réfraction
      					     ! pour chaque modèle d'aérosols.
      DOUBLE PRECISION MI(CTE_NB_AERO)  ! partie imaginaire de l'indice de réfraction.
      					     ! pour chaque modèle d'aérosols.
      
      DOUBLE PRECISION MR_CPT1		!
      DOUBLE PRECISION MR_CPT2		! Valeurs intermediaires
      DOUBLE PRECISION MI_CPT1		!
      DOUBLE PRECISION MI_CPT2		!
      
      DOUBLE PRECISION RM1(CTE_NB_SF)	!
      DOUBLE PRECISION RM2(CTE_NB_SF)	! utilise pour le stockage de
      DOUBLE PRECISION MR1(CTE_NB_RH)   ! valeurs intermediaires 
      DOUBLE PRECISION MR2(CTE_NB_RH)   ! 
      DOUBLE PRECISION MI1(CTE_NB_RH)   ! 
      DOUBLE PRECISION MI2(CTE_NB_RH)   !
	
      INTEGER*4 I,H	   ! Indice
      INTEGER*4 CPT	   ! Compteur 
      INTEGER*4 IFIN	   ! Indicateur de fin d'une chaine de caractères
      
      INTEGER*4 IER        ! code d'erreur =0 si pas d'erreur, =-1 sinon
      
C* Initialisation
C  --------------
       IER=0
       
       DO I=1,CTE_NB_SF
          VARGRANU1I(I)=0.
          VARGRANU2I(I)=0.	  
	   MR(I)=0.
    	   MI(I)=0. 
       ENDDO

       
       ! Initialisation des noms des différents fichiers contenant 
       ! les données du modèle S&F
       
       IFIN=INDEX(DIRFIC,' ')
       IFIN=IFIN-1
       IF (IFIN.LE.0) IFIN=CTE_LENDIR
      
       FICDATASF = DIRFIC(1:IFIN)//'/'//CTE_AER_DATASF
       FIC_SR_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_SR_SF
       FIC_LR_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_LR_SF
       FIC_SU_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_SU_SF
       FIC_LU_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_LU_SF
       FIC_OM_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_OM_SF

C* Ouverture du fichier contenant les paramètres de granulométrie
C* des composants de S&F
C  --------------------------------------------------------------
       OPEN(9,FILE=FICDATASF,STATUS='OLD',ERR=991)
   	  
C* Lecture du log decimal de la variance
C  ------------------------------------- 	 
       READ(9,222,err=993) VARGRANU2I(1),VARGRANU2I(2),VARGRANU2I(3),
     &                     VARGRANU2I(4),VARGRANU2I(5) 

C* Conversion en variance pour LND compatible avec la définition classique
C  ------------------------------------- 	 
       DO I=1,5
          VARGRANU2I(I) = VARGRANU2I(I) * DLOG(10.D+00)
       ENDDO
                    
C* Lecture du rayon modal en fonction de l'humidite relative
C  --------------------------------------------------------- 	 
       READ(9,333,err=993) RH1,RM1(1),RM1(2),
     &                      RM1(3),RM1(4),RM1(5)
        
	CPT=1
	
	! Cas ou RH = 0%
	IF(RH1.EQ.RH) THEN
	
	   DO I=1,CTE_NB_SF
	      VARGRANU1I(I)=RM1(I)
	   ENDDO
	   GOTO 884
	   
	ENDIF

C* Recherche des valeurs de RH encadrant la valeur selectionnee
C  ------------------------------------------------------------
  883   READ(9,333,err=993) RH2,RM2(1),RM2(2),
     &                      RM2(3),RM2(4),RM2(5)
     
        CPT=CPT+1	! Compteur indiquant la position de 
			! la valeur RH2 juste superieure a RH
	
	! Cas ou RH1<RH<=RH2, interpolation des indices de réfraction
	IF((RH1.LT.RH).AND.(RH.LE.RH2)) THEN
	
	   DO I=1,CTE_NB_SF	  
	      CALL SOS_INTERPOL(RM1(I),RM2(I),RH1,RH2,RH,VARGRANU1I(I))
     	   ENDDO
	   GOTO 884
	       
        ELSE	! Reinitialisation des valeurs
       
          RH1=RH2
	  DO I=1,CTE_NB_SF
	     RM1(I)=RM2(I)
     	  ENDDO
	  GOTO 883	! passage a la ligne suivante dans le fichier	       
	  
       ENDIF
       
  884  CONTINUE


C* Ouverture des fichiers contenant l'indice de réfraction pour 
C  chaque constituant des modèles S&F
C  -----------------------------------------------------------
       OPEN(11,FILE=FIC_SR_SF,STATUS='OLD',ERR=992)
       OPEN(12,FILE=FIC_LR_SF,STATUS='OLD',ERR=992)
       OPEN(13,FILE=FIC_SU_SF,STATUS='OLD',ERR=992)
       OPEN(14,FILE=FIC_LU_SF,STATUS='OLD',ERR=992)
       OPEN(15,FILE=FIC_OM_SF,STATUS='OLD',ERR=992)
  
       DO I=1,CTE_NB_SF

C* Determination des indices de réfraction par interpolation lineaire
C  ------------------------------------------------------------------

C* Initialisation avec les premieres valeurs du fichier
C  ----------------------------------------------------
          READ((I+10),555,END=999,err=993)WA1,MR1(1),MI1(1),MR1(2),
     &                             MI1(2),MR1(3),MI1(3),MR1(4),MI1(4),
     &				   MR1(5),MI1(5),MR1(6),MI1(6),MR1(7),
     &				   MI1(7),MR1(8),MI1(8)

C* Pour chaque ligne du fichier, on recherche la longueur d'onde
C  juste superieure a celle selectionnee
C  --------------------------------------------------------------
  885     READ((I+10),555,END=999,err=993)WA2,MR2(1),MI2(1),MR2(2),
     &                             MI2(2),MR2(3),MI2(3),MR2(4),MI2(4),
     &				   MR2(5),MI2(5),MR2(6),MI2(6),MR2(7),
     &				   MI2(7),MR2(8),MI2(8)
     
          ! Cas ou WA1<=WA<=WA2, interpolation des données
	  IF((WA1.LE.WA).AND.(WA.LE.WA2)) THEN
	  
	     ! Cas ou RH = 0%, interpolation de WA mais pas de RH
	     IF(CPT.EQ.1) THEN
	        CALL SOS_INTERPOL(MR1(CPT),MR2(CPT),WA1,WA2,WA,MR(I))
		 CALL SOS_INTERPOL(MI1(CPT),MI2(CPT),WA1,WA2,WA,MI(I))
	     
	     ELSE
	      ! Interpolations pour la partie réelle de l'indice de réfraction	     
		! Interpolation de WA pour RH1(CPT-1) et RH2(CPT) encadrant RH
		 CALL SOS_INTERPOL(MR1(CPT-1),MR2(CPT-1),
     &		                  WA1,WA2,WA,MR_CPT1)
	        CALL SOS_INTERPOL(MR1(CPT),MR2(CPT),WA1,WA2,WA,MR_CPT2)
		! Interpolation de RH
	        CALL SOS_INTERPOL(MR_CPT1,MR_CPT2,RH1,RH2,RH,MR(I))
	      
	      ! Interpolations pour la partie imaginaire
	        CALL SOS_INTERPOL(MI1(CPT-1),MI2(CPT-1),
     &				  WA1,WA2,WA,MI_CPT1)
	        CALL SOS_INTERPOL(MI1(CPT),MI2(CPT),WA1,WA2,WA,MI_CPT2)
	        CALL SOS_INTERPOL(MI_CPT1,MI_CPT2,RH1,RH2,RH,MI(I))
	     ENDIF
		
          ELSE
       
            WA1=WA2	! Reinitilisation 
	     DO H=1,CTE_NB_RH
	        MR1(H)=MR2(H)
	        MI1(H)=MI2(H)
     	     ENDDO
	     GOTO 885	! Passage a la ligne suivante du fichier
	       
          ENDIF
	  
C* Conversion aux formats F5.3 et F8.3	
	  MR(I) =  DNINT( MR(I)*1000.D+00) / 1000.D+00
	  MI(I) = -DNINT(-MI(I)*100000.D+00) / 100000.D+00
	  	    
       ENDDO
       
  999  CONTINUE
  	    
       CLOSE(9)
       CLOSE(11)
       CLOSE(12)
       CLOSE(13)
       CLOSE(14)
       CLOSE(15)

C* Fin nominale 
C-------------------
       GOTO 9999
              
C* Cas d'erreur et retour du status -1 au programme appelant
C----------------------------------------------------------
  991  WRITE(6,*) 
     &     'Error while opening the Shettle&Fenn component  datafile'
       IER=-1
       GOTO 9999
  992  WRITE(6,*) 
     &	  'Error while opening the S&F  refractive indexes datafile'
       IER=-1
       GOTO 9999     
  993  WRITE(6,*) 
     &    'Error while reading the Shettle&Fenn component  datafile'
       IER=-1
       GOTO 9999
  994  WRITE(6,*) 
     &	  'Error while reading the S&F  refractive indexes datafile'
       IER=-1
       GOTO 9999
       
C* Format
C---------   
  222  FORMAT(5(1X,F9.5))
  333  FORMAT(1X,F5.2,5(1X,F9.5))
  555  FORMAT(17(1X,F9.5))

       
 9999  CONTINUE
 
       RETURN
       END	!FIN DE LA PROCEDURE SOS_INIT_PARAMSF




C==============================================================================
C PROCEDURE: SOS_INTERPOL
C ==========
C      Cette procédure permet d'interpoler lineairement en X une fonction f
C      dont on connait les valeurs (X1,Y1) et (X2,Y2).
C
C Description des paramètres
C --------------------------
C	X1 (double) : (E)
C	X2 (double) : (E)
C	Y1 (double) : (E) Valeur f(X1)
C	Y2 (double) : (E) Valeur f(X2)
C	X (double)  : (E) Abscisse pour laquelle on interpole la valeur f(X)
C   Y (double)  : (S) Valeur f(X) interpolee
C 
C
C==============================================================================
      SUBROUTINE SOS_INTERPOL(Y1,Y2,X1,X2,X,Y)

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------

      DOUBLE PRECISION Y	
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Y2
      DOUBLE PRECISION X
      DOUBLE PRECISION X1
      DOUBLE PRECISION X2      
      
C* Calcul du paramètre interpole
C  -----------------------------
      Y = ((Y2-Y1)/(X2-X1))*(X-X2) + Y2

      END	!FIN DE LA PROCEDURE SOS_INTERPOL


            
C==============================================================================
C PROCEDURE: SOS_DECOMPO_LEGENDRE
C ==========
C       Cette procédure calcule les paramètres alpha, beta, gamma delta et zeta
C	de la decomposition de la focntion de phase en polynomes de Legendre.
C	Application eventuelle de la troncature.
C
C       Les fonctions de phase sont stockees dans des tableaux definis jusqu'à
C       CTE_MIE_NBMU_MAX mais sont utilisees jusqu'à MIE_NBMU (angles effectifs).
C
C Description des paramètres
C --------------------------
C	ITRONC (I4) : (E/S) Option de troncature (1 pour la realiser).
C                           Mise a 0 en sortie si la troncature est trop faible.
C	TRACE (logical) : (E)  =vrai si écriture dans le fichier trace
C      MIE_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C	XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double) : (E) Cosinus des angles de Gauss.
C	XHR(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double) : (E) Poids des angles de Gauss.
C	P11(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double)  : (E) Fonction de phase calculée 
C					         à partir des fonctions la théorie de MIE.
C	P12(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double)  : (E) Elément P12 de la matrice de phase.
C	P33(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	(double) : (E) Elément P33 de la matrice de phase.
C	TTT(-CTE_MIE_NBMU_MAX:CTE_MIE_NBM_MAXU) (double) : (S) Copie de la fonction de phase (P11) 
C					 	               avant troncature.
C	COEF_TRONCA (double)	     : (S) Coefficient de troncature.
C	Z1 (double)		     : (S) Valeur de BETA11(0) avant normalisation.
C	ALP(0:CTE_OS_NB_MAX) (double)     : (S) Paramètres alpha(k).
C	BETA11(0:CTE_OS_NB_MAX) (double)  : (S) Paramètres beta11(k).
C	BETA22(0:CTE_OS_NB_MAX) (double)  : (S) Paramètres beta22(k).
C	GAMMA12(0:CTE_OS_NB_MAX) (double) : (S) Paramètres gamma12(k).
C	ZETA(0:CTE_OS_NB_MAX) (double)    : (S) Paramètres zeta(k).
C	DELTA33(0:CTE_OS_NB_MAX) (double)  : (S) Paramètres delta33(k).
C	IER (I4)		     : (S) Code d'erreur =0 si pas d'erreur, sinon =1.	
C
C 
C Description des fichiers utilises
C ---------------------------------
C     Aucun
C
C Constantes utilisees :
C --------------------
C    CTE_MIE_NBMU_MAX : Nombre maximal d'angles positifs.
C    CTE_OS_NB_MAX :  Valeur maximale pour l'ordre limite du développement en 
C                     polynomes de Legendre.
C    SOS_PH_SEUIL_TRONCA : Valeur seuil pour la troncature.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     -Ecriture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C==============================================================================
      SUBROUTINE SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,
     &				          OS_NB,P11,TTT,P12,P22,P33,
     &				          COEF_TRONCA,Z1,
     &				          ALP,BETA11,BETA22,
     &                                GAMMA12,DELTA33,ZETA,IER)

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------
      DOUBLE PRECISION XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)  
                                                  !Cosinus des angles de Gauss.
      DOUBLE PRECISION XHR(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)  
                                                  !Poids des angles de Gauss.
      
      DOUBLE PRECISION P11(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Fonction de phase calculée à partir des fonctions de 
		       !la théorie de MIE.
      DOUBLE PRECISION P12(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P12 de la matrice de phase.
      DOUBLE PRECISION P22(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P22 de la matrice de phase.
      DOUBLE PRECISION P33(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P33 de la matrice de phase.

      DOUBLE PRECISION TTT(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Copie de la fonction de phase (P11) avant troncature.
		       
		
      DOUBLE PRECISION PL(-1:CTE_OS_NB_MAX+1) 
                       !Polynomes de Legendre a l'ordre L calcule pour un angle.
      DOUBLE PRECISION POL(0:CTE_OS_NB_MAX+1) 
                       !Fonction de Legendre d'ordre S=2 et d'ordre L variable
		       !(calcule pour un angle).
		       			
      DOUBLE PRECISION D	!Variable pour le calcul de POL.
      DOUBLE PRECISION E	!Variable pour le calcul de POL.
      		       
      DOUBLE PRECISION X
      DOUBLE PRECISION XRMU	!Valeur utilisee du tableau des angles de Gauss.
				
      DOUBLE PRECISION XX	
      DOUBLE PRECISION XXX

      DOUBLE PRECISION COEF_TRONCA	!Coefficient de troncature.    
      DOUBLE PRECISION AA	!Pente de la droite de troncature.
      DOUBLE PRECISION X1	
      DOUBLE PRECISION X2


      DOUBLE PRECISION CO1	!Variable pour le calcul des ALP(k) et ZETA(k).
      DOUBLE PRECISION CO2	!Variable pour le calcul des ALP(k) et ZETA(k).
      DOUBLE PRECISION CO3	!Variable pour le calcul des ALP(k) et ZETA(k).

      DOUBLE PRECISION SOM1  	!Variable pour le calcul des ALP(k) et ZETA(k).
      DOUBLE PRECISION SOM2	!Variable pour le calcul des ALP(k) et ZETA(k).
      DOUBLE PRECISION SOM3	!Variable pour le calcul des ALP(k) et ZETA(k).
      DOUBLE PRECISION SOM4	!Variable pour le calcul des ALP(k) et ZETA(k).
      
      DOUBLE PRECISION Z1	!Valeur de BETA11(0) avant normalisation.

      DOUBLE PRECISION ALP(0:CTE_OS_NB_MAX)	! Paramètres alpha(k).
      DOUBLE PRECISION BETA11(0:CTE_OS_NB_MAX)	! Paramètres beta11(k).
      DOUBLE PRECISION BETA22(0:CTE_OS_NB_MAX)	! Paramètres beta22(k).
      DOUBLE PRECISION GAMMA12(0:CTE_OS_NB_MAX)	! Paramètres gamma12(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Paramètres zeta(k).
      DOUBLE PRECISION DELTA33(0:CTE_OS_NB_MAX)	! Paramètres delta33(k).
                        
	!   Angles bornes pour la troncature (pour conversion en double precision)
	!-----------------------------------
      DOUBLE PRECISION MU1_TRONCA ! cosinus des angles Theta1 et Theta2 pour la définition
      DOUBLE PRECISION MU2_TRONCA ! du segment tronquant la fonction de phase

      INTEGER*4 MIE_NBMU   ! Nombre d'angles (positifs) effectivement utiles 
      INTEGER*4 ITRONC	   ! Option de troncature (1 pour la realiser).
      INTEGER*4 OS_NB 	   ! Ordre des développements des fonctions de phase 
			   ! en Polynomes de Legendre.      
      INTEGER*4 K          ! Indice de tableau.
      INTEGER*4 KK	   ! Indice du tableau des angles de Gauss pour le mu
      			   ! le plus proche de 0.94 par valeur superieure.			               
      INTEGER*4 NN
      INTEGER*4 MM				
      INTEGER*4 I	   ! Indice.  
      INTEGER*4 J	   ! Indice.      
      INTEGER*4 IER	   ! Code d'erreur =0 si pas d'erreur, =1 sinon.         

      LOGICAL TESS	   !Variable pour le calcul de Troncature. 
      LOGICAL TRACE	   !=vrai si écriture dans le fichier trace 
    



C* Copie de la fonction de phase en intensite avant la troncature.
C-----------------------------------------------------------------
      DO 55 J=-MIE_NBMU,MIE_NBMU
         TTT(J)=P11(J)
  55  CONTINUE 
      
C==============
C* TRONCATURE
C============== 
      IF(ITRONC.EQ.0) GO TO 59

      MU1_TRONCA = CTE_AER_MU1_TRONCA	
      MU2_TRONCA = CTE_AER_MU2_TRONCA


C* Troncature de la fonction de phase en intensite pour les angles 
C* inferieurs a 19.95 degres (mu = 0.94).
C------------------------------------------------------------------------------

C* Determination des angles de Gauss juste inferieurs a 36.87 degres (mu=0.8)
C* et a 19.95 degres (mu = 0.94)
C* Détermination des angles de Gauss juste inférieurs à 36.87 degrés (mu=0.8 : MU1_TRONCA)
C* et à 19.95 degrés (mu = 0.94 : MU2_TRONCA) par exemple typiquement pour les aérosols.
C--------------------------------------------------------------------
      TESS=.TRUE.
      DO 56 J=1,MIE_NBMU
         IF (.NOT.(TESS.AND.(XMU(J).GT.MU1_TRONCA)
     &	     .AND.(XHR(J).NE.0.))) GO TO 56
         K=J-1
         TESS=.FALSE.
  56  CONTINUE

      TESS=.TRUE.
      DO 57 J=1,MIE_NBMU
         IF (.NOT.(TESS.AND.(XMU(J).GT.MU2_TRONCA)
     &	     .AND.(XHR(J).NE.0.)))  GO TO 57
         KK=J-1
         TESS=.FALSE.
  57  CONTINUE

      IF (TRACE) THEN
         WRITE(88,*,err=921) 
     &	    'Truncation applied for a segment defined between :'
         WRITE(88,*,err=921) 'Mu1 :',MU1_TRONCA, 
     &	       ' --> Theta1 = ', DACOS(MU1_TRONCA)*180./INCTE_PI
         WRITE(88,*,err=921) 'Mu1_effective :',XMU(K), 
     &	       ' --> Theta1_eff = ', DACOS(XMU(K))*180./INCTE_PI
	  WRITE(88,*,err=921) 'Mu2 :',MU2_TRONCA, 
     &	       ' --> Theta2 = ', DACOS(MU2_TRONCA)*180./INCTE_PI
         WRITE(88,*,err=921) 'Mu2_effective :',XMU(KK), 
     &	       ' --> Theta2_eff = ', DACOS(XMU(KK))*180./INCTE_PI
	 WRITE(88,*,err=921) ""
      ENDIF

C* Determination de la pente de la droite "log10(p(cos(teta))) = AA * teta + BB"
C* passant par les valeurs de log10 (p) pour les mu de Gauss ci-dessus.
C--------------------------------------------------------------------------
      AA=
     &(DLOG10(P11(KK))-DLOG10(P11(K)))/(DACOS(XMU(KK))-DACOS(XMU(K)))

C* Application de la troncature
C-------------------------------
      X1=DLOG10(P11(KK))
      X2=DACOS(XMU(KK))
      
      DO 58 J=KK+1,MIE_NBMU
         IF(DABS(XMU(J)-1.).LE.1D-08) COEF_TRONCA=X1-AA*X2
         COEF_TRONCA=X1+AA*(DACOS(XMU(J))-X2)
         P11(J)=10**COEF_TRONCA
  58  CONTINUE

  59  CONTINUE


C* Calcul des paramètres BETA11(k) de la decomposition de la fonction
C* de phase (éventuellement tronquée) en polynomes de Legendre.
C------------------------------------------------------------------
      DO 61 J=-MIE_NBMU,MIE_NBMU
         IF(J.EQ.0) GO TO 61
         X=P11(J)*XHR(J)
         XRMU=XMU(J)

C* Calcul des polynomes de Legendre d'ordre K, pour la valeur de mu fixee.
C------------------------------------------------------------------------
         DO 60 K=0,OS_NB
            PL(-1)=0.
            PL(0)=1.
            PL(K+1)=((2*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
            BETA11(K)=BETA11(K)+X*PL(K)
  60     CONTINUE
  61  CONTINUE

      DO 62 K=0,OS_NB
         BETA11(K)=(2*K+1)*BETA11(K)*.5
  62  CONTINUE



C* Calcul du coefficient de troncature
C--------------------------------------
      IF (ITRONC.EQ.1) THEN
         COEF_TRONCA=2*(1-BETA11(0))
      ELSE
         COEF_TRONCA=0.0
      ENDIF

C* Si le coefficient de troncature est inferieur au seuil Seuil_Tronca
C* la troncature est annulee.   
C----------------------------------------------------------------------	
      IF((ITRONC.EQ.1).AND.(COEF_TRONCA.LT.CTE_PH_SEUIL_TRONCA)) THEN
      
         IF (TRACE) THEN
            WRITE(88,*,err=921) 'Truncation coefficient = ',
     &      COEF_TRONCA,' Value too weak : TRUNCATION IS CANCELED'
   
	 ENDIF
	 
C* Recopie de la fonction de phase en intensite avant la troncature
C* et re-initialisation des paramètres de la decomposition de la matrice
C* de phase
C-------------------------------------------------------------------
      	 DO 66 J=-MIE_NBMU,MIE_NBMU
            P11(J)=TTT(J)
  66  	 CONTINUE
   
         ITRONC=0
	 DO 67 K=0,OS_NB
             BETA11(K)=0.D+00
	     BETA22(K)=0.D+00
             GAMMA12(K)=0.D+00
             DELTA33(K)=0.D+00
             ALP(K)=0.D+00
             ZETA(K)=0.D+00
 67      CONTINUE
  
	 GOTO 59 
      ENDIF    !Fin test si troncature trop faible.
  

C* Calcul des paramètres GAMMA12(k), BETA22(K) et DELTA33(k) de la decomposition des 
C* fonctions de phase polarisees P12, P22 et P33 en polynomes ou fonctions de Legendre.
C* Application de l'eventuelle troncature par le rapport de la fonction de phase
C* en intensité tronquée sur celle initiale: P11(J)/TTT(J).
C-----------------------------------------------------------------------------
      DO 11 J=-MIE_NBMU,MIE_NBMU       
	 IF(J.EQ.0)  GO TO 11
       
         XRMU=XMU(J)
	 
         PL(-1)=0.
         PL(0)=1.
	  POL(0)=0.
         POL(1)=0.
         POL(2)=3.*(1.-XRMU**2)/2./DSQRT(6.0D+00)
	 
	  XXX=XHR(J)*P12(J)*P11(J)/TTT(J)
         DO 12 K=2,OS_NB
            D=(2.*K+1.)/DSQRT(1.D+00*(K+3.)*(K-1.))
            E=DSQRT(1.D+00*(K+2.)*(K-2.))/(2.*K+1.)
            POL(K+1)=D*(XRMU*POL(K)-E*POL(K-1))
            GAMMA12(K)=GAMMA12(K)+XXX*POL(K)
  12     CONTINUE
    	
	 
         X=XHR(J)*P22(J)*(P11(J)/TTT(J))	!Parentheses necessaires pour assurer la non-regression stricte
         XX=XHR(J)*P33(J)*P11(J)/TTT(J)
         DO 14 K=0,OS_NB
            PL(K+1)=((2.*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
	     BETA22(K)=BETA22(K)+X*PL(K)
            DELTA33(K)=DELTA33(K)+XX*PL(K)
  14     CONTINUE
  
  11  CONTINUE
   
      DO 15 K=0,OS_NB
         BETA22(K)  = BETA22(K) *(2.*K+1.)*.5
         DELTA33(K) = DELTA33(K)*(2.*K+1.)*.5
         GAMMA12(K) = GAMMA12(K)*(2.*K+1.)*.5
   15 CONTINUE


C* Calcul des paramètres ALP(k) et ZETA(k) en fonctions des BETA22(k) et DELTA33(k).
C----------------------------------------------------------------------------------
      DO 16 I=2,OS_NB
         CO1=4*(2*I+1.)/I/(I-1.)/(I+1.)/(I+2.)
         CO2=I*(I-1.)/((I+1.)*(I+2.))
         CO3=CO2*DELTA33(I)
         CO2=CO2*BETA22(I)
         NN=I*.5
         MM=(I-1)*.5
         SOM1=0.
         SOM2=0.
         SOM3=0.
         SOM4=0.
         DO 17 J=1,NN
            X2=(I-1.)*(I-1.)-3.*(2*J-1.)*(I-J)
            SOM1=SOM1+X2*BETA22(I-2*J)
            SOM2=SOM2+X2*DELTA33(I-2*J)
  17     CONTINUE
  
         DO 18 J=0,MM
            X2=(I-1.)*(I-1.)-3.*J*(2*I-2*J-1.)
            SOM3=SOM3+X2*BETA22(I-2*J-1)
            SOM4=SOM4+X2*DELTA33(I-2*J-1)
  18     CONTINUE
  
         ZETA(I)=CO3-CO1*(SOM2-SOM3)
         ALP(I)=CO2-CO1*(SOM1-SOM4)
  16  CONTINUE

C* Normalisation des paramètres ALP(K), BETA(K), GAMMA(K), DELTA(K) et ZETA(K)
C* de la matrice de phase.
C----------------------------------------------------------------------------
      Z1=BETA11(0)
      DO 20 K=0,OS_NB
         ALP(K)=ALP(K)/Z1
         BETA11(K)=BETA11(K)/Z1
	  BETA22(K)=BETA22(K)/Z1
         GAMMA12(K)=GAMMA12(K)/Z1
         DELTA33(K)=DELTA33(K)/Z1
         ZETA(K)=ZETA(K)/Z1
  20  CONTINUE

C* Fin nominale 
C-------------------
       GOTO 9999

C* Cas d'erreur et retour du status -1 au programme appelant
C----------------------------------------------------------
  921  WRITE(6,*) 'Error while writing GRANU logfile'
       IER=-1
       GOTO 9999

       
 9999  CONTINUE
 
       RETURN
       END	!FIN DE LA PROCEDURE SOS_DECOMPO_LEGENDRE
       
       
       
       

C*==============================================================================
C* PROCEDURE: SOS_GRANU
C* ==========
C      Cette procédure calcule la fonction de phase P(mu) et les fonctions
C      de phase polarisee Q(mu) et U(mu) à partir des calculs de MIE realises
C      par le programme "SOS_MIE", pour une longueur d'onde donnée et pour
C      une distribution en taille des particules donnée (renseignee par igranu).
C
C      Si igranu = 1 => Distribution Log-Normale : 
C                N(r) = exp(-0.5 *(ln(r/rm)/sig)**2) / r
C		 varGranu1 = rm : rayon modal en microns,
C		 varGranu2 = sig : variance.
C
C      Si igranu = 2 => Loi de Junge : N(r) = Nr0 = R0**(-xnu) pour r < R0
C				           N(r) = r**(-xnu) pour r > R0
C				           N(r) = 0 pour r > Rmax
C		 varGranu1 = R0 : rayon minimal en microns,
C		 varGranu2 = xnu : puissance appliquée au rayon (pente de la loi de Junge).
C		 varGranu3 = Rmax : rayon maximal en microns,

C
C       Les fonctions de phase sont stockees dans des tableaux definis jusqu'à
C       CTE_MIE_NBMU_MAX mais sont utilisees jusqu'à MIE_NBMU (angles effectifs).
C
C
C Description des paramètres
C ----------------------------
C       FICMIE(CHARACTER*CTE_LENFIC2)  (E) : Nom du fichier MIE 
C                                 (répertoire + nom fichier + extension).
C       IGRANU  (I4)         (E) : Numéro du modèle de granulométrie à utiliser.
C       VARGRANU1 (double)   (E) : Paramètre 1 de la granulométrie.
C	 VARGRANU2 (double)   (E) : Paramètre 2 de la granulométrie
C	 VARGRANU3 (double)   (E) : Paramètre 2 de la granulométrie
C	 WA  (double)         (E) : Longueur d'onde en microns.
C       MIE_NBMU (I4) 	 (E) : Nombre d'angles (positifs) effectivement utiles
C	 XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double) : (E) Cosinus des angles de Gauss.
C	 TRACE (logical)      (E) : = vrai,si écriture dans le fichier trace
C                                 (No d'unite logique du fichier trace =88).
C	 KMAT1 (double)       (S) : Section efficace d'extinction (micron^2).
C	 KMAT2 (double)       (S) : Section efficace de diffusion (micron^2).
C	 SOMME_NR (double)    (S) : Nombre de particules integre sur la granulométrie.
C       P11 (double)         (S) : Fonction de phase calculée à partir des 
C                                  fonctions de la théorie de MIE. 
C       P12 (double)         (S) : Elément P12 de la matrice de phase.
C	 P33 (double)         (S) : Elément P33 de la matrice de phase.
C       IER  (I4)            (S) : Code d'erreur =0 si pas d'erreur, =-1 sinon.   
C
C
C Autre résultat fourni par la procédure
C ----------------------------------------
C   Informations dans le fichier de trace (optionnel)
C   (No d'unite logique du fichier trace =88 , 
C    fichier cree/ouvert par le programe appelant)
C
C
C
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier de MIE
C       Fichier binaire non formatte, cree par SOS_MIE
C       1 enregistrement d'entete avec
C                Indice de réfraction (partie réelle)
C                Indice de réfraction (partie imaginaire)
C                Paramètre de taille limite
C                Nombre d'angles de Gauss positifs
C
C      1 enregistrement pour chaque valeur du paramètre de taille avec:
C        - Paramètre de taille (réel),
C        - Coefficient d'efficacite en extinction (réel),
C        - Coefficient d'efficacite en diffusion (réel),
C        - Facteur d'anisotropie (double précision),
C        - Tableau de la fonction de phase Imie(mu) (réel),
C        - Tableau de la fonction de phase polarisee Qmie(mu) (réel),
C        - Tableau de la fonction de phase polarisee Umie(mu) (réel).
C
C    --> Les tableaux ont une taille utile limitée à MIE_NBMU : 
C        ils sont lus avec cette limitation.
C
C
C Constantes utilisees :  
C----------------------
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    INCTE_PI : valeur de PI
C    CTE_MIE_NBMU_MAX : Nombre maximal d'angles positifs.
C    CTE_OS_NB_MAX :  Valeur maximale pour l'ordre limite du développement en 
C                     polynomes de Legendre.
C
C Common utilise:
C --------------
C  Aucun
C
C
C Cas d'erreur :
C ------------
C     -Incoherence sur le nombre d'angles de Gauss
C     -Ouverture d'un fichier
C       Erreur si le fichier de MIE ou des angles de Gauss n'existe pas
C     -Lecture/écriture dans un fichier
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C        
C   
C
C
C* MOD:VERSION 2.0: 24/06/2008 
C	   - Les coefficients KMAT1 et KMAT2 sont normalisés à une particule
C            ==> La routine sort des sections efficaces d'extinction et de diffusion
C                (en mic^2) au lieu de coefficient d'extinction et de diffusion (/mic).
C          - Introduction du calcul de rayon effectif et variance effective.
C
C MOD:VERSION 3.0: 10/01/2010
C          - Passage d'une dimension des tableaux de CTE_MIE_NBMU à CTE_MIE_NBMU_MAX
C          - Boucle sur les fonctions de phase jusqu'à MIE_NBMU
C          - Modification de la ligne de lecture  des tableaux de fonctions de phase 
C            par valeur de paramètre de taille Alpha
C            ==> limitation aux bornes utiles -MIE_NBMU:MIE_NBMU de l'enregistrement effectif
C                des fichiers de MIE
C
C MOD:VERSION 4.0: 03/02/2016
C          - La formulation N(r) = exp(-0.5 *(log10(r/rm)/sig)**2) / r
C            pour sig = log10(variance) est abandonnée en faveur 
C            d'une formulation classique du modèle LND:
C			    N(r) = exp(-0.5 *(log(r/rm)/sig)**2) / r	
C	     pour sig = variance.
C
C MOD:VERSION 4.1: 08/10/2020
C          - Ajout de l'argument VARGRANU3 pour définir le rayon Rmax d'une 
C            granulométrie de Junge
C          
C==============================================================================


      SUBROUTINE SOS_GRANU(FICMIE,IGRANU,VARGRANU1,VARGRANU2,
     &                     VARGRANU3,WA,MIE_NBMU,XMU,TRACE,
     &                     KMAT1,KMAT2,SOMME_NR,P11,P12,P33,IER)

      IMPLICIT NONE


C* Definition des variables
C----------------------------------------------------------------------------
      
      CHARACTER*CTE_LENFIC2 FICMIE	! Chemin et nom du fichier MIE
      
      DOUBLE PRECISION VARGRANU1   	! Paramètre 1 de la granulométrie.
      DOUBLE PRECISION VARGRANU2   	! Paramètre 2 de la granulométrie.
      DOUBLE PRECISION VARGRANU3   	! Paramètre 3 de la granulométrie.
      DOUBLE PRECISION RM		!Rayon modal d'une granulométrie LND.
      DOUBLE PRECISION SIG		!Variance d'une granulométrie LND
      DOUBLE PRECISION B		!Variable intermedaire.
      DOUBLE PRECISION R0		!Rayon seuil min d'une granulométrie en loi de Junge
      DOUBLE PRECISION RMAX		!Rayon seuil max d'une granulométrie en loi de Junge
      DOUBLE PRECISION XNU		!Puissance d'une granulométrie en loi de Junge.
      DOUBLE PRECISION NR0		!Nombre de particules pour un rayon < r0 
                                	!d'une loi de Junge.
					
      DOUBLE PRECISION REFF		!Rayon effectif de la granulométrie (microns).
      DOUBLE PRECISION VEFF		!Variance effective de la granulométrie.
      
      DOUBLE PRECISION WA		!Longueur d'onde en microns.

      DOUBLE PRECISION XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)
                                        !Cosinus des angles de Gauss.
      DOUBLE PRECISION TETA     	!Angle de diffusion.
      
      DOUBLE PRECISION KMAT1		! Section efficace d'extinction (micron^2).
      DOUBLE PRECISION KMAT2 		! Section efficace de diffusion (micron^2). 
      
      DOUBLE PRECISION P11(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Fonction de phase calculée à partir des fonctions de 
		       !la théorie de MIE.
      DOUBLE PRECISION P12(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P12 de la matrice de phase.
      DOUBLE PRECISION P33(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Elément P33 de la matrice de phase.
		      
      DOUBLE PRECISION XG	!Facteur d'anisotropie.


      DOUBLE PRECISION ALPHAF	!Valeur maximale du paramètre de taille contenue
      				!dans le fichier des calculs de MIE.

      DOUBLE PRECISION R	!Rayon de particule pour le paramètre de taille 
      				!ALPHA et la longueur d'onde du traitement WA.
      DOUBLE PRECISION PR	!Pas entre deux valeurs successives des rayons de
      				!particules (en microns).
      DOUBLE PRECISION NR	!NR*PR = Nombre de particules de rayon compris
      				!entre R et R + PR (par micron**3).

      DOUBLE PRECISION SOMME_NR	!Integration sur les rayons de NR * PR 
      				!==> nb de particule par unite de volume.
      DOUBLE PRECISION SOMME_R2	!Integration sur les rayons de R^2 * NR * PR.
      DOUBLE PRECISION SOMME_R3	!Integration sur les rayons de R^3 * NR * PR.
      DOUBLE PRECISION SOMME_R4	!Integration sur les rayons de R^4 * NR * PR.
      DOUBLE PRECISION R2NRPR	!Valeur de R*R * NR * PR.
      
      DOUBLE PRECISION RN	!partie réelle de l'indice de réfraction.
      DOUBLE PRECISION RI	!partie imaginaire de l'indice de réfraction.
	
      DOUBLE PRECISION X1		

      
      REAL*4 ALPHA		!Valeur du paramètre de taille.
      REAL*4 PAS		!Pas entre deux valeurs de ALPHA.
      REAL*4 QEXT		!Coefficient d'efficacite en extinction.
      REAL*4 QSCA		!Coefficient d'efficacite en diffusion.
      REAL*4 IMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)         
                                !Fonction de phase : Imie(mu) pour ALPHA.
      REAL*4 QMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)         
                                !Fonction de phase polarisee : Qmie(mu).
      REAL*4 UMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)         
                                !Fonction de phase polarisee : Umie(mu).

      INTEGER*4 MIE_NBMU  	! Nombre d'angles (positifs) effectivement utiles 
      INTEGER*4 IGRANU  	! Numéro du modèle de granulométrie.
      INTEGER*4 J	   	! Indice.
      INTEGER*4 NBMUMIE	   	! Nombre d'angle de GAUSS positifs utilises pour 
      			   	! les calculs de MIE (enregistrement du fichier MIE).			   
      INTEGER*4 IER        	! Code d'erreur =0 si pas d'erreur, =-1 sinon

      LOGICAL TRACE        	! = vrai,si écriture dans le fichier trace.
      

C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0

C* Initialisation des tableaux de fonction de phase
C------------------------------------------------------
      DO 102 J=-MIE_NBMU,MIE_NBMU
         P11(J)=0.
         P12(J)=0.
         P33(J)=0.
  102 CONTINUE

C* Initialisation des coefficient d'extinction et de diffusion
C-------------------------------------------------------------
      KMAT1=0.
      KMAT2=0.

C* Initialisation des cumuls d'integration sur la granulométrie
C-------------------------------------------------------------
      SOMME_NR=0.
      SOMME_R2=0.
      SOMME_R3=0.
      SOMME_R4=0.
                  
C* Ouverture des fichiers 
C------------------------------------------------------
C     Fichier MIE
      OPEN (1,FILE=FICMIE,FORM='UNFORMATTED',
     &      STATUS='OLD',ERR=991)
 
C* Lecture de l'entete du fichier de MIE
C-------------------------------------------      
      READ(1,ERR=992) RN,RI,ALPHAF,NBMUMIE
        
C* Controle paramètre
C--------------------
      IF (NBMUMIE.NE.MIE_NBMU) GOTO 996

C* Lecture des fonctions calculées par la théorie de MIE
C* pour chaque valeur du paramètre de taille ALPHA
C*
C*      Lecture des tableaux de fonctions de phase
C*      qui ont ete enregistres avec des bornes 
C*      limitées au dimensionnement utile -MIE_NBMU:MIE_NBMU
C*       
C  (boucle sur l'etiquette 123,
C   condition d'arret : fin du fichier MIE)
C------------------------------------------------------
      PAS = 0.0001

 123     READ(1,err=9992) ALPHA,QEXT,QSCA,XG,
     &                    (IMIE(J),J=-MIE_NBMU,MIE_NBMU),
     &                    (QMIE(J),J=-MIE_NBMU,MIE_NBMU),
     &                    (UMIE(J),J=-MIE_NBMU,MIE_NBMU)
		
C* Calcul du rayon de particule associe a ALPHA pour la longueur d'onde WA
C--------------------------------------------------------------------------
         R=ALPHA*WA/2.D+00/INCTE_PI
	 
C* Test sur la limite du fichier de MIE
C--------------------------------------
         IF (ALPHA.GE.(ALPHAF-PAS)) GOTO 40
		
C* Estimation du pas pour la valeur de ALPHA suivante
C------------------------------------------------------
         IF (ALPHA.GT.0.10)  PAS=0.001
         IF (ALPHA.GT.1.00 ) PAS=0.01
         IF (ALPHA.GT.10.)   PAS=0.05
         IF (ALPHA.GT.30.)   PAS= 0.10
         IF (ALPHA.GT.100.)  PAS=1.00


C* Expression des différentes granulométries possibles:
C------------------------------------------------------

C*  Distribution log-normale :
C***************************
C*  RM : Rayon modal (microns).
C*  SIG : Variance
C------------------------------
		
         IF (IGRANU.eq.1) THEN
            RM = VARGRANU1
            SIG= VARGRANU2

            B = dlog(R/RM)/SIG
            NR = dexp(-B*B/2.)/(R*SIG*sqrt(2*INCTE_PI))
         ENDIF

C*  Loi de Junge : NR = plateau Rmin**(-XNU) pour r <= Rmin
C****************
C*  R0 : Rayon minimal (microns).
C---------------------------
         IF (IGRANU.eq.2) THEN
            R0  = VARGRANU1
            XNU = VARGRANU2
            RMAX = VARGRANU3

            NR0= (R0)**(-XNU)

	     IF (R.GT.RMAX) GOTO 40 !Sortie de boucle sur taille des particules
            IF (R.LE.R0) THEN
	       NR = NR0
	     ELSE
	       NR = R**(-XNU)   
	     ENDIF            
         ENDIF


C* Calcul des coefficients d'extinction (kMAT1) et de diffusion (KMAT2).
C* Calcul des fonctions de phase en intensite (P11) et en polarisation(P12 et P33).
C-------------------------------------------------------------------------------
         PR=WA*PAS/2./INCTE_PI
         X1=NR*PR*INCTE_PI*R**2
         KMAT1=KMAT1+X1*QEXT
	  KMAT2=KMAT2+QSCA*X1
	
         X1=QSCA*X1
         DO 1 J=-MIE_NBMU,MIE_NBMU
            P11(J)=P11(J)+IMIE(J)*X1
            P12(J)=P12(J)+QMIE(J)*X1
            P33(J)=P33(J)+UMIE(J)*X1
   1     CONTINUE
   
C* Calcul du nombre de particules total
C-------------------------------------------------------------------------------
      SOMME_NR = SOMME_NR + NR*PR
 
C* Calcul des termes pour le rayon effectif et volume moyen
C-------------------------------------------------------------------------------
      R2NRPR = R*R*NR*PR
      SOMME_R2 = SOMME_R2 + R2NRPR
      SOMME_R3 = SOMME_R3 + R2NRPR*R
      SOMME_R4 = SOMME_R4 + R2NRPR*R*R
        
C* Passage a la valeur suivante de ALPHA
C-----------------------------------------
      GOTO 123 	
   40 CONTINUE

C* Calcul du rayon effectif et de la variance effective
C-------------------------------------------------------------------------------      
      IF (TRACE) THEN
         REFF = SOMME_R3 / SOMME_R2
	  VEFF =  1.D+00 +  
     &	         (SOMME_R4 - 2.*REFF*SOMME_R3) / (REFF*REFF * SOMME_R2)
      END IF  

C* Calcul du volume d'une particule moyenne de la granulométrie
C-------------------------------------------------------------------------------    
c      VMOY = (4.D+00*CTE_PI/3.D+00)*SOMME_R3 / SOMME_NR
      
C* Normalisation des fonctions de phase.
C* Copie de la fonction de phase en intensite.
C------------------------------------------------------
      DO 55 J=-MIE_NBMU,MIE_NBMU
         P11(J)=P11(J)/KMAT2
         P12(J)=P12(J)/KMAT2
         P33(J)=P33(J)/KMAT2
  55  CONTINUE
  

C* Normalisation à une particule des coefficients 
C* de diffusion et d'extinction.
C* ==> Conversion en sections efficaces de diffusion
C*     et d'extinction
C------------------------------------------------------  
       KMAT1= KMAT1 / SOMME_NR
       KMAT2= KMAT2 / SOMME_NR

       
C* Ecriture du fichier trace
C-------------------------------------------
      IF (TRACE) THEN
        WRITE(88,*,err=921) ' '
        WRITE(88,*,err=921) '* SUBROUTINE GRANU'
	 WRITE(88,*,err=921) '=================='
        WRITE(88,*,err=921) ' '
        WRITE(88,*,err=921) 'SIMULATION PARAMETERS'
	 WRITE(88,*,err=921) '---------------------'
	 WRITE(88,*,err=921) '  MIE file : ',FICMIE
	 WRITE(88,*,err=921) '  Wavelength in microns : ',WA
	 IF (IGRANU.EQ.1) THEN
	      WRITE(88,*,err=921)
     &	           '  Model 1 : Log-Normale Distribution'
	      WRITE(88,*,err=921) '     rm in microns : ',VARGRANU1
	      WRITE(88,*,err=921) '     sig : ',VARGRANU2
	      WRITE(88,*,err=921) '     sig/ln10 : ',
     &	                                VARGRANU2/DLOG(10.D+00)
	 ELSE
	      WRITE(88,*,err=921) "  Model 2 : Junge's Law "
	      WRITE(88,*,err=921) '     R0 in microns : ',VARGRANU1
	      WRITE(88,*,err=921) '     Slope  : ',VARGRANU2	      	   
	 ENDIF
        WRITE(88,*,err=921) ' '
        WRITE(88,*,err=921) 'RESULTS'
	 WRITE(88,*,err=921) '-------'
	 WRITE(88,*,err=921)
     &	    ' Extinction coefficient - no truncation applied'//
     &	    ' (undefined unit: inverse dimension to a length) :',
     &	    KMAT1*SOMME_NR
	 WRITE(88,*,err=921) 
     &	    ' Scattering coefficient - no truncation applied'//
     &	    ' (undefined unit: inverse dimension to a length) :',
     &     KMAT2*SOMME_NR
        WRITE(88,*,err=921) 
     &	    ' Number of particles (undefined unit: dimension'//
     &	    ' proportional to a length^3) :', SOMME_NR
        WRITE(88,*,err=921)
     &	    ' Extinction cross-section (mic^2) - no trunca. applied : '
     &     ,KMAT1
	 WRITE(88,*,err=921) 
     &	    ' Scattering cross-section (mic^2) - no trunca. applied : ' 
     &     ,KMAT2
	 WRITE(88,*,err=921) 
     &	    ' Effective radius of the size distribution (microns) : ',
     &	    REFF
        WRITE(88,*,err=921) 
     &	    ' Effective variance of the size distribution : ',VEFF
        WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) '  1st column : Scattering angle'
	 WRITE(88,*,err=921) '  2nd column : Original phase function'
        WRITE(88,*,err=921) '                 P11(mu) '
	 WRITE(88,*,err=921) 
     &       '  3rd column : Original polarized phase function'
        WRITE(88,*,err=921) '                 P12(mu) '
	 WRITE(88,*,err=921) 
     &       '  4th column : Original polarized phase function '
        WRITE(88,*,err=921) '                 P33(mu) '
	 WRITE(88,*,err=921) ' '
        DO 44 J=-MIE_NBMU,MIE_NBMU
	    TETA=DACOS(XMU(J))*180./INCTE_PI
	    WRITE(88,10,err=921) TETA,P11(J),P12(J),P33(J)

   44   CONTINUE
   	 WRITE(88,*,err=921) ' '

      END IF  

 


C* Fermeture fichier
C-------------------
      CLOSE(1)
      
C* Fin nominale 
C-------------------
      GOTO 9999

      
C* Cas d'erreur et retour du status -1 au programme appelant
C----------------------------------------------------------
  991 WRITE(6,*) 'Erreur a l''ouverture du fichier MIE'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'Erreur de lecture de la premiere ligne du fichier MIE'
      WRITE(6,*) '  --> Fichier de MIE probablement incomplet'
      WRITE(6,*) '      ==> Controlez l''espace disque'
      IER=-1
      GOTO 9999
 9992 WRITE(6,*) 'Erreur de lecture sur le fichier MIE'
      WRITE(6,*) '  --> Probable tentative de lecture apres fin fichier'
      WRITE(6,*) '  --> Fichier de MIE probablement incomplet'
      WRITE(6,*) '      ==> Controlez l''espace disque'
      IER=-1
      GOTO 9999    
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace GRANU'
      IER=-1
      GOTO 9999      
  996 WRITE(6,*) 'Incoherence sur le nombre d''angles de Gauss'
      write(6,*) 'NBMUMIE=',NBMUMIE,' MIE_NBMU=',MIE_NBMU
      IER=-1  
      GOTO 9999

C* Format
C---------   
   10 FORMAT(1X,F6.2,6(2X,E12.4))
   
         
 9999 CONTINUE
 
  
      RETURN
      END    !FIN DE LA PROCEDURE SOS_GRANU
      
      
      
      
      


C*==============================================================================
C* PROCEDURE: SOS_INTERPO_SPLINT
C* ==========
C      Cette procédure interpole une fonction definie par (XIN, YIN) sur NBIN noeuds
C      aux abscisses XOUT pour NBOUT noeuds (XOUT, YOUT).
C
C      La demarche repose sur l'utilisation des fonctions SPLINE et SPLINT
C         - Les noeuds (XIN, YIN) sont en premier tries par ordre croissant d'abscisse XIN.
C         - La fonction SPLINE calcule les derivees secondes de YIN_TRIE en XIN_TRIE
C           à partir des valeurs YIN_TRIE de la fonction aux abscisses XIN_TRIE et des
C           valeurs de derivee premiere de la fonction au premier et dernier noeud 
C           (min et max de XIN_TRIE).
C	  - La fonction SPLINT calcule les valeurs YOUT aux abscisses XOUT en fonction
C           des valeurs de la fonction et de sa derivee seconde predecemment calcules
C           aux abscisses XIN_TRIE.
C
C Description des paramètres
C ----------------------------
C       NBIN (I4)		(E) : Nombre d'elements des tableaux XIN et YIN
C       XIN(1:NBIN)   (double)	(E) : Tableau des abscisses de la fonction pour lesquels l'ordonnée est connue.
C       YIN(1:NBIN)   (double)	(E) : Tableau des valeurs connues de la fonction (correspond aux abscisses XIN).
C       NBOUT (I4)		(E) : Nombre d'elements des tableaux XOUT et YOUT
C       XOUT(1:NBOUT) (double)	(E) : Tableau des abscisses souhaites de la fonction.
C       YOUT(1:NBOUT) (double)	(S) : Tableau des valeurs estimees de la fonction aux abscisses XOUT.
C       IER  (I4)               (S) : Code d'erreur =0 si pas d'erreur, =-1 sinon.   
C
C Constantes utilisees :  
C----------------------
C  Aucune
C
C Common utilise:
C --------------
C  Aucun
C
C Cas d'erreur :
C ------------
C     - Retour d'erreur de la routine SOS_SPLINT.
C
C     Affichage d'un message à l'écran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C        
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel aux routines: SOS_SPLINE et SOS_SPLINT
C
C==============================================================================

      SUBROUTINE SOS_INTERPO_SPLINT(NBIN,XIN,YIN,NBOUT,XOUT,YOUT,IER)

      IMPLICIT NONE


C* Definition des variables
C----------------------------------------------------------------------------
      INTEGER*4 NBIN  	! Nombre de noeuds (XIN,YIN) de la fonction Y 
      INTEGER*4 NBOUT  	! Nombre de noeuds (XOUT,YOUT) de la fonction Y 
      INTEGER*4 J	! Indice.
      INTEGER*4 K	! Indice.	   
      INTEGER*4 IER     ! Code d'erreur =0 si pas d'erreur, =-1 sinon
      
      DOUBLE PRECISION XIN(NBIN)   !Valeurs des abscisses des noeuds (XIN,YIN) 
      DOUBLE PRECISION YIN(NBIN)   !Valeurs des ordonnées des noeuds (XIN,YIN) 
      DOUBLE PRECISION XOUT(NBOUT)  !Valeurs des abscisses des noeuds (XOUT,YOUT) 
      DOUBLE PRECISION YOUT(NBOUT)  !Valeurs recherchees des ordonnées des noeuds (XOUT,YOUT) 

      DOUBLE PRECISION XIN_TRIE(NBIN)   !Valeurs croissantes des abscisses des noeuds (XIN,YIN) 
      DOUBLE PRECISION YIN_TRIE(NBIN)   !Valeurs des ordonnées des noeuds (XIN,YIN) pour XIN trie par ordre croissant
      DOUBLE PRECISION VAL		!Valeur temporaire
      
      DOUBLE PRECISION DY1	!Derivee d'ordre 1 de YIN_TRIE en abscisse XIN_TRIE(1)
      DOUBLE PRECISION DYN	!Derivee d'ordre 1 de YIN_TRIE en abscisse XIN_TRIE(NBIN)
      DOUBLE PRECISION DSEC_YIN_TRIE(NBIN)  !Valeurs des derivees secondes de la fonction YIN_TRIE 
      					    !aux noeuds XIN_TRIE  
    

C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0

      
C* Trie des tables (XIN, YIN) par ordre croissant de XIN
C------------------------------------------------------
      DO J=1,NBIN
         XIN_TRIE(J)=XIN(J)
	 YIN_TRIE(J)=YIN(J)
      ENDDO
      
      DO J=1,NBIN	     
	 DO K=J+1,NBIN
	    IF (XIN_TRIE(J).GT.XIN_TRIE(K)) THEN		    		    
	       !Permutation des valeurs XIN(J) et XIN(K)
	        VAL=XIN_TRIE(J)
	        XIN_TRIE(J)=XIN_TRIE(K)
	        XIN_TRIE(K)=VAL
	       !Permutation des valeurs YIN(J) et YIN(K)
		 VAL=YIN_TRIE(J)
		 YIN_TRIE(J)=YIN_TRIE(K)
		 YIN_TRIE(K)=VAL   	    	    
	    ENDIF
	 ENDDO 
      ENDDO 
   
      
C* Calcul des derivees secondes de YIN_TRIE en XIN_TRIE
C------------------------------------------------------
      DO J=1,NBIN
         DSEC_YIN_TRIE(J)=0.D+00
      ENDDO

C     Derivees d'ordre 1 de YIN_TRIE en XIN_TRIE(1) et XIN_TRIE(NBIN)
      DY1 = (YIN_TRIE(2) - YIN_TRIE(1)) / (XIN_TRIE(2) - XIN_TRIE(1))
      DYN = (YIN_TRIE(NBIN) - YIN_TRIE(NBIN-1)) / 
     &      (XIN_TRIE(NBIN) - XIN_TRIE(NBIN-1))
      
      CALL SOS_SPLINE(NBIN,XIN_TRIE,YIN_TRIE,DY1,DYN,DSEC_YIN_TRIE)

C* Application pour l'interpolation SPLINT aux abscisses XOUT
C------------------------------------------------------------       
       DO J=1,NBOUT
	 
	 CALL SOS_SPLINT(NBIN,XIN_TRIE,YIN_TRIE,DSEC_YIN_TRIE,
     &	                 XOUT(J),YOUT(J),IER)
     	 
         IF (IER.NE.0) GOTO 901
      ENDDO
            
      
C* Fin nominale 
C-------------------
      GOTO 9999

      
C* Cas d'erreur et retour du status -1 au programme appelant
C----------------------------------------------------------
  901 WRITE(6,*) '  ERROR for SOS_INTERPO_SPLINT'
      GOTO 9999  
     
 9999 CONTINUE
  
      RETURN
      END    !FIN DE LA PROCEDURE SOS_INTERPO_SPLINT     
      



C*==============================================================================
C* PROCEDURE: SOS_SPLINE
C* ==========
C      Cette procédure calcule les valeurs de derivees secondes d'une fonction definie
C      par (X, Y) sur N noeuds (supposes tries par ordre croissant de X)
C      et par ses valeurs de derivees premieres aux noeuds extremes (X(1),DY1) et (X(N),DYN).
C
C      La routine renvoie les valeurs de derivees secondes en chaque point X.
C
C Description des paramètres
C ----------------------------
C       N  (I4)		(E) : Nombre d'elements des tableaux X, Y et DSEC_Y
C       X(1:N) (double)	(E) : Tableau des abscisses de la fonction pour lesquels l'ordonnée est connue.
C       Y(1:N) (double)	(E) : Tableau des valeurs connues de la fonction (correspond aux abscisses X).
C       DY1    (double)	(E) : Derivee premiere de la fonction en X(1).
C       DYN    (double)	(E) : Derivee premiere de la fonction en X(N).
C       DSEC_Y(1:N) (double) (S) : Tableau des valeurs de derivees secondes de la fonction aux abscisses X.
C
C Constantes utilisees :  
C----------------------
C  Aucune    
C
C Common utilise:
C --------------
C  Aucun
C
C Cas d'erreur :
C ------------
C  Aucun         
C
C==============================================================================

      SUBROUTINE SOS_SPLINE(N,X,Y,DY1,DYN,DSEC_Y)

      IMPLICIT NONE

      INTEGER*4 N		 !Nb d'elements des tableaux d'entrée (X,Y) et et de sortie (DSEC_Y)
      DOUBLE PRECISION X(N)	 !Tableau des abscisses de la fonction 
      DOUBLE PRECISION Y(N)	 !Tableau des valeurs de la fonction aux abscisses X
      DOUBLE PRECISION DY1   !Derivee premiere de la fonction en X(1)
      DOUBLE PRECISION DYN	 !Derivee premiere de la fonction en X(N)
      DOUBLE PRECISION DSEC_Y(N) !Tableau des valeurs de derivees secondes de la fonction aux abscisses X
      
      DOUBLE PRECISION U(N)

      DOUBLE PRECISION SIG
      DOUBLE PRECISION P
      DOUBLE PRECISION QN
      DOUBLE PRECISION UN
      
      INTEGER*4 K ! Indice.
      
C* Definition des variables
C----------------------------------------------------------------------------


      IF (DY1.GT..99E30) THEN
        DSEC_Y(1)=0.
        U(1)=0.
      ELSE
        DSEC_Y(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-DY1)
      ENDIF
      DO 11 K=2,N-1
        SIG=(X(K)-X(K-1))/(X(K+1)-X(K-1))
        P=SIG*DSEC_Y(K-1)+2.
        DSEC_Y(K)=(SIG-1.)/P
        U(K)=(6.*((Y(K+1)-Y(K))/(X(K+1)-X(K))-(Y(K)-Y(K-1))
     *      /(X(K)-X(K-1)))/(X(K+1)-X(K-1))-SIG*U(K-1))/P
11    CONTINUE
      IF (DYN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(DYN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      DSEC_Y(N)=(UN-QN*U(N-1))/(QN*DSEC_Y(N-1)+1.)
      DO 12 K=N-1,1,-1
        DSEC_Y(K)=DSEC_Y(K)*DSEC_Y(K+1)+U(K)
12    CONTINUE
      RETURN
      END   !FIN DE LA PROCEDURE SOS_SPLINE   





C*==============================================================================
C* PROCEDURE: SOS_SPLINT
C* ==========
C
C      Cette procédure interpole une fonction en un point XVAL à partir
C      de trois vecteurs de même dimension X(1:N), Y(1:N) et DSEC_Y(1:N) 
C      qui donnent respectivement les abscisses, les ordonnées et les 
C      dérivées secondes de la fonction.
C
C Description des paramètres
C ----------------------------
C       N  (I4)		(E) : Nombre d'elements des tableaux X, Y et DSEC_Y
C       X(1:N) (double)	(E) : Tableau des abscisses de la fonction pour lesquels l'ordonnée est connue.
C       Y(1:N) (double)	(E) : Tableau des valeurs connues de la fonction (correspond aux abscisses X).
C       DSEC_Y(1:N) (double) (E) : Tableau des valeurs de derivees secondes de la fonction aux abscisses X.
C       XVAL (double)	(E) : Valeur de l'abscisse pour laquelle la fonction doit etre interpolee.
C       YVAL (double)	(S) : Valeur de la fonction interpolee en XVAL.
C       IER  (I4)       (S) : Code d'erreur =0 si pas d'erreur, =-1 sinon.   
C
C Constantes utilisees :  
C----------------------
C  Aucune    
C
C Common utilise:
C --------------
C  Aucun
C
C Cas d'erreur :
C ------------
C  Identite des abscisses de la fonction pour les indices estimant les valeurs 
C  qui encadrent l'abscisse souhaitee.       
C
C==============================================================================

      SUBROUTINE SOS_SPLINT(N,X,Y,DSEC_Y,XVAL,YVAL,IER)

      IMPLICIT NONE

      INTEGER*4 N	!Nb d'elements des tableaux d'entrée (X,Y et DSEC_Y)
      INTEGER*4 KLO	!Lower value of index
      INTEGER*4 KHI	!Higher value of index   
      INTEGER*4 K	! Indice
      INTEGER*4 IER     ! Code d'erreur =0 si pas d'erreur, =-1 sinon
      
      DOUBLE PRECISION X(N)	 !Tableau des abscisses de la fonction 
      DOUBLE PRECISION Y(N)	 !Tableau des valeurs de la fonction aux abscisses X
      DOUBLE PRECISION DSEC_Y(N) !Tableau des valeurs de derivees secondes de la fonction aux abscisses X

      DOUBLE PRECISION XVAL	!Valeur de l'abscisse pour laquelle la fonction doit etre interpolee.
      DOUBLE PRECISION YVAL	!Valeur de la fonction interpolee en XVAL.
      
      DOUBLE PRECISION A,B
      DOUBLE PRECISION H


C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0
      
C* Definition des variables
C----------------------------------------------------------------------------
      
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(X(K).GT.XVAL)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=X(KHI)-X(KLO)
      IF (H.NE.0.) THEN
         A=(X(KHI)-XVAL)/H
         B=(XVAL-X(KLO))/H
         YVAL=A*Y(KLO)+B*Y(KHI)+
     *      ((A**3-A)*DSEC_Y(KLO)+(B**3-B)*DSEC_Y(KHI))*(H**2)/6.
      ELSE
         IER=-1
      ENDIF
      
      IF (IER.EQ.-1) GOTO 901
      
C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status -1 au programme appelant
C----------------------------------------------------------
  901 WRITE(6,*) '  ERROR for SPLINT interpolation'
      WRITE(6,*) '  --> Bad X table input'
      GOTO 9999  
     
 9999 CONTINUE      
      RETURN
      END    !FIN DE LA PROCEDURE SOS_SPLINT     
