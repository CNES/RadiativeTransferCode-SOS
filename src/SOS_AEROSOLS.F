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
C* FICHIER: SOS_AEROSOLS.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module gere le calcul des proprietes radiatives de 
C*       differentes granulometrieF d'aerosols.
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et developpement : CS SI, B.Lafrance, C.Ruffel.
C*
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.1: Mise a niveau du code des OS du LOA
C* MOD:VERSION:1.2: Mise a jour des commentaires.
C*
C* MOD:VERSION:2.0: 24/08/2008 
C*       - Introduction d'une modelisation bi-LND.
C*       - Conversion des sorties KMAT1 et KMAT2 de la routine SOS_GRANU
C*         en sections efficaces d'extinction et de diffusion (micron^2)
C*         au lieu de coefficients d'extinction et de diffusion (/mic)
C*         qui n'etaient pas normalises à une particule.
C*       - Mise a jour des commentaires.
C*
C*
C* MOD:VERSION:2.1: 08/09/2008 
C*       - Correction de depassements de zone d'indentation pour compilation f77.
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
C*       - Introduction du parametre SOS_MIE_ALPHAMIN qui precise la valeur minimale
C*         du parametre de taille ALPHAO des calculs de Mie.
C*
C* MOD:VERSION:4.0: 10/01/2010 
C*   - Suppression de la constante SOS_OS_FICANGLE (de SOS.h) 
C*     au profit d'un passage par argument du nom du fichier des angles (produit par SOS_ANGLES).
C*
C*   - Modification de la taille des tableaux de luminance :
C*     ==> Passage de SOS_MIE_NBMU a SOS_MIE_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'a MIE_NBMU egal au nombre effectif d'angles utilises.
C*                                   
C*   - Modification de la taille des tableaux des fonctions de phase :
C*     ==> Passage de SOS_MIE_NBMU a SOS_MIE_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'a MIE_NBMU egal au nombre effectif d'angles utilises.
C*
C*   - Modification de la taille des tableaux de decomposition en polynomes de Legendre :
C*     ==> Passage de SOS_OS_NB a SOS_OS_NB_MAX.
C*     ==> Boucle sur les ordre jusqu'a OS_NB (lu dans le fichier des angles)
C*
C*   - SOS_MIE : Ajout des arguments d'entree : NBMU_MIE (nb angles utiles) et XMU (tableau des angles)
C*
C*   - SOS_NOMFIC_MIE : Ajout des arguments d'entree : NBMU_GAUSS (nb angles de Gauss) 
C*                                                     et FICANGLES_USER (fichier d'angles utilisateurs)
C*
C*   - SOS_GRANU : 
C*         * Ajout de l'argument d'entree : NBMU_MIE (nb angles utiles)
C*         * Modification de la ligne de lecture  des tableaux 
C*           de fonctions de phase par valeur de parametre de taille Alpha
C*           ==> limitation aux bornes utiles -MIE_NBMU:MIE_NBMU
C*
C*   - SOS_DECOMPO_LEGENDRE : 
C*         * Ajout des arguments d'entree : NBMU_MIE (nb angles utiles) et OS_NB
C*         * Boucles limitees a NBMU_MIE au lieu SOS_MIE_NBMU_MAX
C*  
C*   - Suppression de l'utilisation de la variable d'environnement SOS_RACINE_MU
C*     ==> devient SOS_RACINE_FIC car utilisee seulement pour les fichiers WMO et S&F.
C*
C*   - Ajout d'une option de traitement IMOD=4 pour ingerer une matrice de phase externe
C*     --> Fonctions de phase externes et coefficients de diffusion et d'extinction
C*     ==> Introduction des routines : 
C*            	- SOS_INTERPO_SPLINT
C*     		- SOS_SPLINE
C*    		- SOS_SPLINT
C*     ==> Ajout de la constante (SOS.h) :  SOS_MAXNB_ANG_EXT
C*
C* MOD:VERSION:4.1: 24/06/2010 
C*     Les etiquettes de FORMAT en ,X, passent en ,1X, pour eviter des erreurs de compilation.
C*   
C*
C* MOD:VERSION:4.2: 18/12/2015 
C*    - Suppression des formats à la lecture des paramètres d'entrée en faveur 
C*      de READ(PROG_ARG,*,ERR=NNN)
C*      ==> Permet une compilation sous gfortran et évite les mauvaises conversions 
C*          si la valeur d'entrée ne respecte pas strictement le format attendu 
C*          (par exemple, saisie d'un entier N au lieu d'un réel N. 
C*           aboutissant à une valeur interprétée fausse). 
C*
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*      (incluant au décodage des constante) : requis pour compilation gfortran
C*
C*    - Changement de la manière de lire les fichiers formatés (Angles, ...)
C*      pour être compatible avec le compilateur gfortran.
C*
C*    - Modification du format du fichier GRANU résultat pour permettre sa lecture
C*      par le code SOS.F compilé avec gfortran.
C*
C*    - Modification de la récupération du nom du fichier optionnel des 
C*      angles utilisateur (FICANGLES_USER) pour compatibilité avec gfortran/linux
C*
C*    - Modification de la gestion des variances pour les granulométrie LND:
C*      abandon de l'usage d'un SIG = log10(variance) au profit de l'utilisation
C*      d'une formulation classique du modèle LND en fonction de la variance.
C*		N(r) = exp(-0.5 *(log10(r/rm)/sig)**2) / r
C*	devient 
C*		N(r) = exp(-0.5 *(log(r/rm)/sig)**2) / r
C*      ==> Impacte les valeurs des paramètres en entrée de SOS_AEROSOLS et donc
C*          fournis en ligne de commandes.
C*
C*
C* MOD:VERSION:4.3: 06/09/2017
C*     Changes to be compliant with a f77 compilation
C*         - Change a variable type to read the file FICANGLES  
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    SOS_LENDIR : Longueur des noms de repertoires.
C    SOS_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_LENKEYWORD : Taille maximale des Keywords pour le passage des arguments.
C    SOS_PI : Valeur de PI.
C    SOS_MIE_NBMU_MAX : Nombre maximal d'angles positifs.
C    SOS_MIE_ALPHAMIN : Valeur minimale par defaut du parametre de taille 
C                       Alpha des calculs de Mie.
C    SOS_ALPHAMAX_WMO_DL : Valeur limite du parametre de taille pour le calcul
C                          des fichiers du composant Dust-Like de la WMO.
C    SOS_ALPHAMAX_WMO_WS : Idem pour le composant Water Soluble.
C    SOS_ALPHAMAX_WMO_OC : Idem pour le composant Oceanic.
C    SOS_ALPHAMAX_WMO_SO : Idem pour le composant Soot.
C    SOS_ALPHAMAX_SF_SR : Idem pour le composant Small Rural de Shettle & Fenn.
C    SOS_ALPHAMAX_SF_SU : Idem pour le composant Small Urban de Shettle & Fenn.
C    SOS_COEF_NRMAX : Constante pour la determination du parametre de taille max.
C    SOS_OS_NB_MAX :  Valeur maximale pour l'ordre limite du developpement en 
C                     polynomes de Legendre.
C    SOS_PH_SEUIL_TRONCA : Valeur seuil pour la troncature.
C    SOS_AER_DATAWMO : Nom du fichier contenant les parametres WMO (rm,sig,...).
C    SOS_AER_DATASF  : Nom du fichier contenant rayon modal et log variance S&F.
C    SOS_AER_XX_SF   : Nom des fichiers contenant les indices de refraction S&F
C		       pour chaque constituant XX = SR, LR, SU, LU, OM.
C    SOS_MAXNB_ANG_EXT : Taille maximale des tableaux de fonctions de phase externes
C
C Constantes spécifiques :
C    LENLINE  : Longueur maximale d'une ligne de données lue dans un fichier
C----------------------------------------------------------------------------
#include "SOS.h"
#define LENLINE 1000

C==============================================================================
C PROGRAMME: SOS_AEROSOLS
C ==========  
C  Modeles simulees :
C  ----------------
C      Ce programme calcule les proprietes radiatives d'une granulometrie 
C      d'aerosols mono-modale, multi-modale (de la WMO) ou bi-modale (S&F).
C
C      Les constituants de base des modeles de la WMO sont:
C      dust-like, oceanic, water-soluble et soot. 
C      Quatre modeles WMO types peuvent etre simules: continental, maritime, 
C      urbain et une composition definie par l'utilisateur.
C
C      Les constituants de base des modeles de Shettle & Fenn sont:
C      small rural, large rural, small urban, large urban et oceanic. 
C      Cinq modeles S&F types peuvent etre simules: tropospherique, urbain,  
C      maritime et cotier.
C
C      Si imodalite = 0 => prise en compte d'une seule granulometrie d'aerosols
C			   (Log Normale ou loi de Junge).
C
C      Si imodalite = 1 => cas de granulometries multi-modales, utilisation  
C                          des quatre modeles standards WMO.
C
C      Si imodalite = 2 => cas de granulometries bi-modales, utilisation  
C                          des cinq modeles standards S&F.
C
C      Si imodalite = 3 => cas de granulometries LND bi-modales.
C
C      Si imodalite = 4 => cas de l'utilisation d'une matrice de phase externe 
C                         (fonctions de phase en intensite et polarisation).
C
C  Calculs de Mie et application :
C  -----------------------------
C      Les calculs de MIE donnent les proprietes radiatives d'un composant. 
C      Ils sont realises par la procedure "SOS_MIE".
C      Les fichiers de MIE sont stockes sous un repertoire specifique. 
C
C      Le pas sur le parametre de taille alpha est identique est :
C         si        alpha  <= 0.1  ===> pas = 0.0001
C         si 0.1  < alpha  <= 1.0  ===> pas = 0.001
C         si 1.0  < alpha  <= 10   ===> pas = 0.01
C         si 10   < alpha  <= 30   ===> pas = 0.05
C         si 30   < alpha  <= 100  ===> pas = 0.1
C         si 100  < alpha          ===> pas = 1.0
C
C      A partir des calculs de Mie, on determine la fonction de phase 
C      P(mu) et les fonctions de phase polarisee Q(mu) et U(mu), pour une 
C      longueur d'onde donnee et une distribution en taille des particules 
C      donnee pour chaque modele.
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
C      Si l'epaisseur optique des aerosols est nulle, il y a creation 
C      d'un fichier Aerosols avec des valeurs mises a zero.
C
C
C Donnees en entree du programme
C ------------------------------
C     -->Les parametres/arguments du programme
C
C        Les arguments sont passes sous la forme d'une liste de couples (-Keyword Value)
C
C        Definition des mots clef, statut obligatoire / option, format de la valeur associee :
C         -AER.AngFile  ==> Obligatoire (chemin complet au fichier ANGLES_LUM) 
C             FICANGLES (CHARACTER*SOS_LENFIC2) : nom complet du fichier des angles utilises 
C                                                 pour les calculs de luminance et BRDF
C                 (repertoire + nom fichier + extension)
C
C         -AER.AOT ==> Obligatoire 
C            PROG_ARG (CHARACTER*10)  : Epaisseur optique des aerosols (non tronquee).
C                                       (format de décodage : libre --> Float)
C      		       
C         -AER.Wa ==> Obligatoire  (pour cas general AER.AOT > 0)
C            PROG_ARG (CHARACTER*10)  :  Longueur d'onde de simulation en microns.	
C                                       (format de décodage : libre --> Float)
C     		       
C         -AER.Tronca ==> Obligatoire (pour cas general AER.AOT > 0)
C            PROG_ARG (CHARACTER*10)  :  Option de troncature  (1 pour la realiser).
C                                       (format de décodage : libre --> Integer) 
C
C         -AER.Log ==> Obligatoire  (pour cas general AER.AOT > 0)
C            FICLOGGRANU (CHARACTER*SOS_LENFIC2) : Fichier trace du programme Aerosols 
C						  (0 ou chemin complet au fichier Trace)
C                 (repertoire + nom fichier + extension)
C
C         -AER.MieLog ==> Obligatoire  (pour cas general AER.AOT > 0)
C            FICLOGMIE (CHARACTER*SOS_LENFIC2) : Fichier trace MIE  
C						 (0 ou chemin complet au fichier Trace)
C                 (repertoire + nom fichier + extension)	 
C
C         -AER.DirMie ==> Obligatoire   (pour cas general AER.AOT > 0)
C             DIRS   (CHARACTER*SOS_LENDIR)  : repertoire de stockage des fichiers MIE
C	  
C         -AER.ResFile ==> Obligatoire   (pour cas general AER.AOT > 0)
C             FICGRANU (CHARACTER*SOS_LENFIC2) : nom complet du fichier resultat GRANU
C                 (repertoire + nom fichier + extension)
C      		       
C         -AER.Model ==> Obligatoire  (pour cas general AER.AOT > 0)
C            PROG_ARG (CHARACTER*10)  :  Granulometrie mono-modal ou multi-modale.
C                                       (format de décodage : libre --> Integer) 
C      		        0 : Modeles mono-mode
C			1 : Modeles WMO
C			2 : Modeles Shettle & Fenn
C		        3 : Modeles bimodaux de LND
C                       4 : Utilisation de fonctions de phase externes
C      		       
C         -AER.MMD.Mie.AlphaMax ==> Associe a -AER.Model 0
C            PROG_ARG (CHARACTER*10)  : Valeur maximale du parametre de taille 
C                                       des calculs de Mie.	
C                                       (format de décodage : libre --> Float) 
C	  
C         -AER.MMD.Mie.Filename ==> Optionnel : Associe a -AER.Model 0 
C            FICS  (CHARACTER*SOS_LENFIC1) : Nom du fichier MIE utilisateur (sans repertoire)
C						--> FICS="DEFAULT" si pas d'argument -AER.MMD.Mie.Filename
C                 	(nom fichier + extension)						
C	  
C         -AER.MMD.MRwa ==> Associe a -AER.Model 0
C            PROG_ARG (CHARACTER*10)  : Partie reelle de l'indice de refraction pour modele mono-modal
C                                       a la longueur d'onde de simulation des luminances (-AER.Wa)	
C                                       (format de décodage : libre --> Float)
C
C         -AER.MMD.MIwa ==> Associe a -AER.Model 0
C            PROG_ARG (CHARACTER*10)  : Partie imaginaire de l'indice de refraction pour modele mono-modal
C                                       a la longueur d'onde de simulation des luminances (-AER.Wa)	
C                                      (format de décodage : libre --> Float)
C	  
C         -AER.MMD.SDtype ==> Associe a -AER.Model 0
C            PROG_ARG (CHARACTER*10)  : Indice de type de granulometrie pour modele mono-modal
C                                       (format de décodage : libre --> Integer) 
C                          1 : LND (parametres : rayon modal en microns et variance);  
C                          2 : Loi de Junge (parametres : rayon minimal  
C                              en microns et puissance appliquee au rayon).
C
C         -AER.MMD.SDparam1 ==> Associe a -AER.Model 0
C            PROG_ARG (CHARACTER*10)  : 1er parametre de granulometrie pour modele mono-modal
C                                      (format de décodage : libre --> Float)	
C	  
C         -AER.MMD.SDparam2 ==> Associe a -AER.Model 0
C            PROG_ARG (CHARACTER*10)  : 2nd parametre de granulometrie pour modele mono-modal
C                                      (format de décodage : libre --> Float)	
C
C         -AER.WMO.Model ==> Associe a -AER.Model 1
C            PROG_ARG (CHARACTER*10)  : Indice de choix du modele WMO
C                                       (format de décodage : libre --> Integer) 
C      		        1 : Modele WMO Continental
C			2 : Modele WMO Maritime
C			3 : Modele WMO Urban
C			4 : Modele WMO "user definition"
C
C         -AER.WMO.DL ==> Associe a -AER.Model 1 et -AER.WMO.Model 4
C            PROG_ARG (CHARACTER*10)  : Pourcentage volumique C des particules DL (WMO).
C                                      (format de décodage : libre --> Float)
C
C         -AER.WMO.WS ==> Associe a -AER.Model 1 et -AER.WMO.Model 4
C            PROG_ARG (CHARACTER*10)  : Pourcentage volumique C des particules WS (WMO).
C                                      (format de décodage : libre --> Float)
C
C         -AER.WMO.OC ==> Associe a -AER.Model 1 et -AER.WMO.Model 4
C            PROG_ARG (CHARACTER*10)  : Pourcentage volumique C des particules OC (WMO).
C                                      (format de décodage : libre --> Float)
C
C         -AER.WMO.SO ==> Associe a -AER.Model 1 et -AER.WMO.Model 4
C            PROG_ARG (CHARACTER*10)  : Pourcentage volumique C des particules SO (WMO).
C                                      (format de décodage : libre --> Float)

C
C         -AER.SF.Model ==> Associe a -AER.Model 2 
C      	     MODELE   (CHARACTER*10)  : Indice de choix du modele Shettle & Fenn
C                                       (format de décodage : libre --> Integer) 
C      		        1 : Modele SF Tropospheric
C			2 : Modele SF Urban
C			3 : Modele SF Maritime
C			4 : Modele SF Coastal
C
C         -AER.SF.RH ==> Associe a -AER.Model 2
C      	     DA20  (CHARACTER*10)  : Humidite relative pour modele Shettle & Fenn
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.VCdef ==> Associe a -AER.Model 3
C      	     IAMODE_PARAM_BILND (CHARACTER*10) : Indice de type de description 
C                                                du melange bimodal LND. 
C                                       (format de décodage : libre --> Integer) 
C      		        1 : Utilisation des concentrations volumiques de l'utilisateur
C			2 : Utilisation du rapport AOT_coarse / AOT_totale
C
C         -AER.BMD.CoarseVC ==> Associe a -AER.Model 3 et -AER.BMD.VCdef 1
C      	     DACVC (CHARACTER*10)  : Concentration volumique de l'utilisateur pour le mode "LND coarse". 
C                                      (format de décodage : libre --> Float)
C	 
C         -AER.BMD.FineVC ==> Associe a -AER.Model 3 et -AER.BMD.VCdef 1
C      	     DACVF (CHARACTER*10)  : Concentration volumique de l'utilisateur pour le mode "LND fine". 
C                                      (format de décodage : libre --> Float)
C	 
C         -AER.BMD.RAOT ==> Associe a -AER.Model 3 et -AER.BMD.VCdef 2
C      	     DATAUct (CHARACTER*10) : rapport AOT_coarse / AOT_tot pour la longueur d'onde WAREF.
C                                      (format de décodage : libre --> Float)
C	 
C         -AER.Waref ==> Associe a -AER.Model 3 et -AER.BMD.VCdef 2
C      	     DAWAREF (CHARACTER*10) : longueur d'onde de reference WAREF.
C                                      (format de décodage : libre --> Float)
C	 
C         -AER.AOTref ==> Associe a -AER.Model 3 et -AER.BMD.VCdef 2
C      	     DATAWAREF (CHARACTER*10) : AOT totale pour la longueur d'onde WAREF.
C                                      (format de décodage : libre --> Float)
C	 
C         -AER.BMD.CM.MRwa ==> Associe a -AER.Model 3 
C      	     IARN1 (CHARACTER*10) : Partie reelle de l'indice de refraction
C                                   des particules du mode "LND coarse" pour WA de simulation
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.CM.MIwa ==> Associe a -AER.Model 3 
C      	     IAIN1 (CHARACTER*10) : Partie imaginaire de l'indice de refraction
C                                   des particules du mode "LND coarse" pour WA de simulation
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.CM.MRwaref ==> Associe a -AER.Model 3 et -AER.BMD.VCdef 2 
C      	     IARN1 (CHARACTER*10) : Partie reelle de l'indice de refraction
C                                   des particules du mode "LND coarse" pour WAREF
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.CM.MIwaref ==> Associe a -AER.Model 3 et -AER.BMD.VCdef 2 
C      	     IAIN1 (CHARACTER*10) : Partie imaginaire de l'indice de refraction
C                                   des particules du mode "LND coarse" pour WAREF
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.CM.SDradius ==> Associe a -AER.Model 3 
C      	     DARM1 (CHARACTER*10) : Rayon modal de la "LND coarse"
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.CM.SDvar ==> Associe a -AER.Model 3 
C      	     DAVAR1 (CHARACTER*10) : Variance de la "LND coarse"
C                                      (format de décodage : libre --> Float)
C	 
C         -AER.BMD.FM.MRwa ==> Associe a -AER.Model 3 
C      	     IARN2 (CHARACTER*10) : Partie reelle de l'indice de refraction
C                                   des particules du mode "LND fine" pour WA de simulation
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.FM.MIwa ==> Associe a -AER.Model 3 
C      	     IAIN2 (CHARACTER*10) : Partie imaginaire de l'indice de refraction
C                                   des particules du mode "LND fine" pour WA de simulation
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.FM.MRwaref ==> Associe a -AER.Model 3 et -AER.BMD.VCdef 2
C      	     IARN2 (CHARACTER*10) : Partie reelle de l'indice de refraction
C                                   des particules du mode "LND fine" pour WAREF
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.FM.MIwaref ==> Associe a -AER.Model 3 et -AER.BMD.VCdef 2
C      	     IAIN2 (CHARACTER*10) : Partie imaginaire de l'indice de refraction
C                                   des particules du mode "LND fine" pour WAREF
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.FM.SDradius ==> Associe a -AER.Model 3 
C      	     DARM2 (CHARACTER*10) : Rayon modal de la "LND fine"
C                                      (format de décodage : libre --> Float)
C
C         -AER.BMD.FM.SDvar ==> Associe a -AER.Model 3 
C      	     DAVAR2 (CHARACTER*10) : Variance de la "LND fine"
C                                      (format de décodage : libre --> Float)
C 
C         -AER.ExtData ==> Associe a -AER.Model 4 
C            FICEXTDATA  (CHARACTER*SOS_LENFIC2) : 
C                     Nom du fichier de donnees externes de fonctions de phase
C                     et parametres radiatifs (coefficients extinction et diffusion) 
C                    (repertoire + nom fichier + extension)
C                          
C
C    -->Fichiers
C       Le fichier MIE
C       Le fichier Trace Granu
C       Le fichier Trace Mie
C       Le fichier des angles et parametres internes pour les limites des developpements.
C       Le fichier des fonctions de phase externes (entree optionnelle)
C
C Resultats fournis par le programme
C ----------------------------------
C   Fichier resultat GRANU contenant :
C	- coefficient d'extinction et de diffusion, Beta(0)/3
C	- coefficient de troncature
C	- albedo de simple diffusion tronquee
C	- les parametres de decomposition de la fonction de phase en 
C	  polynomes de Legendre.
C
C Variable d'environnement
C ------------------------
C   SOS_RACINE_FIC : Variable contenant le nom du repertoire de stockage 
C		     des fichiers de la WMO et de S&F.
C                    Elle doit etre definie dans le programme de lancement 
C                    par la commande "export SOS_RACINE_FIC=..."
C
C
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier de MIE
C       Fichier binaire non formatte, cree par SOS_MIE
C       1 enregistrement d'entete avec
C                Indice de refraction (partie reelle)
C                Indice de refraction (partie imaginaire)
C                Parametre de taille limite
C                Nombre d'angles de Gauss positifs
C
C      1 enregistrement pour chaque valeur du parametre de taille avec:
C        - Parametre de taille (reel),
C        - Coefficient d'efficacite en extinction (reel),
C        - Coefficient d'efficacite en diffusion (reel),
C        - Facteur d'anisotropie (double precision),
C        - Tableau de la fonction de phase Mie(mu) (reel),
C        - Tableau de la fonction de phase polarisee Qmie(mu) (reel),
C        - Tableau de la fonction de phase polarisee Umie(mu) (reel).
C
C    --> Les tableaux ont une taille utile limitee a MIE_NBMU : 
C        ils sont crees et lus avec cette limitation.
C
C
C    -->Contenu du fichier des angles:
C       Fichier ascii formatte, cree par SOS_ANGLES.
C       Il contient l'information sur les angles a utiliser 
C       et sur les ordres des developpements adaptes : 
C		NB_TOTAL_ANGLES :  nb d'angles total a utiliser
C               NB_GAUSS_ANGLES :  nb d'angles de Gauss
C               ANGLES_USERFILE :  Nom du fichier des angles utilisateur
C               INTERNAL_OS_NB :   Ordre des developpements des fonctions de phase 
C				   en Polynomes de Legendre a utiliser.
C               INDEX COS_ANGLE WEIGHT : lignes donnant le numero d'angle, son cosinus et son poids
C                                        au format I4,1X,2D21.14
C 
C
C    -->Contenu du fichier des fonctions de phase externes:  
C      Fichier ascii formaté, cree par l'utilisateur.
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
C                           F11 		: fonction de phase en intensiteé
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
C Common utilise:
C --------------
C  Aucun
C
C
C Cas d'erreur :
C ------------        
C     -Decodage des parametres incorrect
C     -Erreur dans un sous-programme
C     -Erreur a l'ouverture du fichier Trace
C     -Parametres non valables
C
C     Affichage d'un message a l'ecran, arret du programme et 
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
      PROGRAM SOS_AEROSOLS
      
      IMPLICIT NONE

C Constantes AEROSOLS
C -------------------

      INTEGER SOS_MODALITE	!Valeur max du parametre modalite 
      PARAMETER (SOS_MODALITE=4)

      INTEGER SOS_NB_MODELE	!Nombre de modeles WMO ou Shettle & Fenn possibles.
      PARAMETER (SOS_NB_MODELE=4)

      INTEGER SOS_NB_WMO	!Nombre de types d'aerosols de la WMO utilises:
      PARAMETER (SOS_NB_WMO=4)  !DL,WS,OC,SO.
      
      INTEGER SOS_NB_SF		!Nombre de types d'aerosols du modele Shettle & Fenn 
      PARAMETER (SOS_NB_SF=5)   !utilises:SR,LR,SU,LU,OM.
          
      INTEGER SOS_NB_AERO	!Initialisation avec le plus grand nombre de types 
      PARAMETER (SOS_NB_AERO=5)	!d'aerosols (cas Shettle & Fenn).	
     

C* Definition des variables                         
C*-----------------------------------------------------------------


      DOUBLE PRECISION XMU(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)  
                                                  !Cosinus des angles de Gauss (et utilisateur).
      DOUBLE PRECISION XHR(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)  
                                                  !Poids des angles de Gauss.				

      DOUBLE PRECISION TA	! Epaisseur optique des aerosols (non tronquee).

      DOUBLE PRECISION WA	! Longueur d'onde en microns selectionnee.	
            						  						  
      DOUBLE PRECISION C(SOS_NB_WMO)	 !Proportion volumique des constituants WMO
      DOUBLE PRECISION C_WMO_DL	    	 ! --> Valeur utilisateur pour particule DL
      DOUBLE PRECISION C_WMO_WS	    	 ! --> Valeur utilisateur pour particule WS
      DOUBLE PRECISION C_WMO_OC	    	 ! --> Valeur utilisateur pour particule OC
      DOUBLE PRECISION C_WMO_SO	    	 ! --> Valeur utilisateur pour particule SO
      DOUBLE PRECISION V(SOS_NB_WMO)	 !Concentration volumique des constituants WMO
      DOUBLE PRECISION N(SOS_NB_WMO)	 ! C/V
      DOUBLE PRECISION NTOT		 ! Somme des coefficients de ponderation 
                                         ! des composants aerosols.
      DOUBLE PRECISION ALPHAF_PARTICULE(4) ! Valeur Alpha Max pour chaque
      					   ! composant de la WMO (DL, WS, OC, SO)
					   
      DOUBLE PRECISION Ni(SOS_NB_SF)	 !densite de particules (Shettle & Fenn)					   
      DOUBLE PRECISION RH                ! % d'humidite relative des modeles de Shettle & Fenn

      DOUBLE PRECISION RN	! partie reelle de l'indice de refraction
      DOUBLE PRECISION IN	! partie imaginaire de l'indice de refraction
      DOUBLE PRECISION ALPHAO   ! valeur de alpha minimale des calculs de Mie		
      DOUBLE PRECISION ALPHAF   ! valeur de alpha maximale des calculs de Mie			

      DOUBLE PRECISION SOMME_NR	! nombre de particules : integration sur les rayons de NR * PR.  
      
      DOUBLE PRECISION MR(SOS_NB_AERO)   ! partie reelle de l'indice de refraction
      					 ! pour chaque composant d'aerosols
      DOUBLE PRECISION MI(SOS_NB_AERO)   ! partie imaginaire de l'indice de refraction.
      					 ! pour chaque composant d'aerosols.

      
      DOUBLE PRECISION VARGRANU1    ! Parametre 1 de la granulometrie.
      DOUBLE PRECISION VARGRANU2    ! Parametre 2 de la granulometrie.
      DOUBLE PRECISION VARGRANU1I(SOS_NB_AERO)   ! Parametre 1 de la granulometrie
      						 ! pour chaque modele de constituant.
      DOUBLE PRECISION VARGRANU2I(SOS_NB_AERO)   ! Parametre 2 de la granulometrie
      						 ! pour chaque modele de constituant.      
	      				

	!Parametres pour un modele bimodal de LND
      INTEGER*2 MODE_PARAM_BILND       	!Indice du choix de description du melange de 2 LND.
      DOUBLE PRECISION CVI(2)		!Concentration volumique de la "LND coarse" (1) et 
      				      	!de la "LND fine" (2)
      DOUBLE PRECISION USER_CV_COARSE	!Valeur utilisateur de la concentration volumique de la "LND coarse"
      DOUBLE PRECISION USER_CV_FINE	!Valeur utilisateur de la concentration volumique de la "LND fine"	
      DOUBLE PRECISION RTAUct_WAREF	!Rapport des epaisseurs optique AOT_coarse / AOT_total 
      DOUBLE PRECISION WAREF		!Longueur d'onde du rapport AOT_coarse / AOT_total
      DOUBLE PRECISION TA_WAREF		!Epaisseur optique aerosol totale 
      					!pour la longueur d'onde de reference.
					
      DOUBLE PRECISION RN_LND(2)	!Partie reelle de l'indice de refraction de la 
      					!"LND coarse" (1) et de la "LND fine" (2) pour la longeur d'onde WA.					
      DOUBLE PRECISION BMD_CM_MRWA	! Valeur lue pour RN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MRWA	! Valeur lue pour RN_LND(2) : fine mode    
       
      DOUBLE PRECISION IN_LND(2) 	!Partie imaginaire de l'indice de refraction de la 
      					!"LND coarse" (1) et de la "LND fine" (2) pour la longeur d'onde WA.					
      DOUBLE PRECISION BMD_CM_MIWA	! Valeur lue pour IN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MIWA	! Valeur lue pour IN_LND(2) : fine mode 
      
      DOUBLE PRECISION RN_LND_WAREF(2)	!Partie reelle de l'indice de refraction de la 
      					!"LND coarse" (1) et de la "LND fine" (2) pour la longueur d'onde WAREF.
      DOUBLE PRECISION BMD_CM_MRWAREF	! Valeur lue pour RN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MIWAREF	! Valeur lue pour RN_LND_WAREF(2) : fine mode 
      
      DOUBLE PRECISION IN_LND_WAREF(2) 	!Partie imaginaire de l'indice de refraction de la 
      					!"LND coarse" (1) et de la "LND fine" (2) pour la longueur d'onde WAREF.
      DOUBLE PRECISION BMD_FM_MRWAREF	! Valeur lue pour IN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MIWAREF	! Valeur lue pour IN_LND_WAREF(2) : fine mode      
      
      DOUBLE PRECISION RMODAL_LND(2)  	!Rayon modal des "LND coarse" (1) et "fine" (2).
      DOUBLE PRECISION BMD_CM_RMODAL	! Valeur lue pour RMODAL_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_RMODAL	! Valeur lue pour RMODAL_LND(2) : fine mode  
      
      DOUBLE PRECISION VAR_LND(2)	!Variance des "LND coarse" (1) et "fine" (2).
      DOUBLE PRECISION BMD_CM_VAR	! Valeur lue pour VAR_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_VAR	! Valeur lue pour VAR_LND(2) : fine mode 
          
      
      DOUBLE PRECISION KMAT1c		! Section efficace d'extinction (micron^2) pour le mode 
      					! "coarse" a la longueur d'onde WAREF. Non tronque.
      DOUBLE PRECISION KMAT1f		! Section efficace d'extinction (micron^2) pour le mode 
      					! "fine" a la longueur d'onde WAREF. Non tronque.
     
      DOUBLE PRECISION KMAT1I	! Section efficace d'extinction (micron^2)
      				! pour un constituant I. Non tronque.
      DOUBLE PRECISION KMAT2I 	! Section efficace de diffusion (micron^2)
      				! pour un constituant I.  Non tronque.       
      DOUBLE PRECISION KMAT1	! Section efficace d'extinction (micron^2). Non tronque.
      DOUBLE PRECISION KMAT2 	! Section efficace de diffusion (micron^2). Non tronque.      
      DOUBLE PRECISION A	! Coefficient de troncature.
      DOUBLE PRECISION PIZ	! Albedo de simple diffusion.
      DOUBLE PRECISION PIZTR	! Albedo de simple diffusion apres troncature.
      
      
      DOUBLE PRECISION P11I(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Fonction de phase calculee a partir des fonctions de 
		       !la theorie de MIE pour le constituant I.
      DOUBLE PRECISION P12I(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P12 de la matrice de phase pour le constituant I.
      DOUBLE PRECISION P33I(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P33 de la matrice de phase pour le constituant I.
      
      DOUBLE PRECISION P11(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Fonction de phase du modele simule.
      DOUBLE PRECISION P12(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P12 de la matrice de phase.
      DOUBLE PRECISION P22(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P22 de la matrice de phase.
      DOUBLE PRECISION P33(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P33 de la matrice de phase.
		       
      DOUBLE PRECISION TTT(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Copie de la fonction de phase (P11) avant troncature.
     
      DOUBLE PRECISION Z1	!Valeur de BETA11(0) avant normalisation.
      DOUBLE PRECISION TETA     !Angle de diffusion.
      
                  
      DOUBLE PRECISION ALP(0:SOS_OS_NB_MAX)	! Parametres alpha(k).
      DOUBLE PRECISION BETA11(0:SOS_OS_NB_MAX)	! Parametres beta11(k).
      DOUBLE PRECISION BETA22(0:SOS_OS_NB_MAX)	! Parametres beta22(k).
      DOUBLE PRECISION GAMMA12(0:SOS_OS_NB_MAX)	! Parametres gamma12(k).
      DOUBLE PRECISION ZETA(0:SOS_OS_NB_MAX)	! Parametres zeta(k).
      DOUBLE PRECISION DELTA33(0:SOS_OS_NB_MAX)	! Parametres delta33(k).
     
	    
      DOUBLE PRECISION X, RMAX   ! Parametres pour le calcul de AlphaMax 
      				 ! selon la taille des particules.
               
      DOUBLE PRECISION VAL	! Variable pour la lecture de valeurs inutilisees
      
      DOUBLE PRECISION ANG_EXT(SOS_MAXNB_ANG_EXT)	 ! Table des angles des fonctions de phase externes
      DOUBLE PRECISION MU_EXT(SOS_MAXNB_ANG_EXT)	 ! Table des cosinus des angles des fonctions de phase externes   
      DOUBLE PRECISION F11_EXT(SOS_MAXNB_ANG_EXT)	 ! Fonction de phase P11 de donnees externes
      DOUBLE PRECISION F12_EXT(SOS_MAXNB_ANG_EXT)	 ! Fonction de phase P11 de donnees externes
      DOUBLE PRECISION F22_EXT(SOS_MAXNB_ANG_EXT)	 ! Fonction de phase P11 de donnees externes
      DOUBLE PRECISION F33_EXT(SOS_MAXNB_ANG_EXT)	 ! Fonction de phase P11 de donnees externes
      DOUBLE PRECISION F12sF11_EXT(SOS_MAXNB_ANG_EXT) ! Fonction de phase P12/P11 de donnees externes
      DOUBLE PRECISION F22sF11_EXT(SOS_MAXNB_ANG_EXT) ! Fonction de phase P22/P11 de donnees externes
      DOUBLE PRECISION F33sF11_EXT(SOS_MAXNB_ANG_EXT) ! Fonction de phase P33/P11 de donnees externes   
      
      INTEGER*2 IMOD	   ! Granulometrie mono ou multi-modale.
      INTEGER*2 IMODELE_WMO ! Type du modele WMO.
      INTEGER*2 IMODELE_SF  ! Type du modele Shettle & Fenn.

      INTEGER*4 NBMU_GAUSS  ! Nombre d'angles de Gauss utilises 
      INTEGER*4 MIE_NBMU   ! Nombre d'angles effectifs des simulations
      INTEGER*4 OS_NB 	   ! Ordre des developpements des fonctions de phase 
			   ! en Polynomes de Legendre.      
      INTEGER*4 IGRANU	   ! Numero du modele de granulometrie
      INTEGER*4 ITRONC	   ! Option de troncature  (1 pour la realiser).  
      INTEGER*4 LDIRS	   ! longueur de la chaine de caractere DIRS
      INTEGER*4 LFICS	   ! longueur de la chaine de caractere FICS  
      INTEGER*4 IFIN       ! nb caracteres significatifs d'une chaine                
      INTEGER*4 I	   ! Indice de tableau 
      INTEGER*4 J	   ! Indice de tableau 
      INTEGER*4 K          ! Indice de tableau        
      INTEGER*4 IER        ! code d'erreur =0 si pas d'erreur, =-1 sinon
      INTEGER*4 IARGC      ! Type de la fonction IARGC
      INTEGER*4 NPARAM 	   ! Numero du parametre en entrée (type INTEGER*4 obligatoire
                           ! pour le f77)
      INTEGER*4 NB_FICEXT_ANGLES ! Nombre d'angles du fichier de fonctions de phase externe
      INTEGER*4 NVAL 	   ! Integer pour la lecture de fichier ANGLES
      
      INTEGER*2 NBARGS		! Nombre d'arguments en entrée
      INTEGER*2 DEF_TA	        ! Indicateur d'affectation d'une valeur a TA 
      INTEGER*2 DEF_WA	 	! Indicateur d'affectation d'une valeur a WA 
      INTEGER*2 DEF_ITRONC	! Indicateur d'affectation d'une valeur a ITRONC
      INTEGER*2 DEF_FICLOGGRANU ! Indicateur d'affectation d'une valeur a FICLOGGRANU
      INTEGER*2 DEF_FICLOGMIE 	! Indicateur d'affectation d'une valeur a FICLOGMIE
      INTEGER*2 DEF_FICGRANU	! Indicateur d'affectation d'une valeur a FICGRANU
      INTEGER*2 DEF_IMOD	! Indicateur d'affectation d'une valeur a IMOD
      INTEGER*2 DEF_FICANGLES   ! Indicateur d'affectation d'une valeur a FICANGLES_USER
      INTEGER*2 DEF_ALPHAF	! Indicateur d'affectation d'une valeur a ALPHAF
      INTEGER*2 DEF_FICS	! Indicateur d'affectation d'une valeur a FICS
      INTEGER*2 DEF_DIRS	! Indicateur d'affectation d'une valeur a DIRS
      INTEGER*2 DEF_RN		! Indicateur d'affectation d'une valeur a RN
      INTEGER*2 DEF_IN		! Indicateur d'affectation d'une valeur a IN
      INTEGER*2 DEF_IGRANU	! Indicateur d'affectation d'une valeur a IGRANU
      INTEGER*2 DEF_VARGRANU1	! Indicateur d'affectation d'une valeur a VARGRANU1
      INTEGER*2 DEF_VARGRANU2	! Indicateur d'affectation d'une valeur a VARGRANU2
      INTEGER*2 DEF_IMODELE_WMO ! Indicateur d'affectation d'une valeur a IMODELE_WMO
      INTEGER*2 DEF_WMO_DL 	! Indicateur d'affectation d'une valeur a C_WMO_DL pour cas WMO
      INTEGER*2 DEF_WMO_WS 	! Indicateur d'affectation d'une valeur a C_WMO_WS pour cas WMO
      INTEGER*2 DEF_WMO_OC 	! Indicateur d'affectation d'une valeur a C_WMO_OC pour cas WMO
      INTEGER*2 DEF_WMO_SO 	! Indicateur d'affectation d'une valeur a C_WMO_SO pour cas WMO
      INTEGER*2 DEF_IMODELE_SF  ! Indicateur d'affectation d'une valeur a IMODELE_SF
      INTEGER*2 DEF_RH  	! Indicateur d'affectation d'une valeur a RH pour cas S&F
      INTEGER*2 DEF_MODE_PARAM_BILND ! Indicateur d'affectation d'une valeur a DEF_MODE_PARAM_BILND
      INTEGER*2 DEF_USER_CV_COARSE   ! Indicateur d'affectation d'une valeur a USER_CV_COARSE
      INTEGER*2 DEF_USER_CV_FINE     ! Indicateur d'affectation d'une valeur a USER_CV_FINE
      INTEGER*2 DEF_RTAUct_WAREF     ! Indicateur d'affectation d'une valeur a RTAUct_WAREF
      INTEGER*2 DEF_WAREF     	     ! Indicateur d'affectation d'une valeur a WAREF
      INTEGER*2 DEF_TA_WAREF         ! Indicateur d'affectation d'une valeur a TA_WAREF
      INTEGER*2 DEF_BMD_CM_MRWA		! Indicateur d'affectation d'une valeur a BMD_CM_MRWA 
      INTEGER*2 DEF_BMD_CM_MIWA		! Indicateur d'affectation d'une valeur a BMD_CM_MIWA
      INTEGER*2 DEF_BMD_CM_MRWAREF	! Indicateur d'affectation d'une valeur a BMD_CM_MRWAREF
      INTEGER*2 DEF_BMD_CM_MIWAREF	! Indicateur d'affectation d'une valeur a BMD_CM_MIWAREF
      INTEGER*2 DEF_BMD_CM_RMODAL	! Indicateur d'affectation d'une valeur a BMD_CM_RMODAL
      INTEGER*2 DEF_BMD_CM_VAR		! Indicateur d'affectation d'une valeur a BMD_CM_VAR
      INTEGER*2 DEF_BMD_FM_MRWA		! Indicateur d'affectation d'une valeur a BMD_FM_MRWA
      INTEGER*2 DEF_BMD_FM_MIWA		! Indicateur d'affectation d'une valeur a BMD_FM_MIWA
      INTEGER*2 DEF_BMD_FM_MRWAREF	! Indicateur d'affectation d'une valeur a BMD_FM_MRWAREF
      INTEGER*2 DEF_BMD_FM_MIWAREF	! Indicateur d'affectation d'une valeur a BMD_FM_MIWAREF
      INTEGER*2 DEF_BMD_FM_RMODAL	! Indicateur d'affectation d'une valeur a BMD_FM_RMODAL
      INTEGER*2 DEF_BMD_FM_VAR		! Indicateur d'affectation d'une valeur a BMD_FM_VAR
      INTEGER*2 DEF_FICEXTDATA		! Indicateur d'affectation d'une valeur a FICEXTDATA
       
      CHARACTER*10 PROG_ARG		        ! Arguments du programme
     	   
      CHARACTER*SOS_LENDIR DIRS	                ! chemin d'acces au fichier MIE
      CHARACTER*SOS_LENFIC1 FICS		! nom du fichier MIE 
      CHARACTER*SOS_LENFIC2 FICMIE	        ! chemin et nom du fichier MI
      CHARACTER*SOS_LENFIC2 FICLOGMIE           ! nom du fichier TRACE de MIE  
      CHARACTER*SOS_LENFIC2 FICLOGGRANU        	! nom du fichier TRACE de GRANU
      CHARACTER*SOS_LENFIC2 FICGRANU		! nom du fichier resultat GRANU  
      CHARACTER*SOS_LENFIC2 FICANGLES           ! chemin et nom du fichier des angles effectifs    
      CHARACTER*SOS_LENFIC2 FICANGLES_USER      ! chemin et nom du fichier des angles "utilisateurs"
      CHARACTER*SOS_LENFIC2 FICEXTDATA          ! chemin et nom du fichier de donnees externes de fonctions
      						! de phase et parametres radiatifs (coefficients extinction et
						! diffusion)
      
      CHARACTER*SOS_LENDIR DIRE	                ! chemin d'acces au fichier des parametres des modeles WMO et S&F
      CHARACTER*SOS_LENFIC2 FICDATAWMO          ! chemin et nom du fichier Data_WMO
      
      CHARACTER*15 PARTICULE(SOS_NB_AERO)       ! nom des differentes particules de la WMO
      					        ! ou des modeles S&F	

      CHARACTER*SOS_LENKEYWORD KEYWORD	       	! Mot clef du passage d'arguments.

      CHARACTER*LENLINE LIGNE_TEXTE             ! Ligne de caractères pour décodage de paramètres      
      CHARACTER*1000 CVAL	                ! Chaine de caracteres pour la lecture de valeurs inutilisees


      LOGICAL TRACE        ! = vrai,si ecriture dans le fichier trace
      LOGICAL EX	   ! Verifie l'existence d'un fichier

 
     
      
C* Lecture des parametres/arguments d'entree communs et decodage 
C---------------------------------------------------------------
C --- Affectation de la valeur minimale du parametre de taille 
C     (passage de la constante en double precision)
      ALPHAO = SOS_MIE_ALPHAMIN      
      
      IF ((ALPHAO.LT.1D-4).OR.(ALPHAO.GE.10D-00)) GOTO 1010

C --- Initialisation des indicateurs de definition des parametres (0 pour indefinis)
      DEF_TA=0
      DEF_WA=0
      DEF_ITRONC=0
      DEF_FICLOGGRANU=0
      DEF_FICLOGMIE=0
      DEF_FICGRANU=0
      DEF_IMOD=0
      DEF_FICANGLES=0
      DEF_ALPHAF=0
      DEF_FICS=0
      DEF_DIRS=0
      DEF_RN=0
      DEF_IN=0
      DEF_IGRANU=0
      DEF_VARGRANU1=0
      DEF_VARGRANU2=0
      DEF_IMODELE_WMO=0
      DEF_WMO_DL=0	    
      DEF_WMO_WS=0
      DEF_WMO_OC=0
      DEF_WMO_SO=0
      DEF_IMODELE_SF=0
      DEF_RH=0 
      DEF_MODE_PARAM_BILND=0
      DEF_USER_CV_COARSE=0
      DEF_USER_CV_FINE=0
      DEF_RTAUct_WAREF=0
      DEF_WAREF=0
      DEF_TA_WAREF=0	   
      DEF_BMD_CM_MRWA=0
      DEF_BMD_CM_MIWA=0
      DEF_BMD_CM_MRWAREF=0
      DEF_BMD_CM_MIWAREF=0
      DEF_BMD_CM_RMODAL=0
      DEF_BMD_CM_VAR=0
      DEF_BMD_FM_MRWA=0
      DEF_BMD_FM_MIWA=0
      DEF_BMD_FM_MRWAREF=0
      DEF_BMD_FM_MIWAREF=0
      DEF_BMD_FM_RMODAL=0
      DEF_BMD_FM_VAR=0
      DEF_FICEXTDATA=0

C --- Nombre d'arguments en entree
      NBARGS = IARGC()

C --- Boucle de lecture des couples d'arguments       
      DO NPARAM=1,NBARGS, 2
         CALL GETARG(NPARAM,KEYWORD)
	  	
	 IF (KEYWORD.EQ."-AER.AOT") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Epaisseur optique des aerosols
	     READ(PROG_ARG,*,ERR=990) TA
	     DEF_TA=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.Wa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)     	!Longueur d'onde
	     READ(PROG_ARG,*,ERR=990) WA
	     DEF_WA=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.Tronca") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)       	!Indice de troncature
	     READ(PROG_ARG,*,ERR=990) ITRONC
	     DEF_ITRONC=1
	 ENDIF	 
	 
	 IF (KEYWORD.EQ."-AER.Log") THEN
	     CALL GETARG (NPARAM+1,FICLOGGRANU) !Fichier Trace du programme Aerosols
	     DEF_FICLOGGRANU=1
	 ENDIF	 	  	 

	 IF (KEYWORD.EQ."-AER.MieLog") THEN
	     CALL GETARG (NPARAM+1,FICLOGMIE)	!Fichier Trace des calculs de Mie   
	     DEF_FICLOGMIE=1
	 ENDIF	 	  	 

         IF (KEYWORD.EQ."-AER.DirMie") THEN
	     CALL GETARG (NPARAM+1,DIRS)	!Repertoire de stockage des fichiers de Mie
	     DEF_DIRS=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICGRANU)	!Fichier resultat
	     DEF_FICGRANU=1
	 ENDIF	 	  	
	 	 	  
	 IF (KEYWORD.EQ."-AER.Model") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Type de modele d'aerosols
	     READ(PROG_ARG,*,ERR=990) IMOD
	     DEF_IMOD=1
	 ENDIF	 				   

	 IF (KEYWORD.EQ."-AER.AngFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES)	!Nom du fichier des angles (chemin complet)
	     DEF_FICANGLES=1
	 ENDIF	  

	 IF (KEYWORD.EQ."-AER.MMD.Mie.AlphaMax") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Valeur maximale du parametre de taille pour modele mono-mode
	     READ(PROG_ARG,*,ERR=990) ALPHAF
	     DEF_ALPHAF=1
	 ENDIF	 			  	 	  
			  
	 IF (KEYWORD.EQ."-AER.MMD.Mie.Filename") THEN
	     CALL GETARG (NPARAM+1,FICS)		!Nom du fichier de Mie (si donnee utilisateur)	
	     DEF_FICS=1
	 ENDIF	 				  

	 IF (KEYWORD.EQ."-AER.MMD.MRwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!partie reelle de l'indice de refraction pour modele mono-mode
	     READ(PROG_ARG,*,ERR=990) RN
	     DEF_RN=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.MMD.MIwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!partie imaginaire de l'indice de refraction pour modele mono-mode
	     READ(PROG_ARG,*,ERR=990) IN 
	     DEF_IN=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.SDtype") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!indice de type de granulometrie pour modele mono-mode
	     READ(PROG_ARG,*,ERR=990) IGRANU
	     DEF_IGRANU=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.MMD.SDparam1") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 		!1er parametre de granulometrie pour modele mono-mode
	     READ(PROG_ARG,*,ERR=990) VARGRANU1  
	     DEF_VARGRANU1=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.SDparam2") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!2nd parametre de granulometrie pour modele mono-mode
	     READ(PROG_ARG,*,ERR=990) VARGRANU2  
	     DEF_VARGRANU2=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.Model") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Indice de choix du modele WMO
	     READ(PROG_ARG,*,ERR=990) IMODELE_WMO 
	     DEF_IMODELE_WMO=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.DL") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Pourcentage volumique des particules DL du modele WMO
	     READ(PROG_ARG,*,ERR=990) C_WMO_DL
	     DEF_WMO_DL=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.WMO.WS") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)    	!Pourcentage volumique des particules WS du modele WMO
	     READ(PROG_ARG,*,ERR=990) C_WMO_WS
	     DEF_WMO_WS=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.OC") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   	!Pourcentage volumique des particules OC
	     READ(PROG_ARG,*,ERR=990) C_WMO_OC
	     DEF_WMO_OC=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.SO") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   	!Pourcentage volumique des particules SO
	     READ(PROG_ARG,*,ERR=990) C_WMO_SO
	     DEF_WMO_SO=1
	 ENDIF	
	 
	 IF (KEYWORD.EQ."-AER.SF.Model") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Indice de choix du modele S&F
	     READ(PROG_ARG,*,ERR=990) IMODELE_SF  
	     DEF_IMODELE_SF=1
	 ENDIF 	
	 
	 IF (KEYWORD.EQ."-AER.SF.RH") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Humidite relative pour modele S&F
	     READ(PROG_ARG,*,ERR=990) RH 
	     DEF_RH=1
	 ENDIF 	 
	 	 
	 IF (KEYWORD.EQ."-AER.BMD.VCdef") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Indice de type de description du melange bimodal LND
	     READ(PROG_ARG,*,ERR=990) MODE_PARAM_BILND 
	     DEF_MODE_PARAM_BILND=1
	 ENDIF 	 	 

	 IF (KEYWORD.EQ."-AER.BMD.CoarseVC") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Concentration volumique du mode "coarse"
	     READ(PROG_ARG,*,ERR=990) USER_CV_COARSE
	     DEF_USER_CV_COARSE=1
	 ENDIF 	 

	 IF (KEYWORD.EQ."-AER.BMD.FineVC") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Concentration volumique du mode "fine"
	     READ(PROG_ARG,*,ERR=990) USER_CV_FINE 
	     DEF_USER_CV_FINE=1
	 ENDIF 	 	 

	 IF (KEYWORD.EQ."-AER.BMD.RAOT") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	 !AOT_coarse / AOT_tot pour la longueur d'onde WAREF
	     READ(PROG_ARG,*,ERR=990) RTAUct_WAREF 
	     DEF_RTAUct_WAREF=1
	 ENDIF 	 

	 IF (KEYWORD.EQ."-AER.Waref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	 !longueur d'onde de reference WAREF
	     READ(PROG_ARG,*,ERR=990) WAREF
	     DEF_WAREF=1
	 ENDIF 	 
	 
	 IF (KEYWORD.EQ."-AER.AOTref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !AOT totale pour la longueur d'onde WAREF 
	     READ(PROG_ARG,*,ERR=990) TA_WAREF
	     DEF_TA_WAREF=1
	 ENDIF 	 

	 IF (KEYWORD.EQ."-AER.BMD.CM.MRwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_CM_MRWA
	     DEF_BMD_CM_MRWA=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.CM.MIwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_CM_MIWA
	     DEF_BMD_CM_MIWA=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.CM.MRwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_CM_MRWAREF
	     DEF_BMD_CM_MRWAREF=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.CM.MIwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_CM_MIWAREF
	     DEF_BMD_CM_MIWAREF=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.CM.SDradius") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_CM_RMODAL
	     DEF_BMD_CM_RMODAL=1
	 ENDIF 

	 IF (KEYWORD.EQ."-AER.BMD.CM.SDvar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_CM_VAR
	     DEF_BMD_CM_VAR=1
	 ENDIF 
	 	 
	 IF (KEYWORD.EQ."-AER.BMD.FM.MRwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_FM_MRWA
	     DEF_BMD_FM_MRWA=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.FM.MIwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_FM_MIWA	
	     DEF_BMD_FM_MIWA=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.FM.MRwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_FM_MRWAREF
	     DEF_BMD_FM_MRWAREF=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.FM.MIwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_FM_MIWAREF
	     DEF_BMD_FM_MIWAREF=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.FM.SDradius") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_FM_RMODAL  
	     DEF_BMD_FM_RMODAL=1
	 ENDIF 

	 IF (KEYWORD.EQ."-AER.BMD.FM.SDvar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=990) BMD_FM_VAR	 
	     DEF_BMD_FM_VAR=1
	 ENDIF 	
	 
	 IF (KEYWORD.EQ."-AER.ExtData") THEN
	     CALL GETARG (NPARAM+1,FICEXTDATA) 
	     DEF_FICEXTDATA=1
	 ENDIF 	 	 	 	 	 	 		  	 	 	 	 	 	 	 		  
      ENDDO	!Fin de boucle de lecture des arguments
	        	       	    		
	 
C --- Gestion de la coherence des arguments 
C---- (hors orientation vers choix du modele : gere cas par cas dans les blocs correspondant)

C     Argument TA obligatoire, controle de sa valeur, controle des autres arguments obligatoires  
      IF (DEF_TA.EQ.0) THEN
          GOTO 801
      ELSE   !Cas : valeur de TA definie 
      		
          IF (TA.LT.0.D-00) GOTO 1007  !Sortie en erreur
	  
	  !Cas TA >= 0: Il faut que le nom du fichier resultat soit defini
	  IF (DEF_FICGRANU.EQ.0) GOTO 802
        
	  !Argument FICANGLES obligatoire       
          IF (DEF_FICANGLES.EQ.0) GOTO 800   
	             
	  !Cas TA > 0:
	  IF (TA.GT.0.D-00) THEN
	     
	     !Si TA non nul : il faut que l'option de Trace soit definie
	     IF (DEF_FICLOGGRANU.EQ.0) THEN
                GOTO 803   
	     ELSE !Determination de la requete de suivi des traitements
		IF ((FICLOGGRANU(1:1).EQ.'0').AND.(FICLOGGRANU(2:2).EQ.' ')) THEN
	           TRACE=.FALSE.
                ELSE
                   TRACE=.TRUE. 
		ENDIF
             ENDIF 
	     
	     !Si TA non nul : Argument WA requis, controle de sa valeur 
	     IF (DEF_WA.EQ.0) THEN
                GOTO 804
             ELSE    
                IF ((WA.LT.0.2D-00).OR.(WA.GT.4.0D-00)) GOTO 1008
             ENDIF 
	     
	     !Si TA non nul : Choix du modele d'aerosols requis, controle de sa valeur 
	     IF (DEF_IMOD.EQ.0) THEN
                GOTO 805
             ELSE    
                IF ((IMOD.LT.0).OR.(IMOD.GT.SOS_MODALITE)) GOTO 1006 
             ENDIF 
	     
	     !Si TA non nul : Option de troncature requise, controle de sa valeur 
	     IF (DEF_ITRONC.EQ.0) THEN
                GOTO 806
             ELSE    
                IF ((ITRONC.LT.0).OR.(ITRONC.GT.1)) GOTO 1014
             ENDIF 
	     
	     !Si TA non nul : il faut que l'option de Trace MIE soit definie
	     IF (DEF_FICLOGMIE.EQ.0) THEN
                GOTO 807
             ENDIF 
	     
	     !Si TA non nul : il faut que le repertoire des fichiers de MIE soit defini
	     IF (DEF_DIRS.EQ.0) THEN
                GOTO 808
             ENDIF 	     
	     	     	     
	  ENDIF  !Fin test TA > 0
	  
      ENDIF !Fin test sur definition de TA
 
      	  
C* Ouverture du fichier Trace 
C-----------------------------
      IF (TRACE) THEN
         OPEN (88,FILE=FICLOGGRANU,ERR=911)
      END IF

C* Initialisation de l'indice d'erreur
C* -----------------------------------
      IER=0

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
	  READ(10,*,ERR=931) NVAL, XMU(J), XHR(J)
      ENDDO 
      CLOSE(10)
      
      DO J=1,MIE_NBMU
         XMU(-J) = -XMU(J)
	 XHR(-J) =  XHR(J)
      ENDDO 
                
        
      IF (TRACE) THEN
         WRITE (88,*,ERR=921) "ANGLES UTILISES"
	 WRITE (88,*,ERR=921) "****************"
	 WRITE (88,*,ERR=921) " "
	 WRITE (88,*,ERR=921) "Nombre d'angles au total : ", MIE_NBMU
	 WRITE (88,*,ERR=921) "Nombre d'angles de Gauss : ", NBMU_GAUSS
	 WRITE (88,*,ERR=921) " "
         IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	     WRITE(88,*,ERR=921)  
     &	       "Usage du fichier d'angles complementaires :",
     &         FICANGLES_USER
     	     WRITE(88,*,ERR=921) "Nb d'angles complementaires :",
     &         (MIE_NBMU-NBMU_GAUSS)
	 ELSE
	     WRITE(88,*,ERR=921)  
     &	       "Pas de fichier d'angles complementaires"
	 ENDIF
         WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) "Ordre max NB ", OS_NB
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,err=921) "Liste des angles utilises :"
	 WRITE (88,*,err=921) "   Indice, mu, angle (degres), poids"
	 DO 1000 J=-MIE_NBMU,MIE_NBMU
	    IF (J.EQ.0) GOTO 1000
	    IF (J.LT.0) THEN 
	        WRITE(88,222,ERR=921) J, XMU(J),
     &	               -DACOS(-XMU(J))*180./SOS_PI, 
     &			XHR(J)
            ELSE
	        WRITE(88,222,ERR=921) J, XMU(J),
     &	                DACOS(XMU(J))*180./SOS_PI, 
     &			XHR(J)
	    ENDIF	    
1000     CONTINUE
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) " "  	
      ENDIF !Fichier Trace	 


C* Recuperation du nom du repertoire contenant les fichiers de parametres WMO et Shettle & Fenn
C----------------------------------------------------------------------------------------------
      CALL GETENV('SOS_RACINE_FIC',DIRE)
      IFIN=INDEX(DIRE,' ')
      IFIN=IFIN-1
      IF (IFIN.LE.0) IFIN=SOS_LENDIR
            
C* Initialisation des proportions de constituants des modeles WMO et Shettle & Fenn
C----------------------------------------------------------------------------------      
      DO I=1,SOS_NB_WMO
	  C(I)=0.
	  N(I)=0.
      ENDDO
      
       DO I=1,SOS_NB_SF
	  Ni(I)=0.
      ENDDO

      NTOT=0.
            
C* Initialisation des tableaux de parametres de la matrice de phase
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
                        
C* Calcul des parametres de la matrice de phase des aerosols
C* et de l'albedo de simple diffusion.
C------------------------------------------------------------
      IF (TA.EQ.0.D+00) THEN	!Pas de calcul sur la granulometrie
      
      	PIZ = 0.D+00
	GOTO 1111	!Enregistrement du fichier resultat GRANU
      
      ENDIF

      !Poursuite des traitements pour TA > 0
      !=====================================
      
      	
C* Longueur du nom du repertoire des fichiers de Mie
C---------------------------------------------------
      LDIRS = INDEX(DIRS,' ')
      LDIRS = LDIRS - 1
      IF(LDIRS.LE.0) LDIRS = 100
      

C ===================
C* Calcul mono-modal
C ===================

      IF(IMOD.EQ.0) THEN
	 
C* Controle de la recuperation des parametres specifiques au calcul mono modal
C  ---------------------------------------------------------------------------
	 IF (DEF_ALPHAF.EQ.0) THEN
             GOTO 809
	 ELSE ! Controle de la plage de valeur des parametres de taille 
	     IF (ALPHAO.GT.ALPHAF) GOTO 1009
             IF (ALPHAF.GE.1D5) GOTO 1010      
         ENDIF 
		
	 !Si absence de fichier de Mie utilisateur precalcule
	 IF (DEF_FICS.EQ.0) THEN
             FICS="DEFAULT"     
         ENDIF 

	 IF ((DEF_RN.EQ.0).OR.(DEF_IN.EQ.0)) GOTO 820	 
	 IF (DEF_IGRANU.EQ.0) GOTO 822
	 IF ((DEF_VARGRANU1.EQ.0).OR.(DEF_VARGRANU2.EQ.0)) GOTO 823
	 	
      
C* Controle le signe de l'argument IN
C  -----------------------------------------
	 IF (IN.GT.0.D+00) GOTO 10100
	 
C* Controle la validite de l'argument IGRANU
C  -----------------------------------------
	 IF ((IGRANU.LT.0).OR.(IGRANU.GT.2)) GOTO 1011	 

C* Valeurs de RN et IN forcees au format F5.3 et F8.5
C  -----------------------------------------
	 RN = DNINT(RN*1000.D+00) / 1000.D+00
	 IN = -	DNINT(-IN*100000.D+00) / 100000.D+00 	 

      IF (TRACE) THEN
         WRITE (88,*,ERR=921) "MODELE MONO-MODAL"
	 WRITE (88,*,ERR=921) "*****************"
	 WRITE (88,*,ERR=921) " "
	 WRITE(88,*,ERR=921) 
     & 	 "Indice de refraction (partie reelle) : ",RN
	 WRITE(88,*,ERR=921) 
     & 	 "Indice de refraction (partie imaginaire) : ",IN
	 IF(IGRANU.EQ.1) THEN
	    WRITE(88,*,ERR=921) "Granulometrie : log-normale"
	    WRITE(88,*,ERR=921) "  Rayon modal (microns) : " , 
     &      VARGRANU1 
	    WRITE(88,*,ERR=921) "  Variance (sig) : " , VARGRANU2 
	    WRITE(88,*,ERR=921) "  sig / ln10  : " , 
     &      VARGRANU2  / DLOG(10.D+00)
	 ELSE
	    WRITE(88,*,ERR=921) "Granulometrie : loi de Junge"
	    WRITE(88,*,ERR=921) "  Rayon minimal (microns) : ", 
     &      VARGRANU1
	    WRITE(88,*,ERR=921) "  Puissance : " , VARGRANU2 
	 ENDIF   
	 WRITE (88,*,err=921) "Longueur d'onde en microns : ",WA
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) " "  		
      ENDIF !Fichier Trace	 
      	 
C* Constitution du nom du fichier MIE
C  ----------------------------------
	 IF(FICS.EQ.'DEFAULT') THEN
	    CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &                          RN,IN,ALPHAO,ALPHAF,FICS) 	     
	 ENDIF
	    
	 FICMIE = DIRS(1:LDIRS)//'/'//FICS	 
	 LFICS=INDEX(FICMIE,' ')-1
	 WRITE(6,*) FICMIE(1:LFICS)	
	 
C* Si le fichier MIE existe deja, le programme n'est pas lance
C  -----------------------------------------------------------	 
         INQUIRE(FILE=FICMIE,EXIST=EX)
	 
         IF (EX) THEN  
            WRITE(6,*)'Fichier MIE deja calcule'
	  
	 ELSE   
      
C* Calcul du fichier de MIE
C  ------------------------ 
            CALL SOS_MIE(MIE_NBMU,XMU,XHR,RN,IN,ALPHAO,ALPHAF,
     &	                 FICMIE,FICLOGMIE,IER)
            IF (IER.NE.0) GOTO 993   ! cas d'erreur
	       
         ENDIF
	 
C* Calcul des sections efficaces d'extinction et de diffusion 
C  et des fonctions de phase I, Q, U
C  ----------------------------------------------------------          
         CALL SOS_GRANU(FICMIE,IGRANU,VARGRANU1,VARGRANU2,WA,
     &                  MIE_NBMU,XMU,TRACE,KMAT1,KMAT2,SOMME_NR,
     &                  P11,P12,P33,IER)
         IF (IER.NE.0) GOTO 995   ! cas d'erreur 

C* Cas de particules spheriques : P22(J) = P11(J)
C  -----------------------------------------------------
      DO 1021 J=-MIE_NBMU,MIE_NBMU
         P22(J)=P11(J)
 1021 CONTINUE	 
  
C* Calcul des parametres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				   P11,TTT,P12,P22,P33,A,Z1,
     &				   ALP,BETA11,BETA22,
     &                             GAMMA12,DELTA33,ZETA,IER)  
         IF (IER.NE.0) GOTO 996   ! cas d'erreur

      ENDIF	!Fin test sur Modalite IMOD=0 (Mono-modal)

C ====================================================
C* Calcul multi-modal avec les constituants de la WMO
C ====================================================
	    
      IF(IMOD.EQ.1) THEN
	 
C* Controle de la recuperation des parametres specifiques au calcul WMO
C  --------------------------------------------------------------------       
	 IF (DEF_IMODELE_WMO.EQ.0) GOTO 830
	     
C* Controle la validite du parametre IMODELE
C  -----------------------------------------	    
	 IF (IMODELE_WMO.LT.1) GOTO 1012
	 IF (IMODELE_WMO.GT.SOS_NB_MODELE) GOTO 1012	
	 
C* Initialisation
C  --------------
         ! Nom des differentes particules de la WMO
	 PARTICULE(1)='DUST LIKE'
         PARTICULE(2)='WATER SOLUBLE'
         PARTICULE(3)='OCEANIC'
         PARTICULE(4)='SOOT'

         ALPHAF_PARTICULE(1) = SOS_ALPHAMAX_WMO_DL
	 ALPHAF_PARTICULE(2) = SOS_ALPHAMAX_WMO_WS
	 ALPHAF_PARTICULE(3) = SOS_ALPHAMAX_WMO_OC
	 ALPHAF_PARTICULE(4) = SOS_ALPHAMAX_WMO_SO
	 
	 IGRANU =   1	! distribution Log Normal par defaut
	    
C* Initialisation des parametres de la WMO
C  ---------------------------------------
	 FICDATAWMO=DIRE(1:IFIN)//'/'//SOS_AER_DATAWMO

	 CALL SOS_INIT_PARAMWMO(FICDATAWMO,WA,VARGRANU1I,
     &                          VARGRANU2I,MR,MI,V,IER)
         IF (IER.NE.0) GOTO 994   ! cas d'erreur 

C* Definition de la proportion volumique de chaque constituant
C  pour le modele selectionne
C  -----------------------------------------------------------	 
       	 IF(IMODELE_WMO.EQ.1) THEN	!Modele Continental
	    C(1) = 0.70
	    C(2) = 0.29
	    C(4) = 0.01      
	    
       	 ELSE IF(IMODELE_WMO.EQ.2) THEN 	!Modele Maritime
	    C(2) = 0.05
	    C(3) = 0.95
	 
         ELSE IF(IMODELE_WMO.EQ.3) THEN	!Modele Urbain
	    C(1) = 0.17
	    C(2) = 0.61
	    C(4) = 0.22
	       
         ELSE IF(IMODELE_WMO.EQ.4) THEN	!Modele defini par l'utilisateur
	 
C*          Controle de la recuperation des proportions volumiques 
C*          definies par l'utilisateur
C           ------------------------------------------------------
	    IF ( (DEF_WMO_DL.EQ.0).OR.(DEF_WMO_WS.EQ.0)
     &	          .OR.(DEF_WMO_OC.EQ.0).OR.(DEF_WMO_SO.EQ.0) ) GOTO 825 
	    
	    C(1) = C_WMO_DL    !Pourcentage volumique des particules DL
	    C(2) = C_WMO_WS    !Pourcentage volumique des particules WS
	    C(3) = C_WMO_OC    !Pourcentage volumique des particules OC
	    C(4) = C_WMO_SO    !Pourcentage volumique des particules SO	       	    

         ENDIF	!Fin definition des modeles
	    
         DO 103 I=1,SOS_NB_WMO	    
            N(I) = C(I)/V(I)
	    NTOT = NTOT + N(I)
  103    CONTINUE

      IF (TRACE) THEN
         WRITE (88,*,ERR=921) "MODELE WMO"
	 WRITE (88,*,ERR=921) "**********"
	 WRITE (88,*,ERR=921) " "
	 IF(IMODELE_WMO.EQ.1) THEN
            WRITE(88,*,ERR=921) "Modele CONTINENTAL MOYEN"
         ENDIF
	 IF(IMODELE_WMO.EQ.2) WRITE(88,*,ERR=921) "Modele MARITIME"
	 IF(IMODELE_WMO.EQ.3) WRITE(88,*,ERR=921) "Modele URBAIN"
	 WRITE (88,*,ERR=921) " "
	 WRITE (88,*,ERR=921) "  Coefficients alphaK du mélange :"
	 DO I=1,SOS_NB_WMO
	    WRITE (88,*,ERR=921) "    ",PARTICULE(I)," : ", N(I)/NTOT
         ENDDO
	 WRITE (88,*,ERR=921) " "
	 IF(IMODELE_WMO.EQ.4) THEN
	     WRITE (88,*,ERR=921) "Modele UTILISATEUR "
	     DO I=1,4
	        WRITE (88,*,ERR=921) "    ",PARTICULE(I)," : ",
     &		                     100*C(I)," % volume"
             ENDDO
	 ENDIF 
	 WRITE (88,*,err=921) "Longueur d'onde en microns : ",WA
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) " "  		
      ENDIF !Fichier Trace
      
      

C* Pour chaque type d'aerosols de la WMO
C  -------------------------------------	       
	 DO 105 I=1,SOS_NB_WMO
	    
	    IF(C(I).EQ.0.) GOTO 105
		
C* Constitution du nom du fichier MIE
C  ----------------------------------
	    ALPHAF = ALPHAF_PARTICULE(I)

            CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                        MR(I),MI(I),ALPHAO,ALPHAF,FICS)       
	    FICMIE = DIRS(1:LDIRS)//'/'//FICS
            LFICS=INDEX(FICMIE,' ')-1
	    WRITE(6,*) FICMIE(1:LFICS)	
	    
C* Si le fichier MIE existe deja, le programme n'est pas lance
C  -----------------------------------------------------------		              
            INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	    IF (EX) THEN  
               WRITE(6,*)'Fichier MIE deja calcule'
	  
	    ELSE   
      
C* Calcul du fichier de MIE
C  ------------------------ 
              CALL SOS_MIE(MIE_NBMU,XMU,XHR,MR(I),MI(I),ALPHAO,ALPHAF,
     &                      FICMIE,FICLOGMIE,IER)
    	       IF (IER.NE.0) GOTO 993   ! cas d'erreur
         
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
     &                     VARGRANU2I(I),WA,MIE_NBMU,XMU,TRACE,
     &			   KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
            IF (IER.NE.0) GOTO 995   ! cas d'erreur  
	    

C* Calcul des sections efficaces d'extinction et de diffusion du melange
C  ---------------------------------------------------------------------	     
	    KMAT1 = KMAT1 + (N(I)/NTOT)*KMAT1I
	    KMAT2 = KMAT2 + (N(I)/NTOT)*KMAT2I
	    
	    DO 104 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) + (N(I)/NTOT)*P11I(J)*KMAT2I
	       P12(J) = P12(J) + (N(I)/NTOT)*P12I(J)*KMAT2I
	       P33(J) = P33(J) + (N(I)/NTOT)*P33I(J)*KMAT2I
  104	    CONTINUE  !Fin boucle sur les angles de gauss
    
  105	 CONTINUE	 !Fin boucle sur les constituants de la WMO


C* Normalisation des fonctions de phase du melange
C  -----------------------------------------------	    
	 DO 106 J=-MIE_NBMU,MIE_NBMU   
	    P11(J) = P11(J) / KMAT2
	    P12(J) = P12(J) / KMAT2
	    P33(J) = P33(J) / KMAT2
  106	 CONTINUE

C* Cas de particules spheriques : P22(J) = P11(J)
C  -----------------------------------------------------
      DO 1022 J=-MIE_NBMU,MIE_NBMU
         P22(J)=P11(J)
 1022 CONTINUE	 
 
C* Calcul des parametres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------  
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				   P11,TTT,P12,P22,P33,A,Z1,
     &				   ALP,BETA11,BETA22,
     &                             GAMMA12,DELTA33,ZETA,IER)  
         IF (IER.NE.0) GOTO 996   ! cas d'erreur


      ENDIF	!Fin test sur Modalite IMOD=1 (WMO) 
  
       
C =====================================================================
C* Calcul bi-modal avec les constituants du modele Shettle & Fenn (S&F)
C =====================================================================
      
      IF(IMOD.EQ.2) THEN
      
C* Controle de la recuperation des parametres specifiques au calcul S&F
C  --------------------------------------------------------------------
         IF (DEF_IMODELE_SF.EQ.0) GOTO 840
	 IF (DEF_RH.EQ.0) GOTO 841

C* Controle la validite des parametres IMODELE et RH
C  -------------------------------------------------		 
	 IF (IMODELE_SF.LT.1) GOTO 1012				!Erreur choix modele
	 IF (IMODELE_SF.GT.SOS_NB_MODELE) GOTO 1012  		!Erreur choix modele
	 IF ((RH.LT.0.D-00).OR.(RH.GT.99.D-00)) GOTO 1013       !Erreur sur % humidite
	 
C* Initialisation
C  --------------
         ! Nom des differentes particules du modele S&F
	 PARTICULE(1)='SMALL_RURAL'
         PARTICULE(2)='LARGE RURAL'
         PARTICULE(3)='SMALL_URBAN'
         PARTICULE(4)='LARGE_URBAN'
	 PARTICULE(5)='OCEANIC_MODEL'

	 IGRANU =   1	! distribution Log Normal par defaut
	    
C* Initialisation des parametres du modele S&F
C  -------------------------------------------
	 CALL SOS_INIT_PARAMSF(DIRE,WA,RH,VARGRANU1I,
     &                         VARGRANU2I,MR,MI,IER)
         IF (IER.NE.0) GOTO 997   ! cas d'erreur 

C* Definition de la densite de chaque constituant pour le 
C  modele selectionne
C  -----------------------------------------------------------	 
       	 IF(IMODELE_SF.EQ.1) THEN		!Modele Tropospherique
	    Ni(1) = 1.0     
	    
       	 ELSE IF(IMODELE_SF.EQ.2) THEN 	!Modele Urbain
	    Ni(3) = 0.999875
	    Ni(4) = 0.000125
	 
         ELSE IF(IMODELE_SF.EQ.3) THEN	!Modele Maritime
	    Ni(1) = 0.99
	    Ni(5) = 0.01
	       
         ELSE IF(IMODELE_SF.EQ.4) THEN	!Modele cotier
	    Ni(1) = 0.995
	    Ni(5) = 0.005
	       
         ENDIF	!Fin definition des modeles S&F


      IF (TRACE) THEN
         WRITE (88,*,ERR=921) "MODELE SHETTLE & FENN"
	 WRITE (88,*,ERR=921) "*********************"
	 WRITE (88,*,ERR=921) " "
	 IF(IMODELE_SF.EQ.1) WRITE (88,*,ERR=921) "Modele TROPOSPHERIQUE"
	 IF(IMODELE_SF.EQ.2) WRITE (88,*,ERR=921) "Modele URBAIN"
	 IF(IMODELE_SF.EQ.3) WRITE (88,*,ERR=921) "Modele MARITIME"
	 IF(IMODELE_SF.EQ.4) WRITE (88,*,ERR=921) "Modele COTIER"
	 WRITE (88,*,ERR=921) " "
	 WRITE (88,*,ERR=921) "  Coefficients alphaK du mélange :"
	 DO I=1,SOS_NB_SF
	    WRITE (88,*,ERR=921) "    ",PARTICULE(I)," : ", Ni(I)
         ENDDO
         WRITE (88,*,ERR=921) " "
	 WRITE (88,*,err=921) "Pourcentage d'humidite relative : ",RH
	 WRITE (88,*,err=921) "Longueur d'onde en microns : ",WA 
	 WRITE (88,*,ERR=921) " "   
	 WRITE (88,*,ERR=921) " "		
      ENDIF !Fichier Trace
      
C* Pour chaque type d'aerosols du modele S&F
C  ----------------------------------------     
	 DO 109 I=1,SOS_NB_SF
	    
	    IF(Ni(I).EQ.0.) GOTO 109
		
C* Constitution du nom du fichier MIE
C  ----------------------------------	   
            IF (I.EQ.1) ALPHAF = SOS_ALPHAMAX_SF_SR
	    IF (I.EQ.3) ALPHAF = SOS_ALPHAMAX_SF_SU
	    IF ((I.NE.1).AND.(I.NE.3)) THEN
	        X = VARGRANU2I(I)  
	        RMAX = VARGRANU1I(I) * DEXP(X*X) *
     &	               DEXP(X*DSQRT(-2.D+00*LOG(SOS_COEF_NRMAX)))
     	
	        ALPHAF = REAL( 100 +
     &           100*DINT(2.*SOS_PI*RMAX / (100.*WA)) )
	    ENDIF

            CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                        MR(I),MI(I),ALPHAO,ALPHAF,FICS)       
	    FICMIE = DIRS(1:LDIRS)//'/'//FICS
            LFICS=INDEX(FICMIE,' ')-1
	    WRITE(6,*)FICMIE(1:LFICS)	
	    
C* Si le fichier MIE existe deja, le programme n'est pas lance
C  -----------------------------------------------------------           	              
            INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	    IF (EX) THEN  
               WRITE(6,*)'Fichier MIE deja calcule'
	  
	    ELSE   
      
C* Calcul du fichier de MIE
C  ------------------------ 
               CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &	                    MR(I),MI(I),ALPHAO,ALPHAF,
     &                      FICMIE,FICLOGMIE,IER)
     	       IF (IER.NE.0) GOTO 993   ! cas d'erreur
         
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
		 WRITE(88,*,err=921) 'Longueur d''onde (microns) : ',WA
		 WRITE(88,*,err=921) 'Parametre de taille max : ',ALPHAF
	       ENDIF
	    ENDIF
	       
            CALL SOS_GRANU(FICMIE,IGRANU,VARGRANU1I(I),
     &                     VARGRANU2I(I),WA,MIE_NBMU,XMU,TRACE,
     &			   KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
            IF (IER.NE.0) GOTO 995   ! cas d'erreur  
	     
C* Calcul des sections efficaces d'extinction et de diffusion du melange
C  ---------------------------------------------------------------------
	    KMAT1 = KMAT1 + Ni(I)*KMAT1I
	    KMAT2 = KMAT2 + Ni(I)*KMAT2I
	    
	    DO 108 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) + Ni(I)*P11I(J)*KMAT2I
	       P12(J) = P12(J) + Ni(I)*P12I(J)*KMAT2I
	       P33(J) = P33(J) + Ni(I)*P33I(J)*KMAT2I
  108	    CONTINUE  !Fin boucle sur les angles de gauss
    
  109	 CONTINUE	 !Fin boucle sur les constituants du modele Shettle et Fenn
	    

C* Normalisation des fonctions de phase du melange
C  -----------------------------------------------
	 DO 110 J=-MIE_NBMU,MIE_NBMU   
	    P11(J) = P11(J) / KMAT2
	    P12(J) = P12(J) / KMAT2
	    P33(J) = P33(J) / KMAT2
  110	 CONTINUE

C* Cas de particules spheriques : P22(J) = P11(J)
C  -----------------------------------------------------
      DO 1023 J=-MIE_NBMU,MIE_NBMU
         P22(J)=P11(J)
 1023 CONTINUE	 
 
C* Calcul des parametres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------  
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				   P11,TTT,P12,P22,P33,A,Z1,
     &				   ALP,BETA11,BETA22,
     &                             GAMMA12,DELTA33,ZETA,IER)  
         IF (IER.NE.0) GOTO 996   ! cas d'erreur	 

      	 
      ENDIF	!Fin test sur Modalite IMOD=2 (Shettle&Fenn) 




C =====================================================================
C* Calcul des proprietes radiatives d'une granulometrie LND bimodale
C =====================================================================
      
      IF(IMOD.EQ.3) THEN
	
C*    Controle de la recuperation des parametres specifiques au calcul bimodal LND
C     ----------------------------------------------------------------------------
         IF (DEF_MODE_PARAM_BILND.EQ.0) THEN
	     GOTO 850
         ELSE
	     !Controle la validite du parametre MODE_PARAM_BILND
	     IF (MODE_PARAM_BILND.LT.1) GOTO 1015
             IF (MODE_PARAM_BILND.GT.2) GOTO 1015
         ENDIF 	
       
C*       Recuperation de l'information de la description du melange
C        ----------------------------------------------------------      
C        Utilisation des concentrations volumiques de l'utilisateur
C         --> Controle que les parametres utilisateurs sont renseignes pour l'option
            IF (MODE_PARAM_BILND.EQ.1) THEN 
	         
	     IF (DEF_USER_CV_COARSE.EQ.0) GOTO 851 
	     CVI(1) = USER_CV_COARSE !mode "coarse"
		 
	     IF (DEF_USER_CV_FINE.EQ.0)   GOTO 851
	     CVI(2) = USER_CV_FINE   !mode "fine"
	  
         ENDIF
    
C        Utilisation du rapport AOT_coarse / AOT_totale
C        --> Controle que les parametres sont tous renseignes pour l'option
         IF (MODE_PARAM_BILND.EQ.2) THEN 
	     IF (DEF_RTAUct_WAREF.EQ.0) GOTO 852
	     IF (DEF_WAREF.EQ.0)        GOTO 852
	     IF (DEF_TA_WAREF.EQ.0)     GOTO 852
         ENDIF
   
	    
C*       Recuperation de l'information d'indice et granulometrie des 2 LND
C        -----------------------------------------------------------------         	 
         IF (DEF_BMD_CM_MRWA.EQ.0)     GOTO 853
         IF (DEF_BMD_CM_MIWA.EQ.0)     GOTO 853 
         IF (DEF_BMD_CM_RMODAL.EQ.0)   GOTO 853
         IF (DEF_BMD_CM_VAR.EQ.0)      GOTO 853 
         IF (DEF_BMD_FM_MRWA.EQ.0)     GOTO 854
         IF (DEF_BMD_FM_MIWA.EQ.0)     GOTO 854 
         IF (DEF_BMD_FM_RMODAL.EQ.0)   GOTO 854
         IF (DEF_BMD_FM_VAR.EQ.0)      GOTO 854 
      
         RN_LND(1) = BMD_CM_MRWA  !mode "coarse"
         RN_LND(2) = BMD_FM_MRWA  !mode "fine"
      
         IN_LND(1) = BMD_CM_MIWA  !mode "coarse"
         IN_LND(2) = BMD_FM_MIWA  !mode "fine"
      
         RMODAL_LND(1) = BMD_CM_RMODAL  !mode "coarse"
         RMODAL_LND(2) = BMD_FM_RMODAL  !mode "fine"
      
         VAR_LND(1) = BMD_CM_VAR  !mode "coarse"
         VAR_LND(2) = BMD_FM_VAR  !mode "fine"
      
      
         IF (MODE_PARAM_BILND.EQ.2) THEN 
	    IF (DEF_BMD_CM_MRWAREF.EQ.0) GOTO 855
	    IF (DEF_BMD_CM_MIWAREF.EQ.0) GOTO 855
	    IF (DEF_BMD_FM_MRWAREF.EQ.0) GOTO 855
	    IF (DEF_BMD_FM_MIWAREF.EQ.0) GOTO 855
	  
	    RN_LND_WAREF(1) = BMD_CM_MRWAREF !mode "coarse"
	    RN_LND_WAREF(2) = BMD_FM_MRWAREF !mode "fine"
	  
	    IN_LND_WAREF(1) = BMD_CM_MIWAREF !mode "coarse"
	    IN_LND_WAREF(2) = BMD_FM_MIWAREF !mode "fine"
         ENDIF
 
         

         DO 111 I=1,2
	 
C* Controle le signe de l'argument IN
C  -----------------------------------------
	    IF (IN_LND(I).GT.0.D+00) GOTO 10100	 

C* Valeurs de RN et IN forcees au format F5.3 et F8.5
C  -----------------------------------------
	    RN_LND(I) = DNINT(RN_LND(I)*1000.D+00) / 1000.D+00
	    IN_LND(I) = -DNINT(-IN_LND(I)*100000.D+00) / 100000.D+00 
	    
111      CONTINUE	 
	 



         IF (TRACE) THEN
            WRITE (88,*,ERR=921) "MODELE BIMODAL DE LND"
	    WRITE (88,*,ERR=921) "*********************"
	    WRITE (88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) "Granulometrie : log-normale coarse"
	    WRITE(88,*,ERR=921) "                ------------------"
	    WRITE(88,*,ERR=921) 
     & 	    "   Indice de refraction (partie reelle) : ",RN_LND(1)
	    WRITE(88,*,ERR=921) 
     & 	    "   Indice de refraction (partie imaginaire) : ",IN_LND(1)
	    WRITE(88,*,ERR=921) 
	    WRITE(88,*,ERR=921) "   Rayon modal (microns) : " , 
     &	    RMODAL_LND(1) 
            WRITE(88,*,ERR=921) "   Variance (sig) : " , VAR_LND(1) 
            WRITE(88,*,ERR=921) "   sig /ln10 : " , 
     &	    VAR_LND(1)/DLOG(10.D+00) 
	    WRITE (88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) "Granulometrie : log-normale fine"
	    WRITE(88,*,ERR=921) "                ----------------"
	    WRITE(88,*,ERR=921) 
     & 	    "   Indice de refraction (partie reelle) : ",RN_LND(2)
	    WRITE(88,*,ERR=921) 
     & 	    "   Indice de refraction (partie imaginaire) : ",IN_LND(2)
	    WRITE(88,*,ERR=921) 
	    WRITE(88,*,ERR=921) "   Rayon modal (microns) : " , 
     &	    RMODAL_LND(2) 
	    WRITE(88,*,ERR=921) "   Variance (sig)  : " , VAR_LND(2) 
	    WRITE(88,*,ERR=921) "   sig /ln10 : " , 
     &	    VAR_LND(2)/DLOG(10.D+00) 
	    WRITE(88,*,ERR=921) " "
	 
	    IF (MODE_PARAM_BILND.EQ.1) THEN 
	       WRITE(88,*,ERR=921) 
     &         "Representation de la proportion des composants" 
	       WRITE(88,*,ERR=921) "par les concentrations volumiques :"	    
	       WRITE(88,*,ERR=921) 
     & 	       "   Concentration volumique du mode LND coarse : ",CVI(1)
	       WRITE(88,*,ERR=921) 
     & 	       "   Concentration volumique du mode LND fine : ",CVI(2)
	       WRITE(88,*,ERR=921) " "
	    ENDIF 
	  
	    IF (MODE_PARAM_BILND.EQ.2) THEN 
	       WRITE(88,*,ERR=921) 
     &         "Representation de la proportion des composants" 
	       WRITE(88,*,ERR=921) 
     &        "par le rapport AOT_coarse / AOT_totale :" 
     	       WRITE(88,*,ERR=921) 
     &	       "   Longueur d'onde de reference (microns) : ",WAREF
     	       WRITE(88,*,ERR=921) 
     &	       "   --> Epaisseur optique aerosols totale associee : ",
     &         TA_WAREF	    
	       WRITE(88,*,ERR=921) 
     &         "   --> Rapport des AOT coarse / totale : ",RTAUct_WAREF
               WRITE(88,*,ERR=921) " "
	       WRITE(88,*,ERR=921) 
     &	       "   --> Indice de refraction du mode coarse "
	       WRITE(88,*,ERR=921) 
     & 	       "          - Partie reelle : ",RN_LND_WAREF(1)
	       WRITE(88,*,ERR=921) 
     & 	       "          - Partie imaginaire : ",IN_LND_WAREF(1)
	       WRITE(88,*,ERR=921) 
     &	       "   --> Indice de refraction du mode fin "
	       WRITE(88,*,ERR=921) 
     & 	       "          - Partie reelle : ",RN_LND_WAREF(2)
	       WRITE(88,*,ERR=921) 
     & 	       "          - Partie imaginaire : ",IN_LND_WAREF(2)
	       WRITE(88,*,ERR=921) " "
	    ENDIF  	        	    
	 
         ENDIF !Fichier Trace	
      
      	 
C* Initialisation
C  --------------
 	 IGRANU =   1	! distribution Log Normal par defaut


C* Calcul des coefficients CVI pour une expression du
C  melange en AOT_coarse / AOT_totale
C  --------------------------------------------------
         IF (MODE_PARAM_BILND.EQ.2) THEN
	    !Utilisation du rapport AOT_coarse / AOT_totale
 
	    IF (TRACE) THEN
	       WRITE (88,*,ERR=921) " "
	       WRITE (88,*,ERR=921) 
     &	       "Estimation des coefficients du melange"
	       WRITE (88,*,ERR=921) 
     &	       "======================================"
               WRITE (88,*,ERR=921) "  "
	    ENDIF
	    
C* Pour chacun des deux modes (1 : coarse, 2 : fine)
C  -------------------------------------------------  
            DO 114 I=1,2  
		
C* Constitution du nom du fichier MIE (pour WAREF)
C  -----------------------------------------------	               
	       X = VAR_LND(I) 
	       RMAX = RMODAL_LND(I) * DEXP(X*X) *
     &	              DEXP(X*DSQRT(-2.D+00*LOG(SOS_COEF_NRMAX)))

               ALPHAF = REAL( 100 +
     &	       	 100*DINT(2.*SOS_PI*RMAX / (100.*WAREF)) )  
	              
               CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                           RN_LND_WAREF(I),IN_LND_WAREF(I),
     &	                           ALPHAO,ALPHAF,FICS)       
	       FICMIE = DIRS(1:LDIRS)//'/'//FICS
               LFICS=INDEX(FICMIE,' ')-1
	       WRITE(6,*)FICMIE(1:LFICS)	
	    
C* Si le fichier MIE existe deja, le programme n'est pas lance
C  -----------------------------------------------------------           	              
               INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	       IF (EX) THEN  
                  WRITE(6,*)'Fichier MIE deja calcule'
	  
	       ELSE   
      
C* Calcul du fichier de MIE (pour WAREF)
C  ------------------------------------ 
                  CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &		               RN_LND_WAREF(I),IN_LND_WAREF(I),
     &                         ALPHAO,ALPHAF,
     &                         FICMIE,FICLOGMIE,IER)
     	          IF (IER.NE.0) GOTO 993   ! cas d'erreur
         
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
     &	       '      Longueur d''onde (microns) : ',WAREF
	       WRITE(88,*,err=921) 
     &	       '      Parametre de taille max : ',ALPHAF	  
	       WRITE(88,*,err=921) ' '     
	    ENDIF
	       
            CALL SOS_GRANU(FICMIE,IGRANU,RMODAL_LND(I),
     &                     VAR_LND(I),WAREF,MIE_NBMU,XMU,TRACE,
     &			   KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
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
     &	      '      Estimation : CV mode coarse (non normalisé):',CVI(1)
	       WRITE(88,*,err=921) 
     &	       '      Estimation : CV mode fine (non normalisé):',CVI(2)
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
	    WRITE(88,*,err=921)'      CV mode coarse (normalisé):', CVI(1)
	    WRITE(88,*,err=921)'      CV mode fine   (normalisé):', CVI(2)	  
	    WRITE(88,*,err=921) ' '  
	    WRITE(88,*,err=921) ' '   
	ENDIF

	    
C* Pour chacun des deux modes (1 : coarse, 2 : fine)
C  ------------------------------------------------- 
	IF (TRACE) THEN
	    WRITE (88,*,ERR=921) 
     &	    "Calcul des proprietes radiative du melange"
	    WRITE (88,*,ERR=921) 
     &	    "=========================================="
            WRITE (88,*,ERR=921) "  "
	 ENDIF
	    
	        
	 DO 112 I=1,2
	    	    
	    IF(CVI(I).EQ.0.) GOTO 112
		
C* Constitution du nom du fichier MIE (pour la longueur d'onde de simulation)
C  ----------------------------------	               
	    X = VAR_LND(I)
	    RMAX = RMODAL_LND(I) * DEXP(X*X) *
     &	           DEXP(X*DSQRT(-2.D+00*LOG(SOS_COEF_NRMAX)))
     	
	    ALPHAF = REAL( 100 +
     &        100*DINT(2.*SOS_PI*RMAX / (100.*WA)) )  

            CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     & 	                        RN_LND(I),IN_LND(I),ALPHAO,ALPHAF,FICS)  
	    FICMIE = DIRS(1:LDIRS)//'/'//FICS
            LFICS=INDEX(FICMIE,' ')-1
	    WRITE(6,*)FICMIE(1:LFICS)	
	    
C* Si le fichier MIE existe deja, le programme n'est pas lance
C  -----------------------------------------------------------           	              
            INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	    IF (EX) THEN  
               WRITE(6,*)'Fichier MIE deja calcule'
	  
	    ELSE   
      
C* Calcul du fichier de MIE  (pour la longueur d'onde de simulation)
C  ------------------------ 
               CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &                      RN_LND(I),IN_LND(I),ALPHAO,ALPHAF,
     &                      FICMIE,FICLOGMIE,IER)
     	       IF (IER.NE.0) GOTO 993   ! cas d'erreur
         
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
     &	       '      Longueur d''onde (microns) : ',WA
	       WRITE(88,*,err=921) 
     &	       '      Parametre de taille max : ',ALPHAF	  
	       WRITE(88,*,err=921) ' '     
	    ENDIF
	       
            CALL SOS_GRANU(FICMIE,IGRANU,RMODAL_LND(I),
     &                     VAR_LND(I),WA,MIE_NBMU,XMU,TRACE,
     &			   KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
            IF (IER.NE.0) GOTO 995   ! cas d'erreur  

		 
C* Calcul des coefficients d'extinction et de diffusion du melange
C  ---------------------------------------------------------------
	    KMAT1 = KMAT1 + CVI(I)*KMAT1I
	    KMAT2 = KMAT2 + CVI(I)*KMAT2I
	    
	    DO 113 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) + CVI(I)*P11I(J)*KMAT2I
	       P12(J) = P12(J) + CVI(I)*P12I(J)*KMAT2I
	       P33(J) = P33(J) + CVI(I)*P33I(J)*KMAT2I	       
  113	    CONTINUE  !Fin boucle sur les angles de gauss
    
  112	 CONTINUE	 !Fin boucle sur les mode coarse et fine de la LND bimodale
	    
	 		
		
C* Normalisation des fonctions de phase du melange
C  -----------------------------------------------
	 DO 120 J=-MIE_NBMU,MIE_NBMU   
	    P11(J) = P11(J) / KMAT2
	    P12(J) = P12(J) / KMAT2
	    P33(J) = P33(J) / KMAT2
  120	 CONTINUE

C* Cas de particules spheriques : P22(J) = P11(J)
C  -----------------------------------------------------
      DO 1024 J=-MIE_NBMU,MIE_NBMU
         P22(J)=P11(J)
 1024 CONTINUE	 
 
C* Calcul des parametres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------  
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				   P11,TTT,P12,P22,P33,A,Z1,
     &				   ALP,BETA11,BETA22,
     &                             GAMMA12,DELTA33,ZETA,IER)     
         IF (IER.NE.0) GOTO 996   ! cas d'erreur	 

          	 
      ENDIF	!Fin test sur Modalite IMOD=3 (LND bimodale) 




C ==============================================================================
C* Donnees externes :
C* --> Fonctions de phase externes et coefficients de diffusion et d'extinction
C*
C*     NB : Developpe sur la base des sorties de l'outil Oleg Dubovik (LOA)
C*          qui fournit les propriétés de diffusion primaire de particules
C*          spheroides homogenes orientees de facon aleatoire.
C*          Un parametre precisant le nombre d'angles decrivant les fonctions
C*          a toutefois du etre ajoute.
C*
C*          Bien faire attention au fait qu'on lit un rapport -F12 / F11.
C*          Importance du signe : F12 < 0 pour la diffusion rayleigh.

C ==============================================================================

      IF(IMOD.EQ.4) THEN
	 
C* Controle de la recuperation du nom du fichier de donnees externes
C  -----------------------------------------------------------------
	 IF (DEF_FICEXTDATA.EQ.0) GOTO 860
    	
C* Lecture du fichier de donnees externes
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
	 IF (NB_FICEXT_ANGLES.GT.SOS_MAXNB_ANG_EXT) GOTO 950
	
	 READ(20,*,ERR=941) 

	 DO J=1,NB_FICEXT_ANGLES
		
            READ(20,*,ERR=941,END=942) 
     &	                 ANG_EXT(J), F11_EXT(J), F12sF11_EXT(J),
     &	                 F22sF11_EXT(J), F33sF11_EXT(J)
         ENDDO 
	
	 CLOSE(20)

	 
	IF (TRACE) THEN
            WRITE (88,*,ERR=921) "DONNEES EXTERNES"
	    WRITE (88,*,ERR=921) "****************"
	    WRITE (88,*,ERR=921) " "
	    WRITE (88,*,ERR=921)  "Nom du fichier : ",FICEXTDATA
	    WRITE (88,*,ERR=921) " "  		 
            WRITE (88,*,ERR=921) "Coef extinction : ", KMAT1
	    WRITE (88,*,ERR=921) "Coef diffusion  : ", KMAT2
	    WRITE (88,*,ERR=921) "Nb d'angles : ", NB_FICEXT_ANGLES
	    WRITE (88,*,ERR=921) 
     &	    "  ANGLE     F11         -F12/F11      F22/F11      F33/F11"
	    DO J=1,NB_FICEXT_ANGLES
               WRITE (88,402,ERR=921) ANG_EXT(J), F11_EXT(J), 
     &	       			      F12sF11_EXT(J), F22sF11_EXT(J), 
     &	                              F33sF11_EXT(J)
            ENDDO 
	    WRITE (88,*,ERR=921) " "
         ENDIF !Fichier Trace
	 

C* Passage aux fonctions F12, F22 et F33 pour les mu des angles 
C* decrivant les fonctions externes
C  -----------------------------------------------------------------	 
	 DO J=1,NB_FICEXT_ANGLES
	    MU_EXT(J)  = DCOS(ANG_EXT(J)*SOS_PI/180.D+00)
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
     &	                        2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                          P11(-MIE_NBMU), IER)	
	IF (IER.NE.0) GOTO 945   ! cas d'erreur
	
C*      P12 :
C*      ---	
	CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F12_EXT,
     &	                        2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                          P12(-MIE_NBMU), IER)	
	IF (IER.NE.0) GOTO 945   ! cas d'erreur
	
C*      P22 :
C*      ---	
	CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F22_EXT,
     &	                        2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                          P22(-MIE_NBMU), IER)	
	IF (IER.NE.0) GOTO 945   ! cas d'erreur
		
C*      P33 :
C*      ---	
	CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F33_EXT,
     &	                        2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                          P33(-MIE_NBMU), IER)		
	IF (IER.NE.0) GOTO 945   ! cas d'erreur 
	
	!C* Ajouter module de recomposition pour fichier TRACE
	!C* Entetes de fichiers
	!AJUSTER DEFINITION DE LA LIMITE DE TABLE F11, etc. en dur ==> dans SOS.h
	!Attention signe : -F12sF11
	 
	 IF (TRACE) THEN
	    WRITE (88,*,ERR=921) " Fonctions F11, F12, F22 et F33 "
	    WRITE (88,*,ERR=921) " aux angles du fichier de donnees"
	    WRITE (88,*,ERR=921) " ================================"
	    WRITE (88,*,ERR=921) 
     &	    "  ANGLE     F11          F12          F22          F33"
	    DO J=1,NB_FICEXT_ANGLES
               WRITE (88,402,ERR=921) ANG_EXT(J), F11_EXT(J), 
     &	       			      F12_EXT(J), F22_EXT(J), F33_EXT(J)
            ENDDO 
	    WRITE (88,*,ERR=921) " " 
	    WRITE (88,*,ERR=921) " Fonctions F11, F12, F22 et F33 "
	    WRITE (88,*,ERR=921) " aux angles de Gauss "
	    WRITE (88,*,ERR=921) " ================================" 
	    DO J=-MIE_NBMU, MIE_NBMU
	       TETA=DACOS(XMU(J))*180./SOS_PI
               WRITE (88,402,ERR=921) TETA,P11(J),P12(J),P22(J),P33(J)
            ENDDO 	  		
         ENDIF !Fichier Trace	 
	 
	 	
	
	
	 
C* Calcul des parametres alpha, beta, gamma, delta, zeta
C  avec application ou non de la troncature
C  -----------------------------------------------------
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				   P11,TTT,P12,P22,P33,A,Z1,
     &				   ALP,BETA11,BETA22,
     &                             GAMMA12,DELTA33,ZETA,IER)  
         IF (IER.NE.0) GOTO 996   ! cas d'erreur

      ENDIF	!Fin test sur Modalite IMOD=4 (donnees externes)
      
      
      
      
      
      
C ======================================
C* Calculs communs pour tous les modeles
C ======================================
	 
C* Calcul de l'albedo de simple diffusion
C-------------------------------------------
      PIZ=KMAT2/KMAT1

C* Troncature de l'albedo de simple diffusion des aerosols
C---------------------------------------------------------
      PIZTR = PIZ*(1.-A/2.)/(1.-PIZ*A/2.)      
    
C* Ecriture du fichier trace
C-------------------------------------------
      IF (TRACE) THEN
         WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) '================================'
         WRITE(88,*,err=921) '* ROUTINE AEROSOLS : RESULTATS *'
	 WRITE(88,*,err=921) '================================'
         WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) ' '
	
	 IF(ITRONC.EQ.1) THEN
	      WRITE(88,*,err=921) '  Calcul avec troncature '
	 ELSE
	      WRITE(88,*,err=921) '  Calcul sans troncature '	   
	 ENDIF
	 
         WRITE(88,*,err=921) ' '
            	 
	 IF((A.NE.0.D+00).OR.(IMOD.NE.0)) THEN
	    WRITE(88,*,err=921) '  Coefficient de troncature : ',A
	    WRITE(88,*,err=921) 
     &	     '  Section efficace d''extinction non tronque (micron^2) : '
     &          ,KMAT1
	    WRITE(88,*,err=921) 
     &	     '  Section efficace de diffusion non tronque (micron^2) : ' 
     &          ,KMAT2
	    IF (A.NE.0.D+00) THEN
	       WRITE(88,*,err=921) 
     &		 '  Albedo de simple diffusion tronque : ',PIZTR 
	       WRITE(88,*,err=921) 
     &		 '  Albedo de simple diffusion non tronque : ',
     &             PIZTR / (1.-A/2.+PIZTR*A/2.)      
	       WRITE(88,*,err=921)
     &		 '  Beta11(0) avant normalisation tronque : ',Z1  
               WRITE(88,*,err=921) 
     &	         '  Parametre d''asymetrie tronque: ',BETA11(1)/3.
               WRITE(88,*,err=921) 
     &	        '  Parametre d''asymetrie non tronque: ',
     &          A/2.+(1.-A/2.)*BETA11(1)/3.     
	    ELSE
	       WRITE(88,*,err=921) 
     & 		 '  Albedo de simple diffusion non tronque : ',PIZTR
	       WRITE(88,*,err=921) 
     & 		 '  Beta11(0) avant normalisation non tronque : ',Z1
               WRITE(88,*,err=921) 
     &	         '  Parametre d''asymetrie non tronque: ',BETA11(1)/3.
     	    ENDIF !Fin test si troncature   
            WRITE(88,*,err=921) ' '
	    WRITE(88,*,err=921) '  1ere colonne : angle de diffusion '
	    WRITE(88,*,err=921) 
     &	      '  2eme colonne : fonction de phase tronquee' 	 
	    WRITE(88,*,err=921)
     &        '  3eme colonne : fonction de phase initiale '
	    WRITE(88,*,err=921) 
     &        '  4eme colonne : fonction de phase polarisee tronquee '
            WRITE(88,*,err=921) '		  Q(mu) tronquee'
	    WRITE(88,*,err=921)
     &        '  5eme colonne : fonction de phase polarisee initiale '
            WRITE(88,*,err=921) '		  Q(mu) '
	    WRITE(88,*,err=921)
     &        '  6eme colonne : fonction de phase polarisee tronquee '
            WRITE(88,*,err=921) '                 U(mu) tronquee'
	    WRITE(88,*,err=921) 
     &        '  7eme colonne : fonction de phase polarisee initiale '
            WRITE(88,*,err=921) '                 U(mu) '
	    WRITE(88,*,err=921) ' '	    
            DO 44 J=-MIE_NBMU,MIE_NBMU
	       TETA=DACOS(XMU(J))*180./SOS_PI
	       WRITE(88,10,err=921) TETA,P11(J),TTT(J),
     &                              P12(J)*P11(J)/TTT(J),P12(J),
     &                              P33(J)*P11(J)/TTT(J),P33(J)
   44       CONTINUE
         END IF
    
         WRITE(88,*,err=921) ' '  
	 WRITE(88,35,err=921)'I','ALP(I)','BETA11(I)','GAMMA12(I)',
     &                       'DELTA33(I)','ZETA(I)'  
         DO 45 K=0,OS_NB
	    WRITE(88,25,err=921) K,ALP(K),BETA11(K),GAMMA12(K),DELTA33(K),
     &	                         ZETA(K)
   45    CONTINUE
   	 WRITE(88,*,err=921) ' '
	 	   	   	   	   	   
      END IF

C* Enregistrement des parametres de granulometrie dans le fichier Resultat
C  -----------------------------------------------------------------------
 1111 OPEN(1,FILE=FICGRANU,ERR=913)

      ! Section efficace d'extinction
      WRITE(1,40,ERR=923) KMAT1

      ! Section efficace de diffusion
      WRITE(1,41,ERR=923) KMAT2
     
      ! Facteur d'asymétrie (non tronqué)
      WRITE(1,42,ERR=923) ( A/2.+(1.-A/2.)*BETA11(1)/3. )

      ! Coefficient de troncature      
      WRITE(1,46,ERR=923) A
      
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
      CLOSE(99)		! Fermeture du fichier trace MIE
      CLOSE(88)         ! Fermeture du fichier trace GRANU
      CLOSE(1)		! Fermeture du fichier resultat GRANU
           
C* Fin nominale 
C-------------------
      GOTO 9999
            
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
 800  WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AER.AngFile Value'
      CALL EXIT(1)
      GOTO 9999 
      
  801 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AER.AOT Value (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999  

  802 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AER.ResFile Value'
      CALL EXIT(1)
      GOTO 9999  
      
  803 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AER.Log Value'
      CALL EXIT(1)
      GOTO 9999  
            
  804 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AER.Wa Value (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999        
            
  805 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AER.Model Value (format : I2)'
      CALL EXIT(1)
      GOTO 9999        
            
  806 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AER.Tronca Value (format : I4)'
      CALL EXIT(1)
      GOTO 9999          
           
  807 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AER.MieLog Value'
      CALL EXIT(1)
      GOTO 9999          

  808 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -AER.DirMie Value'
      CALL EXIT(1)
      GOTO 9999          
           
  809 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Mono-modal model requires (-Keyword Value):',
     &           '      -AER.MMD.Mie.AlphaMax Value (format : F9.2)'
      CALL EXIT(1)
      GOTO 9999  
          
  820 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Mono-modal model requires (-Keyword Value):'
      WRITE(6,*) '        -AER.MMD.MRwa Value (format : F5.3)'
      WRITE(6,*) '        -AER.MMD.MIwa Value (format : F8.5)'
      CALL EXIT(1)
      GOTO 9999    
      
  822 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Mono-modal model requires (-Keyword Value):',
     &           '  :     -AER.MMD.SDtype Value (format : I4)'
      CALL EXIT(1)
      GOTO 9999      
      
  823 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Mono-modal model requires (-Keyword Value):'
      WRITE(6,*) '        -AER.MMD.SDparam1 Value (format : F9.5)'
      WRITE(6,*) '        -AER.MMD.SDparam2 Value (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999  
                    
  825 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> User volumetric concentrations for WMO models'
      WRITE(6,*) '      require  (-Keyword Value) definition :'
      WRITE(6,*) '         -AER.WMO.DL Value (format : F9.5)'
      WRITE(6,*) '         -AER.WMO.OC Value (format : F9.5)'
      WRITE(6,*) '         -AER.WMO.WS Value (format : F9.5)'
      WRITE(6,*) '  :      -AER.WMO.SO Value (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999  
                
  830 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> WMO model requires (-Keyword Value):',
     &           '  :     -AER.WMO.Model Value (format : I2)'
      CALL EXIT(1)
      GOTO 9999        	 
	   
  840 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Shettle&Fenn model requires (-Keyword Value):',
     &           '  :     -AER.SF.Model Value (format : I2)'
      CALL EXIT(1)
      GOTO 9999        	 

  841 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Shettle&Fenn model requires (-Keyword Value):',
     &           '  :     -AER.SF.RH Value (format : F5.2)'
      CALL EXIT(1)
      GOTO 9999  	         
	   
  850 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Bimodal LND model requires (-Keyword Value):',
     &           '  :     -AER.BMD.VCdef Value (format : I2)'
      CALL EXIT(1)
      GOTO 9999        	 
	   
  851 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> User volumetric concentration for bimodal LND'
      WRITE(6,*) '      requires (-Keyword Value): '
      WRITE(6,*) '        -AER.BMD.CoarseVC (format : F8.5)'
      WRITE(6,*) '        -AER.BMD.FineVC   (format : F8.5)'
      CALL EXIT(1)
      GOTO 9999       

  852 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Volumetric concentration for bimodal LND'
      WRITE(6,*) '      defined by AOT ratio coarse / total'
      WRITE(6,*) '      requires (-Keyword Value): '
      WRITE(6,*) '        -AER.BMD.RAOT (format : F9.5)'
      WRITE(6,*) '        -AER.Waref    (format : F9.5)'
      WRITE(6,*) '        -AER.AOTref   (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999  
      	 
  853 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Bimodal LND coarse mode for wavelength -AER.Wa'
      WRITE(6,*) '      requires (-Keyword Value): '
      WRITE(6,*) '        -AER.BMD.CM.MRwa (format : F5.3)'
      WRITE(6,*) '        -AER.BMD.CM.MIwa    (format : F8.5)'
      WRITE(6,*) '        -AER.BMD.CM.SDradius   (format : F9.5)'
      WRITE(6,*) '        -AER.BMD.CM.SDvar   (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999  	

  854 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Bimodal LND fine mode for wavelength -AER.Wa'
      WRITE(6,*) '      requires (-Keyword Value): '
      WRITE(6,*) '        -AER.BMD.FM.MRwa (format : F5.3)'
      WRITE(6,*) '        -AER.BMD.FM.MIwa    (format : F8.5)'
      WRITE(6,*) '        -AER.BMD.FM.SDradius   (format : F9.5)'
      WRITE(6,*) '        -AER.BMD.FM.SDvar   (format : F9.5)'
      CALL EXIT(1)
      GOTO 9999  		 

  855 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> Volumetric concentration for bimodal LND'
      WRITE(6,*) '      defined by AOT ratio coarse / total'
      WRITE(6,*) '      requires (-Keyword Value) '
      WRITE(6,*) '      for the refractive indexes of coarse'
      WRITE(6,*) '      and fine modes at wavelength -AER.Waref :'
      WRITE(6,*) '        -AER.BMD.CM.MRwaref (format : F5.3)'
      WRITE(6,*) '        -AER.BMD.CM.MIwaref (format : F8.5)'
      WRITE(6,*) '        -AER.BMD.FM.MRwaref (format : F5.3)'
      WRITE(6,*) '        -AER.BMD.FM.MIwaref (format : F8.5)'
      CALL EXIT(1)
      GOTO 9999	 

  860 WRITE(6,*) '  ERROR on parameters for SOS_AEROSOLS'
      WRITE(6,*) '  --> External phaze functions file -AER.ExtData'
      WRITE(6,*) '      requires (-Keyword Value)  '
      CALL EXIT(1)
      GOTO 9999  		 
      	 
  990 WRITE(6,*) '  ERROR on argument loading for SOS_AEROSOLS'
      WRITE(6,*) '  --> Check required formats.'
      CALL EXIT(1)
      GOTO 9999
  911 WRITE(6,*) '  ERROR on logfile opening for GRANU'
      CALL EXIT(1)
      GOTO 9999
  930 WRITE(6,*) '  ERROR on angles file opening for SOS_AEROSOLS'
      CALL EXIT(1)
      GOTO 9999
  931 WRITE(6,*) '  ERROR on angles file reading for SOS_AEROSOLS'
      CALL EXIT(1)
      GOTO 9999
  940 WRITE(6,*) '  ERROR on external phase functions file opening'
      WRITE(6,*) '  for SOS_AEROSOLS'
      CALL EXIT(1)
      GOTO 9999
  950 WRITE(6,*) '  ERROR : The number of angles for the external phase'
      WRITE(6,*) '  functions description is higher than the defined '
      WRITE(6,*) '  size of internal tables'
      WRITE(6,*) '  --> Check the value of SOS_MAXNB_ANG_EXT in SOS.h'
      WRITE(6,*) '      compared to the value associated to NbAng '
      WRITE(6,*) '      in the external phase functions file'
      CALL EXIT(1)
      GOTO 9999
  941 WRITE(6,*) '  ERROR on external phaze functions file reading'
      WRITE(6,*) '  for SOS_AEROSOLS'
      CALL EXIT(1)
      GOTO 9999 
  942 WRITE(6,*) '  ERROR on external phaze functions file reading'
      WRITE(6,*) '  for SOS_AEROSOLS : End of file encountered'
      WRITE(6,*) '  --> Check if the value associated to NbAng is '
      WRITE(6,*) '      correct, in the external phase functions file'
      CALL EXIT(1)
      GOTO 9999       
  945 WRITE(6,*) '  ERROR while interpolating external phaze functions'
      WRITE(6,*) '  for SOS_AEROSOLS'
      WRITE(6,*) '  Check that the abscissas of the phaze functions'
      WRITE(6,*) '  are correct : do not use twice same values'
      CALL EXIT(1)
      GOTO 9999             
  913 WRITE(6,*) '  ERROR on GRANU result file opening'
      CALL EXIT(1)
      GOTO 9999      
  921 WRITE(6,*) '  ERROR on logfile writing for GRANU'
      CALL EXIT(1)
      GOTO 9999
  922 WRITE(6,*) '  ERROR on Gauss angles file reading'
      CALL EXIT(1)
      GOTO 9999
  923 WRITE(6,*) '  ERROR on GRANU result file writing'
      CALL EXIT(1)
      GOTO 9999      
  993 WRITE(6,*) '  ERROR on subroutine SOS_MIE'
      CALL EXIT(1)
      GOTO 9999
  994 WRITE(6,*) '  ERROR on subroutine SOS_INIT_PARAMWMO'
      CALL EXIT(1)
      GOTO 9999
  995 WRITE(6,*) '  ERROR on subroutine SOS_GRANU'
      CALL EXIT(1)
      GOTO 9999 
  996 WRITE(6,*) '  ERROR on subroutine SOS_DECOMPO_LEGENDRE'
      CALL EXIT(1)
      GOTO 9999 
  997 WRITE(6,*) '  ERROR on subroutine SOS_INIT_PARAMSF'
      CALL EXIT(1)
      GOTO 9999                  
 1006 WRITE(6,*) '  ERROR on -AER.Model value : out of available range.'
      CALL EXIT(1)
      GOTO 9999     
 1007 WRITE(6,*) '  ERROR on -AER.AOT value : wrong negative value.'
      CALL EXIT(1)
      GOTO 9999                    
 1008 WRITE(6,*) '  ERROR on -AER.Wa value : '
      WRITE(6,*) '  --> wavelenth has to be defined in ]0.2, 4.0[ mic.'
      CALL EXIT(1)
      GOTO 9999
 1009 WRITE(6,*) '  ERROR : SOS_MIE_ALPHAMIN > AER.MMD.Mie.AlphaMax'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
      CALL EXIT(1) 
      GOTO 9999
 1010 WRITE(6,*) '  ERROR : SOS_MIE_ALPHAMIN F6.4'
      WRITE(6,*) '          AER.MMD.Mie.AlphaMax F9.2'
      WRITE(6,*) '  --> check SOS.h and argument value & compile again'
      CALL EXIT(1)
      GOTO 9999  
10100 WRITE(6,*) '  ERROR : --> imaginary part of refractive'
      WRITE(6,*) '              indexes have to be negative'
      CALL EXIT(1)
      GOTO 9999          
 1011 WRITE(6,*) '  ERROR on -AER.MMD.SDtype value  : '
      WRITE(6,*) '  --> out of available range for mono-modal model.'
      CALL EXIT(1)
      GOTO 9999 
 1012 WRITE(6,*) '  ERROR on -AER.WMO.Model or -AER.SF.Model value  : '
      WRITE(6,*) '  --> out of available range for WMO '
      WRITE(6,*) '      or Shettle&Fenn model option.'     
      CALL EXIT(1)
      GOTO 9999
 1013 WRITE(6,*) '  ERROR on -AER.SF.RH : Relative Humidity (0 to 99%)'
      CALL EXIT(1)
      GOTO 9999 
 1014 WRITE(6,*) '  ERROR on -AER.Tronca : Troncature option  (0 or 1)'
      CALL EXIT(1)
      GOTO 9999       
 1015 WRITE(6,*) '  ERROR on -AER.BMD.VDdef : Option for volumetric '
      WRITE(6,*) '  concentration definition for bimodal LND (1 or 2)'
      CALL EXIT(1)
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
  
 9999 STOP
      END   	!FIN DU PROGRAMME SOS_AEROSOLS

 
 
 
C==============================================================================
C PROCEDURE: SOS_NOM_FICMIE
C ==========
C      Cette procedure genere le nom du fichier de Mie en fonction de :
C      	- l'indice de refraction : partie reelle RN
C	- l'indice de refraction : partie imaginaire IN
C	- le parametre de taille minimal ALPHAO
C	- le parametre de taille maximal ALPHAF
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
C      On convertit RN en une chaine de caracteres CRN(1:4) 
C      Pour cela on passe du Format F5.3 * 1000 => A4
C           Le premier caractere CRN(1) est le chiffre unite de RN,
C           Le second caractere CRN(2) est la premiere decimale de RN,
C           etc.
C
C
C Description des parametres entree / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*SOS_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	RN (double)	: (E) partie reelle de l'indice de refraction
C	IN (double)	: (E) partie imaginaire de l'indice de refraction
C       ALPHAO (double) : (E) parametre de taille minimal
C       ALPHAF (double) : (E) parametre de taille maximal
C	FICS (CHARACTER*SOS_LENFIC1) : (S) nom du fichier de Mie (sans arborescence)	
C
C Constante utilise:
C ------------------
C    SOS_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
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
         
      DOUBLE PRECISION RN	! partie reelle de l'indice de refraction.
      DOUBLE PRECISION IN	! partie imaginaire de l'indice de refraction.
      DOUBLE PRECISION ALPHAO   ! valeur de alpha a partir de laquelle 
      				! on commence le calcul de Mie.				
      DOUBLE PRECISION ALPHAF   ! valeur de alpha pour laquelle on arrete le calcul.

      CHARACTER*SOS_LENFIC2 FICANGLES_USER  ! chemin et nom du fichier des angles "utilisateurs"				
      CHARACTER*SOS_LENFIC1 FICS	    ! nom du fichier MIE 	
          
      CHARACTER*4 CRN		! Chaines de caracteres
      CHARACTER*6 CIN 
      CHARACTER*5 CALPHAO   
      CHARACTER*7 CALPHAF
      CHARACTER*4 CNBMU
      
      INTEGER*2 I
      INTEGER*4 IDEB, IDEB_CH, IFIN, IFIN_CH  ! position de debut et fin d'une chaine de caracteres  
      INTEGER*4 NBMU_GAUSS  ! Nombre d'angles de Gauss utilises 
      			      	 
	    
C --------- Conversion en chaines de caracteres des parametres d'entree	   
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
      
C -------- Construction de la chaine de caracteres pour l'indice RN	    
      FICS(4:4) = CRN(1:1)     ! Chiffre des unités
      FICS(5:5) = '.'	       ! Virgule
      FICS(6:8) = CRN(2:4)     ! Trois premieres decimales   
      FICS(9:9) = '-'
      
C -------- Construction de la chaine de caracteres pour l'indice IN	    
      FICS(10:10) = CIN(1:1)     ! Chiffre des unités
      FICS(11:11) = '.'	         ! Virgule
      FICS(12:16) = CIN(2:6)     ! Cinq premieres decimales   	    
      FICS(17:17) = '-'	
      
C -------- Construction de la chaine de caracteres pour AlphaO	    
      FICS(18:18) = CALPHAO(1:1)     ! Chiffre des unités
      FICS(19:19) = '.'	             ! Virgule
      FICS(20:23) = CALPHAO(2:5)     ! Quatre premieres decimales   	    
      FICS(24:24) = '-'	

C -------- Construction de la chaine de caracteres pour AlphaF	    
      FICS(25:29) = CALPHAF(1:5)     ! Chiffres des unités
      FICS(30:30) = '.'	             ! Virgule
      FICS(31:32) = CALPHAF(6:7)     ! Deux premieres decimales   	    
      FICS(33:33) = '-'	            

C -------- Construction de la chaine de caracteres pour MUgauss
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
  
       
C -------- Construction de la chaine de caracteres pour NomFicUser
C -------- Structuration pour avoir NBMU_GAUSS a gauche (sans "blancs")	            


      IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN

C       Determination de la position du nom de fichier sans son arborescence          
	  IFIN_CH=INDEX(FICANGLES_USER,' ')
          IFIN_CH=IFIN_CH-1
	  IF (IFIN_CH.LE.0) IFIN_CH=SOS_LENDIR
	  
	  IDEB_CH=0
	  	  
10	  IDEB=INDEX(FICANGLES_USER(IDEB_CH:IFIN_CH),'/')
	  IDEB_CH=IDEB_CH+IDEB
          IF (IDEB.NE.0) GOTO 10	  

C       Ajout du nom de fichier d'angles "utilisateur" au nom du fichier de MIE    
          IFIN=INDEX(FICS,' ')  
	  FICS(IFIN:IFIN) = '-'
	  
	  IFIN=IFIN+1
	  FICS(IFIN:IFIN+IFIN_CH-IDEB_CH+1) 
     &      = FICANGLES_USER(IDEB_CH:IFIN_CH)
	  
      ENDIF 
     	
	
      END !SOS_NOM_FICMIE
 
 
 
 
C==============================================================================
C PROCEDURE: SOS_INIT_PARAMWMO
C ==========
C      Cette procedure initialise l'ensemble des parametres de chaque
C      constituant de la WMO pour la longueur d'onde de simulation :
C      	- le rayon modal
C	- la variance (log10)
C	- concentration volumique
C	- partie reelle de l'indice de refraction
C	- partie imaginaire
C
C      Interpolation lineaire des indices pour la longueur d'onde d'entree.
C
C      Force les indices resultats au format F5.3 pour la partie reelle
C      et F8.5 pour la partie imaginaire (valeur negative).
C
C Description des parametres entree / sortie
C ------------------------------------------
C	FICDATAWMO (CHARACTER*SOS_LENFIC2)	: (E) chemin et nom du fichier Data_WMO
C	WA (double)	: (E) longueur d'onde selectionnee
C       VARGRANU1I(SOS_NB_WMO) (double) : (S) rayon modal de chaque constituant
C       VARGRANU2I(SOS_NB_WMO) (double) : (S) log decimal de la variance de 
C					  chaque constituant
C	MR(SOS_NB_WMO) (double)	: (S) partie reelle de l'indice de refraction
C       MI(SOS_NB_WMO) (double) : (S) partie imaginaire de l'indice de refraction
C	V(SOS_NB_WMO) (double)	: (S) concentration volumique micron3/cm3
C	IER (I4)	: (S) code d'erreur = 0 si pas d'erreur, sinon =1	
C
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier des donnees de la WMO
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
C                            longueur d'onde, partie reelle et imaginaire
C			     de l'indice de refraction de chaque constituant.   
C
C
C Constantes utilisees :
C --------------------
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------    
C     -Ouverture d'un fichier
C       Erreur si le fichier des donnees de la WMO n'existe pas
C     -Lecture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel a la routine: SOS_INTERPOL
C
C
C  MOD:VERSION 2.0: 18/12/2015
C      Conversion de la variance du fichier de données WMO (en log décimal)
C      en variance conforme à la définition classique des LND: VAR_WMO * LOG(10)
C
C==============================================================================
      SUBROUTINE SOS_INIT_PARAMWMO(FICDATAWMO,WA,VARGRANU1I,
     &				   VARGRANU2I,MR,MI,V,IER)

      IMPLICIT NONE
      
C Constantes AEROSOLS
C -------------------

      INTEGER SOS_NB_WMO	!Nombre de types d'aerosols de la WMO utilises: DL,WS,OC,SO
      PARAMETER (SOS_NB_WMO=4)

      INTEGER SOS_NB_AERO	!Dimension des tableaux resultats 
      PARAMETER (SOS_NB_AERO=5)	!(4 composants WMO, 5 composants Shettle & Fenn)	
      
            
C* Definition des variables       
C*-----------------------------------------

      CHARACTER*SOS_LENFIC2 FICDATAWMO	! chemin et nom du fichier Data_WMO
      
      DOUBLE PRECISION WA	! Longueur d'onde en microns selectionnee.	
      DOUBLE PRECISION WA1	! Longueur d'onde en microns lue dans le
      DOUBLE PRECISION WA2	! fichier Data_WMO.
      
      DOUBLE PRECISION VARGRANU1I(SOS_NB_AERO)   ! Parametre 1 de la granulometrie
      						 ! pour chaque modele de constituant.
      DOUBLE PRECISION VARGRANU2I(SOS_NB_AERO)   ! Parametre 2 de la granulometrie
      						 ! pour chaque modele de constituant.
      
      DOUBLE PRECISION MR(SOS_NB_AERO)   ! partie reelle de l'indice de refraction
      					 ! pour chaque modele d'aerosols.
      DOUBLE PRECISION MI(SOS_NB_AERO)	! partie imaginaire de l'indice de refraction.
      					! pour chaque modele d'aerosols.
      DOUBLE PRECISION MR1(SOS_NB_WMO)	!
      DOUBLE PRECISION MI1(SOS_NB_WMO)  ! utilise pour le stockage de 
      DOUBLE PRECISION MR2(SOS_NB_WMO)  ! valeurs intermediaires
      DOUBLE PRECISION MI2(SOS_NB_WMO)  !
      
      DOUBLE PRECISION V(SOS_NB_WMO)	!Concentration volumique des constituants
      
      INTEGER*4 I	   !Indice.       
      INTEGER*4 IER        ! code d'erreur =0 si pas d'erreur, =-1 sinon
      
C* Initialisation
C  --------------
       IER=0	 
       
       DO I=1,SOS_NB_WMO
          VARGRANU1I(I)=0.
          VARGRANU2I(I)=0.
	  MR(I)=0.
    	  MI(I)=0.
	  V(I)=0.
       ENDDO

C* Ouverture du fichier contenant les donnees de la WMO
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

C* Conversion en variance pour LND compatible avec la définition de OSOAA
C  ------------------------------------- 	 
       DO I=1,4
          VARGRANU2I(I) = VARGRANU2I(I) * DLOG(10.D+00)
       ENDDO

C* Lecture de la concentration volumique
C  ------------------------------------- 	 
       READ(8,444,err=992) V(1),V(2),V(3),V(4) 

C* Determination des indices de refraction par interpolation lineaire
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
       
	  DO I=1,SOS_NB_WMO	  
	     CALL SOS_INTERPOL(MR1(I),MR2(I),WA1,WA2,WA,MR(I))
	     CALL SOS_INTERPOL(MI1(I),MI2(I),WA1,WA2,WA,MI(I))
	     
C******** Conversion aux formats F5.3 et F8.3	   
	     MR(I) =  DNINT( MR(I)*1000.D+00) / 1000.D+00
	     MI(I) = -DNINT(-MI(I)*100000.D+00) / 100000.D+00	    
	  ENDDO
	  GOTO 881
	       
       ELSE
       	 
          WA1=WA2	! Reinitialisation: ligne1 recoit ligne2
	  DO I=1,SOS_NB_WMO
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
C      Cette procedure initialise l'ensemble des parametres de chaque
C      constituant des modeles de Shettle et Fenn pour la longueur d'onde 
C      et l'humidite relative selectionnees :
C      	- le rayon modal
C	- la variance
C	- partie reelle de l'indice de refraction
C	- partie imaginaire
C
C Description des parametres
C --------------------------
C	DIRE (CHARACTER*SOS_LENDIR)  : (E) chemin d'acces au fichier de donnees
C	WA (double)	      : (E) longueur d'onde selectionne
C	RH (double)	      : (E) humidite relative
C       VARGRANU1I(SOS_NB_AERO) (double) : (S) rayon modal de chaque constituant
C       VARGRANU2I(SOS_NB_AERO) (double) : (S) log decimal de la variance de 
C					  chaque constituant
C	MR(SOS_NB_AERO) (double) : (S) partie reelle de l'indice de refraction
C       MI(SOS_NB_AERO) (double) : (S) partie imaginaire de l'indice de refraction
C	IER (I4) : (S) code d'erreur = 0 si pas d'erreur, sinon =1	
C
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier des parametres de granulometrie des composants
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
C      pour chaque ligne: longeur d'onde, partie reelle et partie imaginaire 
C			  de l'indice de refraction pour 8 humidites 	
C			  relatives : 0%, 50%, 70%, 80%, 90%, 95%, 98%, 99%.
C
C
C Constantes utilisees :
C --------------------
C    SOS_LENDIR : Longueur des noms de repertoires.
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_AER_DATASF  : Nom du fichier contenant rayon modal et log variance S&F.
C    SOS_AER_XX_SF   : Nom des fichiers contenant les indices de refraction S&F
C		       pour chaque constituant XX = SR, LR, SU, LU, OM.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------    
C     -Ouverture d'un fichier
C       Erreur si l'un des fichiers de donnees n'existe pas
C     -Lecture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel a la routine: SOS_INTERPOL
C
C
C  MOD:VERSION 2.0: 18/12/2015
C      Conversion de la variance du fichier de données S&F (en log décimal)
C      en variance conforme à la définition classique des LND: VAR_WMO * LOG(10)
C
C==============================================================================
      SUBROUTINE SOS_INIT_PARAMSF(DIRE,WA,RH,VARGRANU1I,
     &                         	  VARGRANU2I,MR,MI,IER)

      IMPLICIT NONE
      
C Constantes AEROSOLS
C -------------------

      INTEGER SOS_NB_SF	        !Nombre des composants des modeles Shettle & Fenn
      PARAMETER (SOS_NB_SF=5)

      INTEGER SOS_NB_RH		!Nombre de valeurs d'humidite relative
      PARAMETER (SOS_NB_RH=8)
 
      INTEGER SOS_NB_AERO	!Dimension des tableaux resultats 
      PARAMETER (SOS_NB_AERO=5)	!(4 composants WMO, 5 composants Shettle & Fenn)	
      
            
C* Definition des variables       
C*-----------------------------------------

      CHARACTER*SOS_LENDIR DIRE	! Chemin d'acces aux fichiers S&F
      CHARACTER*SOS_LENFIC2 FICDATASF	! Nom du fichier Data_SF
      CHARACTER*SOS_LENFIC2 FIC_SR_SF 	! Fichier IRefrac_SR
      CHARACTER*SOS_LENFIC2 FIC_LR_SF 	! Fichier IRefrac_LR
      CHARACTER*SOS_LENFIC2 FIC_SU_SF 	! Fichier IRefrac_SU
      CHARACTER*SOS_LENFIC2 FIC_LU_SF 	! Fichier IRefrac_LU
      CHARACTER*SOS_LENFIC2 FIC_OM_SF 	! Fichier IRefrac_OM
      
      DOUBLE PRECISION WA	! Longueur d'onde en microns selectionnee.	
      DOUBLE PRECISION WA1	! Longueur d'onde en microns lue dans les fichiers
      DOUBLE PRECISION WA2	
      
      DOUBLE PRECISION RH	! Humidite relative (%)
      DOUBLE PRECISION RH1	! Humidite relative lue dans les fichiers
      DOUBLE PRECISION RH2	!
      
      DOUBLE PRECISION VARGRANU1I(SOS_NB_AERO)  ! Parametre 1 de la granulometrie
      						! pour chaque modele de constituant.
      DOUBLE PRECISION VARGRANU2I(SOS_NB_AERO)  ! Parametre 2 de la granulometrie
      						! pour chaque modele de constituant.
      
      DOUBLE PRECISION MR(SOS_NB_AERO)  ! partie reelle de l'indice de refraction
      					! pour chaque modele d'aerosols.
      DOUBLE PRECISION MI(SOS_NB_AERO)	! partie imaginaire de l'indice de refraction.
      					! pour chaque modele d'aerosols.
      
      DOUBLE PRECISION MR_CPT1		!
      DOUBLE PRECISION MR_CPT2		! Valeurs intermediaires
      DOUBLE PRECISION MI_CPT1		!
      DOUBLE PRECISION MI_CPT2		!
      
      DOUBLE PRECISION RM1(SOS_NB_SF)	!
      DOUBLE PRECISION RM2(SOS_NB_SF)	! utilise pour le stockage de
      DOUBLE PRECISION MR1(SOS_NB_RH)   ! valeurs intermediaires 
      DOUBLE PRECISION MR2(SOS_NB_RH)   ! 
      DOUBLE PRECISION MI1(SOS_NB_RH)   ! 
      DOUBLE PRECISION MI2(SOS_NB_RH)   !
	
      INTEGER*4 I,H	   ! Indice
      INTEGER*4 CPT	   ! Compteur 
      INTEGER*4 IFIN	   ! Indicateur de fin d'une chaine de caracteres
      
      INTEGER*4 IER        ! code d'erreur =0 si pas d'erreur, =-1 sinon
      
C* Initialisation
C  --------------
       IER=0
       
       DO I=1,SOS_NB_SF
          VARGRANU1I(I)=0.
          VARGRANU2I(I)=0.	  
	  MR(I)=0.
    	  MI(I)=0. 
       ENDDO

       
       ! Initialisation des noms des differents fichiers contenant 
       ! les donnees du modele S&F
       
       IFIN=INDEX(DIRE,' ')
       IFIN=IFIN-1
       IF (IFIN.LE.0) IFIN=SOS_LENDIR
      
       FICDATASF = DIRE(1:IFIN)//'/'//SOS_AER_DATASF
       FIC_SR_SF = DIRE(1:IFIN)//'/'//SOS_AER_SR_SF
       FIC_LR_SF = DIRE(1:IFIN)//'/'//SOS_AER_LR_SF
       FIC_SU_SF = DIRE(1:IFIN)//'/'//SOS_AER_SU_SF
       FIC_LU_SF = DIRE(1:IFIN)//'/'//SOS_AER_LU_SF
       FIC_OM_SF = DIRE(1:IFIN)//'/'//SOS_AER_OM_SF

C* Ouverture du fichier contenant les parametres de granulometrie
C* des composants de S&F
C  --------------------------------------------------------------
       OPEN(9,FILE=FICDATASF,STATUS='OLD',ERR=991)
   	  
C* Lecture du log decimal de la variance
C  ------------------------------------- 	 
       READ(9,222,err=993) VARGRANU2I(1),VARGRANU2I(2),VARGRANU2I(3),
     &                     VARGRANU2I(4),VARGRANU2I(5) 

C* Conversion en variance pour LND compatible avec la définition de OSOAA
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
	
	   DO I=1,SOS_NB_SF
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
	
	! Cas ou RH1<RH<=RH2, interpolation des indices de refraction
	IF((RH1.LT.RH).AND.(RH.LE.RH2)) THEN
	
	   DO I=1,SOS_NB_SF	  
	      CALL SOS_INTERPOL(RM1(I),RM2(I),RH1,RH2,RH,VARGRANU1I(I))
     	   ENDDO
	   GOTO 884
	       
        ELSE	! Reinitialisation des valeurs
       
          RH1=RH2
	  DO I=1,SOS_NB_SF
	     RM1(I)=RM2(I)
     	  ENDDO
	  GOTO 883	! passage a la ligne suivante dans le fichier	       
	  
       ENDIF
       
  884  CONTINUE


C* Ouverture des fichiers contenant l'indice de refraction pour 
C  chaque constituant des modeles S&F
C  -----------------------------------------------------------
       OPEN(11,FILE=FIC_SR_SF,STATUS='OLD',ERR=992)
       OPEN(12,FILE=FIC_LR_SF,STATUS='OLD',ERR=992)
       OPEN(13,FILE=FIC_SU_SF,STATUS='OLD',ERR=992)
       OPEN(14,FILE=FIC_LU_SF,STATUS='OLD',ERR=992)
       OPEN(15,FILE=FIC_OM_SF,STATUS='OLD',ERR=992)
  
       DO I=1,SOS_NB_SF

C* Determination des indices de refraction par interpolation lineaire
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
     
          ! Cas ou WA1<=WA<=WA2, interpolation des donnees
	  IF((WA1.LE.WA).AND.(WA.LE.WA2)) THEN
	  
	     ! Cas ou RH = 0%, interpolation de WA mais pas de RH
	     IF(CPT.EQ.1) THEN
	        CALL SOS_INTERPOL(MR1(CPT),MR2(CPT),WA1,WA2,WA,MR(I))
		CALL SOS_INTERPOL(MI1(CPT),MI2(CPT),WA1,WA2,WA,MI(I))
	     
	     ELSE
	      ! Interpolations pour la partie reelle de l'indice de refraction	     
		! Interpolation de WA pour RH1(CPT-1) et RH2(CPT) encadrant RH
		CALL SOS_INTERPOL(MR1(CPT-1),MR2(CPT-1),
     &                            WA1,WA2,WA,MR_CPT1)
	        CALL SOS_INTERPOL(MR1(CPT),MR2(CPT),WA1,WA2,WA,MR_CPT2)
		! Interpolation de RH
	        CALL SOS_INTERPOL(MR_CPT1,MR_CPT2,RH1,RH2,RH,MR(I))
	      
	      ! Interpolations pour la partie imaginaire
	        CALL SOS_INTERPOL(MI1(CPT-1),MI2(CPT-1),
     &                            WA1,WA2,WA,MI_CPT1)
	        CALL SOS_INTERPOL(MI1(CPT),MI2(CPT),WA1,WA2,WA,MI_CPT2)
	        CALL SOS_INTERPOL(MI_CPT1,MI_CPT2,RH1,RH2,RH,MI(I))
	     ENDIF
		
          ELSE
       
             WA1=WA2	! Reinitilisation 
	     DO H=1,SOS_NB_RH
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
C      Cette procedure permet d'interpoler lineairement en X une fonction f
C      dont on connait les valeurs (X1,Y1) et (XOUT,Y2).
C
C Description des parametres
C --------------------------
C	X1 (double) : (E)
C	X2 (double) : (E)
C	Y1 (double) : (E) Valeur f(X1)
C	Y2 (double) : (E) Valeur f(X2)
C	X (double)  : (E) Abscisse pour laquelle on interpole la valeur f(X)
C       Y (double)  : (S) Valeur f(X) interpolee
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
      
C* Calcul du parametre interpole
C  -----------------------------
      Y = ((Y2-Y1)/(X2-X1))*(X-X2) + Y2

      END	!FIN DE LA PROCEDURE SOS_INTERPOL


            
C==============================================================================
C PROCEDURE: SOS_DECOMPO_LEGENDRE
C ==========
C       Cette procedure calcule les parametres alpha, beta, gamma delta et zeta
C	de la decomposition de la focntion de phase en polynomes de Legendre.
C	Application eventuelle de la troncature.
C
C       Les fonctions de phase sont stockees dans des tableaux definis jusqu'a
C       SOS_MIE_NBMU_MAX mais sont utilisees jusqu'a MIE_NBMU (angles effectifs).
C
C Description des parametres
C --------------------------
C	ITRONC (I4) : (E/S) Option de troncature (1 pour la realiser).
C                           Mise a 0 en sortie si la troncature est trop faible.
C	TRACE (logical) : (E)  =vrai si ecriture dans le fichier trace
C       MIE_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C	XMU(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX) (double) : (E) Cosinus des angles de Gauss.
C	XHR(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX) (double) : (E) Poids des angles de Gauss.
C	P11(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX) (double)  : (E) Fonction de phase calculee 
C					         a partir des fonctions la theorie de MIE.
C	P12(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX) (double)  : (E) Element P12 de la matrice de phase.
C	P33(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	(double) : (E) Element P33 de la matrice de phase.
C	TTT(-SOS_MIE_NBMU_MAX:SOS_MIE_NBM_MAXU) (double) : (S) Copie de la fonction de phase (P11) 
C					 	               avant troncature.
C	A (double)		     : (S) Coefficient de troncature.
C	Z1 (double)		     : (S) Valeur de BETA11(0) avant normalisation.
C	ALP(0:SOS_OS_NB_MAX) (double)     : (S) Parametres alpha(k).
C	BETA11(0:SOS_OS_NB_MAX) (double)  : (S) Parametres beta11(k).
C	BETA22(0:SOS_OS_NB_MAX) (double)  : (S) Parametres beta22(k).
C	GAMMA12(0:SOS_OS_NB_MAX) (double) : (S) Parametres gamma12(k).
C	ZETA(0:SOS_OS_NB_MAX) (double)    : (S) Parametres zeta(k).
C	DELTA33(0:SOS_OS_NB_MAX) (double)  : (S) Parametres delta33(k).
C	IER (I4)		     : (S) Code d'erreur =0 si pas d'erreur, sinon =1.	
C
C 
C Description des fichiers utilises
C ---------------------------------
C     Aucun
C
C Constantes utilisees :
C --------------------
C    SOS_MIE_NBMU_MAX : Nombre maximal d'angles positifs.
C    SOS_OS_NB_MAX :  Valeur maximale pour l'ordre limite du developpement en 
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
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C
C==============================================================================
      SUBROUTINE SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,
     &				      OS_NB,P11,TTT,P12,P22,P33,A,Z1,
     &				      ALP,BETA11,BETA22,
     &                                GAMMA12,DELTA33,ZETA,IER)

      IMPLICIT NONE

C* Definition des variables       
C*-----------------------------------------
      DOUBLE PRECISION XMU(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)  
                                                  !Cosinus des angles de Gauss.
      DOUBLE PRECISION XHR(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)  
                                                  !Poids des angles de Gauss.
      
      DOUBLE PRECISION P11(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Fonction de phase calculee a partir des fonctions de 
		       !la theorie de MIE.
      DOUBLE PRECISION P12(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P12 de la matrice de phase.
      DOUBLE PRECISION P22(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P22 de la matrice de phase.
      DOUBLE PRECISION P33(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P33 de la matrice de phase.

      DOUBLE PRECISION TTT(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Copie de la fonction de phase (P11) avant troncature.
		       
		
      DOUBLE PRECISION PL(-1:SOS_OS_NB_MAX+1) 
                       !Polynomes de Legendre a l'ordre L calcule pour un angle.
      DOUBLE PRECISION POL(0:SOS_OS_NB_MAX+1) 
                       !Fonction de Legendre d'ordre S=2 et d'ordre L variable
		       !(calcule pour un angle).
		       			
      DOUBLE PRECISION D	!Variable pour le calcul de POL.
      DOUBLE PRECISION E	!Variable pour le calcul de POL.
      		       
      DOUBLE PRECISION X
      DOUBLE PRECISION XRMU	!Valeur utilisee du tableau des angles de Gauss.
				
      DOUBLE PRECISION XX	
      DOUBLE PRECISION XXX

      DOUBLE PRECISION A	!Coefficient de troncature.    
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

      DOUBLE PRECISION ALP(0:SOS_OS_NB_MAX)	! Parametres alpha(k).
      DOUBLE PRECISION BETA11(0:SOS_OS_NB_MAX)	! Parametres beta11(k).
      DOUBLE PRECISION BETA22(0:SOS_OS_NB_MAX)	! Parametres beta22(k).
      DOUBLE PRECISION GAMMA12(0:SOS_OS_NB_MAX)	! Parametres gamma12(k).
      DOUBLE PRECISION ZETA(0:SOS_OS_NB_MAX)	! Parametres zeta(k).
      DOUBLE PRECISION DELTA33(0:SOS_OS_NB_MAX)	! Parametres delta33(k).
                        

      INTEGER*4 MIE_NBMU   ! Nombre d'angles (positifs) effectivement utiles 
      INTEGER*4 ITRONC	   ! Option de troncature (1 pour la realiser).
      INTEGER*4 OS_NB 	   ! Ordre des developpements des fonctions de phase 
			   ! en Polynomes de Legendre.      
      INTEGER*4 K          ! Indice de tableau.
      INTEGER*4 KK	   ! Indice du tableau des angles de Gauss pour le mu
      			   ! le plus proche de 0.94 par valeur superieure.			               
      INTEGER*4	NN
      INTEGER*4	MM				
      INTEGER*4 I	   ! Indice.  
      INTEGER*4 J	   ! Indice.      
      INTEGER*4 IER	   ! Code d'erreur =0 si pas d'erreur, =1 sinon.         

      LOGICAL TESS	   !Variable pour le calcul de Troncature. 
      LOGICAL TRACE	   !=vrai si ecriture dans le fichier trace 
      
      
C* Copie de la fonction de phase en intensite avant la troncature.
C-----------------------------------------------------------------
      DO 55 J=-MIE_NBMU,MIE_NBMU
         TTT(J)=P11(J)
  55  CONTINUE 
      
C==============
C* TRONCATURE
C============== 

      IF(ITRONC.EQ.0) GO TO 59

C* Troncature de la fonction de phase en intensite pour les angles 
C* inferieurs a 19.95 degres (mu = 0.94).
C------------------------------------------------------------------------------

C* Determination des angles de Gauss juste inferieurs a 36.87 degres (mu=0.8)
C* et a 19.95 degres (mu = 0.94)
C--------------------------------------------------------------------
      TESS=.TRUE.
      DO 56 J=1,MIE_NBMU
         IF (.NOT.(TESS.AND.(XMU(J).GT.0.8))) GO TO 56
         K=J-1
         TESS=.FALSE.
  56  CONTINUE
  
      TESS=.TRUE.
      DO 57 J=1,MIE_NBMU
         IF (.NOT.(TESS.AND.(XMU(J).GT.0.94))) GO TO 57
         KK=J-1
         TESS=.FALSE.
  57  CONTINUE

C* Determination de la pente de la droite "log10(p(cos(teta))) = AA * teta + BB"
C* passant par les valeurs de log10 (p) pour les mu de Gauss ci-dessus.
C--------------------------------------------------------------------------
      AA=(DLOG10(P11(KK))-DLOG10(P11(K)))/(DACOS(XMU(KK))-DACOS(XMU(K)))

C* Application de la troncature
C-------------------------------
      X1=DLOG10(P11(KK))
      X2=DACOS(XMU(KK))
      
      DO 58 J=KK+1,MIE_NBMU
         IF(DABS(XMU(J)-1.).LE.1D-08) A=X1-AA*X2
         A=X1+AA*(DACOS(XMU(J))-X2)
         P11(J)=10**A
  58  CONTINUE

  59  CONTINUE


C* Calcul des parametres BETA11(k) de la decomposition de la fonction
C* de phase (eventuellement tronquee) en polynomes de Legendre.
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
         A=2*(1-BETA11(0))
      ELSE
         A=0.0
      ENDIF

C* Si le coefficient de troncature est inferieur au seuil Seuil_Tronca
C* la troncature est annulee.   
C----------------------------------------------------------------------	
      IF((ITRONC.EQ.1).AND.(A.LT.SOS_PH_SEUIL_TRONCA)) THEN
      
         IF (TRACE) THEN
            WRITE(88,*,err=921) 'Coefficient de troncature A = ',A,
     &            ' Valeur trop faible : SUPPRESSION de la TRONCATURE'
   
	 ENDIF
	 
C* Recopie de la fonction de phase en intensite avant la troncature
C* et re-initialisation des parametres de la decomposition de la matrice
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
  

C* Calcul des parametres GAMMA12(k), BETA22(K) et DELTA33(k) de la decomposition des 
C* fonctions de phase polarisees P12, P22 et P33 en polynomes ou fonctions de Legendre.
C* Application de l'eventuelle troncature par le rapport de la fonction de phase
C* en intensite tronquee sur celle initiale: P11(J)/TTT(J).
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


C* Calcul des parametres ALP(k) et ZETA(k) en fonctions des BETA22(k) et DELTA33(k).
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

C* Normalisation des parametres ALP(K), BETA(K), GAMMA(K), DELTA(K) et ZETA(K)
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
C      Cette procedure calcule la fonction de phase P(mu) et les fonctions
C      de phase polarisee Q(mu) et U(mu) a partir des calculs de MIE realises
C      par le programme "SOS_MIE", pour une longueur d'onde donnee et pour
C      une distribution en taille des particules donnee (renseignee par igranu).
C
C      Si igranu = 1 => Distribution Log-Normale : 
C                N(r) = exp(-0.5 *(ln(r/rm)/sig)**2) / r
C		 varGranu1 = rm : rayon modal en microns,
C		 varGranu1 = sig : variance.
C
C      Si igranu = 2 => Loi de Junge : N(r) = Nr0 = R0**(-xnu) pour r < R0
C				       N(r) = r**(-xnu) pour r > R0
C		 varGranu1 = R0 : rayon minimal en microns,
C		 varGranu1 = xnu : puissance appliquee au rayon.
C
C
C       Les fonctions de phase sont stockees dans des tableaux definis jusqu'a
C       SOS_MIE_NBMU_MAX mais sont utilisees jusqu'a MIE_NBMU (angles effectifs).
C
C
C Description des parametres
C ----------------------------
C       FICMIE(CHARACTER*SOS_LENFIC2)  (E) : Nom du fichier MIE 
C                                 (repertoire + nom fichier + extension).
C       IGRANU  (I4)         (E) : Numero du modele de granulometrie a utiliser.
C       VARGRANU1 (double)   (E) : Parametre 1 de la granulometrie.
C	VARGRANU2 (double)   (E) : Parametre 2 de la granulometrie
C	WA  (double)         (E) : Longueur d'onde en microns.
C       MIE_NBMU (I4) 	     (E) : Nombre d'angles (positifs) effectivement utiles
C	XMU(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX) (double) : (E) Cosinus des angles de Gauss.
C	TRACE (logical)      (E) : = vrai,si ecriture dans le fichier trace
C                                 (No d'unite logique du fichier trace =88).
C	KMAT1 (double)       (S) : Section efficace d'extinction (micron^2).
C	KMAT2 (double)       (S) : Section efficace de diffusion (micron^2).
C	SOMME_NR (double)    (S) : Nombre de particules integre sur la granulometrie.
C       P11 (double)          (S) : Fonction de phase calculee a partir des 
C                                  fonctions de la theorie de MIE. 
C       P12 (double)          (S) : Element P12 de la matrice de phase.
C	P33 (double)          (S) : Element P33 de la matrice de phase.
C       IER  (I4)            (S) : Code d'erreur =0 si pas d'erreur, =-1 sinon.   
C
C
C Autre resultat fourni par la procedure
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
C                Indice de refraction (partie reelle)
C                Indice de refraction (partie imaginaire)
C                Parametre de taille limite
C                Nombre d'angles de Gauss positifs
C
C      1 enregistrement pour chaque valeur du parametre de taille avec:
C        - Parametre de taille (reel),
C        - Coefficient d'efficacite en extinction (reel),
C        - Coefficient d'efficacite en diffusion (reel),
C        - Facteur d'anisotropie (double precision),
C        - Tableau de la fonction de phase Imie(mu) (reel),
C        - Tableau de la fonction de phase polarisee Qmie(mu) (reel),
C        - Tableau de la fonction de phase polarisee Umie(mu) (reel).
C
C    --> Les tableaux ont une taille utile limitee a MIE_NBMU : 
C        ils sont lus avec cette limitation.
C
C
C Constantes utilisees :  
C----------------------
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_PI : valeur de PI
C    SOS_MIE_NBMU_MAX : Nombre maximal d'angles positifs.
C    SOS_OS_NB_MAX :  Valeur maximale pour l'ordre limite du developpement en 
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
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
C     retour du status -1 au programme appelant
C        
C   
C
C
C* MOD:VERSION 2.0: 24/06/2008 
C	   - Les coefficients KMAT1 et KMAT2 sont normalises à une particule
C            ==> La routine sort des sections efficaces d'extinction et de diffusion
C                (en mic^2) au lieu de coefficient d'extinction et de diffusion (/mic).
C          - Introduction du calcul de rayon effectif et variance effective.
C
C MOD:VERSION 3.0: 10/01/2010
C          - Passage d'une dimension des tableaux de SOS_MIE_NBMU a SOS_MIE_NBMU_MAX
C          - Boucle sur les fonctions de phase jusqu'a MIE_NBMU
C          - Modification de la ligne de lecture  des tableaux de fonctions de phase 
C            par valeur de parametre de taille Alpha
C            ==> limitation aux bornes utiles -MIE_NBMU:MIE_NBMU de l'enregistrement effectif
C                des fichiers de MIE
C
C MOD:VERSION 4.0: 18/12/2015
C          - La formulation N(r) = exp(-0.5 *(log10(r/rm)/sig)**2) / r
C            pour sig = log10(variance) est abandonnée en faveur 
C            d'une formulation classique du modèle LND:
C			    N(r) = exp(-0.5 *(log(r/rm)/sig)**2) / r	
C	     pour sig = variance.
C
C==============================================================================


      SUBROUTINE SOS_GRANU(FICMIE,IGRANU,VARGRANU1,
     &                     VARGRANU2,WA,MIE_NBMU,XMU,TRACE,
     &                     KMAT1,KMAT2,SOMME_NR,P11,P12,P33,IER)

      IMPLICIT NONE


C* Definition des variables
C----------------------------------------------------------------------------
      
      CHARACTER*SOS_LENFIC2 FICMIE	!Chemin et nom du fichier MIE
      
      DOUBLE PRECISION VARGRANU1   	!Parametre 1 de la granulometrie.
      DOUBLE PRECISION VARGRANU2   	!Parametre 2 de la granulometrie.
      DOUBLE PRECISION RM		!Rayon modal d'une granulometrie LND.
      DOUBLE PRECISION SIG		!Variance d'une granulometrie LND
      DOUBLE PRECISION B		!Variable intermedaire.
      DOUBLE PRECISION R0		!Rayon seuil d'une granulometrie en loi de Junge
      DOUBLE PRECISION XNU		!Puissance d'une granulometrie en loi de Junge.
      DOUBLE PRECISION NR0		!Nombre de particules pour un rayon < r0 
                                	!d'une loi de Junge.
					
      DOUBLE PRECISION REFF		!Rayon effectif de la granulometrie (microns).
      DOUBLE PRECISION VEFF		!Variance effective de la granulometrie.
      
      DOUBLE PRECISION WA		!Longueur d'onde en microns.

      DOUBLE PRECISION XMU(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)
                                        !Cosinus des angles de Gauss.
      DOUBLE PRECISION TETA     	!Angle de diffusion.
      
      DOUBLE PRECISION KMAT1		! Section efficace d'extinction (micron^2).
      DOUBLE PRECISION KMAT2 		! Section efficace de diffusion (micron^2). 
      
      DOUBLE PRECISION P11(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Fonction de phase calculee a partir des fonctions de 
		       !la theorie de MIE.
      DOUBLE PRECISION P12(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P12 de la matrice de phase.
      DOUBLE PRECISION P33(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)	
                       !Element P33 de la matrice de phase.
		      
      DOUBLE PRECISION XG	!Facteur d'anisotropie.


      DOUBLE PRECISION ALPHAF	!Valeur maximale du parametre de taille contenue
      				!dans le fichier des calculs de MIE.

      DOUBLE PRECISION R	!Rayon de particule pour le parametre de taille 
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
      
      DOUBLE PRECISION RN	!partie reelle de l'indice de refraction.
      DOUBLE PRECISION RI	!partie imaginaire de l'indice de refraction.
	
      DOUBLE PRECISION X1		

      
      REAL*4 ALPHA		!Valeur du parametre de taille.
      REAL*4 PAS		!Pas entre deux valeurs de ALPHA.
      REAL*4 QEXT		!Coefficient d'efficacite en extinction.
      REAL*4 QSCA		!Coefficient d'efficacite en diffusion.
      REAL*4 IMIE(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)         
                                !Fonction de phase : Imie(mu) pour ALPHA.
      REAL*4 QMIE(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)         
                                !Fonction de phase polarisee : Qmie(mu).
      REAL*4 UMIE(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)         
                                !Fonction de phase polarisee : Umie(mu).

      INTEGER*4 MIE_NBMU  	! Nombre d'angles (positifs) effectivement utiles 
      INTEGER*4 IGRANU  	! Numero du modele de granulometrie.
      INTEGER*4 J	   	! Indice.
      INTEGER*4 NBMUMIE	   	! Nombre d'angle de GAUSS positifs utilises pour 
      			   	! les calculs de MIE (enregistrement du fichier MIE).			   
      INTEGER*4 IER        	! Code d'erreur =0 si pas d'erreur, =-1 sinon

      LOGICAL TRACE        	! = vrai,si ecriture dans le fichier trace.
      

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

C* Initialisation des cumuls d'integration sur la granulometrie
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
        
C* Controle parametre
C--------------------
      IF (NBMUMIE.NE.MIE_NBMU) GOTO 996

C* Lecture des fonctions calculees par la theorie de MIE
C* pour chaque valeur du parametre de taille ALPHA
C*
C*      Lecture des tableaux de fonctions de phase
C*      qui ont ete enregistres avec des bornes 
C*      limitees au dimensionnement utile -MIE_NBMU:MIE_NBMU
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
         R=ALPHA*WA/2.D+00/SOS_PI
	 
C* Test sur la limite du fichier de MIE
C--------------------------------------
         IF (ALPHA.GE.(ALPHAF-PAS)) GOTO 40
		
C* Estimation du pas pour la valeur de ALPHA suivante
C------------------------------------------------------
         IF (ALPHA.GT.0.10) PAS=0.001
         IF (ALPHA.GT.1.00 ) PAS=0.01
         IF (ALPHA.GT.10.) PAS=0.05
         IF (ALPHA.GT.30.) PAS= 0.10
         IF (ALPHA.GT.100.) PAS=1.00


C* Expression des differentes granulometries possibles:
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
            NR = dexp(-B*B/2.)/(R*SIG*sqrt(2*SOS_PI))
         ENDIF

C*  Loi de Junge :
C****************
C*  R0 : Rayon minimal (microns).
C---------------------------
         IF (IGRANU.eq.2) THEN
            R0  = VARGRANU1
            XNU = VARGRANU2

            NR0= (R0)**(-XNU)

            IF (R.GT.R0) THEN
              NR = R**(-XNU)
            ELSE
              NR = NR0
            ENDIF	
         ENDIF


C* Calcul des coefficients d'extinction (kMAT1) et de diffusion (KMAT2).
C* Calcul des fonctions de phase en intensite (P11) et en polarisation(P12 et P33).
C-------------------------------------------------------------------------------
         PR=WA*PAS/2./SOS_PI
         X1=NR*PR*SOS_PI*R**2
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
 
C* Calcul des termes pour le rayon effectif.
C-------------------------------------------------------------------------------
      IF (TRACE) THEN 
	 R2NRPR = R*R*NR*PR
	 SOMME_R2 = SOMME_R2 + R2NRPR
	 SOMME_R3 = SOMME_R3 + R2NRPR*R
	 SOMME_R4 = SOMME_R4 + R2NRPR*R*R
      END IF  
C* Passage a la valeur suivante de ALPHA
C-----------------------------------------
      GOTO 123 	
  999 CONTINUE
   40 CONTINUE

C* Calcul du rayon effectif et de la variance effective
C-------------------------------------------------------------------------------      
      IF (TRACE) THEN
         REFF = SOMME_R3 / SOMME_R2
	 VEFF =  1. +  
     &     ( (SOMME_R4 - 2.*REFF*SOMME_R3) / (REFF*REFF * SOMME_R2) )
      END IF  


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
         WRITE(88,*,err=921) '* ROUTINE GRANU'
	 WRITE(88,*,err=921) '==============='
         WRITE(88,*,err=921) ' '
         WRITE(88,*,err=921) 'PARAMETRES DE SIMULATIONS'
	 WRITE(88,*,err=921) '------------------------- '
	 WRITE(88,*,err=921) '  Fichier MIE : ',FICMIE
	 WRITE(88,*,err=921) '  Longueur d''onde en microns : ',WA
	 IF (IGRANU.EQ.1) THEN
	      WRITE(88,*,err=921)
     &	           '  Modele 1 : Distribution Log-Normale'
	      WRITE(88,*,err=921) '     rm en microns : ',VARGRANU1
	      WRITE(88,*,err=921) '     sig : ',VARGRANU2
	      WRITE(88,*,err=921) '     sig/ln10 : ',
     &	                                VARGRANU2/DLOG(10.D+00)
	 ELSE
	      WRITE(88,*,err=921) '  Modele 2 : Loi de Junge'
	      WRITE(88,*,err=921) '	R0 en microns : ',VARGRANU1
	      WRITE(88,*,err=921) '     nu  : ',VARGRANU2	      	   
	 ENDIF
         WRITE(88,*,err=921) ' '
         WRITE(88,*,err=921) 'RESULTATS'
	 WRITE(88,*,err=921) '--------- '
	 WRITE(88,*,err=921)
     &	     '  Coefficient d''extinction non tronque (par micron) : '
     &	      ,KMAT1*SOMME_NR
	 WRITE(88,*,err=921) 
     &	     '  Coefficient de diffusion non tronque (par micron) : '
     &       ,KMAT2*SOMME_NR
         WRITE(88,*,err=921) 
     &	     '  Nombre de particules (par unite de volume micron^3) : '
     &       ,SOMME_NR
         WRITE(88,*,err=921)
     &	     '  Section efficace d''extinction non tronque (micron^2) : '
     &       ,KMAT1
	 WRITE(88,*,err=921) 
     &	     '  Section efficace de diffusion non tronque (micron^2) : ' 
     &       ,KMAT2
	 WRITE(88,*,err=921) 
     &	     '  Rayon effectif de la granulometrie (microns) : ',REFF
         WRITE(88,*,err=921) 
     &	     '  Variance effective de la granulometrie : ',VEFF
         WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) '  1ere colonne : angle de diffusion '
	 WRITE(88,*,err=921) 
     &       '  2eme colonne : fonction de phase initiale '
         WRITE(88,*,err=921) '                 P11(mu) '
	 WRITE(88,*,err=921) 
     &       '  3eme colonne : fonction de phase polarisee initiale '
         WRITE(88,*,err=921) '                 P12(mu) '
	 WRITE(88,*,err=921) 
     &       '  4eme colonne : fonction de phase polarisee initiale '
         WRITE(88,*,err=921) '                 P33(mu) '
	 WRITE(88,*,err=921) ' '
         DO 44 J=-MIE_NBMU,MIE_NBMU
	    TETA=DACOS(XMU(J))*180./SOS_PI
	    WRITE(88,10,err=921) TETA,P11(J),P12(J),P33(J)

   44    CONTINUE
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
  921 WRITE(6,*) 'Erreur d''ecriture sur le fichier trace GRANU'
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
C      Cette procedure interpole une fonction definie par (XIN, YIN) sur NBIN noeuds
C      aux abscisses XOUT pour NBOUT noeuds (XOUT, YOUT).
C
C      La demarche repose sur l'utilisation des fonctions SPLINE et SPLINT
C         - Les noeuds (XIN, YIN) sont en premier tries par ordre croissant d'abscisse XIN.
C         - La fonction SPLINE calcule les derivees secondes de YIN_TRIE en XIN_TRIE
C           a partir des valeurs YIN_TRIE de la fonction aux abscisses XIN_TRIE et des
C           valeurs de derivee premiere de la fonction au premier et dernier noeud 
C           (min et max de XIN_TRIE).
C	  - La fonction SPLINT calcule les valeurs YOUT aux abscisses XOUT en fonction
C           des valeurs de la fonction et de sa derivee seconde predecemment calcules
C           aux abscisses XIN_TRIE.
C
C Description des parametres
C ----------------------------
C       NBIN (I4)		(E) : Nombre d'elements des tableaux XIN et YIN
C       XIN(1:NBIN)   (double)	(E) : Tableau des abscisses de la fonction pour lesquels l'ordonnee est connue.
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
C     Affichage d'un message a l'ecran, la routine interrompt ses calculs et 
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
      DOUBLE PRECISION YIN(NBIN)   !Valeurs des ordonnees des noeuds (XIN,YIN) 
      DOUBLE PRECISION XOUT(NBOUT)  !Valeurs des abscisses des noeuds (XOUT,YOUT) 
      DOUBLE PRECISION YOUT(NBOUT)  !Valeurs recherchees des ordonnees des noeuds (XOUT,YOUT) 

      DOUBLE PRECISION XIN_TRIE(NBIN)   !Valeurs croissantes des abscisses des noeuds (XIN,YIN) 
      DOUBLE PRECISION YIN_TRIE(NBIN)   !Valeurs des ordonnees des noeuds (XIN,YIN) pour XIN trie par ordre croissant
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
C      Cette procedure calcule les valeurs de derivees secondes d'une fonction definie
C      par (X, Y) sur N noeuds (supposes tries par ordre croissant de X)
C      et par ses valeurs de derivees premieres aux noeuds extremes (X(1),DY1) et (X(N),DYN).
C
C      La routine renvoie les valeurs de derivees secondes en chaque point X.
C
C Description des parametres
C ----------------------------
C       N  (I4)		(E) : Nombre d'elements des tableaux X, Y et DSEC_Y
C       X(1:N) (double)	(E) : Tableau des abscisses de la fonction pour lesquels l'ordonnee est connue.
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

      INTEGER*4 N		 !Nb d'elements des tableaux d'entree (X,Y) et et de sortie (DSEC_Y)
      DOUBLE PRECISION X(N)	 !Tableau des abscisses de la fonction 
      DOUBLE PRECISION Y(N)	 !Tableau des valeurs de la fonction aux abscisses X
      DOUBLE PRECISION DY1       !Derivee premiere de la fonction en X(1)
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
C      Cette procedure interpole une fonction en un point XVAL a partir
C      de trois vecteurs de même dimension X(1:N), Y(1:N) et DSEC_Y(1:N) 
C      qui donnent respectivement les abscisses, les ordonnées et les 
C      dérivées secondes de la fonction.
C
C Description des parametres
C ----------------------------
C       N  (I4)		(E) : Nombre d'elements des tableaux X, Y et DSEC_Y
C       X(1:N) (double)	(E) : Tableau des abscisses de la fonction pour lesquels l'ordonnee est connue.
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

      INTEGER*4 N	!Nb d'elements des tableaux d'entree (X,Y et DSEC_Y)
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
