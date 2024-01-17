C******************************************************************************
C* FICHIER: SOS.h
C* PROJET: Ordres successifs de diffusion
C* ROLE: Définition des constantes des programmes du code des OS.
C*
C* AUTEUR: B.Lafrance ( CS )
C* DATE: 30/01/2002
C*
C* MOD:VERSION:3.0
C*
C* MOD:VERSION:3.1 22/02/2016
C*    - Introduction de la constante CTE_SOLAR_DISC_SOLID_ANGLE
C*    - Introduction de la constante CTE_REMOVE_PREVIOUS_BIN_FILE
C*
C* MOD:VERSION:3.2 10/11/2020
C*    - Introduction de la constante CTE_SEUIL_NUM
C*    - Ajustement des valeurs des dimensions des données CKD : 
C*      CTE_CKD_NWVL_MAX, CTE_CKD_NAI_MAX, CTE_CKD_NT_MAX, CTE_CKD_NP_MAX
C*    - Introduction de la constante CTE_HT_STD_PSURF
C*    - Introduction de la constante CTE_MAX_NB_MODE_MIXTURE
C*
C* MOD:VERSION:3.3 31/08/2022
C*    - Introduction de la constante CTE_THRESHOLD_Q_U_NULL
C*    - Introduction des constantes 

C* MOD:VERSION:3.4 15/06/2023
C*    - Ajustement des constantes de dimensionnement des fichiers CKD: 
C*      CTE_CKD_NUMAX (passe de 25000 à 27500 cmm1) et CTE_CKD_NB_NU_PER_FILE (passe de 225 à 50)
C*      CTE_CKD_NWVL_MAX (passe de 225 à 50), CTE_CKD_NAI_MAX(passe de 20 à 5)
C*      CTE_CKD_NT_MAX  (passe de 8 à 9) et CTE_CKD_NP_MAX  (passe de 25 à 31)
C*    - Ajout du paramètre CTE_CKD_NC_MAX qui précise le nombre max de concentrations en vapeur d'eau

C* MOD:VERSION:3.5 27/07/2023
C*    - Ajout de constantes pour la définition du profil:
C*      CTE_TOA_FIRST_LAYER_OPT_THICKNESS : Epaisseur optique au premier niveau du profil
C*      CTE_OS_NT_MIN : Nombre minimal de couches
C*      CTE_DELTA_Z : Pas en altitude pour la recherche de l'altitude du premier niveau
C*    - Le paramètre CTE_SEUIL_TAUABS est renommé CTE_THRESHOLD_TAUABS 
C*
C*    - Ajout du paramètre CTE_GAP_TOLER_SUM_RATES qui définit l'erreur maximale tolérée sur l'écart à 1 
C*      de la somme des taux d'épaisseur optique d'un mélange d'aérosols
C*
C* MOD:VERSION:3.6 18/09/2023
C*    - Conversion de CTE_TOA_ALT de 210. en 120.D+00 pour compilation sous gcc/10.2.0
C*
C#######################################################################
C    1. Constantes communes                        
C#######################################################################

C Longueur des chaines de caracteres
C-----------------------------------
C Repertoires
#define CTE_LENDIR 350

C Nom des fichiers sans arborescence
#define CTE_LENFIC1 150

C Nom complet des fichiers (avec arborescence)
#define CTE_LENFIC2 500

C Taille maximale des Keywords pour le passage des arguments
#define CTE_LENKEYWORD 30

C Length of command system chains
#define CTE_LENCOM 400


C Domaine de validite (en microns)
C---------------------------------
#define CTE_WAMIN 0.364
#define CTE_WAMAX 4.0

C ID for values not defined by the user
C--------------------------------------
C --> Integer type :
#define CTE_NOT_DEFINED_VALUE_INT -999
C --> Double type :
#define CTE_NOT_DEFINED_VALUE_DBLE -999.D+00

C Fichiers par défaut
C--------------------
#define CTE_DEFAULT_FICANGLES_RES_LUM  "SOS_UsedAngles.txt"
#define CTE_DEFAULT_FICANGLES_RES_MIE  "Aer_UsedAngles.txt"
#define CTE_DEFAULT_FICGRANU  "Aerosols.txt"
#define CTE_DEFAULT_RESBIN   "SOS_Result.bin"
#define CTE_DEFAULT_RESUP  "SOS_Up.txt"
#define CTE_DEFAULT_RESDOWN "SOS_Down.txt"
      
C#######################################################################
C    2. Constantes spécifiques aux routines du programme SOS_AEROSOLS                         
C#######################################################################

C Dimension des tableaux --> cste spécifique au programme SOS_MIE
C Doit être supérieur à 2*AlphaMax+20 (AlphaMax: param de taille max).
C----------------------------------------------------------------
#define CTE_MIE_DIM 10000

C Taille maximale des tableaux de fonctions de phase externes
C --> cste spécifique au programme SOS_AEROSOLS
C----------------------------------------------------------------
#define CTE_MAXNB_ANG_EXT 200


C Default value for the maximum radius of aerosols particles
C for a Junge model of size distribution
C    --> Specific constant to the file SOS_AEROSOLS.F
C----------------------------------------------------------------
#define CTE_DEFAULT_AER_JUNGE_RMAX 50.


C Valeur minimale du paramètre de taille Alpha des calculs de MIE.
C (doit être comprise entre 0.0001 et 10.0 
C  NB: impérativement en double précision)
C    --> cste spécifique au programme SOS_AEROSOLS
C----------------------------------------------------------------
#define CTE_MIE_ALPHAMIN 0.0001D+00

C Valeur limite du paramètre de taille pour le calcul des fichiers 
C de MIE des constituants des modèles WMO et Shettle & Fenn
C    --> cste spécifique au programme SOS_AEROSOLS
C----------------------------------------------------------------
#define CTE_ALPHAMAX_WMO_DL   4000.
#define CTE_ALPHAMAX_WMO_WS   50.
#define CTE_ALPHAMAX_WMO_OC   800.
#define CTE_ALPHAMAX_WMO_SO   10.
#define CTE_ALPHAMAX_SF_SR    70.
#define CTE_ALPHAMAX_SF_SU    90.

C Valeur minimale du rapport de granulométrie n(r) / Nmax 
C pour la détermination du paramètre de taille limite
C utile aux calculs de Mie.
C    --> cste spécifique au programme SOS_AEROSOLS
C----------------------------------------------------------------
#define CTE_COEF_NRMAX 0.0001


C Nom du fichier contenant le rayon modal, le log de la variance,
C la concentration volumique et les valeurs des indices de réfraction
C (partie reelle et partie imaginaire) en fonction de la 
C longueur d'onde.
C    --> cste spécifique au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_AER_DATAWMO		"Data_WMO_cor_2015_12_16"

C Nom du fichier contenant le log de la variance et les valeurs 
C du rayon modal en fonction de l'humidite relative.
C    --> cste spécifique au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_AER_DATASF		"Data_SF_cor_2015_12_16"

C Nom des fichiers contenant les valeurs des indices de réfraction
C (partie reelle et partie imaginaire) en fonction de la 
C longueur d'onde et de l'humidite relative.
C    --> cstes spécifiques au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_AER_SR_SF		"IRefrac_SR_cor_2015_12_16"
#define CTE_AER_LR_SF		"IRefrac_LR"
#define CTE_AER_SU_SF		"IRefrac_SU_cor_2015_12_16"
#define CTE_AER_LU_SF		"IRefrac_LU_cor_2015_12_16"
#define CTE_AER_OM_SF		"IRefrac_OM_cor_2015_12_16"

C Angles used to define the angular range for the phaze function 
C linearization while applying a truncation.
C    --> Specific constant to the files SOS_AEROSOLS.F
C---------------------------------------------------------------
#define CTE_AER_MU1_TRONCA 0.8
#define CTE_AER_MU2_TRONCA 0.94

C Valeur seuil pour la troncature
C    --> cstes spécifiques au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_PH_SEUIL_TRONCA 0.1


C Nombre maximal de modes pour un mélange d'aérosols
C    --> cstes spécifiques au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_MAX_NB_MODE_MIXTURE 20

C Seuil de tolérance d'écart à 1 de la somme des taux d'épaisseur 
C optique d'un mélange d'aérosols
C    --> cstes spécifiques au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_GAP_TOLER_SUM_RATES 0.000001

      
C#######################################################################
C    3. Constantes spécifiques à la définition du profil atmosphérique  
C#######################################################################

C Pression atmosphérique standard
#define CTE_HT_STD_PSURF 1013.0

C Altitude max (en km) du profil atmosphérique (niveau TOA)
C       --> cste spécifique aux programmes SOS, SOS_PROFIL et SOS_OS
C-------------------------------------------------------------------
#define CTE_TOA_ALT 120.D+00

C Nombre de couches maximum du profil atmosphérique.
C       --> cste spécifique aux programmes SOS, SOS_PROFIL et SOS_OS  
C-------------------------------------------------------------------
#define CTE_OS_NT 600

C Epaisseur optique minimale d'une couche du profil atmosphérique 
C (hors premier niveau).
C       --> cste spécifique aux programmes SOS, SOS_PROFIL et SOS_OS  
C-------------------------------------------------------------------
#define CTE_TCOUCHE 0.005

C Epaisseur optique de la première couche du profil atmosphérique.
C       --> cste spécifique au programme SOS_PROFIL   
C-------------------------------------------------------------------
#define CTE_TOA_FIRST_LAYER_OPT_THICKNESS 0.0002

C Pas (en km) pour la recherche de l'altitude du premier niveau
C       --> cste spécifique au programme SOS_PROFIL   
C-------------------------------------------------------------------
#define CTE_DELTA_Z 0.05

C Seuil (en km) pour la comparaison des altitudes lors de la définition 
C du profil 
C       --> cste spécifique au programme SOS_PROFIL   
C-------------------------------------------------------------------
#define CTE_THRESHOLD_DZ 0.001

C Nombre de couches minimum du profil atmosphérique.
C       --> cste spécifique aux programmes SOS, SOS_PROFIL et SOS_OS
C-------------------------------------------------------------------
#define CTE_OS_NT_MIN 100

C Nombre minimal de sous-couches moléculaires pour un profil avec
C positionnement variable des aerosols entre deux altitudes.
C       --> cste spécifique au programme SOS_PROFIL
C-------------------------------------------------------------------
#define CTE_PROFIL_MIN_NBC 3

C Définition de l'épaisseur de la couche de transition (km)
C    --> cste spécifique au programme SOS_PROFIL
C-------------------------------------------------------------------
#define CTE_DZTRANSI	0.010

C Nombre de corps absorbants (il y a 8 espèces à prendre en compte, 
C incluant les raies d'absorption et/ou les continuums) :
C --> H2O CO2 O3 N2O CO CH4 O2 NO2
C----------------------------------------------------------------------------
#define CTE_NBABS  8

C Nombre de niveaux sur les profils des gaz d'absorption 
C-------------------------------------------------------
#define CTE_ABS_NBLEV  50

C Nombre de colonnes des fichiers de profils des gaz d'absorption 
C----------------------------------------------------------------
#define CTE_ABS_NBCOL  13

C Définition des bandes d'absorption à considérer (1 pour prendre en compte, 0 sinon)
C Absorbants: H2O CO2 O3 N2O CO CH4 O2 NO2
C Les absorptions du NO2 sont uniquement définis par continuums.
C Les autres gaz contiennent des raies et continuums.
C-------------------------------------------------------------------------------
#define CTE_ABS_H2O  1
#define CTE_ABS_CO2  1
#define CTE_ABS_O3   1
#define CTE_ABS_N2O  1
#define CTE_ABS_CO   1
#define CTE_ABS_CH4  1
#define CTE_ABS_O2   1
#define CTE_CONT_NO2 1

C Dimensionnement maximal des fichiers de coefficients CKD
C-------------------------------------------------------------------------------
C CTE_CKD_NWVL_MAX : Nombre maximal de longueurs d'onde dans les tables de coeff CKD
C CTE_CKD_NAI_MAX  : Nombre maximal d'exponentiels des développements 
C CTE_CKD_NT_MAX   : Nombre maximal de températures
C CTE_CKD_NP_MAX   : Nombre maximal de pressions
C CTE_CKD_NC_MAX   : Nombre maximal de concentrations H2O (spécifique à la vapeur d'eau)
C-------------------------------------------------------------------------------
#define CTE_CKD_NWVL_MAX  50
#define CTE_CKD_NAI_MAX  5
#define CTE_CKD_NT_MAX  9
#define CTE_CKD_NP_MAX  31
#define CTE_CKD_NC_MAX  12

C Définition du découpage spectral des fichiers de coefficients CKD
C-------------------------------------------------------------------------------
C CTE_CKD_NUMAX   : Nombre d'onde maximal couvert par les fichiers de coeff CKD(en cm-1)
C CTE_CKD_NUMIN   : Nombre d'onde minimal couvert par les fichiers de coeff CKD (en cm-1)
C CTE_CKD_NB_NU_PER_FILE : Nombre de longueurs d'onde par fichier CKD (doit être identique pour toutes les résolutions spectrales)
C-----------------------------------------------
#define CTE_CKD_NUMAX  27500
#define CTE_CKD_NUMIN   2500
#define CTE_CKD_NB_NU_PER_FILE  50


C Epaisseur optique d'absorption maximale (applicable aux cas extrêmes si transmission nulle)
C-------------------------------------------------------------------------------
#define CTE_TAUABS_MAX   999.

C Seuil sur l'épaisseur optique d'absorption pour optimisation du profil atmosphérique
C ------------------------------------------------------------------------------------
#define CTE_THRESHOLD_TAUABS 1.5


C#######################################################################
C    4. Constantes spécifiques aux programmes de calculs des matrices
C       de reflexion de surface
C#######################################################################

C Facteur pour la comparaison de GMIN et GMAX
C    --> cste spécifique au programme SOS_GLITTER
C--------------------------------------------------------------
#define CTE_PH_TEST 10000


C Nombre d'angles azimutaux (2**CTE_PH_NQ)
C    --> cste spécifique aux programmes
C        SOS_GLITTER, SOS_ROUJEAN et SOS_NADAL
C--------------------------------------------------------------
#define CTE_PH_NU 1024


C Exposant de la puissance de 2 de CTE_PH_NU (NU=2**NQ)
C    --> cste spécifique au programme SOS_GLITTER
C--------------------------------------------------------------
#define CTE_PH_NQ 10


C Valeur seuil pour l'arret de la decomposition en series de 
C Fourier de la fonction (F21 Nadal) / (F21 Fresnel)
C    --> cste spécifique au programme SOS_NADAL
C--------------------------------------------------------------
#define CTE_SEUIL_SF_NADAL 0.001


C Valeur seuil pour l'arret de la decomposition en series de 
C Fourier de la BRDF de Roujean
C    --> cste spécifique au programme SOS_ROUJEAN
C--------------------------------------------------------------
#define CTE_SEUIL_SF_ROUJEAN 0.001


C Angle zenithal solaire limite de calcul de la BRDF de Roujean
C (en degres)
C       --> cste spécifique aux programmes SOS_ROUJEAN 
C           et SOS_NOM_FIC_SURFACE
C--------------------------------------------------------------
#define CTE_TETAS_LIM_ROUJEAN 60


C Angle zenithal de visee limite de calcul de la BRDF de Roujean
C (en degres)
C       --> cste spécifique aux programmes SOS_ROUJEAN 
C           et SOS_NOM_FIC_SURFACE
C--------------------------------------------------------------
#define CTE_TETAV_LIM_ROUJEAN 60


C Valeur seuil pour tests d'erreurs d'arrondies numériques
C    --> cste spécifique au programme SOS_ROUJEAN
C--------------------------------------------------------------
#define CTE_SEUIL_NUM 1.D-10




C#######################################################################
C    5. Constantes spécifiques aux routines du programme SOS.
C#######################################################################

C Factor of molecular depolarization in the air
C    --> cste spécifique au programme SOS_OS
C--------------------------------------------------------------
#define CTE_MDF 0.0279

C Ordre minimal du developpement en series de Fourier
C    --> cste spécifique au programme SOS_OS
C------------------------------------------------------------
#define CTE_OS_IBOR 0

C Nombre maximal d'ordres de diffusion par défaut
C    --> cste spécifique au programme SOS_OS
C------------------------------------------------------------
#define CTE_DEFAULT_IGMAX 100


C Seuil de convergence en serie geometrique
C    --> cste spécifique au programme SOS_OS
C------------------------------------------------------------
#define CTE_PH_SEUIL_CV_SG  0.00001

C Seuils d'arret des diffusions multiples
C    --> cste spécifiques au programme SOS_OS
C------------------------------------------------------------
#define CTE_PH_SEUIL_SUMDIF  0.00001
#define CTE_PH_SEUIL_VALDIF  1.0D-50

C Seuil pour l'arret de la decomposition en series de Fourier
C    --> cste spécifique au programme SOS_OS
C------------------------------------------------------------
#define CTE_PH_SEUIL_SF  0.00001



C Valeur seuil pour le calcul des angles de rotation
C    --> cste spécifique au fichier SOS_TRPHI (routine SOS_ANGLE)
C--------------------------------------------------------------
#define CTE_SEUIL_Z 0.0001

C Valeur seuil pour le calcul des elements M21 et M31 de la matrice
C de reflexion exprimee dans un repere lie au plan meridien
C    --> cste spécifique au fichier SOS_TRPHI (routine SOS_MATRIC)
C--------------------------------------------------------------
#define CTE_SEUIL_X 0.00001

C Threshold value under which Q or U is fixed to be null
C    --> Specific constant to the file SOS_TRPHI.F
C--------------------------------------------------------------
#define CTE_THRESHOLD_Q_U_NULL 1.D-15



C Angle solide du disque solaire (sr) pour 
C la distance Terre - Soleil moyenne
C    --> cste spécifique au fichier SOS_TRPHI
C--------------------------------------------------------------
#define CTE_SOLAR_DISC_SOLID_ANGLE 6.8D-05

C Option de suppression ou non du fichier -SOS.ResBin s'il existe déjà
C d'une simulation précédente.
C    0 : pas de suppression --> Le code s'arrete avec un message d'arret.
C    1 : le code supprime le fichier avec un message d'avertissement.
#define CTE_REMOVE_PREVIOUS_BIN_FILE 1


C#######################################################################
C    6. Constantes spécifiques à la définition des angles et ordres
C       des developpements
C#######################################################################

C ===================================================
C 6.1 Constantes pour le dimensionnement des tableaux
C     en angles et ordres limites des developpements
C     en series de Fourier et polynomes de Legendre.
C ===================================================

C    Nombre maximal d'angles (positifs) pour le dimensionnement 
C    des tableaux de fonctions de phase.
C    Conseil : * Doit au moins valoir la valeur par défaut CTE_DEFAULT_NBMU_MIE
C                + le nombre max d'angles utilisateur CTE_NBMAX_USER_ANGLES 
C              * A ajuster, si l'usage commun est d'utiliser plus d'angles de 
C                Gauss que la valeur par défaut.
C              * Mais, ne pas définir une valeur abusivement grande pour ne pas 
C                reduire les temps de traitement.
C    --> cste spécifique aux programmes SOS_AEROSOLS, SOS_MIE 
C        et SOS_ANGLES
C--------------------------------------------------------------
#define CTE_MIE_NBMU_MAX 100		

C    Nombre maximal d'angles (positifs) pour le dimensionnement des 
C    tableaux de champs de luminance et de reflexion de surface
C    Conseil : * Doit au moins valoir la valeur par défaut CTE_DEFAULT_NBMU_LUM
C                + le nombre max d'angles utilisateur CTE_NBMAX_USER_ANGLES 
C              * A ajuster, si l'usage commun est d'utiliser plus d'angles de 
C                Gauss que la valeur par défaut.
C              * Mais, ne pas définir une valeur abusivement grande pour ne pas 
C                reduire les temps de traitement.
C    --> cste spécifique aux programmes :
C        SOS_SURFACE, SOS_GLITTER, SOS_ROUJEAN, SOS_RONDEAUX_BREON
C        SOS_NADAL, SOS, SOS_OS, SOS_TRPHI et SOS_ANGLES
C-----------------------------------------------------------------------
#define CTE_OS_NBMU_MAX 80    

C    Ordre maximal du developpement en series de Fourier
C    et ordre limite du developpement en polynomes de Legendre
C    pour le dimensionnement des tableaux
C    Conseil : Prendre CTE_OS_NB_MAX >= 2*CTE_MIE_NBMU_MAX
C       --> cste spécifique aux programmes 
C           SOS, SOS_OS, SOS_AEROSOLS, SOS_SURFACE, SOS_ROUJEAN et SOS_NADAL
C---------------------------------------------------------------------------
#define CTE_OS_NB_MAX 200

C    Ordre maximal du developpement en polynomes de Legendre
C    des elements de la matrice de Fresnel
C    pour le dimensionnement des tableaux
C    Conseil : Prendre  CTE_OS_NS_MAX = 2*CTE_OS_NBMU_MAX
C       --> cste spécifique au programme SOS_SURFACE
C---------------------------------------------------------------------------
#define CTE_OS_NS_MAX 136

C    Ordre maximal du developpement en series de Fourier
C    de la fonction G pour le dimensionnement des tableaux
C    Conseil : Prendre  CTE_OS_NM_MAX = CTE_OS_NB_MAX + CTE_OS_NS_MAX
C       --> cste spécifique aux programmes
C           SOS_SURFACE, SOS_GLITTER et SOS_RONDEAUX_BREON
C---------------------------------------------------------------------------
#define CTE_OS_NM_MAX 336


C ===================================================
C 6.2 Nombre d'angles a utiliser et ordres maximaux
C     des developpements en absence d'une définition
C     des angles par l'utilisateur.
C ===================================================

C    Nombre d'angles de Gauss positifs par défaut pour les calculs aerosols
C    --> cste spécifique au programme SOS_ANGLES 
C--------------------------------------------------------------
#define CTE_DEFAULT_NBMU_MIE 40

C    Nombre d'angles de Gauss positifs par défaut pour les calculs de 
C    champs de luminance
C    --> cste spécifique au programme SOS_ANGLES 
C--------------------------------------------------------------
#define CTE_DEFAULT_NBMU_LUM 24

C    Valeur par défaut de l'ordre limite du developpement en polynomes de 
C    Legendre et series de Fourier du champ de luminance :
C    Generalement : CTE_DEFAULT_OS_NB = 2*CTE_DEFAULT_NBMU_MIE
C       --> cste spécifique au programme SOS_ANGLES
C---------------------------------------------------------------------------
#define CTE_DEFAULT_OS_NB 80

C    Valeur par défaut de l'ordre limite du developpement en polynomes de 
C    Legendre des elements de la matrice de Fresnel
C    Generalement : CTE_DEFAULT_OS_NS = 2*CTE_DEFAULT_NBMU_LUM
C       --> cste spécifique au programme SOS_ANGLES
C---------------------------------------------------------------------------
#define CTE_DEFAULT_OS_NS 48

C    Valeur par défaut de l'ordre limite du developpement en series de 
C    Fourier de la fonction G.
C    Il faut que CTE_DEFAULT_OS_NM >= CTE_DEFAULT_OS_NB + CTE_DEFAULT_OS_NS
C       --> cste spécifique au programme SOS_ANGLES
C---------------------------------------------------------------------------
#define CTE_DEFAULT_OS_NM 128


C ===================================================
C 6.3 Limitation du nombres d'angles pour le calcul 
C     des fichiers définissant les angles a utiliser
C ===================================================

C    Nombre maximal d'angles "utilisateurs" 
C    complementaires des angles de Gauss
C       --> cste spécifique aux programmes SOS_ANGLES
C--------------------------------------------------------------
#define CTE_NBMAX_USER_ANGLES 20

C    Nombre maximal d'angles positifs pour les calculs des angles 
C    utilises par les codes (angles de Gauss + angles utilisateurs 
C    + 1 = angle zenithal solaire )
C    Conseil : correspond au max de CTE_MIE_NBMU_MAX et CTE_OS_NBMU_MAX
C       --> cste spécifique aux programmes SOS_ANGLES
C--------------------------------------------------------------
#define CTE_NBANGLES_MAX 100

C    Seuil entre l'ecart d'un cosinus(X) et cosinus(TETAS) pour 
C    considerer que TETAS = X 
C       --> cste spécifique aux programmes SOS_ANGLES
C--------------------------------------------------------------
#define CTE_SEUIL_ECART_MUS 0.00001
