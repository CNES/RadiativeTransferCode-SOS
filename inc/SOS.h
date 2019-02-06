C******************************************************************************
C* FICHIER: SOS.h
C* PROJET: Ordres successifs de diffusion
C* ROLE: Definition des constantes des programmes du code des OS.
C*
C* AUTEUR: B.Lafrance ( CS )
C* DATE: 30/01/2002
C*
C* MOD:VERSION:3.0
C* 
C******************************************************************************

C#######################################################################
C    1. Constantes communes                        
C#######################################################################

C Longueur des chaines de caracteres
C-----------------------------------
C Repertoires
#define SOS_LENDIR 150

C Nom des fichiers sans arborescence
#define SOS_LENFIC1 150

C Nom complet des fichiers (avec arborescence)
#define SOS_LENFIC2 300

C Taille maximale des Keywords pour le passage des arguments
#define SOS_LENKEYWORD 30


C Valeur de PI
C-------------
#define SOS_PI 3.141592653585



C#######################################################################
C    2. Constantes specifiques aux routines du programme SOS_AEROSOLS                         
C#######################################################################

C Dimension des tableaux --> cste specifique au programme SOS_MIE
C----------------------------------------------------------------
#define SOS_MIE_DIM 10000

C Taille maximale des tableaux de fonctions de phase externes
C --> cste specifique au programme SOS_AEROSOLS
C----------------------------------------------------------------
#define SOS_MAXNB_ANG_EXT 200


C Valeur minimale du parametre de taille Alpha des calculs de MIE.
C (doit etre comprise entre 0.0001 et 10.0 
C  NB: imperativement en double precision)
C    --> cste specifique au programme SOS_AEROSOLS
C----------------------------------------------------------------
#define SOS_MIE_ALPHAMIN 0.0001D+00


C Valeur limite du parametre de taille pour le calcul des fichiers 
C de MIE des constituants des modeles WMO et Shettle & Fenn
C    --> cste specifique au programme SOS_AEROSOLS
C----------------------------------------------------------------
#define SOS_ALPHAMAX_WMO_DL   4000.
#define SOS_ALPHAMAX_WMO_WS   50.
#define SOS_ALPHAMAX_WMO_OC   800.
#define SOS_ALPHAMAX_WMO_SO   10.

#define SOS_ALPHAMAX_SF_SR    70.
#define SOS_ALPHAMAX_SF_SU    90.

C Valeur minimale du rapport de granulometrie n(r) / Nmax 
C pour la determination du parametre de taille limite
C utile aux calculs de Mie.
C    --> cste specifique au programme SOS_AEROSOLS
C----------------------------------------------------------------
#define SOS_COEF_NRMAX 0.0001


C Nom du fichier contenant le rayon modal, le log de la variance,
C la concentration volumique et les valeurs des indices de refraction
C (partie reelle et partie imaginaire) en fonction de la 
C longueur d'onde.
C    --> cste specifique au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define SOS_AER_DATAWMO		"Data_WMO_cor_2015_12_16"

C Nom du fichier contenant le log de la variance et les valeurs 
C du rayon modal en fonction de l'humidite relative.
C    --> cste specifique au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define SOS_AER_DATASF		"Data_SF_cor_2015_12_16"

C Nom des fichiers contenant les valeurs des indices de refraction
C (partie reelle et partie imaginaire) en fonction de la 
C longueur d'onde et de l'humidite relative.
C    --> cstes specifiques au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define SOS_AER_SR_SF		"IRefrac_SR_cor_2015_12_16"
#define SOS_AER_LR_SF		"IRefrac_LR"
#define SOS_AER_SU_SF		"IRefrac_SU_cor_2015_12_16"
#define SOS_AER_LU_SF		"IRefrac_LU_cor_2015_12_16"
#define SOS_AER_OM_SF		"IRefrac_OM_cor_2015_12_16"

C Valeur seuil pour la troncature
C    --> cstes specifiques au programme SOS_AEROSOLS
C---------------------------------------------------------------
#define SOS_PH_SEUIL_TRONCA 0.1



C#######################################################################
C    3. Constantes specifiques a la definition du profil atmospherique  
C#######################################################################

C Nombre de couches du profil atmospherique.
C       --> cste specifique aux programmes SOS, SOS_PROFIL et SOS_OS
C-------------------------------------------------------------------
#define SOS_OS_NT 100


C Nombre minimal de sous-couches moleculaires pour un profil avec
C positionnement variable des aerosols entre deux altitudes.
C       --> cste specifique au programme SOS_PROFIL
C-------------------------------------------------------------------
#define SOS_PROFIL_MIN_NBC 3


C Definition de l'epaisseur de la couche de transition (km)
C    --> cste specifique au programme SOS_PROFIL
C-------------------------------------------------------------------
#define SOS_DZTRANSI	0.010



C#######################################################################
C    4. Constantes specifiques aux programmes de calculs des matrices
C       de reflexion de surface
C#######################################################################

C Facteur pour la comparaison de GMIN et GMAX
C    --> cste specifique au programme SOS_GLITTER
C--------------------------------------------------------------
#define SOS_PH_TEST 10000


C Nombre d'angles azimutaux (2**SOS_PH_NQ)
C    --> cste specifique aux programmes
C        SOS_GLITTER, SOS_ROUJEAN et SOS_NADAL
C--------------------------------------------------------------
#define SOS_PH_NU 1024


C Exposant de la puissance de 2 de SOS_PH_NU (NU=2**NQ)
C    --> cste specifique au programme SOS_GLITTER
C--------------------------------------------------------------
#define SOS_PH_NQ 10


C Valeur seuil pour l'arret de la decomposition en series de 
C Fourier de la fonction (F21 Nadal) / (F21 Fresnel)
C    --> cste specifique au programme SOS_NADAL
C--------------------------------------------------------------
#define SOS_SEUIL_SF_NADAL 0.001


C Valeur seuil pour l'arret de la decomposition en series de 
C Fourier de la BRDF de Roujean
C    --> cste specifique au programme SOS_ROUJEAN
C--------------------------------------------------------------
#define SOS_SEUIL_SF_ROUJEAN 0.001


C Angle zenithal solaire limite de calcul de la BRDF de Roujean
C (en degres)
C       --> cste specifique aux programmes SOS_ROUJEAN 
C           et SOS_NOM_FIC_SURFACE
C--------------------------------------------------------------
#define SOS_TETAS_LIM_ROUJEAN 60


C Angle zenithal de visee limite de calcul de la BRDF de Roujean
C (en degres)
C       --> cste specifique aux programmes SOS_ROUJEAN 
C           et SOS_NOM_FIC_SURFACE
C--------------------------------------------------------------
#define SOS_TETAV_LIM_ROUJEAN 60





C#######################################################################
C    5. Constantes specifiques aux routines du programme SOS.
C#######################################################################

C Ordre minimal du developpement en series de Fourier
C    --> cste specifique au programme SOS_OS
C------------------------------------------------------------
#define SOS_OS_IBOR 0

C Seuil de convergence en serie geometrique
C    --> cste specifique au programme SOS_OS
C------------------------------------------------------------
#define SOS_PH_SEUIL_CV_SG  0.01

C Seuil d'arret des diffusions multiples
C    --> cste specifique au programme SOS_OS
C------------------------------------------------------------
#define SOS_PH_SEUIL_SUMDIF  0.00001

C Seuil pour l'arret de la decomposition en series de Fourier
C    --> cste specifique au programme SOS_OS
C------------------------------------------------------------
#define SOS_PH_SEUIL_SF  0.0001



C Valeur seuil pour le calcul des angles de rotation
C    --> cste specifique au fichier SOS_TRPHI (routine S0S_ANGLE)
C--------------------------------------------------------------
#define SOS_SEUIL_Z 0.0001

C Valeur seuil pour le calcul des elements M21 et M31 de la matrice
C de reflexion exprimee dans un repere lie au plan meridien
C    --> cste specifique au fichier SOS_TRPHI (routine S0S_MATRIC)
C--------------------------------------------------------------
#define SOS_SEUIL_X 0.00001


C Angle solide du disque solaire (sr) pour 
C la distance Terre - Soleil moyenne
C    --> cste specifique au fichier SOS_TRPHI
C--------------------------------------------------------------
#define SOS_SOLAR_DISC_SOLID_ANGLE 6.8D-05



C#######################################################################
C    6. Constantes specifiques a la definition des angles et ordres
C       des developpements
C#######################################################################

C ===================================================
C 6.1 Constantes pour le dimensionnement des tableaux
C     en angles et ordres limites des developpements
C     en series de Fourier et polynomes de Legendre.
C ===================================================

C    Nombre maximal d'angles (positifs) pour le dimensionnement 
C    des tableaux de fonctions de phase.
C    Conseil : * Doit au moins valoir la valeur par defaut SOS_DEFAULT_NBMU_MIE
C                + le nombre max d'angles utilisateur SOS_NBMAX_USER_ANGLES 
C              * A ajuster, si l'usage commun est d'utiliser plus d'angles de 
C                Gauss que la valeur par defaut.
C              * Mais, ne pas definir une valeur abusivement grande pour ne pas 
C                reduire les temps de traitement.
C    --> cste specifique aux programmes SOS_AEROSOLS, SOS_MIE 
C        et SOS_ANGLES
C--------------------------------------------------------------
#define SOS_MIE_NBMU_MAX 100		

C    Nombre maximal d'angles (positifs) pour le dimensionnement des 
C    tableaux de champs de luminance et de reflexion de surface
C    Conseil : * Doit au moins valoir la valeur par defaut SOS_DEFAULT_NBMU_LUM
C                + le nombre max d'angles utilisateur SOS_NBMAX_USER_ANGLES 
C              * A ajuster, si l'usage commun est d'utiliser plus d'angles de 
C                Gauss que la valeur par defaut.
C              * Mais, ne pas definir une valeur abusivement grande pour ne pas 
C                reduire les temps de traitement.
C    --> cste specifique aux programmes :
C        SOS_SURFACE, SOS_GLITTER, SOS_ROUJEAN, SOS_RONDEAUX_BREON
C        SOS_NADAL, SOS, SOS_OS, SOS_TRPHI et SOS_ANGLES
C-----------------------------------------------------------------------
#define SOS_OS_NBMU_MAX 68     

C    Ordre maximal du developpement en series de Fourier
C    et ordre limite du developpement en polynomes de Legendre
C    pour le dimensionnement des tableaux
C    Conseil : Prendre SOS_OS_NB_MAX >= 2*SOS_MIE_NBMU_MAX
C       --> cste specifique aux programmes 
C           SOS, SOS_OS, SOS_AEROSOLS, SOS_SURFACE, SOS_ROUJEAN et SOS_NADAL
C---------------------------------------------------------------------------
#define SOS_OS_NB_MAX 200

C    Ordre maximal du developpement en polynomes de Legendre
C    des elements de la matrice de Fresnel
C    pour le dimensionnement des tableaux
C    Conseil : Prendre  SOS_OS_NS_MAX = 2*SOS_OS_NBMU_MAX
C       --> cste specifique au programme SOS_SURFACE
C---------------------------------------------------------------------------
#define SOS_OS_NS_MAX 136

C    Ordre maximal du developpement en series de Fourier
C    de la fonction G pour le dimensionnement des tableaux
C    Conseil : Prendre  SOS_OS_NM_MAX = SOS_OS_NB_MAX + SOS_OS_NS_MAX
C       --> cste specifique aux programmes
C           SOS_SURFACE, SOS_GLITTER et SOS_RONDEAUX_BREON
C---------------------------------------------------------------------------
#define SOS_OS_NM_MAX 336


C ===================================================
C 6.2 Nombre d'angles a utiliser et ordres maximaux
C     des developpements en absence d'une definition
C     des angles par l'utilisateur.
C ===================================================

C    Nombre d'angles de Gauss positifs par defaut pour les calculs aerosols
C    --> cste specifique au programme SOS_ANGLES 
C--------------------------------------------------------------
#define SOS_DEFAULT_NBMU_MIE 40

C    Nombre d'angles de Gauss positifs par defaut pour les calculs de 
C    champs de luminance
C    --> cste specifique au programme SOS_ANGLES 
C--------------------------------------------------------------
#define SOS_DEFAULT_NBMU_LUM 24

C    Valeur par defaut de l'ordre limite du developpement en polynomes de 
C    Legendre et series de Fourier du champ de luminance :
C    Generalement : SOS_DEFAULT_OS_NB = 2*SOS_DEFAULT_NBMU_MIE
C       --> cste specifique au programme SOS_ANGLES
C---------------------------------------------------------------------------
#define SOS_DEFAULT_OS_NB 80

C    Valeur par defaut de l'ordre limite du developpement en polynomes de 
C    Legendre des elements de la matrice de Fresnel
C    Generalement : SOS_DEFAULT_OS_NS = 2*SOS_DEFAULT_NBMU_LUM
C       --> cste specifique au programme SOS_ANGLES
C---------------------------------------------------------------------------
#define SOS_DEFAULT_OS_NS 48

C    Valeur par defaut de l'ordre limite du developpement en series de 
C    Fourier de la fonction G.
C    Il faut que SOS_DEFAULT_OS_NM >= SOS_DEFAULT_OS_NB + SOS_DEFAULT_OS_NS
C       --> cste specifique au programme SOS_ANGLES
C---------------------------------------------------------------------------
#define SOS_DEFAULT_OS_NM 128


C ===================================================
C 6.3 Limitation du nombres d'angles pour le calcul 
C     des fichiers definissant les angles a utiliser
C ===================================================

C    Nombre maximal d'angles "utilisateurs" 
C    complementaires des angles de Gauss
C       --> cste specifique aux programmes SOS_ANGLES
C--------------------------------------------------------------
#define SOS_NBMAX_USER_ANGLES 20

C    Nombre maximal d'angles positifs pour les calculs des angles 
C    utilises par les codes (angles de Gauss + angles utilisateurs 
C    + 1 = angle zenithal solaire )
C    Conseil : correspond au max de SOS_MIE_NBMU_MAX et SOS_OS_NBMU_MAX
C       --> cste specifique aux programmes SOS_ANGLES
C--------------------------------------------------------------
#define SOS_NBANGLES_MAX 100

C    Seuil entre l'ecart d'un cosinus(X) et cosinus(TETAS) pour 
C    considerer que TETAS = X 
C       --> cste specifique aux programmes SOS_ANGLES
C--------------------------------------------------------------
#define SOS_SEUIL_ECART_MUS 0.00001
