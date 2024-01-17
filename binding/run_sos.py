#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
author: Thierry Erudel
organization: CS GROUP France
date: 05/09/2022
version: 0.1

"""

#Version 1.1 : 09/08/2023
# --> Ajout de coef_tronca en paramètre de sortie


import os
from pdb import set_trace

import numpy as np

import sos

SOS_NOT_DEFINED_VALUE_INT = -999
SOS_NOT_DEFINED_VALUE_DBLE = -999.0
SOS_DEFAULT_AER_JUNGE_RMAX = 50.

SOS_DEFAULT_FICANGRESLUM = "SOS_UsedAngles.txt"
SOS_DEFAULT_FICANGRESMIE = "Aer_UsedAngles.txt"
SOS_DEFAULT_FICGRANU = "Aerosols.txt"
SOS_DEFAULT_RESBIN = "SOS_Result.bin"
SOS_DEFAULT_RESUP = "SOS_Up.txt"
SOS_DEFAULT_RESDOWN = "SOS_Down.txt"

# Initialisation des paramètres optionnels avec leur valeur par défaut
# Initialisation des valeurs par défaut des paramètres pour les angles

# Nombre d'angles de Gauss pour les calculs de luminance : non défini
NBMU_GAUSS_LUM = SOS_NOT_DEFINED_VALUE_INT

# Nombre d'angles de Gauss pour les calculs de Mie : non défini
NBMU_GAUSS_MIE = SOS_NOT_DEFINED_VALUE_INT

# Absence de fichier d'angles "utilisateur" pour les calculs de luminance  :
# non défini
FICANGLES_USER_LUM = "NO_USER_ANGLES"

# Absence de fichier d'angles "utilisateur" pour les calculs de Mie : non
# défini
FICANGLES_USER_MIE = "NO_USER_ANGLES"

# Pas de fichier trace de la définition des angles
FICANGLOG = "NO_LOG_FILE"

# Fichiers résultats
FICANGLES_RES_LUM = SOS_DEFAULT_FICANGRESLUM
FICANGLES_RES_MIE = SOS_DEFAULT_FICANGRESMIE

# Initialisation des valeurs par défaut des paramètres pour les aérosols

# Paramètres pour le model mono-modal
RN = SOS_NOT_DEFINED_VALUE_DBLE
IN = SOS_NOT_DEFINED_VALUE_DBLE
RN_WAREF = SOS_NOT_DEFINED_VALUE_DBLE
IN_WAREF = SOS_NOT_DEFINED_VALUE_DBLE
IGRANU = SOS_NOT_DEFINED_VALUE_INT
LND_RADIUS_MMD_AER = SOS_NOT_DEFINED_VALUE_DBLE
LND_LNVAR_MMD_AER = SOS_NOT_DEFINED_VALUE_DBLE
JD_SLOPE_MMD_AER = SOS_NOT_DEFINED_VALUE_DBLE
JD_RMIN_MMD_AER = SOS_NOT_DEFINED_VALUE_DBLE
JD_RMAX_MMD_AER = SOS_DEFAULT_AER_JUNGE_RMAX  # Rmax par défaut
VARGRANU1_MMD_AER = SOS_NOT_DEFINED_VALUE_DBLE
VARGRANU2_MMD_AER = SOS_NOT_DEFINED_VALUE_DBLE
VARGRANU3_MMD_AER = SOS_NOT_DEFINED_VALUE_DBLE

# Paramètres model WMO
C_WMO_DL = SOS_NOT_DEFINED_VALUE_DBLE
C_WMO_WS = SOS_NOT_DEFINED_VALUE_DBLE
C_WMO_OC = SOS_NOT_DEFINED_VALUE_DBLE
C_WMO_SO = SOS_NOT_DEFINED_VALUE_DBLE
IModele_WMO = SOS_NOT_DEFINED_VALUE_INT

# Paramètres S&F
IModele_SF = SOS_NOT_DEFINED_VALUE_INT
RH = SOS_NOT_DEFINED_VALUE_DBLE

# Paramètres bimode LND
MODE_PARAM_BILND = SOS_NOT_DEFINED_VALUE_INT
USER_CV_COARSE = SOS_NOT_DEFINED_VALUE_DBLE
USER_CV_FINE = SOS_NOT_DEFINED_VALUE_DBLE
RTAUct_WAREF = SOS_NOT_DEFINED_VALUE_DBLE
BMD_CM_MRWA = SOS_NOT_DEFINED_VALUE_DBLE
BMD_CM_MIWA = SOS_NOT_DEFINED_VALUE_DBLE
BMD_CM_MRWAREF = SOS_NOT_DEFINED_VALUE_DBLE
BMD_CM_MIWAREF = SOS_NOT_DEFINED_VALUE_DBLE
BMD_CM_RMODAL = SOS_NOT_DEFINED_VALUE_DBLE
BMD_CM_VAR = SOS_NOT_DEFINED_VALUE_DBLE
BMD_FM_MRWA = SOS_NOT_DEFINED_VALUE_DBLE
BMD_FM_MIWA = SOS_NOT_DEFINED_VALUE_DBLE
BMD_FM_MRWAREF = SOS_NOT_DEFINED_VALUE_DBLE
BMD_FM_MIWAREF = SOS_NOT_DEFINED_VALUE_DBLE
BMD_FM_RMODAL = SOS_NOT_DEFINED_VALUE_DBLE
BMD_FM_VAR = SOS_NOT_DEFINED_VALUE_DBLE

# Paramètres donnéns utilisateur
FICEXTDATA_AER = "NO_USER_AEROSOLS_PHAZE_FCT"

# Paramètres mélange de modèles par utilisateur
FICMIXTURE_AER = "NO_USER_AEROSOLS_MIXTURE"

# Nom du fichier d'aerosols utilisateur
FICUSER_AER = "NO_USER_AEROSOLS"

# Nom par défaut du fichier des propriétés radiatives des aérosols
FICGRANU = SOS_DEFAULT_FICGRANU

# Pas de fichier trace des calculs de Mie pour les aérosols
FICMIE_LOG = "NO_LOG_FILE"

# Pas de fichier trace de la définition des matrices de phase aérosols
FICGRANU_LOG = "NO_LOG_FILE"

# Troncature de la pointe avant de la fonction de phase aérosols
ITRONC_AER = 1

# Initialisation des valeurs par défaut des paramètres de surface

# Pas de fichier trace de la définition des matrices de phase aérosols
FICSURF_LOG = "NO_LOG_FILE"

# Indice de réfraction
SURF_IND = SOS_NOT_DEFINED_VALUE_DBLE

# Paramètres glitter
WIND = SOS_NOT_DEFINED_VALUE_DBLE

# Paramètres Roujean
K0_ROUJEAN = SOS_NOT_DEFINED_VALUE_DBLE
K1_ROUJEAN = SOS_NOT_DEFINED_VALUE_DBLE
K2_ROUJEAN = SOS_NOT_DEFINED_VALUE_DBLE

# Paramètres Nadal
ALPHA_NADAL = SOS_NOT_DEFINED_VALUE_DBLE
BETA_NADAL = SOS_NOT_DEFINED_VALUE_DBLE

# Paramètres Maignan
COEF_C_MAIGNAN = SOS_NOT_DEFINED_VALUE_DBLE

# Fichier résultat (fichier utilise dans SOS, qui remplace celui produit pas
# SOS_surface)
FICSURF = "DEFAULT"

# Répertoire de sortie
DIR_ROUJ = "NO_ROUJEAN_DIR"

# Répertoire de sortie
SURF_OUTDIR = "NO_OUTPUT_DIR"

# Initialisation des valeurs par défaut des paramètres de profil

# Epaisseur optique moléculaire
TR = SOS_NOT_DEFINED_VALUE_DBLE

# Echelles de hauteur
HA = SOS_NOT_DEFINED_VALUE_DBLE
HR = SOS_NOT_DEFINED_VALUE_DBLE

# Pression
PSURF = SOS_NOT_DEFINED_VALUE_DBLE

# Contenus en gaz absorbants
H2O = SOS_NOT_DEFINED_VALUE_DBLE
O3 = SOS_NOT_DEFINED_VALUE_DBLE
CO2 = SOS_NOT_DEFINED_VALUE_DBLE
CH4 = SOS_NOT_DEFINED_VALUE_DBLE

# Altitudes min et max du profil aérosol
ZMIN = SOS_NOT_DEFINED_VALUE_DBLE
ZMAX = SOS_NOT_DEFINED_VALUE_DBLE

# Fichier de profil de gaz absorbants
FICABSPROFIL = "NO_USER_ABS_PROFILE_FILE"

# Pas de fichier trace
FICPROFIL_LOG = "NO_LOG_FILE"

# Initialisation des valeurs par défaut des paramètres SOS
# Pas de fichier trace
FICSOS_LOG = "NO_LOG_FILE"

# Option de coupure de la polarisation : optionnelle
IPOLAR = 1

# Niveau de sortie des simulations de luminance : optionnel
ZOUT = -1.0

# Fichiers de sortie sur angles utilisateurs
FICUP_ANGLES_USER = "NO_OUTPUT"
FICDOWN_ANGLES_USER = "NO_OUTPUT"

# Fichier de sortie des transmissions
FICTRANS = "NO_OUTPUT"

# Fichier de sortie des flux
FICFLUX = "NO_OUTPUT"

# Fichiers de résultats
FICSOS_RES_BIN = SOS_DEFAULT_RESBIN
FICUP = SOS_DEFAULT_RESUP
FICDOWN = SOS_DEFAULT_RESDOWN

# Définition des azimut
PHIOS = SOS_NOT_DEFINED_VALUE_DBLE
PAS_PHI = SOS_NOT_DEFINED_VALUE_INT

# Init des épaisseurs optiques pour aggrégation
TTOT_VRAI = 0.
TTOT_TRONC = 0.
TAUOUT = 0.

def gen_hdr_sos_output(sosView,updown,zalt):

    sep='#-----------------------------------------------------------------------------------------------------\n'
    azimutTxt=''
    if sosView == 1 :
        azimutTxt= 'THE AZIMUTH ANGLE AND '
	
    if updown == 1 :
        hdr = '#UPWARD RADIANCE FIELD VERSUS %sVIEWING ZENITH ANGLE\n' %(azimutTxt)
    else :
        hdr = '#DOWNWARD RADIANCE FIELD VERSUS %sVIEWING ZENITH ANGLE\n' %(azimutTxt)

    if sosView == 1 :
        txt = '# (RELATIVE AZIMUTH AND ALTITUDE ARE FIXED)\n'
    else :
        txt = '# (ALTITUDE IS FIXED)'
    hdr = hdr + txt
    hdr = hdr + sep
    
    if sosView == 1 :
        hdr = hdr + '# Relative azimuth (degrees) :\n#\n'
	
    hdr = hdr + '#      Relative azimuth convention :\n'
    
    if sosView == 1 & updown == 2 :
        hdr = hdr + '#        180 deg <-> Viewing direction and Sun in the same half-plane\n'
        hdr = hdr + '#          0 deg <-> Viewing direction and Sun in opposite half-planes with respect to the zenith direction\n#\n'              
    else :
        hdr = hdr + '#        180 deg <-> Satellite and Sun in the same half-plane\n'
        hdr = hdr + '#          0 deg <-> Satellite and Sun in opposite half-planes with respect to the zenith direction\n#\n'
	
    hdr = hdr + '# Value of the selected altitude for the output (km) : %s\n#\n' %(zalt)
    
    hdr = hdr + '# Columns parameters :\n'
    if sosView == 2 :  
        hdr = hdr + '#   PHI     :  Relative azimuth Angle (in degrees)\n'
    hdr = hdr + '#   VZA     :  Viewing Zenith Angle (in degrees)\n'
    hdr = hdr + '#   SCA_ANG :  Scattering angle (in degrees)\n'
    hdr = hdr + '#   I       :  Stokes parameter I at output altitude z (in sr-1)\n'
    hdr = hdr + '#              normalised to the extraterrestrial solar irradiance (PI * L(z) / Esun)\n'
    hdr = hdr + '#              normalised to the extraterrestrial solar irradiance\n'
    hdr = hdr + '#   Q       :  Stokes parameter Q at output altitude z (in sr-1)\n'
    hdr = hdr + '#              normalised to the extraterrestrial solar irradiance \n'    
    hdr = hdr + '#   U       :  Stokes parameter U at output altitude z (in sr-1)\n'
    hdr = hdr + '#              normalised to the extraterrestrial solar irradiance \n'
    hdr = hdr + '#   POL_ANG :  Polarization angle (in degrees). Note: if undefined the value is -999.00\n'
    hdr = hdr + '#   POL_RATE:  Degree of polarization (in %)\n'
    hdr = hdr + '#   IPOL    :  Polarized intensity at level z (in sr-1)\n'
    hdr = hdr + '#              normalised to the extraterrestrial solar irradiance (PI * Lpol(z) / Esun)\n'
    hdr = hdr + sep
         
    azimutTxt=''	 
    if sosView == 2 :
        azimutTxt= 'PHI      '	 
    hdr = hdr + '#   %sVZA     SCA_ANG        I              Q              U       POL_ANG  POL_RATE    IPOL\n' %(azimutTxt)	
    if sosView == 2 :
        azimutTxt= '(degrees) '	
    hdr = hdr + '#%s(degrees) (degrees)  (no unit)      (no unit)      (no unit)   (degrees) (pcts)  (no unit)\n' %(azimutTxt)
     
    return hdr

def gen_sos_output(repOut,sosView,updown,zalt,nblum,pas_phi,phi,vza,sca_ang,i_out,q_out,u_out,pol_ang_out,pol_rate_out,l_pol_out):

    if updown == 1 :
        nom_fic_out='SOS_Up.txt'
    else :
        nom_fic_out='SOS_Down.txt'
	
    fic = open(os.path.join(repOut,nom_fic_out),'w')
    
    #Ecriture de l'entete
    hdr = gen_hdr_sos_output(sosView,updown,zalt)
    fic.write(hdr)
    
    if sosView == 1 :
        
        #On ecrit d'abord le demi-plan PHI + 180°
        for iteta in range(nblum-1,-1,-1) :
            txt='  %7.2f %7.2f  %13.6e  %13.6e  %13.6e  %7.2f %7.2f %13.6e\n' \
            %(-vza[iteta],sca_ang[0,iteta],i_out[0,iteta],q_out[0,iteta],u_out[0,iteta],pol_ang_out[0,iteta],pol_rate_out[0,iteta],l_pol_out[0,iteta])
            fic.write(txt)
	    
        #Puis le demi-plan PHI	
        for iteta in range(nblum) :
            txt='  %7.2f %7.2f  %13.6e  %13.6e  %13.6e  %7.2f %7.2f %13.6e\n' \
            %(vza[iteta],sca_ang[1,iteta],i_out[1,iteta],q_out[1,iteta],u_out[1,iteta],pol_ang_out[1,iteta],pol_rate_out[1,iteta],l_pol_out[1,iteta])
            fic.write(txt)

    else : 
    
        nbphi = ceil(360./pas_phi)
        for iphi in range(nbphi) :	    	    
            for iteta in range(nblum) :
                txt=' %7.2f %7.2f %7.2f  %13.6e  %13.6e  %13.6e  %7.2f %7.2f %13.6e\n)' \
                %(phi[iphi],vza[iteta],sca_ang[iphi,iteta],i_out[iphi,iteta],q_out[iphi,iteta],\
                u_out[iphi,iteta],pol_ang_out[iphi,iteta],pol_rate_out[iphi,iteta],l_pol_out[iphi,iteta])
                fic.write(txt)

    fic.close()
        	
def set_sos_params(dict_sos, trace=True):
    resroot = dict_sos['-SOS_Main.ResRoot']
    ficmain_log = dict_sos['-SOS_Main.Log']
    wa_simu = dict_sos['-SOS_Main.Wa']
    nbmu_gauss_lum = dict_sos['-ANG.Rad.NbGauss']
    ficangles_user_lum = dict_sos['-ANG.Rad.UserAngFile']
    tetas = dict_sos['-ANG.Thetas']
    ficangles_res_lum = dict_sos['-ANG.Rad.ResFile']
    nbmu_gauss_mie = dict_sos['-ANG.Aer.NbGauss']
    ficangles_user_mie = dict_sos['-ANG.Aer.UserAngFile']
    ficangles_res_mie = dict_sos['-ANG.Aer.ResFile']
    ficanglog = dict_sos['-ANG.Log']
    waref_aot = dict_sos['-AER.Waref']
    aot_ref = dict_sos['-AER.AOTref']
    itronc_aer = dict_sos['-AER.Tronca']
    ficgranu_log = dict_sos['-AER.Log']
    ficmie_log = dict_sos['-AER.MieLog']
    dir_mie = dict_sos['-AER.DirMie']
    ficgranu = dict_sos['-AER.ResFile']
    imod_aer = dict_sos['-AER.Model']
    rn_wa = dict_sos['-AER.MMD.MRwa']
    in_wa = dict_sos['-AER.MMD.MIwa']
    rn_waref = dict_sos['-AER.MMD.MRwaref']
    in_waref = dict_sos['-AER.MMD.MIwaref']
    igranu = dict_sos['-AER.MMD.SDtype']
    lnd_radius_mmd_aer = dict_sos['-AER.MMD.LNDradius']
    lnd_lnvar_mmd_aer = dict_sos['-AER.MMD.LNDvar']
    jd_slope_mmd_aer = dict_sos['-AER.MMD.JD.slope']
    jd_rmin_mmd_aer = dict_sos['-AER.MMD.JD.rmin']
    jd_rmax_mmd_aer = dict_sos['-AER.MMD.JD.rmax']
    imodel_wmo = dict_sos['-AER.WMO.Model']
    c_wmo_dl = dict_sos['-AER.WMO.DL']
    c_wmo_ws = dict_sos['-AER.WMO.WS']
    c_wmo_oc = dict_sos['-AER.WMO.OC']
    c_wmo_so = dict_sos['-AER.WMO.SO']
    imodele_sf = dict_sos['-AER.SF.Model']
    rh = dict_sos['-AER.SF.RH']
    mode_param_bilnd = dict_sos['-AER.BMD.VCdef']
    user_cv_coarse = dict_sos['-AER.BMD.CoarseVC']
    user_cv_fine = dict_sos['-AER.BMD.FineVC']
    rtauct_waref = dict_sos['-AER.BMD.RAOT']
    bmd_cm_mrwa = dict_sos['-AER.BMD.CM.MRwa']
    bmd_cm_miwa = dict_sos['-AER.BMD.CM.MIwa']
    bmd_cm_mrwaref = dict_sos['-AER.BMD.CM.MRwaref']
    bmd_cm_miwaref = dict_sos['-AER.BMD.CM.MIwaref']
    bmd_cm_rmodal = dict_sos['-AER.BMD.CM.SDradius']
    bmd_cm_var = dict_sos['-AER.BMD.CM.SDvar']
    bmd_fm_mrwa = dict_sos['-AER.BMD.FM.MRwa']
    bmd_fm_miwa = dict_sos['-AER.BMD.FM.MIwa']
    bmd_fm_mrwaref = dict_sos['-AER.BMD.FM.MRwaref']
    bmd_fm_miwaref = dict_sos['-AER.BMD.FM.MIwaref']
    bmd_fm_rmodal = dict_sos['-AER.BMD.FM.SDradius']
    bmd_fm_var = dict_sos['-AER.BMD.FM.SDvar']
    ficextdata_aer = dict_sos['-AER.ExtData']
    ficmixture_aer = dict_sos['-AER.DefMixture']
    ficuser_aer = dict_sos['-AER.UserFile']
    ficprofil_log = dict_sos['-AP.Log']
    tr = dict_sos['-AP.MOT']
    hr = dict_sos['-AP.HR']
    ha = dict_sos['-AP.AerHS.HA']
    iprofil = dict_sos['-AP.AerProfile.Type']
    zmin = dict_sos['-AP.AerLayer.Zmin']
    zmax = dict_sos['-AP.AerLayer.Zmax']
    psurf = dict_sos['-AP.Psurf']
    h2o = dict_sos['-AP.H2O']
    o3 = dict_sos['-AP.O3']
    co2 = dict_sos['-AP.CO2']
    ch4 = dict_sos['-AP.CH4']
    absprofil = dict_sos['-AP.AbsProfile.Type']
    ficabsprofil = dict_sos['-AP.AbsProfile.UserFile']
    nustep = dict_sos['-AP.SpectralResol']
    isurf = dict_sos['-SURF.Type']
    dir_surf = dict_sos['-SURF.Dir']
    ficsurf_log = dict_sos['-SURF.Log']
    surf_ind = dict_sos['-SURF.Ind']
    wind = dict_sos['-SURF.Glitter.Wind']
    k0_roujean = dict_sos['-SURF.Roujean.K0']
    k1_roujean = dict_sos['-SURF.Roujean.K1']
    k2_roujean = dict_sos['-SURF.Roujean.K2']
    alpha_nada = dict_sos['-SURF.Nadal.Alpha']
    beta_nada = dict_sos['-SURF.Nadal.Beta']
    coef_c_maignan = dict_sos['-SURF.Maignan.C']
    rho = dict_sos['-SURF.Alb']
    ficsurf = dict_sos['-SURF.File']
    ficsos_log = dict_sos['-SOS.Log']
    ficsos_res_bin = dict_sos['-SOS.ResBin']
    # ficup = dict_sos['-SOS.ResFileUp']
    # ficdown = dict_sos['-SOS.ResFileDown']
    # ficup_angles_user = dict_sos['-SOS.ResFileUp.UserAng']
    # ficdown_angles_user = dict_sos['-SOS.ResFileDown.UserAng']
    fictrans = dict_sos['-SOS.Trans']
    ficflux = dict_sos['-SOS.Flux']
    zout = dict_sos['-SOS.OutputAlt']
    igmax = dict_sos['-SOS.IGmax']
    ipolar = dict_sos['-SOS.Ipolar']
    itrphi = dict_sos['-SOS.View']
    phios = dict_sos['-SOS.View.Phi']
    pas_phi = dict_sos['-SOS.View.Dphi']
    imode_ckd_calcul = dict_sos['-SOS.AbsModeCKD']
    ier = 0  # initialisation de l'indicateur d'erreur de traitement pour la
    # procédure
    trace = trace

    return resroot, ficmain_log, wa_simu, nbmu_gauss_lum, \
           ficangles_user_lum, tetas, ficangles_res_lum, nbmu_gauss_mie, \
           ficangles_user_mie, ficangles_res_mie, ficanglog, waref_aot, \
           aot_ref, itronc_aer, ficgranu_log, ficmie_log, dir_mie, ficgranu,\
           imod_aer, rn_wa, in_wa, rn_waref, in_waref, igranu, \
           lnd_radius_mmd_aer, lnd_lnvar_mmd_aer, jd_slope_mmd_aer, \
           jd_rmin_mmd_aer, jd_rmax_mmd_aer, imodel_wmo, c_wmo_dl, c_wmo_ws,\
           c_wmo_oc, c_wmo_so, imodele_sf, rh, mode_param_bilnd, \
           user_cv_coarse, user_cv_fine, rtauct_waref, bmd_cm_mrwa, \
           bmd_cm_miwa, bmd_cm_mrwaref, bmd_cm_miwaref, bmd_cm_rmodal, \
           bmd_cm_var, bmd_fm_mrwa, bmd_fm_miwa, bmd_fm_mrwaref, \
           bmd_fm_miwaref, bmd_fm_rmodal, bmd_fm_var, ficextdata_aer, \
           ficmixture_aer, ficuser_aer, ficprofil_log, tr, hr, ha, iprofil, \
           zmin, zmax, psurf, h2o, o3, co2, ch4, absprofil, ficabsprofil, \
           nustep, isurf, dir_surf, ficsurf_log, surf_ind, wind, k0_roujean, \
           k1_roujean, k2_roujean, alpha_nada, beta_nada, coef_c_maignan, \
           rho, ficsurf, ficsos_log, ficsos_res_bin, fictrans, ficflux,\
           zout, \
           igmax, ipolar, itrphi, phios, pas_phi, imode_ckd_calcul, \
           ier, trace


if __name__ == '__main__':
    
    
    ROOT="./"
    REP_RES = os.path.join(ROOT,'RESULTS_TEST_BINDING/')
    print(REP_RES)
    REP_MIE = os.path.join(REP_RES,'MIE')
    REP_SURF = os.path.join(REP_RES,'SURFACE')
    REP_DATA = "./TestFiles"
    FIC_ANG_LUM = os.path.join(REP_DATA,'ficAngRad.txt')
    print(FIC_ANG_LUM)
    FIC_ANG_AER = os.path.join(REP_DATA,'ficAngAer.txt')
    print(FIC_ANG_AER)
    
    
    dict_sos_all_parameters = {
        "-SOS_Main.ResRoot": "",
        "-SOS_Main.Log": "SOS_Main.Log",
        "-SOS_Main.Wa": SOS_NOT_DEFINED_VALUE_DBLE,
        "-ANG.Rad.NbGauss": NBMU_GAUSS_LUM,
        "-ANG.Rad.UserAngFile": FICANGLES_USER_LUM,
        "-ANG.Thetas": SOS_NOT_DEFINED_VALUE_DBLE,
        "-ANG.Rad.ResFile": FICANGLES_RES_LUM,
        "-ANG.Aer.NbGauss": NBMU_GAUSS_MIE,
        "-ANG.Aer.UserAngFile": FICANGLES_USER_MIE,
        "-ANG.Aer.ResFile": FICANGLES_RES_MIE,
        "-ANG.Log": "Angles.Log",
        "-AER.Waref": SOS_NOT_DEFINED_VALUE_DBLE,
        "-AER.AOTref": SOS_NOT_DEFINED_VALUE_DBLE,
        "-AER.Tronca": ITRONC_AER,
        "-AER.Log": FICGRANU_LOG,
        "-AER.MieLog": FICMIE_LOG,
        "-AER.DirMie": "",
        "-AER.ResFile": FICGRANU,
        "-AER.Model": SOS_NOT_DEFINED_VALUE_INT,
        "-AER.MMD.MRwa": RN,
        "-AER.MMD.MIwa": IN,
        "-AER.MMD.MRwaref": RN_WAREF,
        "-AER.MMD.MIwaref": IN_WAREF,
        "-AER.MMD.SDtype": IGRANU,
        "-AER.MMD.LNDradius": LND_RADIUS_MMD_AER,
        "-AER.MMD.LNDvar": LND_LNVAR_MMD_AER,
        "-AER.MMD.JD.slope": JD_SLOPE_MMD_AER,
        "-AER.MMD.JD.rmin": JD_RMIN_MMD_AER,
        "-AER.MMD.JD.rmax": JD_RMAX_MMD_AER,
        "-AER.WMO.Model": IModele_WMO,
        "-AER.WMO.DL": C_WMO_DL,
        "-AER.WMO.WS": C_WMO_WS,
        "-AER.WMO.OC": C_WMO_OC,
        "-AER.WMO.SO": C_WMO_SO,
        "-AER.SF.Model": IModele_SF,
        "-AER.SF.RH": RH,
        "-AER.BMD.VCdef": MODE_PARAM_BILND,
        "-AER.BMD.CoarseVC": USER_CV_COARSE,
        "-AER.BMD.FineVC": USER_CV_FINE,
        "-AER.BMD.RAOT": RTAUct_WAREF,
        "-AER.BMD.CM.MRwa": BMD_CM_MRWA,
        "-AER.BMD.CM.MIwa": BMD_CM_MIWA,
        "-AER.BMD.CM.MRwaref": BMD_CM_MRWAREF,
        "-AER.BMD.CM.MIwaref": BMD_CM_MIWAREF,
        "-AER.BMD.CM.SDradius": BMD_CM_RMODAL,
        "-AER.BMD.CM.SDvar": BMD_CM_VAR,
        "-AER.BMD.FM.MRwa": BMD_FM_MRWA,
        "-AER.BMD.FM.MIwa": BMD_FM_MIWA,
        "-AER.BMD.FM.MRwaref": BMD_FM_MRWAREF,
        "-AER.BMD.FM.MIwaref": BMD_FM_MIWAREF,
        "-AER.BMD.FM.SDradius": BMD_FM_RMODAL,
        "-AER.BMD.FM.SDvar": BMD_FM_VAR,
        "-AER.ExtData": FICEXTDATA_AER,
        "-AER.DefMixture": FICMIXTURE_AER,
        "-AER.UserFile": FICUSER_AER,
        "-AP.Log": "Aerosols.Log",
        "-AP.MOT": TR,
        "-AP.HR": HR,
        "-AP.AerHS.HA": HA,
        "-AP.AerProfile.Type": 1,
        "-AP.AerLayer.Zmin": ZMIN,
        "-AP.AerLayer.Zmax": ZMAX,
        "-AP.Psurf": PSURF,
        "-AP.H2O": H2O,
        "-AP.O3": O3,
        "-AP.CO2": CO2,
        "-AP.CH4": CH4,
        "-AP.AbsProfile.Type": SOS_NOT_DEFINED_VALUE_INT,
        "-AP.AbsProfile.UserFile": FICABSPROFIL,
        "-AP.SpectralResol": SOS_NOT_DEFINED_VALUE_INT,
        "-SURF.Type": 0,
        "-SURF.Dir": "",
        "-SURF.Log": FICSURF_LOG,
        "-SURF.Ind": SURF_IND,
        "-SURF.Glitter.Wind": WIND,
        "-SURF.Roujean.K0": K0_ROUJEAN,
        "-SURF.Roujean.K1": K1_ROUJEAN,
        "-SURF.Roujean.K2": K2_ROUJEAN,
        "-SURF.Nadal.Alpha": ALPHA_NADAL,
        "-SURF.Nadal.Beta": BETA_NADAL,
        "-SURF.Maignan.C": COEF_C_MAIGNAN,
        "-SURF.Alb": SOS_NOT_DEFINED_VALUE_DBLE,
        "-SURF.File": FICSURF,
        "-SOS.Log": "SOS.Log",
        "-SOS.ResBin": SOS_DEFAULT_RESBIN,
        "-SOS.ResFileUp": FICUP,
        "-SOS.ResFileDown": FICDOWN,
        # "-SOS.ResFileUp.UserAng": FICUP_ANGLES_USER,
        # "-SOS.ResFileDown.UserAng": FICDOWN_ANGLES_USER,
        "-SOS.Trans": FICTRANS,
        "-SOS.Flux": 'FicFlux.txt',
        "-SOS.OutputAlt": ZOUT,
        "-SOS.IGmax": 200,
        "-SOS.Ipolar": IPOLAR,
        "-SOS.View": 2,
        "-SOS.View.Phi": PHIOS,
        "-SOS.View.Dphi": PAS_PHI,
        "-SOS.MDF": 0.0279,
        "-SOS.AbsModeCKD": 1
    }

    dict_sos_user = {
        '-SOS_Main.Wa': 0.440,
        '-SOS_Main.ResRoot': REP_RES,
        '-SOS_Main.Log': 'SOS_Main.Log',
        '-AER.DirMie': REP_MIE,
        '-AP.AbsProfile.Type': 2,
        '-AP.SpectralResol': 10.0,
        '-SOS.AbsModeCKD': 2,
        '-ANG.Thetas':  40.00,
        '-ANG.Rad.NbGauss':  24,
        '-ANG.Aer.NbGauss':  40,
        '-ANG.Log': 'Angles.Log',
        '-AP.MOT': 0.230,
        '-AP.AerProfile.Type': 1,
        '-AP.HR':  8,
        '-AP.AerHS.HA': 2,
        '-AP.H2O': 2.906,
        '-AP.O3': 335,
        '-AP.Psurf': 1013,
        '-AER.Waref':  0.55,
        '-AER.AOTref': 0.300,
        '-AER.ResFile': 'Aerosols.txt',
        '-AER.Log': 'Aerosols.Log',
        '-AER.MieLog': 'Mie.Log',
        '-SURF.Alb':  0.000000,
        '-AER.Tronca':  1,
        '-AER.Model': 1,
        '-AER.WMO.Model': 1,
        '-SURF.Dir': REP_SURF,
        '-SURF.Type': 1,
        '-SURF.Glitter.Wind': 2.0,
        ' -SURF.Alb': 0.02,
        '-SURF.Ind': 1.33,
        '-SURF.Log': 'Surface.Log',
        '-SOS.View': 1,
        '-SOS.View.Phi': 35,
        '-SOS.IGmax':  30,
        # '-SOS.ResFileUp': 'SOS_Up.txt',
        # '-SOS.ResFileDown': 'SOS_Down.txt',
        '-SOS.ResBin': 'SOS_Result.bin',
        '-SOS.Log': 'SOS.Log',
        '-SOS.Flux': 'FicFlux.txt',
        '-SOS.Trans': 'SOS_Transm.txt',
    }

    dict_sos_updated = dict_sos_all_parameters.copy()
    for k, v in dict_sos_user.items():
        if k in dict_sos_all_parameters.keys():
            dict_sos_updated[k] = v

    RESROOT, FICMAIN_LOG, WA_SIMU, NBMU_GAUSS_LUM, FICANGLES_USER_LUM, \
    TETAS, FICANGLES_RES_LUM, NBMU_GAUSS_MIE, FICANGLES_USER_MIE, \
    FICANGLES_RES_MIE, FICANGLOG, WAREF_AOT, AOT_REF, ITRONC_AER, \
    FICGRANU_LOG, FICMIE_LOG, DIR_MIE, FICGRANU, IMOD_AER, RN, IN, RN_WAREF,\
    IN_WAREF, IGRANU, LND_RADIUS_MMD_AER, LND_LNVAR_MMD_AER, \
    JD_SLOPE_MMD_AER, JD_RMIN_MMD_AER, JD_RMAX_MMD_AER, IModele_WMO, \
    C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO, IModele_SF, RH, \
    MODE_PARAM_BILND, USER_CV_COARSE, USER_CV_FINE, RTAUct_WAREF, \
    BMD_CM_MRWA, BMD_CM_MIWA, BMD_CM_MRWAREF, BMD_CM_MIWAREF, BMD_CM_RMODAL,\
    BMD_CM_VAR, BMD_FM_MRWA, BMD_FM_MIWA, BMD_FM_MRWAREF, BMD_FM_MIWAREF, \
    BMD_FM_RMODAL, BMD_FM_VAR, FICEXTDATA_AER, FICMIXTURE_AER, FICUSER_AER, \
    FICPROFIL_LOG, TR, HR, HA, IPROFIL, ZMIN, ZMAX, PSURF, H2O, O3, CO2, \
    CH4, ABSPROFIL, FICABSPROFIL, NUSTEP, ISURF, DIR_SURF, FICSURF_LOG, \
    SURF_IND, WIND, K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN, ALPHA_NADAL, \
    BETA_NADAL, COEF_C_MAIGNAN, RHO, FICSURF, FICSOS_LOG, FICSOS_RES_BIN, \
    FICTRANS, FICFLUX, \
    ZOUT, IGMAX, IPOLAR, ITRPHI, PHIOS, PAS_PHI, \
    IMODE_CKD_CALCUL, IER, TRACE = set_sos_params(
        dict_sos_updated, trace=True)


#    nblum, ind_angout, phi, vza, sca_ang_up, i_up, q_up, u_up, pol_ang_up, pol_rate_up, l_pol_up, \
#    sca_ang_down, i_down, q_down, u_down, pol_ang_down, pol_rate_down, l_pol_down, \
#    flux_dir_down, flux_diff_down, flux_tot_down, flux_diff_up, \
#    coef_tronca = \
    nblum, ind_angout, phi, vza, sca_ang_up, i_up, q_up, u_up, pol_ang_up, pol_rate_up, l_pol_up, \
    sca_ang_down, i_down, q_down, u_down, pol_ang_down, pol_rate_down, l_pol_down, \
    flux_dir_down, flux_diff_down, flux_tot_down, flux_diff_up, \
    coef_tronca = \
        sos.sos_proc(resroot=RESROOT, ficmain_log=FICMAIN_LOG,
                          wa_simu=WA_SIMU, nbmu_gauss_lum=NBMU_GAUSS_LUM,
                          ficangles_user_lum=FICANGLES_USER_LUM,
                          tetas=TETAS, ficangles_res_lum=FICANGLES_RES_LUM,
                          nbmu_gauss_mie=NBMU_GAUSS_MIE,
                          ficangles_user_mie=FICANGLES_USER_MIE,
                          ficangles_res_mie=FICANGLES_RES_MIE,
                          ficanglog=FICANGLOG,
                          waref_aot=WAREF_AOT, aot_ref=AOT_REF,
                          itronc_aer=ITRONC_AER,
                          ficgranu_log=FICGRANU_LOG, ficmie_log=FICMIE_LOG,
                          dir_mie=DIR_MIE,
                          ficgranu=FICGRANU, imod_aer=IMOD_AER, rn_wa=RN,
                          in_wa=IN, rn_waref=RN_WAREF, in_waref=IN_WAREF,
                          igranu=IGRANU,
                          lnd_radius_mmd_aer=LND_RADIUS_MMD_AER,
                          lnd_lnvar_mmd_aer=LND_LNVAR_MMD_AER,
                          jd_slope_mmd_aer=JD_SLOPE_MMD_AER,
                          jd_rmin_mmd_aer=JD_RMIN_MMD_AER,
                          jd_rmax_mmd_aer=JD_RMAX_MMD_AER,
                          imodele_wmo=IModele_WMO, c_wmo_dl=C_WMO_DL,
                          c_wmo_ws=C_WMO_WS, c_wmo_oc=C_WMO_OC,
                          c_wmo_so=C_WMO_SO, imodele_sf=IModele_SF, rh=RH,
                          mode_param_bilnd=MODE_PARAM_BILND,
                          user_cv_coarse=USER_CV_COARSE,
                          user_cv_fine=USER_CV_FINE, rtauct_waref=RTAUct_WAREF,
                          bmd_cm_mrwa=BMD_CM_MRWA, bmd_cm_miwa=BMD_CM_MIWA,
                          bmd_cm_mrwaref=BMD_CM_MRWAREF,
                          bmd_cm_miwaref=BMD_CM_MIWAREF,
                          bmd_cm_rmodal=BMD_CM_RMODAL,
                          bmd_cm_var=BMD_CM_VAR, bmd_fm_mrwa=BMD_FM_MRWA,
                          bmd_fm_miwa=BMD_FM_MIWA,
                          bmd_fm_mrwaref=BMD_FM_MRWAREF,
                          bmd_fm_miwaref=BMD_FM_MIWAREF,
                          bmd_fm_rmodal=BMD_FM_RMODAL,
                          bmd_fm_var=BMD_FM_VAR,
                          ficextdata_aer=FICEXTDATA_AER,
                          ficmixture_aer=FICMIXTURE_AER,
                          ficuser_aer=FICUSER_AER,
                          ficprofil_log=FICPROFIL_LOG, tr=TR, hr=HR, ha=HA,
                          iprofil=IPROFIL, zmin=ZMIN, zmax=ZMAX,
                          psurf=PSURF, h2o=H2O, o3=O3, co2=CO2, ch4=CH4,
                          absprofil=ABSPROFIL, ficabsprofil=FICABSPROFIL,
                          nustep=NUSTEP, isurf=ISURF, dir_surf=DIR_SURF,
                          ficsurf_log=FICSURF_LOG, surf_ind=SURF_IND,
                          wind=WIND, k0_roujean=K0_ROUJEAN,
                          k1_roujean=K1_ROUJEAN, k2_roujean=K2_ROUJEAN,
                          alpha_nadal=ALPHA_NADAL, beta_nadal=BETA_NADAL,
                          coef_c_maignan=COEF_C_MAIGNAN, rho=RHO,
                          ficsurf=FICSURF, ficsos_log=FICSOS_LOG,
                          ficsos_res_bin=FICSOS_RES_BIN, fictrans=FICTRANS,
                          ficflux=FICFLUX,
                          zout=ZOUT, igmax=IGMAX,
                          ipolar=IPOLAR, itrphi=ITRPHI, phios=PHIOS,
                          pas_phi=PAS_PHI, imode_ckd_calcul=IMODE_CKD_CALCUL,
                          ier=IER, trace=TRACE)
#    set_trace()
#    a=1

print(i_up.shape)

print("Plan retrodiff")
print(i_up[0,:])

print("Plan diff avant")
print(i_up[1,:])

########Generation du fichier SOS_Up.txt
repOut=os.path.join(REP_RES,'SOS')
gen_sos_output(repOut,ITRPHI,1,ZMAX,nblum,PAS_PHI,phi,vza,sca_ang_up,i_up,q_up,u_up,pol_ang_up,pol_rate_up,l_pol_up)

########Generation du fichier SOS_Down.txt
repOut=os.path.join(REP_RES,'SOS')
gen_sos_output(repOut,ITRPHI,2,ZMAX,nblum,PAS_PHI,phi,vza,sca_ang_down,i_down,q_down,u_down,pol_ang_down,pol_rate_down,l_pol_down)

########Flux
print("FLux direct descendant : %s" %(flux_dir_down))
print("FLux diffus descendant : %s" %(flux_diff_down))
print("FLux total descendant : %s" %(flux_tot_down))
print("FLux diffus montant : %s" %(flux_diff_up))
print("coef troncature fonction de phase : %s" %(coef_tronca))

