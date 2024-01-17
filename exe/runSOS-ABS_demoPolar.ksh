#!/bin/ksh
#################################################################
# FILE: runOS-ABS_demo.ksh
#
# PROJECT: Successive Orders of Scattering 
#
# RULE: Example of command file for a SOS simulation with absorption
#       -> Define the simulation parameters (-keyword -value).
#       -> Call the SOS_ABS_MAIN executable.
#    
# This shell is available for the SOS-ABS code version 5.1
#                                 ------------------------ 
#
#    
# In order to use this shell, the environment variable SOS_ABS_ROOT
# has to be previously defined :                       ------------
#       SOS_ABS_ROOT: complete access path to the software directory 
#
#
# AUTHOR: Bruno LAFRANCE ( CS Group)
# DATE: 2022 / 06 / 24
#


####################################################################
#                                                                  #
#  Environment variables                                           #
#                                                                  #
####################################################################

##----------------------------------------------------------
# SOS_ABS_ROOT : the path to the software directory is supposed to be defined
##----------------------------------------------------------

##-----------------------------------------------------
# dirRESULT : access path to the results storage directory 
##-----------------------------------------------------
export dirRESULT=${SOS_ABS_ROOT}/DEMO_RESULTS

##-----------------------------------------------------
# Storage directory of Mie aerosols files 
##-----------------------------------------------------
export dirMIE=$dirRESULT/MIE

##-----------------------------------------------------
# Storage directory of Surface BRDF/BPDF files 
##-----------------------------------------------------
export dirSURF=$dirRESULT/BRDF_BPDF


####################################################################
#                                                                  #
#  Parameters of simulation                                        #
#                                                                  #
####################################################################
#
# Example 1 :
# ----------
#
# Simulation of radiance over sea at 910 nm (water vapour absorption):
#
#     The molecular optical thickness is calculated from a surface pressure of 1013 mb
#
#     The aerosols correspond to the Maritime WMO Model.
#     The Aerosol optical thickness is 0.3 at 550 nm.
#
#     Aerosols and molecules profiles of optical thicknesses are defined by 
#     an exponential decrease with scale heights : 2 and 8 km, respectively.
#
#     The gas absorption is calculated using the tropical profile for a 10 cm-1 spectral resolution of CKD coefficients,
#     with the fine CKD modelling of the absorption / scattering coupling.
#
#     Wind velocity of 2 m/s
#     Water / atmosphere refractive index of 1.34
#     Sea albedo is 0.00 (at 910 nm)
#
#     The solar zenithal angle is 35° .
#     The radiance is calculated in the solar principal plane (relative azimuth angle is null).
#     40 Gauss angles are used to describe the radiance field.

	 
   $SOS_ABS_ROOT/exe/SOS_ABS_MAIN.exe -SOS_Main.Wa 0.910 -SOS_Main.ResRoot ${dirRESULT} -SOS_Main.Log SOS_Main_Demo.Log \
          -ANG.Rad.NbGauss 40  -ANG.Thetas 35. -SOS.View 2 -SOS.View.Dphi 30  \
	  -SOS.ResFileUp SOS_Up_DemoPolar.txt  -SOS.ResFileDown SOS_Down_DemoPolar.txt\
          -AP.Log Profile_Demo.Log \
	  -AP.Psurf 1013  -AP.AerProfile.Type 1 -AP.HR 8.0 -AP.AerHS.HA 2.0 	\
	  -AP.AbsProfile.Type 1 -AP.SpectralResol 10. -SOS.AbsModeCKD 1 \
	  -AER.DirMie ${dirMIE} -AER.Model 1 -AER.WMO.Model 2  -AER.Waref 0.550  -AER.AOTref 0.3 \
	  -AER.ResFile Aerosols_Demo.txt -AER.Log Aerosols_Demo.Log  \
          -SURF.Dir ${dirSURF}  -SURF.Type 1 -SURF.Alb 0.00 -SURF.Ind 1.34 -SURF.Glitter.Wind 2.0  	 


#
