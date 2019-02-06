#!/bin/ksh
#################################################################
# FILE: SOS.ksh
#
# PROJECT: Successive Orders of Scattering.
#
# RULE: Command file for a SOS simulation.
#       -> Definition of simulation parameters (-keyword -value).
#       -> Call of the main_SOS.ksh shell.
#    
# This shell is available for the SOS code version 5.0,  
# not for previous ones.          -------------------- 
#
#    
# In order to use this shell, the environment variable RACINE
# has to be previously defined :                       ------
#       RACINE: complete access path to the software directory 
#
#
# AUTHOR: Bruno Lafrance ( CS )
# DATE: 2009 / 11 / 12
#
# MOD:VERSION:4.0: main adjustement for a parameter definition by
#                  couples : -Keyword Value
#
# MOD:VERSION:4.1: :  2015 / 12 / 18
#      - Correction of the comments for the -AP.HR parameter.
#
#      - Deleting the non relevant specified formats 
#        for the values of parameters. 
#
#      - Clarification for LND size distributions.
#        Now, they are based on the common definition :
#            N(r) = exp(-0.5 *(ln(r/rm)/sig)**2) / r
#          
#        (   rather than previously 
#            N(r) = exp(-0.5 *(log10(r/rm)/sig_old)**2) / r  )
#
#        This impacts the input value for the parameter sig
#        which must be changed by a factor ln(10) : 
#        sig = sig_old x ln(10)
#        
#        Concerned parameters are:
#           -AER.MMD.SDparam2 :  Standard deviation of the LND 
#                                for case -AER.MMD.SDtype 1
#           -AER.BMD.CM.SDvar : Standard deviation of the 
#        			"LND coarse mode"
#           -AER.BMD.FM.SDvar : Standard deviation of the 
#        			"LND fine mode"
#
#################################################################



####################################################################
#                                                                  #
#  Environment variables                                           #
#                                                                  #
####################################################################

##-----------------------------------------------------
# SOS_RACINE : access path to the software directory  
# --> SOS_RACINE/exe contains the compilated codes 
##-----------------------------------------------------
export SOS_RACINE=$RACINE

##-----------------------------------------------------
# SOS_RESULT : access path to the results storage directory 
##-----------------------------------------------------
export SOS_RESULT=$RACINE/SOS_TEST

##----------------------------------------------------------
# SOS_RACINE_FIC : access path to WMO and Shettle&Fenn files 
##----------------------------------------------------------
export SOS_RACINE_FIC=$RACINE/fic

##-----------------------------------------------------
# Storage directory of BRDF and BPDF files 
##-----------------------------------------------------

  ## Storage directory of reflection files to Cox & Munk's Sun glint model
  ##----------------------------------------------------------------------
  export dirSUNGLINT=$SOS_RESULT/SURFACE/GLITTER

  ## Storage directory of reflection files to Roujean's BRDF model
  ##--------------------------------------------------------------
  export dirROUJEAN=$SOS_RESULT/SURFACE/ROUJEAN

  ## Storage directory of reflection files to Rondeaux & Herman's BPDF model
  ##-------------------------------------------------------------------------
  export dirRH=$SOS_RESULT/SURFACE/RH

  ## Storage directory of reflection files to Breon's BPDF model
  ##------------------------------------------------------------
  export dirBREON=$SOS_RESULT/SURFACE/BREON

  ## Storage directory of reflection files to Nadal's BPDF model
  ##------------------------------------------------------------
  export dirNADAL=$SOS_RESULT/SURFACE/NADAL


##-----------------------------------------------------
# Storage directory of Mie aerosols files 
##-----------------------------------------------------
  export dirMIE=$SOS_RESULT/MIE

##-----------------------------------------------------
# Storage directory of log files 
##-----------------------------------------------------
  export dirLOG=$SOS_RESULT/LOG

##----------------------------------------------------------------------------------------------
# Storage directory of simulation results (profiles, aerosols parameters, radiance simulations)
##----------------------------------------------------------------------------------------------
  export dirRESULTS=$SOS_RESULT/SOS



####################################################################
#                                                                  #
#  Parametres de simulation                                        #
#                                                                  #
####################################################################


##################################################################
##################################################################
###            CALL OF THE main_SOS.ksh SHELL
##################################################################
##################################################################
#
# List of all available -Keyword parameter
# (see the end of this script to a complete description of -Keyword parameters)
#
#   Angles calculation parameters :
#   --------------------------------  
#        -ANG.Thetas	\
#        -ANG.Rad.NbGauss -ANG.Rad.UserAngFile -ANG.Rad.ResFile \
#        -ANG.Aer.NbGauss -ANG.Aer.UserAngFile -ANG.Aer.ResFile \
#        -ANG.Log
#
#   Radiance calculation parameters :
#   --------------------------------
#        -SOS.Wa \
#	 -SOS.View  -SOS.View.Phi  -SOS.View.Dphi  \
#        -SOS.OutputLevel   -SOS.IGmax  -SOS.Ipolar -SOS.MDF \
#        -SOS.ResBin  -SOS.ResFileUp  -SOS.ResFileDown  -SOS.Log  \
#        -SOS.ResFileUp.UserAng  -SOS.ResFileDown.UserAng \
#        -SOS.Trans -SOS.Config
#
#   Atmospheric Profile parameters :
#   --------------------------------
#        -AP.ResFile  -AP.UserFile -AP.Log  \
#	 -AP.MOT   \
#        -AP.Type  -AP.HR  -AP.AerHS.HA  -AP.AerLayer.Zmin  -AP.AerLayer.Zmax 	\
#
#   Aerosols parameters :
#   ---------------------
#	 -AER.Waref  -AER.AOTref  \
#	 -AER.ResFile -AER.UserFile -AER.Log  -AER.MieLog  \
#	 -AER.Tronca  -AER.Model  \
#     Aerosols parameters for mono-modal models :
#        -AER.MMD.Mie.Filename   -AER.MMD.Mie.AlphaMax -AER.MMD.MRwa -AER.MMD.MIwa  -AER.MMD.MRwaref -AER.MMD.MIwaref  \
#	 -AER.MMD.SDtype  -AER.MMD.SDparam1  -AER.MMD.SDparam2 	\
#     Aerosols parameters for WMO models :
#        -AER.WMO.Model -AER.WMO.DL -AER.WMO.WS -AER.WMO.OC  -AER.WMO.SO  \
#     Aerosols parameters for Shettle&Fenn models :
#	 -AER.SF.Model  -AER.SF.RH  \
#     Aerosols parameters for LND bi-modal models :
#        -AER.BMD.VCdef -AER.BMD.CoarseVC  -AER.BMD.FineVC  -AER.BMD.RAOT  \
#	 -AER.BMD.CM.MRwa -AER.BMD.CM.MIwa -AER.BMD.CM.MRwaref -AER.BMD.CM.MIwaref -AER.BMD.CM.SDradius -AER.BMD.CM.SDvar 	\
#        -AER.BMD.FM.MRwa -AER.BMD.FM.MIwa -AER.BMD.FM.MRwaref -AER.BMD.FM.MIwaref -AER.BMD.FM.SDradius -AER.BMD.FM.SDvar 	\
#     Aerosols parameters for external data (phaze functions, scattering and extinction coefficients) :
#	 -AER.ExtData
#
#   Surface parameters :
#   --------------------
#         -SURF.Log  -SURF.File  \
#         -SURF.Type -SURF.Alb  -SURF.Ind  \
#         -SURF.Glitter.Wind \
#         -SURF.Roujean.K0 -SURF.Roujean.K1  -SURF.Roujean.K2  \
#         -SURF.Nadal.Alpha -SURF.Nadal.Beta 
#
#
#
#   In order to perform a simulation  :
#       $RACINE/exe/main_SOS.ksh -Keyword Value   -Keyword Value -Keyword Value \
#                                -Keyword Value   -Keyword Value -Keyword Value
#
#
# Example 1 :
# ----------
#
# Simulation of radiance over sea for :
#
#     The radiance simulation is applied at 440 nm
#     --> The molecular optical thickness is 0.230
#
#     Wind velocity of 2 m/s
#     Water / atmosphere refractive index of 1.33
#     Sea albedo is 0.02 (at 440 nm)
#
#     Aerosols corresponding to the Maritime WMO Model
#     Aerosol optical thickness is 0.3 at 550 nm.
#
#     Aerosols and molecules profiles of optical thicknesses
#     are defined by exponential decrease wirh scale heights : 2 and 8 km.
#
#     The solar zenithal angle is 32.48° which correspond to the angle number 9 in file MU24 (cos = 0.84358826162439).
#     The radiance is calculated in the solar principal plane (relative azimuth angle is null).



     $RACINE/exe/main_SOS.ksh -SOS.Wa 0.440	\
         -ANG.Rad.NbGauss 24    -ANG.Rad.ResFile ${dirRESULTS}/SOS_UsedAngles.txt \
	 -ANG.Aer.NbGauss 40    -ANG.Aer.ResFile ${dirRESULTS}/AER_UsedAngles.txt \
	 -ANG.Log ${dirLOG}/Angles.Log \
         -ANG.Thetas 32.48 -SOS.View 1 -SOS.View.Phi 0.  \
         -SOS.IGmax 30 \
	 -SOS.ResFileUp ${dirRESULTS}/SOS_Up.txt -SOS.ResFileDown ${dirRESULTS}/SOS_Down.txt \
         -SOS.ResBin ${dirRESULTS}/SOS_Result.bin  -SOS.Log ${dirLOG}/SOS.Log \
	 -SOS.Config ${dirRESULTS}/SOS_config.txt \
	 -SOS.Trans ${dirRESULTS}/SOS_transm.txt \
         -AP.ResFile ${dirRESULTS}/Profile.txt -AP.Log ${dirLOG}/Profile.Log \
	 -AP.MOT 0.230 -SOS.MDF 0.0279 \
         -AP.Type 1 -AP.HR 8.0 -AP.AerHS.HA 2.0 	\
	 -AER.Waref 0.550 -AER.AOTref 0.300 \
	 -AER.ResFile ${dirRESULTS}/Aerosols.txt -AER.Log ${dirLOG}/Aerosols.Log -AER.MieLog 0 \
	 -AER.Tronca 1 -AER.Model 1 \
         -AER.WMO.Model 2 \
         -SURF.Log ${dirLOG}/Surface.Log  -SURF.File DEFAULT \
         -SURF.Type 1 -SURF.Alb 0.02 -SURF.Ind 1.33 \
         -SURF.Glitter.Wind 2.0 
	 
	 

#-----------------------------------------------------------------
#  LIST OF POSSIBLE ARGUMENTS 
#-----------------------------------------------------------------
#    Use of couples (-Keyword Value) for arguments definition
#
##    Keyword  : -ANG.Thetas  
#    ----------------------   
#       Status : Required    
#       Value :  Solar zenithal angle (degrees) 
#       Value format : Encoding as float (10 characters max).		
#    
#    Keyword  : -ANG.Rad.NbGauss  
#    ---------------------------   
#       Status : Optional (if not defined, the default value SOS_DEFAULT_NBMU_LUM is used : see SOS.h )       
#       Value :  Number of Gauss angles to be used for radiance computations 
#       Value format : Encoding as integer (4 characters max)
#
#    Keyword  : -ANG.Rad.UserAngFile  
#    -------------------------------   
#       Status : Optional      
#       Value :  Filename of the complementary list of user's angles to complete the 
#                ANG.Rad.NbGauss angles (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -ANG.Rad.ResFile  
#    ---------------------------   
#       Status : Required     
#       Value :  Filename of list of angles used to BRDF/BPDF and radiance computations. 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -ANG.Aer.NbGauss  
#    ---------------------------   
#       Status : Optional (if not defined, the default value SOS_DEFAULT_NBMU_MIE is used : see SOS.h )       
#       Value :  Number of Gauss angles to be used for phaze matrix computations 
#       Value format : Encoding as integer (4 characters max)
#
#    Keyword  : -ANG.Aer.UserAngFile 
#    -------------------------------   
#       Status : Optional      
#       Value :  Filename of the complementary list of user's angles to complete the 
#                ANG.Aer.NbGauss angles (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -ANG.Aer.ResFile  
#    ----------------------------   
#       Status : Required     
#       Value :  Filename of list of angles used to the matrix phase function computations. 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -ANG.Log
#    -------------------   
#       Status : Required 
#       Value :  log filename for SOS computations (complete path).
#      		 Assign 0 to avoid generating the log file. 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted
#
#    Keyword  : -SOS.Wa  
#    ------------------  
#       Status : Required      
#       Value :  Wavelength of radiance calculation (microns).
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -SOS.MDF  
#    -------------------   
#       Status : Required      
#       Value :  Factor of molecular depolarization.
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -SOS.View  
#    --------------------   
#       Status : Required      
#       Value :  Index for output type
#                  1 : radiance for viewing directions in a constante azimutal plane. 
#                  2 : radiance for viewing directions covering half the space, with a constant step on azimut.  
#       Value format : Encoding as integer (10 characters max).   	 	 	
#		 	 	 
#    Keyword  : -SOS.View.Phi  
#    ------------------------   
#       Status : Associated to -SOS.View 1      
#       Value :  Relative azimut angle (degrees).  
#       Value format : Encoding as float (10 characters max).   	 	 	
#
#    Keyword  : -SOS.View.Dphi  
#    -------------------------   
#       Status : Associated to -SOS.View 2      
#       Value :  Constant step on azimut (degrees). Integer value.          
#
#    Keyword  : -SOS.OutputLevel  
#    ---------------------------   
#       Status : Optionnal   
#       Value :  Profile level for radiance results.
#                  If not defined : -1 is the default value attributed
#                  -1 : standard results 
#			--> Upward radiance at TOA
#                       --> Downward radiance at surface level
#                  n  : results for the level n of the atmospheric profile
#			(0 <= n <= SOS_OS_NT)
#			--> Upward and downward radiance are given for the
#			    level N of the profile.
#       Value format : Encoding as integer (4 characters max).
#			
#    Keyword  : -SOS.Ipolar 
#    ----------------------   
#       Status : Optionnal   
#       Value :  Profile level for radiance results.
#                  If not defined : 1 is the default value attributed
#                  0 : radiance calculation without taking account of polarization. 
#                  1 : radiance calculation with polarization processes.
#       Value format : Encoding as integer (10 characters max).
#
#    Keyword  : -SOS.IGmax  
#    ---------------------   
#       Status : Required      
#       Value :  Maximal order of interaction (scattering & surface reflexion). 
#       Value format : Encoding as integer (10 characters max).   	 	 	
#       
#    Keyword  : -SOS.ResBin  
#    ----------------------   
#       Status : Required      
#       Value :  Filename of the binary file resulting from SOS computations (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted     
#
#    Keyword  : -SOS.ResFileUp  
#    -------------------------   
#       Status : Required      
#       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
#                --> upward radiance.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted     
#
#    Keyword  : -SOS.ResFileDown  
#    ---------------------------   
#       Status : Required      
#       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
#                --> downward radiance.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted     
#
#    Keyword  : -SOS.ResFileUp.UserAng 
#    ---------------------------------  
#       Status : Optional  (requires the use of -SOS.UserAng.File)       
#       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
#                giving radiance values only for user's angles 
#                --> upward radiance.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted     
#
#    Keyword  : -SOS.ResFileDown.UserAng 
#    -----------------------------------  
#       Status : Optional  (requires the use of -SOS.UserAng.File)       
#       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
#                giving radiance values only for user's angles 
#                --> downward radiance.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted     
#
#    Keyword  : -SOS.Log
#    -------------------   
#       Status : Required 
#       Value :  log filename for SOS computations (complete path).
#      		 Assign 0 to avoid generating the log file.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted
#
#    Keyword  : -SOS.Config
#    ----------------------   
#       Status : Optional 
#       Value :  Filename for a copy of the applied SOS.h parameters file (complete path).
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted
#
#    Keyword  : -SOS.Trans
#    ---------------------  
#       Status : Optional 
#       Value :  Filename of the file given information about the simulated 
#                atmospheric transmissions and spherical albedo (complete path).
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted
#
#    Keyword  : -AP.UserFile  
#    -----------------------   
#       Status : Optional : has to be defined if AP.ResFile is not defined     
#       Value :  Filename of user's profile file (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -AP.ResFile  
#    ----------------------   
#       Status : Optional : has to be defined if AP.UserFile is not defined         
#       Value :  Filename of the result SOS_PROFIL computations (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted     
#
#    Keyword  : -AP.Log
#    ------------------   
#       Status : Required 
#       Value :  log filename for SOS_PROFIL computations (complete path).
#      		 Assign 0 to avoid generating the log file.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted     
#
#    Keyword  : -AP.MOT  
#    ------------------   
#       Status : Required      
#       Value :  Molecular optical thickness for the wavelength of radiance simulation
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AP.HR  
#    ------------------   
#       Status : Associated to -AP.MOT >= 0.0001      
#       Value :  Molecular heigth scale (km).
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AP.Type  
#    -------------------   
#       Status : Required      
#       Value :  Type of aerosol profile
#      		        1 : Aerosol profile described by an exponential decrease using a heigth scale
#			2 : Homogeneous aerosol layer between two altitudes.
#       Value format : Encoding as integer (10 characters max).
#
#    Keyword  : -AP.AerHS.HA  
#    -----------------------   
#       Status : Associated to -AP.Type 1      
#       Value :  Aerosol heigth scale (km).
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AP.AerLayer.Zmin  
#    ----------------------------   
#       Status : Associated to -AP.Type 2      
#       Value :  Minimal altitude of the aerosol layer (km).
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AP.AerLayer.Zmax  
#    ----------------------------   
#       Status : Associated to -AP.Type 2      
#       Value :  Maximal altitude of the aerosol layer (km).
#       Value format : Encoding as float (10 characters max).
#       	  
#    Keyword  : -AER.AOTref  
#    ----------------------   
#       Status : Required   
#       Value :  Aerosol optical thickness for the reference wavelength.
#       Value format : Encoding as float (10 characters max).       
#       
#    Keyword  : -AER.Waref  
#    ---------------------   
#       Status : Required (for common case with AER.AOTref > 0)      
#       Value :  Wavelength (microns) for reference aerosol optical thickness.
#                --> real value, without truncature applied.
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.UserFile  
#    ------------------------   
#       Status : Optional : has to be defined if AER.ResFile is not defined     
#       Value :  Filename of user's aerosol file (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -AER.ResFile  
#    -----------------------   
#       Status : Optional : has to be defined if AER.UserFile is not defined         
#       Value :  Filename of the result SOS_AEROSOLS computations (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted     
#
#    Keyword  : -AER.Log
#    -------------------   
#       Status : Required (for common case with AER.AOTref > 0)     
#       Value :  log filename for SOS_AEROSOLS computations (complete path).
#      		 Assign 0 to avoid generating the log file.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted     
#       
#    Keyword  : -AER.MieLog  
#    ----------------------   
#       Status : Required      
#       Value :  log filename for MIE computations (complete path). 
#      		 Assign 0 to avoid generating the log file.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted       
#
#    Keyword  : -AER.Tronca 
#    ----------------------   
#       Status : Required (for common case with AER.AOTref > 0)     
#       Value :  Option for aerosol phaze function troncature (1 to apply a troncature)
#       Value format : Encoding as integer (10 characters max).
#
#    Keyword  : -AER.Model 
#    ----------------------   
#       Status : Required (for common case with AER.AOTref > 0)     
#       Value :  Type of size distribution (mono or multi-mode)
#      		        0 : Mono-modal models.
#			1 : WMO models.
#			2 : Shettle & Fenn models.
#		        3 : LND bi-modal models.
#       Value format : Encoding as integer (10 characters max).
#
#    Keyword  : -AER.MMD.Mie.Filename
#    --------------------------------   
#       Status : Optionnal in association to -AER.Model 0 (for common case with AER.AOTref > 0)        
#       Value :  Filename of the results for MIE computations (local name in directory -AER.DirMie). 
#                                                              ----------
#                If this argument is not defined, the code will assign DEFAULT in order to 
#                automatically define the Mie filename.
#       Value format : Encoding as string (SOS_LENFIC1 characters max : see SOS.h).
#
#    Keyword  : -AER.MMD.Mie.AlphaMax  
#    ----------------------------   
#       Status : Associated to -AER.Model 0 (for common case with AER.AOTref > 0)     
#       Value :  Maximal value of the size parameter for Mie calculations.
#       Value format : Encoding as float (10 characters max).
#       
#
#    Keyword  : -AER.MMD.MRwa  
#    ------------------------   
#       Status : Associated to -AER.Model 0 (for common case with AER.AOTref > 0),
#                has to be defined if AER.Waref different of SOS.Wa.    
#       Value :  Real part of the refractive index at the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.MMD.MIwa  
#    ------------------------ 
#       Status : Associated to -AER.Model 0 (for common case with AER.AOTref > 0),
#                has to be defined if AER.Waref different of SOS.Wa.         
#       Value :  Imaginary part of the refractive index for the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F8.5 (10 characters max).
#	  
#    Keyword  : -AER.MMD.MRwaref  
#    ---------------------------   
#       Status : Associated to -AER.Model 0 (for common case with AER.AOTref > 0).
#       Value :  Real part of the refractive index at the reference wavelength of aerosol properties calculation.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.MMD.MIwaref  
#    --------------------------- 
#       Status : Associated to -AER.Model 0 (for common case with AER.AOTref > 0).         
#       Value :  Imaginary part of the refractive index for the reference wavelength of aerosol properties calculation.
#       Value format : Encoding as float F8.5 (10 characters max).
#	                
#    Keyword  : -AER.MMD.SDtype 
#    --------------------------   
#       Status : Associated to -AER.Model 0 (for common case with AER.AOTref > 0)  
#       Value :  Type of mono-modal size distribution 
#			1 : LND.
#			2 : Junge's law.
#       Value format : Encoding as integer (10 characters max).			
#	    
#    Keyword  : -AER.MMD.SDparam1  
#    ---------------------------- 
#       Status : Associated to -AER.Model 0 (for common case with AER.AOTref > 0)     
#       Value :  Depending on the -AER.MMD.SDtype value.
#			-AER.MMD.SDtype 1 : modal radius of the LND (microns).
#			-AER.MMD.SDtype 2 : minimal radius of the Junge's law (microns).
#       Value format : Encoding as float (10 characters max).	  
#	  
#    Keyword  : -AER.MMD.SDparam2  
#    ---------------------------- 
#       Status : Associated to -AER.Model 0 (for common case with AER.AOTref > 0)     
#       Value :  Depending on the -AER.MMD.SDtype value.
#			-AER.MMD.SDtype 1 : Standard deviation of the LND.
#			-AER.MMD.SDtype 2 : power applied to the radius of the Junge's law.
#       Value format : Encoding as float (10 characters max).	  
#	  
#    Keyword  : -AER.WMO.Model  
#    ------------------------- 
#       Status : Associated to -AER.Model 1 (for common case with AER.AOTref > 0)     
#       Value :  Type of WMO model.
#      		        1 : Continental WMO model.
#			2 : Maritime WMO model.
#			3 : Urban WMO model.
#		        4 : WMO model by user definition.
#       Value format : Encoding as integer (10 characters max).	  
#
#    Keyword  : -AER.WMO.DL  
#    ---------------------- 
#       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOTref > 0)     
#       Value :  Volumetric concentration (%) for WMO Dust-Like particles.
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.WMO.WS  
#    ---------------------- 
#       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOTref > 0)     
#       Value :  Volumetric concentration (%) for WMO Water Soluble particles.
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.WMO.OC  
#    ---------------------- 
#       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOTref > 0)     
#       Value :  Volumetric concentration (%) for WMO OCeanic particles.
#       Value format : Encoding as float  (10 characters max).
#
#    Keyword  : -AER.WMO.SO  
#    ---------------------- 
#       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOTref > 0)     
#       Value :  Volumetric concentration (%) for WMO SOot particles.
#       Value format : Encoding as float (10 characters max).	  
#	  
#    Keyword  : -AER.SF.Model  
#    ------------------------ 
#       Status : Associated to -AER.Model 2 (for common case with AER.AOTref > 0)     
#       Value :  Type of Shettle & Fenn model.
#      		        1 : Tropospheric S&F model.
#			2 : Urban S&F model.
#			3 : Maritime S&F model.
#		        4 : Coastal S&F model.
#       Value format : Encoding as integer (10 characters max).	  
#
#    Keyword  : -AER.SF.HR  
#    --------------------- 
#       Status : Associated to -AER.Model 2 (for common case with AER.AOTref > 0)     
#       Value :  Relative humidity (%) for Shettle & Fenn model.
#       Value format : Encoding as float F5.2 (10 characters max).      
#	  
#    Keyword  : -AER.BMD.VCdef 
#    ---------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOTref > 0)  
#       Value :  Type of selected method for bi-modal LND volumetric concentration definition 
#			1 : Use of user volumetric concentrations.
#			2 : Use of user ratio AOT_coarse / AOT_total.
#       Value format : Encoding as integer (10 characters max).			
#
#    Keyword  : -AER.BMD.CoarseVC  
#    ---------------------------- 
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 1 (for common case with AER.AOTref > 0)     
#       Value :  User volumetric concentration for the "LND coarse mode".
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.BMD.FineVC  
#    -------------------------- 
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 1 (for common case with AER.AOTref > 0)     
#       Value :  User volumetric concentration for the "LND fine mode".
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.BMD.RAOT  
#    ------------------------ 
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOTref > 0)     
#       Value :  User value of the ration AOT_coarse / AOT_total for the aerosol reference wavelength.
#       Value format : Encoding as float 10 characters max).
#	         
#    Keyword  : -AER.BMD.CM.MRwa  
#    ---------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOTref > 0),
#                has to be defined if AER.Waref different of SOS.Wa.             
#       Value :  Real part of the refractive index for the "LND coarse mode" 
#                at the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.MIwa  
#    ---------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOTref > 0),
#                has to be defined if AER.Waref different of SOS.Wa.             
#       Value :  Imaginary part of the refractive index for the "LND coarse mode" 
#                at the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F8.5 (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.MRwaref  
#    ------------------------------   
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOTref > 0).    
#       Value :  Real part of the refractive index for the "LND coarse mode" 
#                at the aerosol reference wavelength.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.MIwaref  
#    -----------------------------   
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOTref > 0).    
#       Value :  Imaginary part of the refractive index for the "LND coarse mode" 
#                at the aerosol reference wavelength.
#       Value format : Encoding as float F8.5 (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.SDradius  
#    -------------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOTref > 0).    
#       Value :  Modal radius of the "LND coarse mode" (microns).
#       Value format : Encoding as float (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.SDvar  
#    ----------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOTref > 0).    
#       Value :  Standard deviation of the "LND coarse mode".
#       Value format : Encoding as float (10 characters max).
#	  	 	  
#    Keyword  : -AER.BMD.FM.MRwa  
#    ---------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0),
#                has to be defined if AER.Waref different of SOS.Wa.                 
#       Value :  Real part of the refractive index for the "LND fine mode" 
#                at the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.BMD.FM.MIwa  
#    ---------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0),
#                has to be defined if AER.Waref different of SOS.Wa.                 
#       Value :  Imaginary part of the refractive index for the "LND fine mode" 
#                at the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F8.5 (10 characters max).
#	  
#    Keyword  : -AER.BMD.FM.MRwaref  
#    ------------------------------   
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
#       Value :  Real part of the refractive index for the "LND fine mode" 
#                at the aerosol reference wavelength.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.BMD.FM.MIwaref  
#    -----------------------------   
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
#       Value :  Imaginary part of the refractive index for the "LND fine mode" 
#                at the aerosol reference wavelength.
#       Value format : Encoding as float F8.5 (10 characters max).
#	  
#    Keyword  : -AER.BMD.FM.SDradius  
#    -------------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
#       Value :  Modal radius of the "LND fine mode" (microns).
#       Value format : Encoding as float (10 characters max).
#	  
#    Keyword  : -AER.BMD.FM.SDvar  
#    ----------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
#       Value :  Standard deviation of the "LND fine mode".
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.ExtData
#    -----------------------   
#       Status : Associated to -AER.Model 4 (for common case with AER.AOT > 0)       
#       Value :  Filename (complete path) of user's external phaze functions 
#                and radiative parameters (extinction and scattering coefficients)
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#       Special requirement : the reference aerosol wavelength (-AER.Waref) and 
#                             the radiance simulation wavelength (-SOS.Wa) have to be equal.
#
#    Keyword  : -SURF.Log  
#    --------------------   
#       Status : Required      
#       Value :  log filename for SOS_SURFACE computations (complete path). 
#      		 Assign 0 to avoid generating the log file.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted
#	 	
#    Keyword  : -SURF.File  
#    ---------------------   
#       Status : Associated to  -SURF.Type 1 3 4 5 or 6     
#       Value :  User filename for the surface reflexion matrix (complete path).
#                --> Glitter or others BRDF/BPDF. 
#                DEFAULT if one uses a file automatically generated by the SOS code
#                and stored in the directory SURF.Dir.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -SURF.Type  
#    ---------------------   
#       Status : Required      
#       Value :  Kind of surface 
#      		        0 : lambertian surface for an albedo RHO.
#			1 : lambertian surface + glitter (roughness water).
#			2 : lambertian surface + Fresnel's reflexion (plane water).
#		        3 : lambertian surface + Roujean's BRDF.
#			4 : lambertian surface + Roujean's BRDF + Rondeaux's BPDF.
#			5 : lambertian surface + Roujean's BRDF + Breon's BPDF .
#			6 : lambertian surface + Roujean's BRDF + Nadal's BPDF.
#       Value format : Encoding as integer (10 characters max).
#
#    Keyword  : -SURF.Alb 
#    --------------------   
#       Status : Required      
#       Value :  Surface albedo for the wavelength of radiance simulation
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -SURF.Ind  
#    --------------------   
#       Status : Associated to  -SURF.Type 1 2 4 5 or 6     
#       Value :  Surface / atmosphere refractive index. 
#       Value format : Encoding as float F5.3 (10 characters max).
#	 		 
#    Keyword  : -SURF.Glitter.Wind  
#    -----------------------------   
#       Status : Associated to  -SURF.Type 1    
#       Value :  Wind velocity (m/s). 
#       Value format : Encoding as float F4.1 (10 characters max).
#       
#    Keyword  : -SURF.Roujean.K0  
#    ---------------------------   
#       Status : Associated to  -SURF.Type 3 4 5 or 6      
#       Value :  Parameter K0 of Roujean's BRDF model. 
#       Value format : Encoding as float F7.3 (10 characters max).
#
#    Keyword  : -SURF.Roujean.K1  
#    --------------------------- 
#       Status : Associated to  -SURF.Type 3 4 5 or 6      
#       Value :  Parameter K1 of Roujean's BRDF model. 
#       Value format : Encoding as float F7.3 (10 characters max).
#
#    Keyword  : -SURF.Roujean.K2  
#    ---------------------------   
#       Status : Associated to  -SURF.Type 3 4 5 or 6      
#       Value :  Parameter K2 of Roujean's BRDF model. 
#       Value format : Encoding as float F7.3 (10 characters max).
#
#    Keyword  : -SURF.Nadal.Alpha 
#    ----------------------------   
#       Status : Associated to  -SURF.Type 6     
#       Value :  Parameter Alpha of Nadal's BPDF model. 
#       Value format : Encoding as float F6.4 (10 characters max).
#
#    Keyword  : -SURF.Nadal.Beta  
#    ---------------------------   
#       Status : Associated to  -SURF.Type 6      
#       Value :  Parameter Beta of Nadal's BPDF model. 
#       Value format : Encoding as float F4.1 (10 characters max).
#
