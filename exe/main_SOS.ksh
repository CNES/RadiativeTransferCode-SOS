#!/bin/ksh
#################################################################
# FILE: main_SOS.ksh
#
# PROJECT: Successive Orders of Scattering.
#
# RULE: Command file for a complete execution of a SOS simulation :
#       --> Surface BRDF and BPDF calculations
#       --> Aerosol optical thickness and phaze matrix calculations 
#       --> Atmospheric profile of optical thicknesses calculations 
#       --> Radiance simulation by successive order of scaterring calculations 

#    
# This shell is available for the SOS code version 5.0,  
# not for previous ones.          -------------------- 
#
#    
# In order to use this shell, the environment variable SOS_RACINE
# has to be previously defined :                       ----------
#       SOS_RACINE: complete access path to the software directory 
#
#
# AUTHOR: Bruno Lafrance ( CS )
# DATE: 2009 / 11 / 12
#
#
#
## MOD:VERSION:5.0:  2009 / 11 / 12
#      - main adjustement for a parameter definition by
#        couples : -Keyword Value
# 
## MOD:VERSION:5.1:  2010 / 01 / 15
#      - introduction of new keywords 
#      - Test if the value associated to a keyword is not missing 
#      - The configuration file used to the compilation 
#        is now copied in -SOS.Config (if this one is defined)
#
## MOD:VERSION:5.2:  2015 / 12 / 09
#      - Adjustment for the reading of kmat1_wa and kmat1_waref
#        with respect to the new format of Aerosols files.
#      - Adding information while an external aerosols phaze 
#        function is used.
#      - Deleting the non relevant specified formats 
#        for the values of parameters. 
#      - Correction of an occuring error if -SURF.File is missing
#        while it is not used by the code
#      - Correction of an occuring error if -AER.Model is missing
#        while the user requires using his aerosols file AER.UserFile
#      - Adding the control that the wavelength for radiance 
#        simulations (-SOS.Wa) equals the reference wavelength 
#        (-AER.Waref) for the AOT definition,
#        while the user requires using his profile file AP.UserFile
#################################################################


####################################################################
#                                                                  #
#  Creating the directory structure for the storage of results     #
#  and intermediary files 				           #
#                                                                  #
####################################################################

test ! -d ${dirMIE}      && mkdir -p ${dirMIE}
test ! -d ${dirLOG}      && mkdir -p ${dirLOG}
test ! -d ${dirRESULTS}  && mkdir -p ${dirRESULTS}

test ! -d ${dirSUNGLINT} && mkdir -p ${dirSUNGLINT}
test ! -d ${dirROUJEAN}  && mkdir -p ${dirROUJEAN}
test ! -d ${dirRH}       && mkdir -p ${dirRH}
test ! -d ${dirBREON}    && mkdir -p ${dirBREON}
test ! -d ${dirNADAL}    && mkdir -p ${dirNADAL}
     


##################################################################
##################################################################
###                   Function : SOS_ANGLES
##################################################################
##################################################################
#
#-----------------------------------------------------------------
#  RULE 
#-----------------------------------------------------------------
#    Define the angles to be used to :
#      - Phaze function simulations
#      - BRDF/BPDF & radiance simulations
#
#    Define the internal parameters for SOS_AEROSOLS, SOS_SURFACE  
#    and SOS computations : OS_NB (fonctions de phase), 
#                           OS_NS (ordre max des developpements des
#                                  fonctions de Fresnel et des luminances)
#   			    OS_NM (ordre maximal du developpement en series de Fourier  
#                                  de la fonction G des calculs de matrice de reflexion
#
#-----------------------------------------------------------------
#  ARGUMENTS 
#-----------------------------------------------------------------
#    Use of couples (-Keyword Value) for arguments definition
#
#    Keyword  : -ANG.Rad.NbGauss  
#    ---------------------------   
#       Status : Optional (if not defined, the default value SOS_DEFAULT_NBMU_LUM is used : see SOS.h )       
#       Value :  Number of Gauss angles to be used for radiance computations 
#       Value format : Encoding as integer (10 characters max)
#
#    Keyword  : -ANG.Rad.UserAngFile  
#    -------------------------------   
#       Status : Optional      
#       Value :  Filename of the complementary list of user's angles to complete the 
#                ANG.Rad.NbGauss angles (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -ANG.Rad.Thetas  
#    --------------------------   
#       Status : Required     
#       Solar zenithal angle (degrees) 
#       Value format : Encoding as float (10 characters max).
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
#       Value format : Encoding as integer (10 characters max)
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
#-----------------------------------------------------------------
# Results provided by SOS_ANGLES
#-----------------------------------------------------------------
#   - List of angles to be used, including internal parameters of maximal orders 
#     for Legendre and Fourier's developments.
#   - Log file ANGLES (optional)
#
#------------------------------------------------------------------
#
## Launching ANGLES program
##--------------------------
function Call_def_angles {

   ## Launching SOS_ANGLES program
   ##------------------------------------
   echo "-->Running SOS_ANGLES"   
   $SOS_RACINE/exe/SOS_ANGLES.exe ${Arg_ANGLES}
     
   ## If an error occurs during PROFIL computation : the OS code is stopped
   ##------------------------------------------------------------
   if [ "-$?" != "-0" ]
   then
	 echo "ERROR for SOS_ANGLES process"
         echo "STOP SOS code"
         exit 1
   else
         echo "SOS_ANGLES OK"   
   fi
   
}  #End Call_def_angles



##################################################################
##################################################################
###                   Function : SOS_SURFACE
##################################################################
##################################################################
#
#-----------------------------------------------------------------
#  RULE 
#-----------------------------------------------------------------
#    Running SURFACE BRDF/BPDF file computations
#
#-----------------------------------------------------------------
#  ARGUMENTS 
#-----------------------------------------------------------------
#    Use of couples (-Keyword Value) for arguments definition
#
#    Keyword  : -SURF.AngFile 
#    ------------------------   
#       Status : Required      
#       Value :  Filename of angles to be used in SURFACE computations 
#                (complete path)
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
#    Keyword  : -SURF.Log  
#    --------------------   
#       Status : Required      
#       Value :  log filename for SOS_SURFACE computations (complete path). 
#      		 Assign 0 to avoid generating the log file.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted
#
#    Keyword  : -SURF.Dir  
#    --------------------   
#       Status : Required      
#       Value :  Storage directory for BRDF/BPDF files producted by SOS_SURFACE computations (complete path). 
#       Value format : Encoding as string (SOS_LENDIR characters max : see SOS.h).
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
#    Keyword  : -SURF.DirBRDFtoBPDFcalculation  
#    -----------------------------------------   
#       Status : Associated to  -SURF.Type 3 4 5 or 6     
#       Value :  Storage directory for Roujean's BRDF files, associated to BPDF models (complete path).
#       Value format : Encoding as string (SOS_LENDIR characters max: see SOS.h).
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
#
#-----------------------------------------------------------------
# Results provided by SOS_SURFACE
#-----------------------------------------------------------------
#   - SURFACE BRDF/BPDF file (Glitter, Roujean, ...)
#        ==> Already existing files are not generated again.
#   - Log file SURFACE (optional)
#   - Glitter computation creates 3 temporary files which are deleted
#     at the end of computations (files: RES_GSF, RES_FRESNEL et RES_MAT_REFLEX)
#
#------------------------------------------------------------------
#
## Launching SURFACE program
##--------------------------
function Call_simu_surface {

   repSURFACE=DEFAULT	#Initialization in case of: typeSurf = 0 ou 2.
   repBRDFpourBPDF=$dirROUJEAN	#Storage directory for Roujean's BRDF files
				#(also usefull for Rondeaux, Breon or Nadal's BPDF calculation).
				
   if [ $typeSurf -ne 0 ] && [ $typeSurf -ne 2 ] ; then

       if [ $typeSurf -eq 1 ] ; then
          repSURFACE=$dirSUNGLINT ;
       else 
          if [ $typeSurf -eq 3 ] ; then
             repSURFACE=$dirROUJEAN  ;
          else 
             if [ $typeSurf -eq 4 ] ; then
                repSURFACE=$dirRH  ;
             else 
                if [ $typeSurf -eq 5 ] ; then
                   repSURFACE=$dirBREON ;
                else 
       	           if [ $typeSurf -eq 6 ] ; then
                      repSURFACE=$dirNADAL ;
                   fi ;
       	        fi;
     	     fi;
          fi;
       fi;     	 
 
 				 
       if [ "-$ficSURFACE" = "-DEFAULT" ]
       then
         echo "-->Running SOS_SURFACE"
	 
	 Arg_SURFACE=${Arg_SURFACE}" -SURF.Dir $repSURFACE"
	 Arg_SURFACE=${Arg_SURFACE}" -SURF.DirBRDFtoBPDFcalculation $repBRDFpourBPDF"	 
	 
	 $SOS_RACINE/exe/SOS_SURFACE.exe $Arg_SURFACE

         ## If an error occurs during SURFACE computation : the OS code is stopped
         ##------------------------------------------------------------
         if [ "-$?" != "-0" ]
         then
            echo "ERROR for SOS_SURFACE process"
            echo "Generated BRDF / BPDF files have to be deleted"
            echo "STOP SOS code"
            exit 1
         else        
            echo "SOS_SURFACE OK"   ;
         fi
       fi;    #End of test  ficSURFACE = DEFAULT
    
   fi	#End of test typeSurf = 0 or 2

}  #End Call_simu_surface





##################################################################
##################################################################
###                   Function : SOS_AEROSOLS
##################################################################
##################################################################
#-----------------------------------------------------------------
#  RULE 
#-----------------------------------------------------------------
#    Running AEROSOL file computations
#     --> Calculation of scattering and extinction coefficient.
#     --> Calculation of asymetric factor.
#     --> Calculation of phaze function truncature coefficient.
#     --> Calculation of single scattering albedo.
#     --> Calculation of phaze matrix coefficients.
#
#-----------------------------------------------------------------
#  ARGUMENTS 
#-----------------------------------------------------------------
#    Use of couples (-Keyword Value) for arguments definition
#
#
#    Keyword  : -AER.AngFile 
#    -----------------------   
#       Status : Required      
#       Value :  Filename of angles to be used in AEROSOLS computations 
#                (complete path)
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -AER.AOT  
#    -------------------   
#       Status : Required      
#       Value :  Aerosol optical thickness for the wavelength of aerosol properties calculation
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.Wa  
#    ------------------   
#       Status : Required (for common case with AER.AOT > 0)     
#       Value :  Wavelength of aerosol properties calculation.
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.Tronca 
#    ----------------------   
#       Status : Required (for common case with AER.AOT > 0)     
#       Value :  Option for aerosol phaze function troncature (1 to apply a troncature)
#       Value format : Encoding as integer (4 characters max).
#
#    Keyword  : -AER.Log
#    -------------------   
#       Status : Required (for common case with AER.AOT > 0)     
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
#    Keyword  : -AER.DirMie
#    ----------------------   
#       Status : Required      
#       Value :  Storage directory for MIE files producted by SOS_AEROSOLS computations (complete path). 
#       Value format : Encoding as string (SOS_LENDIR characters max : see SOS.h).
#
#    Keyword  : -AER.ResFile
#    -----------------------   
#       Status : Required      
#       Value :  Filename of the result SOS_AEROSOLS computations or user filename (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -AER.Model 
#    ----------------------   
#       Status : Required (for common case with AER.AOT > 0)     
#       Value :  Type of size distribution (mono or multi-mode)
#      		        0 : Mono-modal models.
#			1 : WMO models.
#			2 : Shettle & Fenn models.
#		        3 : LND bi-modal models.
#                       4 : Use of external phase functions
#       Value format : Encoding as integer (10 characters max).
#
#    Keyword  : -MMD.Mie.AlphaMax  
#    ----------------------------   
#       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
#       Value :  Maximal value of the size parameter for Mie calculations.
#       Value format : Encoding as float F9.2 (10 characters max).
#
#    Keyword  : -AER.MMD.Mie.Filename
#    --------------------------------   
#       Status : Optionnal in association to -AER.Model 0 (for common case with AER.AOT > 0)        
#       Value :  Filename of the results for MIE computations (local name in directory -AER.DirMie). 
#                                                              ----------
#                If this argument is not defined, the code will assign DEFAULT in order to 
#                automatically define the Mie filename.
#       Value format : Encoding as string (SOS_LENFIC1 characters max : see SOS.h).
#					
#    Keyword  : -AER.MMD.MRwa  
#    ------------------------ 
#       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
#       Value :  Real part of the refractive index at the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.MMD.MIwa  
#    ------------------------ 
#       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
#       Value :  Imaginary part of the refractive index for the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F8.5 (10 characters max).
#	  
#    Keyword  : -AER.MMD.SDtype #    Keyword  : -AER.ResFile
#    -----------------------   
#       Status : Required      
#       Value :  Filename of the result SOS_AEROSOLS computations or user filename (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    --------------------------   
#       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)  
#       Value :  Type of mono-modal size distribution 
#			1 : LND.
#			2 : Junge's law.
#       Value format : Encoding as integer (4 characters max).			
#	    
#    Keyword  : -AER.MMD.SDparam1  
#    ---------------------------- 
#       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
#       Value :  Depending on the -AER.MMD.SDtype value.
#			-AER.MMD.SDtype 1 : modal radius of the LND (microns).
#			-AER.MMD.SDtype 2 : minimal radius of the Junge's law (microns).
#       Value format : Encoding as float (10 characters max).	  
#	  
#    Keyword  : -AER.MMD.SDparam2  
#    ---------------------------- 
#       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
#       Value :  Depending on the -AER.MMD.SDtype value.
#			-AER.MMD.SDtype 1 : log10 of the LND variance (microns).
#			-AER.MMD.SDtype 2 : power applied to the radius of the Junge's law (microns).
#       Value format : Encoding as float (10 characters max).	  
#	  
#    Keyword  : -AER.WMO.Model  
#    ------------------------- 
#       Status : Associated to -AER.Model 1 (for common case with AER.AOT > 0)     
#       Value :  Type of WMO model.
#      		        1 : Continental WMO model.
#			2 : Maritime WMO model.
#			3 : Urban WMO model.
#		        4 : WMO model by user definition.
#       Value format : Encoding as integer (10 characters max).	  
#
#    Keyword  : -AER.WMO.DL  
#    ---------------------- 
#       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
#       Value :  Volumetric concentration (%) for WMO Dust-Like particles.
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.WMO.WS  
#    ---------------------- 
#       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
#       Value :  Volumetric concentration (%) for WMO Water Soluble particles.
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.WMO.OC  
#    ---------------------- 
#       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
#       Value :  Volumetric concentration (%) for WMO OCeanic particles.
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.WMO.SO  
#    ---------------------- 
#       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
#       Value :  Volumetric concentration (%) for WMO SOot particles.
#       Value format : Encoding as float (10 characters max).	  
#	  
#    Keyword  : -AER.SF.Model  
#    ------------------------ 
#       Status : Associated to -AER.Model 2 (for common case with AER.AOT > 0)     
#       Value :  Type of Shettle & Fenn model.
#      		        1 : Tropospheric S&F model.
#			2 : Urban S&F model.
#			3 : Maritime S&F model.
#		        4 : Coastal S&F model.
#       Value format : Encoding as integer (10 characters max).	  
#
#    Keyword  : -AER.SF.HR  
#    --------------------- 
#       Status : Associated to -AER.Model 2 (for common case with AER.AOT > 0)     
#       Value :  Relative humidity (%) for Shettle & Fenn model.
#       Value format : Encoding as float F5.2 (10 characters max).
#	  
#    Keyword  : -AER.BMD.VCdef 
#    ---------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)  
#       Value :  Type of selected method for bi-modal LND volumetric concentration definition 
#			1 : Use of user volumetric concentrations.
#			2 : Use of user ratio AOT_coarse / AOT_total.
#       Value format : Encoding as integer (10 characters max).			
#
#    Keyword  : -AER.BMD.CoarseVC  
#    ---------------------------- 
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 1 (for common case with AER.AOT > 0)     
#       Value :  User volumetric concentration for the "LND coarse mode".
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.BMD.FineVC  
#    -------------------------- 
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 1 (for common case with AER.AOT > 0)     
#       Value :  User volumetric concentration for the "LND fine mode".
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -AER.BMD.RAOT  
#    ------------------------ 
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)     
#       Value :  User value of the ration AOT_coarse / AOT_total for the aerosol reference wavelength.
#       Value format : Encoding as float (10 characters max).
#	  
#    Keyword  : -AER.Waref  
#    ---------------------   
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
#       Value :  Reference wavelength for the aerosol optical thickness (microns).
#       Value format : Encoding as float (10 characters max).
#	  
#    Keyword  : -AER.AOTref  
#    ----------------------   
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
#       Value :  Aerosol optical thickness for the reference wavelength.
#       Value format : Encoding as float (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.MRwa  
#    ---------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
#       Value :  Real part of the refractive index for the "LND coarse mode" 
#                at the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.MIwa  
#    ---------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
#       Value :  Imaginary part of the refractive index for the "LND coarse mode" 
#                at the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F8.5 (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.MRwaref  
#    ------------------------------   
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
#       Value :  Real part of the refractive index for the "LND coarse mode" 
#                at the aerosol reference wavelength.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.MIwaref  
#    -----------------------------   
#       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
#       Value :  Imaginary part of the refractive index for the "LND coarse mode" 
#                at the aerosol reference wavelength.
#       Value format : Encoding as float F8.5 (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.SDradius  
#    -------------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
#       Value :  Modal radius of the "LND coarse mode" (microns).
#       Value format : Encoding as float (10 characters max).
#	  
#    Keyword  : -AER.BMD.CM.SDvar  
#    ----------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
#       Value :  Log10 of the "LND coarse mode" variance.
#       Value format : Encoding as float (10 characters max).
#	  
#    Keyword  : -AER.BMD.FM.MRwa  
#    ---------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
#       Value :  Real part of the refractive index for the "LND fine mode" 
#                at the wavelength of aerosol properties calculation.
#       Value format : Encoding as float F5.3 (10 characters max).
#	  
#    Keyword  : -AER.BMD.FM.MIwa  
#    ---------------------------   
#       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
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
#       Value :  Log10 of the "LND fine mode" variance.
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
#-----------------------------------------------------------------
# Results provided by SOS_AEROSOLS
#-----------------------------------------------------------------
#   - MIE files
#        ==> Already existing files are not generated again.
#   - Log file of MIE calculations (optionnal)
#   - Log file of phaze matrix calculations (optionnal)
#   - Aerosols file result : aerosols radiative properties
#        ==> Already existing files are overwritted
#
#------------------------------------------------------------------
#
## Launching AEROSOLS program
##---------------------------
function Call_simu_aer {

   if [ "-$ficAERgranu" = "-DEFAULT" ]
   then

      ## Launching SOS_AEROSOLS program
      ##------------------------------------
      echo "-->Running SOS_AEROSOLS"

      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.DirMie ${dirMIE}"
	
      $SOS_RACINE/exe/SOS_AEROSOLS.exe ${Arg_AEROSOLS_Simu}
	      
      ## If an error occurs during AEROSOLS computation : the OS code is stopped
      ##------------------------------------------------------------
      if [ "-$?" != "-0" ]
      then
            echo "ERROR for SOS_AEROSOLS process"
            echo "Generated MIE files have to be deleted"
            echo "STOP SOS code"
            exit 1
      else
            echo "SOS_AEROSOLS OK"   
      fi
   fi

}  #End Call_simu_aer


          


##################################################################
##################################################################
###                   Function : SOS_PROFIL
##################################################################
##################################################################
#-----------------------------------------------------------------
#  RULE 
#-----------------------------------------------------------------
#    Running PROFILE file computations
#     --> Calculation of molecular and aerosol optical thicknesses 
#         for each level of an atmospheric profile.
#     --> Proportion of molecules and aerosols for each level.
#
#-----------------------------------------------------------------
#  ARGUMENTS 
#-----------------------------------------------------------------
#    Use of couples (-Keyword Value) for arguments definition
#
#    Keyword  : -AP.ResFile  
#    ----------------------   
#       Status : Required      
#       Value :  Filename of the result SOS_PROFIL computations (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
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
#    Keyword  : -AP.AOT  
#    ------------------   
#       Status : Required      
#       Value :  Aerosol optical thickness for the wavelength of radiance simulation
#                --> real value, without truncature applied.
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
#       Status : Associated to -AP.Type 1 and -AP.AOT >= 0.0001     
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
#
#-----------------------------------------------------------------
# Results provided by SOS_PROFIL
#-----------------------------------------------------------------
#   - Log file of PROFIL calculations (optionnal)
#   - Profile file result : Optical thickness profile and proportion of molecules and aerosols
#                           (without aerosol truncature adjustment)
#        ==> Already existing files are overwritted
#------------------------------------------------------------------
#
## Launching PROFIL program
##-------------------------
function Call_simu_profil {
   if [ "-$ficPROFIL" = "-DEFAULT" ]
   then 

      ## Launching SOS_PROFIL program
      ##------------------------------------
      echo "-->Running SOS_PROFIL"
      
      Arg_PROFIL=${Arg_PROFIL}" -AP.AOT ${AOT}"
      
      $SOS_RACINE/exe/SOS_PROFIL.exe ${Arg_PROFIL}
     
      ## If an error occurs during PROFIL computation : the OS code is stopped
      ##------------------------------------------------------------
      if [ "-$?" != "-0" ]
      then
            echo "ERROR for SOS_PROFIL process"
            echo "STOP SOS code"
            exit 1
      else
            echo "SOS_PROFIL OK"   
      fi
   fi
}  #End Call_simu_profil



##################################################################
##################################################################
###                   Function : SOS_OS
##################################################################
##################################################################
#-----------------------------------------------------------------
#  RULE 
#-----------------------------------------------------------------
#    Running procedure to
#     --> read the aerosols radiative parameters
#     --> read the atmospheric profile of molecule 
#         and aerosol optical thicknesses.
#     --> resolve the radiative transfer equation by
#         successive order of interaction (scattering/surface reflexion).
#     --> calculate the Stokes parameters from their decomposition as 
#         Fourier series..
#
#-----------------------------------------------------------------
#  ARGUMENTS 
#-----------------------------------------------------------------
#    Use of couples (-Keyword Value) for arguments definition
#
#    Keyword  : -SOS.AngFile 
#    -----------------------   
#       Status : Required      
#       Value :  Filename of angles to be used in SOS computations 
#                (complete path)
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -SOS.AerFile  
#    -----------------------   
#       Status : Required      
#       Value :  Filename of the result SOS_AEROSOLS computations (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -SOS.ProfileFile  
#    ---------------------------   
#       Status : Required      
#       Value :  Filename of the result SOS_PROFIL computations (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -SOS.ResBin  
#    ----------------------   
#       Status : Required      
#       Value :  Filename of the binary file resulting from SOS computations (complete path). 
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -SOS.ResFileUp  
#    -------------------------   
#       Status : Required      
#       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
#                --> upward radiance.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#
#    Keyword  : -SOS.ResFileDown  
#    ---------------------------   
#       Status : Required      
#       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
#                --> downward radiance.
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
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
#    Keyword  : -SOS.MDF  
#    -------------------   
#       Status : Required      
#       Value :  Factor of molecular depolarization.
#       Value format : Encoding as float (10 characters max).
#
#    Keyword  : -SOS.Trans
#    ---------------------  
#       Status : Optional 
#       Value :  Filename of the file given information about the simulated 
#                atmospheric transmissions (complete path).
#       Value format : Encoding as string (SOS_LENFIC2 characters max : see SOS.h).
#        ==> Already existing file will be overwritted
#
#    Keyword  : -SURF.Alb 
#    --------------------   
#       Status : Required      
#       Value :  Surface albedo for the wavelength of radiance simulation
#       Value format : Encoding as float (10 characters max).
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
#    Keyword  : -SURF.Ind  
#    --------------------   
#       Status : Associated to  -SURF.Type 1 2 4 5 or 6     
#       Value :  Surface / atmosphere refractive index. 
#       Value format : Encoding as float F5.3 (10 characters max).
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
#    Keyword  : -SURF.Dir  
#    --------------------   
#       Status : Associated to  -SURF.Type 1 3 4 5 or 6  and -SURF.File DEFAULT   
#       Value :  Storage directory for BRDF/BPDF files automatically 
#                generated by the SOS code  --> Glitter or others BRDF/BPDF.               
#       Value format : Encoding as string (SOS_LENDIR characters max : see SOS.h).
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
#       Value format : Encoding as integer (4 characters max).
#
#    Keyword  : -SOS.IGmax  
#    ---------------------   
#       Status : Required      
#       Value :  Maximal order of interaction (scattering & surface reflexion). 
#       Value format : Encoding as integer (10 characters max).   	 	 	
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
#       Value format : Encoding as integer (10 characters max).   	 	 	
#
#
#-----------------------------------------------------------------
# Results provided by SOS_OS
#-----------------------------------------------------------------
#   - Log file of SOS calculations (optionnal)
#   - Radiance files result : binary one
#                             ascii upward radiance
#                             ascii downward radiance
#        ==> Already existing files are overwritted
#------------------------------------------------------------------
#
## Launching SOS program
##----------------------
function Call_simu_os { 
   echo "-->Running SOS"
  
   $SOS_RACINE/exe/SOS.exe ${Arg_SOS}						

   ## If an error occurs during SOS computation : the OS code is stopped
   ##------------------------------------------------------------
   if [ "-$?" != "-0" ]
   then
      echo "ERROR for SOS process"
      echo "STOP SOS code"
      exit 1
   else
      echo "SOS OK"   
   fi


}  #End Call_simu_os




####################################################################
#                                                                  #
#  Main program                                                    #
#                                                                  #
####################################################################
#-----------------------------------------------------------------
#  RULE 
#-----------------------------------------------------------------
#    main program wich
#     --> loads the user's parameters for a radiance simulation.
#     --> checks the status for specific arguments required,
#     --> call the function to prepare the angles to be used by simulations,
#     --> call the function to simulate the surface radiative properties,
#     --> prepare the aerosols parameters and call the function to 
#         simulate the aerosols radiative properties for the wavelength of
#         radiance simulation,
#     --> call the function to simulate the molecules and aerosols profile,
#     --> call the function to simulate the diffuse radiance.
#
#-----------------------------------------------------------------
#  ARGUMENTS 
#-----------------------------------------------------------------
#    Use of couples (-Keyword Value) for arguments definition  
#		
#    Keyword  : -ANG.Thetas  
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
#       Value :  Constant step on azimut (degrees). 
#       Value format : Encoding as integer (10 characters max).   	 	      
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
#       Value format : Encoding as integer (10 characters max).
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
#                atmospheric transmissions (complete path).
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
#			4 : Use of external phase functions
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
#    --------------------------------   
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
#			-AER.MMD.SDtype 1 : log10 of the LND variance (microns).
#			-AER.MMD.SDtype 2 : power applied to the radius of the Junge's law (microns).
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
#       Value format : Encoding as float (10 characters max).
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
#       Value format : Encoding as float (10 characters max).
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
#       Value :  Log10 of the "LND coarse mode" variance.
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
#       Value :  Log10 of the "LND fine mode" variance.
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
#-----------------------------------------------------------------
# Results provided by the main program are results from
# functions SOS_ANGLES, SOS_SURFACE, SOS_AEROSOLS, SOS_PROFIL and SOS 
#-----------------------------------------------------------------
#   SOS_ANGLES :
#      - File of angles used to BRDF/BPDF and radiance computations.
#      - File of angles used to the matrix phase function computations.
#
#   SOS_SURFACE :
#      - SURFACE BRDF/BPDF file (Glitter, Roujean, ...)
#        ==> Already existing files are not generated again.
#      - Log file SURFACE (optional)
#      - Glitter computation creates 3 temporary files which are deleted
#        at the end of computations (files: RES_GSF, RES_FRESNEL et RES_MAT_REFLEX)
#
#   SOS_AEROSOLS :
#      - MIE files
#        ==> Already existing files are not generated again.
#      - Log file of MIE calculations (optionnal)
#      - Log file of phaze matrix calculations (optionnal)
#      - Aerosols file result : aerosols radiative properties
#        ==> Already existing files are overwritted
#
#   SOS_PROFIL :
#      - Log file of PROFIL calculations (optionnal)
#      - Profile file result : Optical thickness profile and proportion of molecules and aerosols
#                           (without aerosol truncature adjustment)
#        ==> Already existing files are overwritted
# 
#   SOS :
#      - Log file of SOS calculations (optionnal)
#      - Radiance files result : binary one
#                                ascii upward radiance
#                                ascii downward radiance
#        ==> Already existing files are overwritted
#                                                                    
#########################################################################

## Loading of simulation parameters
############################################
Arg_ANGLES=" "
Arg_SURFACE=" "
Arg_AEROSOLS=" "
Arg_PROFIL=" "
Arg_SOS=" "

#Number of arguments
#-------------------
nbarg=$#

#Initialization of status for specific arguments
#required for the script main_SOS.ksh
#-----------------------------------------------
SOS_Wa_status=0
SOS_Config_status=0
AP_UserFile_status=0
AP_ResFile_status=0
AER_Waref_status=0
AER_AOTref_status=0
AER_UserFile_status=0
AER_ResFile_status=0
AER_Model_status=0
AER_MMD_MRwa_status=0
AER_MMD_MIwa_status=0
AER_MMD_MRwaref_status=0
AER_MMD_MIwaref_status=0
AER_BMD_CM_MRwaref_status=0
AER_BMD_CM_MIwaref_status=0
AER_BMD_FM_MRwaref_status=0
AER_BMD_FM_MIwaref_status=0
AER_ExtData_status=0
SURF_Type_status=0
SURF_File_status=0


#Arguments loading
#-----------------
err=0
numArg=0
while [ ${numArg} -lt ${nbarg} ] && [ ${err} -eq 0 ]
do
    #Keyword=$1  
    #Associated value=$2    

    #Get Keyword and Associated value the the current argument number
    case $1 in 
       -ANG.Thetas )           Arg_ANGLES=${Arg_ANGLES}" -ANG.Rad.Thetas $2" ;;   
 
       -ANG.Rad.NbGauss )      Arg_ANGLES=${Arg_ANGLES}" $1 $2";;			              

       -ANG.Rad.UserAngFile )  Arg_ANGLES=${Arg_ANGLES}" $1 $2";;
       
       -ANG.Rad.ResFile )      Arg_ANGLES=${Arg_ANGLES}" $1 $2"
                               Arg_SOS=${Arg_SOS}" -SOS.AngFile $2"
			       Arg_SURFACE=${Arg_SURFACE}" -SURF.AngFile $2";;
       
       -ANG.Aer.NbGauss )      Arg_ANGLES=${Arg_ANGLES}" $1 $2";;			              

       -ANG.Aer.UserAngFile )  Arg_ANGLES=${Arg_ANGLES}" $1 $2";;
       
       -ANG.Aer.ResFile )      Arg_ANGLES=${Arg_ANGLES}" $1 $2"
                               Arg_AEROSOLS=${Arg_AEROSOLS}" -AER.AngFile $2";;
			       
       -ANG.Log )              Arg_ANGLES=${Arg_ANGLES}" $1 $2";;
     
       -SOS.Wa )               waRadianceSimu=$2 
                               SOS_Wa_status=1 ;;		#Specific case

       -SOS.MDF )              Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.View )             Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.View.Phi )         Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.View.Dphi )        Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.OutputLevel )      Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.IGmax )            Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.Ipolar )           Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.ResBin )           Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.ResFileUp )        Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.ResFileDown )      Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.ResFileUp.UserAng )   Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.ResFileDown.UserAng ) Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SOS.Log )              Arg_SOS=${Arg_SOS}" $1 $2";;
      
       -SOS.Config )           Fic_Copy_Config=$2
                               SOS_Config_status=1;; #Specific case for the current main script

       -SOS.Trans ) 	       Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -AP.UserFile )          AP_UserFile=$2
       			       Arg_PROFIL=${Arg_PROFIL}" $1 $2"
			       AP_UserFile_status=1;; #Specific case for SOS_PROFIL and SOS_OS
       
       -AP.ResFile )           AP_ResFile=$2
       			       Arg_PROFIL=${Arg_PROFIL}" $1 $2"
			       AP_ResFile_status=1;; #Specific case for SOS_PROFIL and SOS_OS
       
       -AP.Log )               Arg_PROFIL=${Arg_PROFIL}" $1 $2";;
              	      
       -AP.MOT )               Arg_PROFIL=${Arg_PROFIL}" $1 $2";;
           
       -AP.HR )                Arg_PROFIL=${Arg_PROFIL}" $1 $2";;

       -AP.Type )              Arg_PROFIL=${Arg_PROFIL}" $1 $2";;
       
       -AP.AerHS.HA )          Arg_PROFIL=${Arg_PROFIL}" $1 $2";;
       
       -AP.AerLayer.Zmin )     Arg_PROFIL=${Arg_PROFIL}" $1 $2";;
       
       -AP.AerLayer.Zmax )     Arg_PROFIL=${Arg_PROFIL}" $1 $2";;
       
       -AER.Waref )            waref=$2
                               Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2"
                               AER_Waref_status=1;;   #Specific case
       
       -AER.AOTref )           AOTref=$2
       			       Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2"
                               AER_AOTref_status=1;;  #Specific case

       -AER.UserFile )         AER_UserFile=$2
       			       Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2"
			       AER_UserFile_status=1;;
       
       -AER.ResFile )          AER_ResFile=$2
                               Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2"
			       AER_ResFile_status=1;;  #Specific case 
       
       -AER.Log )              ficTraceAER=$2
       			       Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;

       -AER.MieLog )           Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.Tronca )           Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.Model )            AerModel=$2
                               AER_Model_status=1
                               Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
                  
       -AER.MMD.Mie.Filename ) Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.MMD.Mie.AlphaMax ) Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.MMD.MRwa )         mr_wa=$2
                               AER_MMD_MRwa_status=1;;  #Specific case
       
       -AER.MMD.MIwa )         mi_wa=$2
                               AER_MMD_MIwa_status=1;;  #Specific case
       
       -AER.MMD.MRwaref )      mr_waref=$2
                               AER_MMD_MRwaref_status=1;;  #Specific case
       
       -AER.MMD.MIwaref )      mi_waref=$2
                               AER_MMD_MIwaref_status=1;;  #Specific case
       

       
       -AER.MMD.SDtype )       Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.MMD.SDparam1 )     Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.MMD.SDparam2 )     Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.WMO.Model )        Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.WMO.DL )           Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.WMO.WS )           Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.WMO.OC )           Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.WMO.SO )           Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.SF.Model )         Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.SF.RH )            Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.BMD.VCdef )        Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2"
                               BMD_VCdef=$2;; #Specific case if reference and radiance 
			       				#simulation wavelengths are equal
       
       -AER.BMD.CoarseVC )     Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.BMD.FineVC )       Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.BMD.RAOT )         Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.BMD.CM.MRwa )      mr_LNDc_wa=$2;; #Specific case : directed to Arg_AEROSOLS depending 
                                               #    on reference or radiance simulation wavelength
       
       -AER.BMD.CM.MIwa )      mi_LNDc_wa=$2;; #Specific case : directed to Arg_AEROSOLS depending 
                                               #    on reference or radiance simulation wavelength
       
       -AER.BMD.CM.MRwaref )   mr_LNDc_waref=$2
                               AER_BMD_CM_MRwaref_status=1
                               Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";; #Specific case
       
       -AER.BMD.CM.MIwaref )   mi_LNDc_waref=$2
                               AER_BMD_CM_MIwaref_status=1
       			       Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";; #Specific case
       
       -AER.BMD.CM.SDradius )  Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";; 
       
       -AER.BMD.CM.SDvar )     Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
     
       
       -AER.BMD.FM.MRwa )      mr_LNDf_wa=$2;; #Specific case : directed to Arg_AEROSOLS depending 
                                               #    on reference or radiance simulation wavelength
       
       -AER.BMD.FM.MIwa )      mi_LNDf_wa=$2;; #Specific case : directed to Arg_AEROSOLS depending 
                                               #    on reference or radiance simulation wavelength
       
       -AER.BMD.FM.MRwaref )   mr_LNDf_waref=$2
                               AER_BMD_FM_MRwaref_status=1
                               Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";; #Specific case
       
       -AER.BMD.FM.MIwaref )   mi_LNDf_waref=$2
                               AER_BMD_FM_MIwaref_status=1
       			       Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";; #Specific case
       
       -AER.BMD.FM.SDradius )  Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";; 
       
       -AER.BMD.FM.SDvar )     Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
       
       -AER.ExtData )          AER_ExtData_status=1
                               Arg_AEROSOLS=${Arg_AEROSOLS}" $1 $2";;
             
       -SURF.Log )             Arg_SURFACE=${Arg_SURFACE}" $1 $2";;
          
       -SURF.File )            ficSURFACE=$2
       			       SURF_File_status=1
                               Arg_SOS=${Arg_SOS}" $1 $2";; #Specific case
       
    
       -SURF.Type )            typeSurf=$2
                               SURF_Type_status=1
                               Arg_SURFACE=${Arg_SURFACE}" $1 $2"
                               Arg_SOS=${Arg_SOS}" $1 $2";; #Specific case
       
       -SURF.Alb )             Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SURF.Ind )             Arg_SURFACE=${Arg_SURFACE}" $1 $2"
                               Arg_SOS=${Arg_SOS}" $1 $2";;
              
       -SURF.Glitter.Wind )    Arg_SURFACE=${Arg_SURFACE}" $1 $2"
                               Arg_SOS=${Arg_SOS}" $1 $2";;     
      
       -SURF.Roujean.K0 )      Arg_SURFACE=${Arg_SURFACE}" $1 $2"
                               Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SURF.Roujean.K1 )      Arg_SURFACE=${Arg_SURFACE}" $1 $2"
                               Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SURF.Roujean.K2 )      Arg_SURFACE=${Arg_SURFACE}" $1 $2"
                               Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SURF.Nadal.Alpha )     Arg_SURFACE=${Arg_SURFACE}" $1 $2"
                               Arg_SOS=${Arg_SOS}" $1 $2";;
       
       -SURF.Nadal.Beta )      Arg_SURFACE=${Arg_SURFACE}" $1 $2"
                               Arg_SOS=${Arg_SOS}" $1 $2";;
    esac
    
    #Test if the "value" attribute is not missing (and then correspond to a keyword or null value)
    case $2 in 
           
       -ANG.Thetas )           err=1;;     	
       -ANG.Rad.NbGauss )      err=1;;			              
       -ANG.Rad.UserAngFile )  err=1;;
       -ANG.Rad.ResFile )      err=1;;
       -ANG.Aer.NbGauss )      err=1;;			              
       -ANG.Aer.UserAngFile )  err=1;;   
       -ANG.Aer.ResFile )      err=1;;
       -ANG.Log )              err=1;;  
       -SOS.Wa )               err=1;; 
       -SOS.MDF )              err=1;;      
       -SOS.View )             err=1;;       
       -SOS.View.Phi )         err=1;;       
       -SOS.View.Dphi )        err=1;;       
       -SOS.OutputLevel )      err=1;;      
       -SOS.IGmax )            err=1;;      
       -SOS.Ipolar )           err=1;;       
       -SOS.ResBin )           err=1;;       
       -SOS.ResFileUp )        err=1;;      
       -SOS.ResFileDown )      err=1;;
       -SOS.ResFileUp.UserAng )   err=1;;      
       -SOS.ResFileDown.UserAng ) err=1;;         
       -SOS.Log )              err=1;;      
       -SOS.Config )           err=1;; 
       -SOS.Trans ) 	       err=1;;      
       -AP.UserFile )          err=1;; 
       -AP.ResFile )           err=1;; 
       -AP.Log )               err=1;;              	      
       -AP.MOT )               err=1;;           
       -AP.HR )                err=1;;
       -AP.Type )              err=1;;       
       -AP.AerHS.HA )          err=1;;       
       -AP.AerLayer.Zmin )     err=1;;       
       -AP.AerLayer.Zmax )     err=1;;       
       -AER.Waref )            err=1;;
       -AER.AOTref )           err=1;;
       -AER.UserFile )         err=1;;       
       -AER.ResFile )          err=1;;
       -AER.Log )              err=1;;
       -AER.MieLog )           err=1;;       
       -AER.Tronca )           err=1;;       
       -AER.Model )            err=1;;
       -AER.MMD.Mie.Filename ) err=1;;
       -AER.MMD.Mie.AlphaMax ) err=1;;
       -AER.MMD.MRwa )         err=1;;
       -AER.MMD.MIwa )         err=1;;
       -AER.MMD.MRwaref )      err=1;; 
       -AER.MMD.MIwaref )      err=1;;
       -AER.MMD.SDtype )       err=1;;
       -AER.MMD.SDparam1 )     err=1;;
       -AER.MMD.SDparam2 )     err=1;;
       -AER.WMO.Model )        err=1;;
       -AER.WMO.DL )           err=1;;
       -AER.WMO.WS )           err=1;;
       -AER.WMO.OC )           err=1;;
       -AER.WMO.SO )           err=1;;
       -AER.SF.Model )         err=1;;
       -AER.SF.RH )            err=1;;
       -AER.BMD.VCdef )        err=1;;
       -AER.BMD.CoarseVC )     err=1;;
       -AER.BMD.FineVC )       err=1;;
       -AER.BMD.RAOT )         err=1;;
       -AER.BMD.CM.MRwa )      err=1;; 
       -AER.BMD.CM.MIwa )      err=1;; 
       -AER.BMD.CM.MRwaref )   err=1;; 
       -AER.BMD.CM.MIwaref )   err=1;; 
       -AER.BMD.CM.SDradius )  err=1;; 
       -AER.BMD.CM.SDvar )     err=1;;
       -AER.BMD.FM.MRwa )      err=1;; 
       -AER.BMD.FM.MIwa )      err=1;; 
       -AER.BMD.FM.MRwaref )   err=1;; 
       -AER.BMD.FM.MIwaref )   err=1;; 
       -AER.BMD.FM.SDradius )  err=1;; 
       -AER.BMD.FM.SDvar )     err=1;;   
       -AER.ExtData )          err=1;;   
       -SURF.Log )             err=1;;
       -SURF.File )            err=1;;     
       -SURF.Type )            err=1;; 
       -SURF.Alb )             err=1;;
       -SURF.Ind )             err=1;;
       -SURF.Glitter.Wind )    err=1;;  
       -SURF.Roujean.K0 )      err=1;;
       -SURF.Roujean.K1 )      err=1;;
       -SURF.Roujean.K2 )      err=1;;
       -SURF.Nadal.Alpha )     err=1;;
       -SURF.Nadal.Beta )      err=1;;   
       ""   )                  err=1;;   
    esac  #Test if the Value attribute is not missing
    
    if [ ${err} -ne 0 ] 
    then  
       echo "ERROR ON PARAMETERS: " 
       echo "   The Value associated to $1 is probably missing"
       exit 1
    fi

    #Next arguments $1 and $2 
    shift 
    shift
    
    #Incrementation of current argument number
    numArg=$((${numArg}+2))
    
done #Loop on arguments



#Check the status for specific arguments required
#------------------------------------------------
#waRadianceSimu
if [ "${SOS_Wa_status}" = "0" ]
then
    echo "ERROR ON PARAMETERS: Required parameter missing :" 
    echo "   The parameter -SOS.Wa  has to be defined !"
    exit 1
fi


# --- Check for the choice of profile data file: user or computed one
if [ ${AP_UserFile_status} -eq 1 ] 
then #UserFile is defined
    if [ ${AP_ResFile_status} -eq 1 ] 
    then #ResFile is defined
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameters -AP.UserFile and -AP.ResFile cannot be defined simultaneously !"
       exit 1
    else #Only UserFile is defined
       ficPROFIL=USERFILE
    fi

    flagWaDiff=`echo ${waref} ${waRadianceSimu} | awk '{if ($1 == $2) print 0 ; else print 1;}'`
    if [ ${flagWaDiff} -eq 1 ] 
    then #Wavelength for aerosol definition and for radiance simulation are different : 
         #the AOT in the user profile file cannot be adjusted to another wavelength
       echo "ERROR ON PARAMETERS: " 
       echo "   When using a user profile (-AP.UserFile) "
       echo "   the reference aerosol wavelength (-AER.Waref)"
       echo "   and the radiance simulation wavelength (-SOS.Wa) have to be equal"
       exit 1
    fi

else #UserFile is not defined
    if [ ${AP_ResFile_status} -eq 1 ] 
    then #Only ResFile is defined
       ficPROFIL=DEFAULT
    else #None UserFile or ResFile is defined
       echo "ERROR ON PARAMETERS: " 
       echo "   One of theses parameters -AP.UserFile or -AP.ResFile has to be defined !"
       exit 1
    fi
fi



# --- Check for the choice of aerosol data file: user or computed one
if [ ${AER_UserFile_status} -eq 1 ] 
then #UserFile is defined
    if [ ${AER_ResFile_status} -eq 1 ] 
    then #ResFile is defined
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameters -AER.UserFile and -AER.ResFile cannot be defined simultaneously !"
       exit 1
    else #Only UserFile is defined
       ficAERgranu=USERFILE
    fi
else #UserFile is not defined
    if [ ${AER_ResFile_status} -eq 1 ] 
    then #Only ResFile is defined
       ficAERgranu=DEFAULT
    else #None UserFile or ResFile is defined
       echo "ERROR ON PARAMETERS: " 
       echo "   One of theses parameters -AER.UserFile or -AER.ResFile has to be defined !"
       exit 1
    fi
fi

# --- Check the wavelength for aerosols reference is the same than for radiance simulation
# --- when using an external file for aerosol phaze functions
if [ ${AER_ExtData_status} -eq 1 ] 
then #An external phaze function is used
    flagWaDiff=`echo ${waref} ${waRadianceSimu} | awk '{if ($1 == $2) print 0 ; else print 1;}'`
    if [ ${flagWaDiff} -eq 1 ] 
    then #Wavelength for aerosol definition and for radiance simulation are different : 
         #no information for the radiance simulation wavelength
       echo "ERROR ON PARAMETERS: " 
       echo "   When using an external aerosol phaze function (-AER.ExtData) "
       echo "   the reference aerosol wavelength (-AER.Waref)"
       echo "   and the radiance simulation wavelength (-SOS.Wa) have to be equal"
       exit 1
    fi
fi


#AOTref
if [ "${AER_AOTref_status}" = "0" ] 	
then
    echo "ERROR ON PARAMETERS: Required parameter missing :" 
    echo "   The parameter -AER.AOTref  has to be defined !"
    exit 1
fi


#Check required parameters for main_SOS.ksh process

#--  flagAOTref : 0 if AOTref=0, 1 if AOTref not null
flagAOTref=`echo $AOTref | awk '{if ($1 == 0.) print 0 ; else print 1;}'`   
#--  flagWaDiff : 0 if waref = waRadianceSimu, 1 if wavelengths are different
flagWaDiff=`echo ${waref} ${waRadianceSimu} | awk '{if ($1 == $2) print 0 ; else print 1;}'` 

#waref         
if [ ${flagAOTref} -eq 1 ] && [ "${AER_Waref_status}" = "0" ]
then
    echo "ERROR ON PARAMETERS: " 
    echo "   The parameter -AER.Waref  has to be defined when AER.AOTref is not null !"
    exit 1
fi

#AerModel
if [ ${flagAOTref} -eq 1 ] && [ "${AER_Model_status}" = "0" ] && [ "${AER_UserFile_status}" = "0" ]
then
    echo "ERROR ON PARAMETERS: " 
    echo "   The parameter -AER.Model  has to be defined when AER.AOTref is not null !"
    exit 1
fi

if [ ${flagAOTref} -eq 1 ] && [ "${AerModel}" = "0" ]
then

    #mr_wa
    if [ "${AER_MMD_MRwa_status}" = "0" ] 
    then   
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameter -AER.MMD.Mie.MRwa  has to be defined when AER.AOTref is not null and AER.Model=0 !"
       exit 1
    fi

    #mi_wa
    if [ "${AER_MMD_MIwa_status}" = "0" ] 
    then   
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameter -AER.MMD.Mie.MIwa  has to be defined when AER.AOTref is not null and AER.Model=0 !"
       exit 1
    fi
        
    #mr_waref
    if [ "${AER_MMD_MRwaref_status}" = "0" ] && [ ${flagWaDiff} -eq 1 ]
    then   
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameter -AER_MMD_MRwaref  has to be defined when AER.AOTref is not null, AER.Model=0 and AER.Waref different of SOS.Wa !"
       exit 1
    fi
    
    #mi_waref
    if [ "${AER_MMD_MIwaref_status}" = "0" ] && [ ${flagWaDiff} -eq 1 ]
    then   
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameter -AER_MMD_MIwaref  has to be defined when AER.AOTref is not null, AER.Model=0 and AER.Waref different of SOS.Wa !"
       exit 1
    fi    
fi

if [ ${flagAOTref} -eq 1 ] && [ "${AerModel}" = "3" ]
then
    #mr_LNDc_waref
    if [ "${AER_BMD_CM_MRwaref_status}" = "0" ] && [ ${flagWaDiff} -eq 1 ]
    then   
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameter -AER.BMD.CM.MRwaref  has to be defined when AER.AOTref is not null, AER.Model=3 and AER.Waref different of SOS.Wa !"
       exit 1
    fi
    
    #mi_LNDc_waref
    if [ "${AER_BMD_CM_MIwaref_status}" = "0" ] && [ ${flagWaDiff} -eq 1 ]
    then   
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameter -AER.BMD.CM.MIwaref  has to be defined when AER.AOTref is not null, AER.Model=3 and AER.Waref different of SOS.Wa !"
       exit 1
    fi

    #mr_LNDf_waref
    if [ "${AER_BMD_FM_MRwaref_status}" = "0" ] && [ ${flagWaDiff} -eq 1 ]
    then   
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameter -AER.BMD.FM.MRwaref  has to be defined when AER.AOTref is not null, AER.Model=3 and AER.Waref different of SOS.Wa !"
       exit 1
    fi
    
    #mi_LNDf_waref
    if [ "${AER_BMD_FM_MIwaref_status}" = "0" ] && [ ${flagWaDiff} -eq 1 ]
    then   
       echo "ERROR ON PARAMETERS: " 
       echo "   The parameter -AER.BMD.FM.MIwaref  has to be defined when AER.AOTref is not null, AER.Model=3 and AER.Waref different of SOS.Wa !"
       exit 1
    fi        
fi


#typeSurf
if [ "${SURF_Type_status}" = "0" ] 
then
    echo "ERROR ON PARAMETERS: Required parameter missing :" 
    echo "   The parameter -SURF.Type  has to be defined !"
    exit 1
fi

#ficSURFACE
if [ "${SURF_File_status}" = "0" ] && [ "${typeSurf}" != "0" ] && [ "${typeSurf}" != "2" ]
then
    echo "ERROR ON PARAMETERS: Required parameter missing :" 
    echo "   The parameter -SURF.File  has to be defined !"
    exit 1
fi



 
#Running calculations
#--------------------

### 1) Definition of angles to be used for simulations
       Call_def_angles

### 2) Simulation of surface radiative properties
       Call_simu_surface
  
### 3) Simulation of aerosols radiative properties       
       flagAOT=`echo $AOTref | awk '{if ($1 == 0.) print 0 ; else print 1;}'`       
       if [ ${flagAOT} -eq 1 ] && [ "-$ficAERgranu" = "-DEFAULT" ]  ; then
	   
	   #-----------------------------------------------------------------------
	   # Case 1 : AOTref not null and automatic calculation (not userfile used)
	   #-----------------------------------------------------------------------
	   echo "--> Aerosols case 1 : AOTref not null and automatic calculation (not userfile used)"

	   if [ ${AER_ExtData_status} -eq 1 ] ; then #An external phaze function is used
		echo "--> Use of an external aerosols phaze function" ;
           fi

	   # Simulation for the reference wavelength 
	   # ==> Calculation of the extinction coefficient for the reference wavelength
	   #                                                       ********************
	   #     ====> Generating : File $AER_ResFile
	      Arg_AEROSOLS_Simu=${Arg_AEROSOLS}" -AER.Wa ${waref}"
	      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.AOT ${AOTref}"
	      
	      #Parameters for SOS_AEROSOLS program
	      # if wa != waref ==> there are parameters defined for waref
	      # if wa  = waref ==> there is not parameter defined for waref : parameters are only defined for wa
	      #Test if wa != waref
	      flag=`echo ${waref} ${waRadianceSimu} | awk '{if ($1 == $2) print 0 ; else print 1;}'`  
	      if [ ${flag} -eq 1 ] ; then  #wa != waref
	      
	          if [ ${AerModel} -eq 0 ] ; then
	              Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.MMD.MRwa ${mr_waref}"
		      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.MMD.MIwa ${mi_waref}"	          
	          fi
	      
	          if [ ${AerModel} -eq 3 ] ; then
	              Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.CM.MRwa ${mr_LNDc_waref}"
		      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.CM.MIwa ${mi_LNDc_waref}"
		      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.FM.MRwa ${mr_LNDf_waref}"
		      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.FM.MIwa ${mi_LNDf_waref}"		  		  
	          fi
		  
	      else #wa = waref
	      
	          if [ ${AerModel} -eq 0 ] ; then
	              Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.MMD.MRwa ${mr_wa}"
		      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.MMD.MIwa ${mi_wa}"	          
	          fi
	      
	          if [ ${AerModel} -eq 3 ] ; then
	              Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.CM.MRwa ${mr_LNDc_wa}"
		      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.CM.MIwa ${mi_LNDc_wa}"
		      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.FM.MRwa ${mr_LNDf_wa}"
		      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.FM.MIwa ${mi_LNDf_wa}"	
		      
		      if [ ${BMD_VCdef} -eq 2 ] ; then  #Specific case for volumetric concentrations 
		      						    #defined by the ratio of AOTcoarse / AOTtot 
								    #for the reference wavelength :
								    #Only one definition of refractive indexes for the
								    #radiance simulation wavelength but values for the 
								    #reference wavelength are expected by SOS_AEROSOLS.
		         Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.CM.MRwaref ${mr_LNDc_wa}"
		         Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.CM.MIwaref ${mi_LNDc_wa}"
		         Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.FM.MRwaref ${mr_LNDf_wa}"
		         Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.FM.MIwaref ${mi_LNDf_wa}"
		      fi	  		  
	          fi
		  
	      fi #End test if wa != waref
	      
	      # --> Call SOS_AEROSOLS program : creating the file ${AER_ResFile}
	      Call_simu_aer 
	      
	      
	      #Test if wa != waref
	      flag=`echo ${waref} ${waRadianceSimu} | awk '{if ($1 == $2) print 0 ; else print 1;}'`  
	      if [ ${flag} -eq 1 ] ; then
	      
	         #---------------------------------------------------------------------
	         # Case 1a : wa != waref 
	         #---------------------------------------------------------------------
		 echo "--> Aerosols case 1a : wa != waref" 
	         # Storage of the extinction coefficient for the reference wavelength
		 kmat1_waref=`head -1 ${AER_ResFile} | cut -d":" -f2 | awk '{print $1}'`
		
	         # Simulation for the wavelength of radiance computation
	         # ==> Calculation of the extinction coefficient for the wavelength of radiance computation
	         #                                                       **********************************
	         # ==> Calculation of the phaze matrix and et troncature coefficient.
	         #     ====> Generating : File $AER_ResFile
	         
		 Arg_AEROSOLS_Simu=${Arg_AEROSOLS}" -AER.Wa ${waRadianceSimu}"
	         Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.AOT 0.1"
		 	#Not null AOT to allow a call of the function SOS_AEROSOL
	                #=> Calculation of the extinction coefficient 
	       		     #The real AOT value for the wavelength of radiance computation
			     #will be deduced from this extinction coefficient and this one 
			     #calculated for the reference wavelength
	         if [ ${AerModel} -eq 0 ] ; then
	             Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.MMD.MRwa ${mr_wa}"
		     Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.MMD.MIwa ${mi_wa}"	          
	         fi
	      
	         if [ ${AerModel} -eq 3 ] ; then
	             Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.CM.MRwa ${mr_LNDc_wa}"
		     Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.CM.MIwa ${mi_LNDc_wa}"
		     Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.FM.MRwa ${mr_LNDf_wa}"
		     Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.BMD.FM.MIwa ${mi_LNDf_wa}"		  
	         fi
		 
	   
	         # Storage of the log file generated from the reference wavelength
		 if [ "-=${ficTraceAER}" != "-0" ] ; then
		   cp ${ficTraceAER} ./.tempo1.txt ;
		 fi
		 
	         ficAERgranu=DEFAULT
		 # --> Call SOS_AEROSOLS program : creating the file ${AER_ResFile}
	         Call_simu_aer 
	   
	         # Storage of the extinction coefficient for the wavelength of radiance computation
	         kmat1_wa=`head -1 ${AER_ResFile} | cut -d":" -f2 | awk '{print $1}'`	     
	         
	         # Calculation of the AOT for the wavelength of radiance computation
	         AOT=`echo $AOTref $kmat1_wa $kmat1_waref |  awk '{aot = $1*$2/$3; printf("%6.4f",aot)}'`
		 
	         # Paste of log files : 
		 #   first for aerosol calculation to the reference wavelength
		 #   first for aerosol calculation to the wavelength of radiance computation
		 if [ "-=${ficTraceAER}" != "-0" ] ; then
		   cp ${ficTraceAER} ./.tempo2.txt 
		   echo " ################################################ " >  ${ficTraceAER}
		   echo " # CALCULS POUR LA LONGUEUR D'ONDE DE REFERENCE # " >> ${ficTraceAER}
		   echo " ################################################ " >>  ${ficTraceAER}
		   echo " "          >> ${ficTraceAER}
		   cat ./.tempo1.txt >> ${ficTraceAER}
		   echo " "          >> ${ficTraceAER}
		   echo " "          >> ${ficTraceAER}
		   echo " "          >> ${ficTraceAER}
		   echo " ######################################################################### " >> ${ficTraceAER}
		   echo " # CALCULS POUR LA LONGUEUR D'ONDE DE SIMULATION DU CHAMP DE RAYONNEMENT # " >> ${ficTraceAER}
		   echo " ######################################################################### " >> ${ficTraceAER}
		   echo " "          >> ${ficTraceAER}
		   cat ./.tempo2.txt >> ${ficTraceAER}
		   echo " "          >> ${ficTraceAER}
		   echo " Estimation de l'epaisseur optique aérosols pour la longueur d'onde de simulation" >> ${ficTraceAER}
		   echo " --------------------------------------------------------------------------------" >> ${ficTraceAER}
		   echo "    - Longueur d'onde de reference (microns) : " ${waref}			    >> ${ficTraceAER}
		   echo "    - Section efficace d'extinction pour la longueur d'onde de reference (micron^2) : " ${kmat1_waref} >> ${ficTraceAER}
		   echo "    - Epaisseur optique aérosols pour la longueur d'onde de reference : " ${AOTref}		>> ${ficTraceAER}
		   echo "    - Longueur d'onde de simulation du champ de rayonnement (microns) : " ${waRadianceSimu}	>> ${ficTraceAER}
		   echo "    - Section efficace d'extinction pour la longueur d'onde de simulation (micron^2) : " ${kmat1_wa}	>> ${ficTraceAER}
		   echo "         ---> Epaisseur optique aérosols pour la longueur d'onde de simulation : " ${AOT}		>> ${ficTraceAER}
		   \rm  ./.tempo1.txt  ./.tempo2.txt ;
		 fi #Log file case
		 
	      else #case : wa = waref
	         echo "--> Aerosols case 1b : wa = waref" 
		 # Storage of the AOT value (in case of wa = waref)
	         AOT=${AOTref}
	      fi #test if wa != waref	
	       
       else #End of Case 1
       
           #---------------------------------------------------------------------
	   # Case 2 : AOTref = 0 and/or use of an aerosol userfile 
	   #---------------------------------------------------------------------           
           if [ ${flagAOT} -eq 1 ] ; then  #AOTref not nul and use of an aerosol userfile 
	      echo "--> Aerosols case 2 : AOTref not nul and use of an aerosol userfile "
	      #Test if wa != waref
	      flag=`echo ${waref} ${waRadianceSimu} | awk '{if ($1 == $2) print 0 ; else print 1;}'`  
	      if [ ${flag} -eq 1 ] ; then
	           echo "--> Aerosols case 2a : wa != waref"
		   echo " "
	           echo "   USE CASE ERROR :  (not allowed) "
		   echo "      - Use of an aerosol userfile : " ${ficAERgranu}
		   echo "      - The wavelength of radiance computation (${waRadianceSimu} microns) is different of the reference one (${waref} microns)"
		   echo "      ==> It is not possible to estimate the aerosol optical thickness for the wavelength of radiance computation !"
		   echo " "
		   echo "   THE PROGRAM IS STOPPED"
		   echo " "
		   echo "   NB : In order to use a precalculated userfile of the aerosol radiative properties,"
		   echo "        it is necessary to define the same wavelength for AOT reference and radiance computation."
		   echo " "
		   exit 1;
	      else
	         echo "--> Aerosols case 2b : wa = waref"
	         AOT=${AOTref}
	      fi
	      
	   else 
	      #---------------------------------------------------------------------
	      # Case : AOTref = 0 ==> AOT = 0  for the wavelength of radiance computation.
	      # ==> The data of the result aerosol file are forced to be zero 
	      # (the userfile is neglected : a result file has be defined)
	      #---------------------------------------------------------------------
	      echo "--> Aerosols case 3 : AOTref = 0"
	      Arg_AEROSOLS_Simu=${Arg_AEROSOLS}" -AER.Wa ${waRadianceSimu}"
	      Arg_AEROSOLS_Simu=${Arg_AEROSOLS_Simu}" -AER.AOT 0.0"
	      AOT=0.0
	      ficAERgranu=DEFAULT
	      Call_simu_aer;
	   fi
	   
       fi  #End of else bloc for test on 'AOTref not null and not userfile is used'	         
       
       	  
### 4) Simulation of atmospheric profile (molecules / aerosols)
       Call_simu_profil
       
       
### 5) SOS simulation of diffuse radiance for the defined wavelength of radiance computation
       #Definition of aerosol file
       if [ "-$ficAERgranu" = "-DEFAULT" ]  ; then
          Arg_SOS=${Arg_SOS}" -SOS.AerFile ${AER_ResFile}";
       else
          Arg_SOS=${Arg_SOS}" -SOS.AerFile ${AER_UserFile}";
       fi
       
       #Definition of profile file
       if [ "-$ficPROFIL" = "-DEFAULT" ] ; then
          Arg_SOS=${Arg_SOS}" -SOS.ProfileFile ${AP_ResFile}";
       else
          Arg_SOS=${Arg_SOS}" -SOS.ProfileFile ${AP_UserFile}";
       fi
  
       #Definition of storage directory for BRDF/BPDF files
       Arg_SOS=${Arg_SOS}" -SURF.Dir $repSURFACE"
       Call_simu_os
              
	      
### 6) Copy of the parameters file used to the compilation of executables
if [ "${SOS_Config_status}" = "1" ]
then
    echo "-->Copy of file : $SOS_RACINE/exe/config.txt to ${Fic_Copy_Config}"
    cp $SOS_RACINE/exe/config.txt ${Fic_Copy_Config}
fi

