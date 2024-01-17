# Successive Orders of Scattering method coupled with the gaseous ABSorptions (SOS-ABS V5.1 - November 9, 2023) 
## Overview
The SOS-ABS code (Successive Orders of Scattering method coupled with the gaseous ABSorptions) is a 1D plan parallel radiative transfer code simulating the polarized radiance of the {Earth surface  atmosphere} system, including the coupling effect of gaseous absorption. It works for the reflective spectral domain, from 2500 to 27500 $cm^{-1}$ (0.364 to 4 $\mu m$), under clear sky conditions.

It is an heritage of the original Successive Orders of Scattering method from the LOA laboratory [[Deuze et al, 1989](https://doi.org/10.1016/0022-4073(89)90118-0) ; [Lenoble et al., 2007](https://doi.org/10.1016/j.jqsrt.2007.03.010)]

The SOS-ABS model includes the initial features and improved methods. SOS-ABS allows simulating:
- **Atmospheric profiles**:  
For the characterization of the atmosphere, the user can define the molecular and aerosol optical thickness.  
The user can also define the amount of gas among ozone (O3), water vapour (H2O), carbon dioxide (CO2) and methane (CH4). The dioxygen (O2), nitrous oxide (N2O), nitrogen dioxide (NO2) and carbon monoxide (CO) are also considered. The vertical distribution is based on typical profiles (Tropical, Middle Latitude Summer or Winter, South Arctic Summer or Winter, US Standard), or given by a user's atmospheric profile.
- **Aerosol models**:  
A wide variety of aerosol size distribution is available in the SOS-ABS model to make the optical properties of the atmosphere highly close to real-world conditions. Log-normal (LND) or Junge mono-modal size distributions, bimodal LND, pre-calculated WMO [WMO, 1986] or Shettle & Fenn models [Shettle and Fenn, 1979] can be used. It is also possible to simulate a combination of various mono-modal modes, as defined by the user, which is compliant with the use of [CAMS](https://doi.org/10.21957/84ya94mls) data.  
Users can also use their own aerosol phase function and radiative properties.
- **Surface reflection models**:  
The reflection of the direct solar beam and of the diffuse light on the surface is considered for both land and water surfaces.  
Over land, a Lambertian target can be simulated as well as the Roujean's BRDF model [[Roujean et al., 1992](https://doi.org/10.1029/92JD01411)]. The polarization of the surface can be considered too with the Rondeaux’s model [[Rondeaux and Herman, 1991](https://doi.org/10.1016/0034-4257(91)90072-E)], the Bréon’s model [[Breon et al., 1995](https://doi.org/10.1109/TGRS.1995.8746030)] or the Maignan's model [[Maignan et al., 2009](https://doi.org/10.1016/j.rse.2009.07.022)].  
Over a sea surface, the air / sea interface can be modelled either for a flat surface or by considering the sea roughness defined by the wind speed and the correlated waves [[Cox and Munk, 1954](https://doi.org/10.1364/JOSA.44.000838)].

The user can define specific angles for which output simulated radiance are required.

By default, SOS-ABS provides the upward radiance field at TOA and the downward radiance at ground level. However, it is also possible to get the upward and downward radiance, in term of intensity and polarized light, for any given altitude in the atmosphere.
## Documentation

- [User Manual](https://github.com/CNES/RadiativeTransferCode-SOS/blob/master/doc/SOS-ABS-User_Manual_V1.0.pdf)
- [Example of SOS-ABS simulation from command lines](https://github.com/CNES/RadiativeTransferCode-SOS/blob/exe/runSOS-ABS_demo.ksh)


## FAQ
To report an issue or interact with the development team, please raise an issue on this repository. 
## Contributors
This project has been supported by CNES.

The approach to modelling the gaseous transmissions was established in collaboration with the LOA laboratory (Pr. Phillipe Dubuisson and Pr. Jerome Riedi). 

The CKD coefficients were calculated by the HYGEOS company (Dr. Mathieu Compiègne).

The specification and implementation for SOS-ABS of the absorption and scattering coupling effect were made by CS GROUP France (Dr. Lafrance Bruno and Dr. Xavier Lenot).

## References 
 
This code can be referenced as :  
- J. Lenoble, M. Herman, J.L. Deuze, B. Lafrance, R. Santer, D. Tanre, 
_A successive order of scattering code for solving the vector equation of transfer in the earth's atmosphere with aerosols_, Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 107, Issue 3, 2007. [doi:10.1016/j.jqsrt.2007.03.010](https://doi.org/10.1016/j.jqsrt.2007.03.010).

## License 
Copyright 2023, Centre National d'Etudes Spatiales (CNES)

SOS-ABS code is licensed under the GNU General Public License v3.0.  
SOS-ABS is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
