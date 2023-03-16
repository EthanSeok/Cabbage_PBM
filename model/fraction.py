""" Fractionation of radiation, LAI, Vcmax
    Input  : doy, latitude, Irrad(W m-2), LAI, Vcmax, hour, pressure
    Output : Ibeam, Idiffuse, Isunlit, Ishaded, LAIsunlit, LAIshaded, VcmaxSunlit, VcmaxShaded
             (Ibeam, Idiffuse, Isunlit, Ishaded : PAR umol m-2 s-1)
"""
##### Constants
KELVIN   = 273.15       # Kelvin temperature
R        = 8.314        # ideal gas constant
SBC      = 5.6697e-8    # Stefan-Boltzmann constant   W m-2 K-4
PSC      = 6.66e-4      # psychrometer constant
Cp       = 29.3         # specific heat of air    J mol-1 C-1
LAMBDA   = 44000.       # Latent heat of vaporization at 25 C J mol-1import numpy as np

### Constant or parameters
SC      = 1361.     # Solar conatant   W m-2
conv    = 4.57      # conversion from W m-2 to umol m-2 s-1
Po      = 101.3     # atomopheric pressure at sea level  kPa 

### Parameters of Chinese cabbage
Nb      = 25         # residual leaf N content (abour 0.5% N) mmol m-2  from de Pury(1997) 
No      = 137        # leaf N content of canopy top         mmol m-2
Nl      = 120        # leaf N content per unit leaf area    mmol m-2
kn      = 0.713      # coefficient for leaf N allocation             
# Vcmax = 110        # rubisco capacity by calculation      umol m-2 s-1 from de Pury(1997) 
# Xn = Vcmax/(Nl-Nb) # ratio of rubisco capacity to leaf N 

### clumping factor
clump   = 0.5

############################################################################################################
###  Irradiance(W m-2), LAI, Vcmax fractionation of sunlit and shaded parts of canopy
###  from de Pury (1997)
############################################################################################################

import numpy as np

class Fractionation():
    
    def __init__(self, latitude, press=100.0):
        latitude = float(latitude)
        self.lat = np.radians(latitude)                                  # lat = radians of latitude
        self.P    = press # pressure(kPa) about 101.3                    # atomopheric pressure   kPa
        
        # Output (Irrad)
        self.It = 0.0        # PAR total
        self.Ib = 0.0        # beam fraction of PAR
        self.Id = 0.0        # diffuse fractio of PAR
        self.Icsun = 0.0     # PAR of sunlit leaf
        self.Icsh  = 0.0     # PAR of shaded leaf
        # Output (lai)
        self.laiSun = 0.0    # LAI for sunlit leaf
        self.laiSh  = 0.0    # lai for shaded leaf
        # Output (rubisco)
        self.VcmaxSun = 0.0  # rubisco for sunlit leaf
        self.VcmaxSh  = 0.0  # rubisco for shaded leaf
                
    def radFraction(self, doy, hour, PPFD, LAI):   # from de Pury(1997) method
        ##################################################################################################
        ## parameters for radiation fractionation from de Pury(1997)
        a       = 0.72    # atmospheric transmission coefficients of PAR  
        fa      = 0.426   # forward scattering coeff of PAR
        rhocd   = 0.036   # canopy reflection coeff for diffuse PAR
        rhoh    = 0.04    # reflection coeff of a canopy with horizontal leaves
        rhol    = 0.10    # reflection coeff for PAR
        taul    = 0.05    # leaf transmissivity to PAR
        sigma   = 0.15    # leaf scattering coeff for PAR (rhol + taul)
        kd      = 0.78 * clump   # diffuse PAR extinction coeff.             0.5 <- clumping
        kdprime = 0.719   # diffuse & scattered diffuse OAR ext coeff  
        
        ## calc values
        decl = -0.4093 * np.cos(2 * np.pi * (doy + 10) / 365)  # sun declination   rad
        ha   = np.pi / 12 * (hour - 12)                        # hour angle   rad
        sin_a = np.sin(decl)*np.sin(self.lat)
        cos_b = np.cos(decl)*np.cos(self.lat)
        incl = np.arccos(sin_a + cos_b * np.cos(ha))           # sun inclination   rad
        sunhgt = max(0.05, np.pi / 2 - incl)                   # solar height   rad
        kb = 0.5/np.sin(sunhgt) * clump                          # beam radiation extinction coeff  0.5 <- clumping
        kbprime = 0.46/np.sin(sunhgt)                          # beam + scattered beam PAR ext coeff
        rhocb = 1-np.exp(-2*rhoh*kb/(1+kb))                         # canopy reflec coeff for beam PAR
        m = self.P / Po /np.sin(sunhgt)                        # optical air mass
        fd =(1-a**m)/(1+a**m * (1/fa-1))                            # fraction of diffuse irradiance
        
        ###################################################################################################
        ## Calculation of PAP fractionation
        ## For rad fractionation, not consider scattering because it will consider in gasexchange step
        lai = LAI
        It = PPFD          
        Id = It * fd                       # diffuse fraction of irradiance PAR umol m-2 s-1
        Ib = It - Id                       # beam fraction of irradiance PAR umol m-2 s-1
        ## for total leaves
        Icbs = Ib*(1-rhocb)*(1-np.exp(-kbprime*lai))       # canopy absorbed beam irrad.
        Icd = Id*(1-rhocd)*(1-np.exp(-kdprime*lai))        # canopy absorbed diffuse irrad.
        Ic = Icbs + Icd                                         # canopy absorbed total irrad.
        ## for sunlit leaves
        Icdb = Ib*(1-sigma)*(1-np.exp(-kb*lai))       # absorbed direct beam by sunlit leaves
        Icdf = Id*(1-rhocd)*(1-np.exp(-(kdprime+kb)*lai)) * \
            kdprime/(kdprime+kb)                                # absorbed diffse beam by sunlit leaves
        Icsc = Ib*((1-rhocb)*(1-np.exp(-(kbprime+kb)*lai))*kbprime/(kbprime+kb)- \
            (1-sigma)*(1-np.exp(-2*kb*lai))/2)             # absorbed scattered beam by sunlit leaves
        Icsun = Icdb + Icdf + Icsc
        ## for shaded leaves
        Icshdf = Id*(1-rhocd)*(1-np.exp(-kdprime*lai)-(1-np.exp(-(kdprime+kb)*lai))* \
            kdprime/(kdprime+kb))                              # absorbed diffuse beam by shaded leaves
        Icshsc = Ib*((1-rhocb)*(1-np.exp(-kbprime*lai)-(1-np.exp(-(kbprime+kb)*lai))* \
            kbprime/(kbprime+kb))-(1-sigma)*(1-np.exp(-kb*lai)-(1-np.exp(-2*kb*lai))/2))
                                                               # absorbed scattered beam by shaded leaves
        Icsh = Icshdf + Icshsc
        self.It = It         # PAR total (PAR umol m-2 s-1)
        self.Ib = Ib         # beam fraction of PAR
        self.Id = Id         # diffuse fractio of PAR
        self.Icsun = Icsun   # PAR of sunlit leaf
        self.Icsh  = Icsh    # PAR of shaded leaf

    ####################################################################################################
    ## LAI fractionation
    def laiFraction(self, doy, hour, LAI):
        lai = LAI
        decl = -0.4093 * np.cos(2 * np.pi * (doy + 10) / 365)  # sun declination   rad
        ha   = np.pi / 12 * (hour - 12)                        # hour angle   rad
        sin_a = np.sin(decl)*np.sin(self.lat)
        cos_b = np.cos(decl)*np.cos(self.lat)
        incl = np.arccos(sin_a + cos_b * np.cos(ha))           # sun inclination   rad
        sunhgt = max(0.05, np.pi / 2 - incl)                   # solar height   rad
        kb = 0.5/np.sin(sunhgt)  * clump                              # beam radiation extinction coeff  0.5 <- clumping
        laiSun = (1 - np.exp(-kb * lai)) / kb
        laiSh  = lai - laiSun
        self.laiSun = laiSun    # LAI for sunlit leaf
        self.laiSh  = laiSh     # lai for shaded leaf
        
    ####################################################################################################
    ## Rubisco fractionation
    def rubFraction(self, doy, hour, LAI, Vcmax=110):
        lai = LAI
        Xn = Vcmax/(Nl-Nb)
        Nc = lai * ((No-Nb)*(1-np.exp(-kn))/kn + Nb)           # leaf N profile of canopy mmol m-2
        decl = -0.4093 * np.cos(2 * np.pi * (doy + 10) / 365)  # sun declination   rad
        ha   = np.pi / 12 * (hour - 12)                        # hour angle   rad
        sin_a = np.sin(decl)*np.sin(self.lat)
        cos_b = np.cos(decl)*np.cos(self.lat)
        incl = np.arccos(sin_a + cos_b * np.cos(ha))           # sun inclination   rad
        sunhgt = max(0.05, np.pi / 2 - incl)                   # solar height   rad
        kb = 0.5/np.sin(sunhgt)  * clump                             # beam radiation extinction coeff  0.5 <- clumping

        Vcmaxtot = lai * Xn * (No-Nb) * (1-np.exp(-kn)) / kn   # Vcmax of whole canopy umol m-2 -s
        Vcsun = lai * Xn * (No-Nb) * (1-np.exp(-kn-kb*lai)) / (kn + kb*lai)
                                                               # Vcmax of sunlit leaves umol m-2 -s
        Vcsh  = Vcmaxtot - Vcsun                               # Vcmax of shaded leaves umol m-2 -s  
        
        self.VcmaxSun = Vcsun     # rubisco for sunlit leaf
        self.VcmaxSh  = Vcsh      # rubisco for shaded leaf      