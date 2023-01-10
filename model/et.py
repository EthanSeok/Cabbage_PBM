"""
ET class
 - calculate hourly et, ets, etc, h, hs, hc (all W/m2), etmm, etsmm, etcmm (mm/hr/m2 groud)
 - input : latitude, doy, hour, Irrad(W/m2), temp, wind(m/s), RH, lai, vwc01(m3/m3), hgt(m), width(m)
 - output : et, ets, etc, etmm, etsmm, etcmm, h, hs, hc
"""

import numpy as np

# parameters
keddy      = 2.0
kwind      = 2.0
rstA1      = 11.255
rstA2      = 0.0051
refhgt     = 2.0              # reference height (2 m)
porosity   = 0.42
poredist   = 0.18
vwcsat     = 0.40

class ET():
    
    def __init__(self): 
        self.lat       = 0.0
        self.doy       = 0.0
        self.hour      = 0.0
        self.Irrad     = 0.0
        self.temp      = 0.0
        self.wind      = 0.0
        self.RH        = 0.0
        self.lai       = 0.0
        self.vwc01     = 0.0
        self.height    = 0.0
        self.leafwidth = 0.0
        
        # calculation results
        self.et        = 0.0
        self.ets       = 0.0
        self.etc       = 0.0
        self.etmm      = 0.0  # mm/hr/m2 ground
        self.etsmm     = 0.0  # mm/hr/m2 ground
        self.etcmm     = 0.0  # mm/hr/m2 ground 
        self.h         = 0.0
        self.hc        = 0.0
        self.hs        = 0.0
        self.Tleaf     = 0.0

        
    ## set input value
    def setValue(self, latitude, doy, hour, Irrad, temp, wind, RH, lai, vwc01=0.15, hgt=0.6, width=0.15):
        self.lat       = np.radians(latitude)
        self.doy       = doy
        self.hour      = hour
        self.Irrad     = Irrad
        self.temp      = temp
        self.wind      = wind
        self.RH        = RH
        self.lai       = max(lai, 0.001)
        self.vwc01     = vwc01           # volumetric water content of surface soil layer (0-1)
        self.height    = hgt             # plant height (m)
        self.leafwidth = width           # leaf width (m)
        
            
    # solar elevation (radians)
    def elev(self):
        dDecl = -0.4093 * np.cos(2 * np.pi * (self.doy + 10) / 365)
        dA    = np.sin(dDecl) * np.sin(self.lat)
        dB    = np.cos(dDecl) * np.cos(self.lat)
        dHa   = np.pi * (self.hour - 12.0) / 12.0
        elevation = np.arcsin(dA + dB * np.cos(dHa))
        return elevation
    
    # kdr calculation
    def kdr(self):
        dElev = self.elev()
        kdr = max(0.5 / np.sin(dElev), 0.0001)
        return kdr
    
    # slope of saturated vapor pressure with temperature(mbar/K)
    def slope(self):
        dA = 25029.4 * np.exp(17.269 * self.temp / (self.temp + 237.3))
        dB = (self.temp + 237.3) ** 2
        return (dA / dB)
    
    ## Resistances
    # Friction velocity by height
    def FrictionVelocity(self):
        height = self.height
        wind = self.wind
        dZ = 0.64 * height
        dR = 0.13 * height
        fricv = 0.4 * wind / np.log((height - dZ) / dR)
        return fricv
    
    ### Calculation aerodynamic resistance (s/m)
    #   rsa and raa calculation
    #   rsa (aerodynamic resistance between mean canopy flow and soil)
    #   raa (aerodynamic resistance between mean canopy flow and reference height) 
    def AeroResistance(self):
        height = self.height
        dZ = 0.64 * height
        dR = 0.13 * height
        dFv = self.FrictionVelocity()
        dKh = 0.4 * dFv * height
        dP  = np.exp(-1* keddy * 0.004 / height) - np.exp (-1*keddy * (dZ + dR) / height)
        rsa = dP * height * np.exp(keddy) / (keddy * dKh)
        dP  = np.exp(keddy * (1.0 - (dZ + dR) / height)) - 1.0
        dP  = dP * (height / (keddy * dKh))
        raa = np.log((refhgt - dZ)/(height - dZ)) / (0.4 * dFv) + dP
        return (rsa, raa)
    
    ## Calculatin boundary layer resistance (s/m)
    # rca calculation
    #    rca (aerodynamic resistance between mean canopy flow and canopy)
    def BoundLayerResistance(self):
        height    = self.height
        leafwidth = self.leafwidth
        lai    = self.lai
        dZ = 0.64 * height
        dR = 0.13 * height
        dFv  = self.FrictionVelocity()
        dUh  = (dFv / 0.4) * np.log((height - dZ) / dR)
        dRca = 0.012 * lai * (1.0 - np.exp(-1 * kwind / 2.0))
        rca  = kwind / (dRca * np.sqrt(dUh / leafwidth))
        return rca   
    
    ### Calculation canopy resistance (s/m)
    # rcs calculation
    #    rcs (canopy resistance)
    def CanopyRestance(self):
        lai    = self.lai
        totrad = self.Irrad
        totrad = max(totrad, 0.01)             # avoid 0 dividing
        dPar = 0.5 * totrad   
        dRcs = (rstA1 + dPar) / (rstA2 * dPar) 
        LAICR = 0.5 * 4.0                      # critical lai
        if lai > LAICR:
            dLai = LAICR
        else:
            dLai = lai
        rcs = dRcs / dLai
        return rcs    
 
    ### Calculation soil resistance (s/m)
    # rss calculation
    #    rss (soil resistance)
    def SoilResistance(self):
        vwc01   = self.vwc01
        dRssdry = (2.0 * 0.02) / (porosity * 0.0000247)
        dExp    = np.exp(-1*(1.0 / poredist) * (vwc01 / vwcsat))
        rss     = dRssdry * dExp
        return rss
 
    ### Net energy available soil and canopy (W/m2)
    def EnergySupply(self):
        acanopy = 0
        asoil   = 0
        lai    = self.lai
        dRn    = self.Irrad 
        kdr    = self.kdr()
        if dRn > 0:
            scateff = np.sqrt(0.5)
            dRsn = dRn * np.exp(-1*kdr * scateff * lai)
            acanopy = dRn - dRsn
            asoil   = dRn  -  self.G(dRsn) - acanopy
        return (acanopy, asoil)
    
    # Soil heat flux (W/m2)
    #   th (local solar time), rsn (solar radiation reaching the soil)
    def G(self, rsn):
        dSolInc = np.pi / 2.0 - self.elev()
        Grd     = 0.35 * np.cos(dSolInc) * rsn
        return Grd

    
    # vapor pressure deficit(vpd) at reference height (mbar)
    def vpdRef(self):
        satvp = 6.1078 * np.exp(17.269 * self.temp / (self.temp + 237.3))
        hourvp = satvp * self.RH / 100.0
        return (satvp - hourvp)
    
    # vapor pressure deficit(vpd) at mean canopy flow(mbar)
    #   a (total energy available), et (total ET), temp (temperature)
    #   raa (aerodynamic resistance between mean canopy flow and reference height)
    def vpdMcf(self, a, et, raa):
        dDelta = self.slope()
        dA     = self.vpdRef()
        dB     = (raa / 1221.09) * (dDelta * a - (dDelta + 0.658) * et)
        return (dA + dB)
    
    
    # Sensible heaf for soil (W/m2)
    #   rsa (aerodynamic resistance between mean canopy flow and soil)
    #   rss (soil resistance), asoil (energy available to soil)
    def Hs(self, raa, rsa, rss, a, asoil, et):
        dA = 0.658 * asoil * (rss + rsa) - 1221.09 * self.vpdMcf(a, et, raa)
        dB = self.slope() * rsa + 0.658 * (rss + rsa)
        self.y = dA / dB
        return (dA / dB)
    
    # Sensible heat for crop (W/m2)
    #   rca (aerodynamic resistance between mean canopy flow and crop)
    #   rcs (canopy resistance), acanopy (energy available to canopy)
    def Hc(self, raa, rca, rcs, a, acanopy, et):
        dA = 0.658 * acanopy * (rcs + rca) - 1221.09 * self.vpdMcf(a, et, raa)
        dB = self.slope() * rca + 0.658 * (rcs + rca)
        self.x = dA / dB
        return (dA / dB)
    
    # Latent heat for soil (W/m2)
    #   rsa (aerodynamic resistance between mean canopy flow and soil)
    #   rss (soil resistance), asoil (energy available to soil)
    def ETs(self, raa, rsa, rss, a, asoil, et):
        dA = self.slope() * asoil + 1221.09 * self.vpdMcf(a, et, raa) / rsa
        dB = self.slope() + 0.658 * (rss + rsa) / rsa
        return (dA / dB)
    
    # Latent heat for crop (W/m2)
    #   rca (aerodynamic resistance between mean canopy flow and canopy)
    #   rcs (canopy resistance), asoil (energy available to soil)
    def ETc(self, raa, rca, rcs, a, acanopy, et):
        dA = self.slope() * acanopy + 1221.09 * self.vpdMcf(a, et, raa) / rca
        dB = self.slope() + 0.658 * (rcs + rca) / rca
        return (dA / dB)
    
    ### Canopy temperature (C)
    #     h (sensible heaf total), hc (sensible heat canopy)
    def calcCanopyTemp(self, h, hc):
        dRsa, dRaa = self.AeroResistance()
        dT0  = ((h / 1221.09) * dRaa) + self.temp
        dRca = self.BoundLayerResistance()
        dTf  = ((hc / 1221.09) * dRca) + dT0
        self.Tleaf = dTf
    
    ### Latent heat and Sensible heat fluxes (W/m2) at time th
    #      et, ets, etc (latent heat for total, soil and canopy)
    #      h, hs, hc (sensible heat for total, soil and canopy) 
    def calcFlux(self):
        # hourly meteorology
        dWind = self.wind
        dTemp = self.temp
        dTotrad = self.Irrad
        
        # resistances
        dRsa, dRaa = self.AeroResistance()
        dRca = self.BoundLayerResistance()
        dRcs = self.CanopyRestance()
        dRss = self.SoilResistance()
        
        # energy availables to soil and plant
        dAc, dAs = self.EnergySupply()
        
        # fluxes
        dDelta = self.slope()
        dRa    = (dDelta + 0.658) * dRaa
        dRc    = (dDelta + 0.658) * dRca + 0.658 * dRcs
        dRs    = (dDelta + 0.658) * dRsa + 0.658 * dRss
        dCc    = 1.0 + (dRc * dRa /(dRs * (dRc + dRa)))
        dCc    = 1.0 / dCc
        dCs    = 1.0 + (dRs * dRa /(dRc * (dRs + dRa)))
        dCs    = 1.0 / dCs
        dA     = dAs + dAc
        dD     = self.vpdRef()
        
        dP     = dDelta * dA + (1221.09*dD - dDelta*dRca*dAs) / (dRaa + dRca)
        dQ     = dDelta + 0.658 * (1.0 + dRcs / (dRaa + dRca))
        dPMc   = dP / dQ
        
        dP     = dDelta * dA + (1221.09*dD - dDelta*dRsa*dAc) / (dRaa + dRsa)
        dQ     = dDelta + 0.658 * (1.0 + dRss / (dRaa + dRsa))
        dPMs   = dP / dQ
        
        et     = dCc * dPMc + dCs * dPMs  
        ets    = self.ETs(dRaa, dRsa, dRss, dA, dAs, et)  #raa, rsa, rss, a, asoil, et  
        etc    = et - ets
        #etc    = dCc * dPMc
        #ets    = et - etc
        hs     = self.Hs(dRaa, dRsa, dRss, dA, dAs, et)   #raa, rsa, rss, a, asoil, et
        hc     = self.Hc(dRaa, dRca, dRcs, dA, dAc, et) 
        h      = hs + hc
        
        self.calcCanopyTemp(h, hc)
        self.et  = et
        self.ets = ets
        self.etc = etc
        self.etmm   = et  / 2454000 * 3600      # convert W/m2 gr-> mm/hr/m2 gr  (W=J/s)
        self.etsmm  = ets / 2454000 * 3600     # convert J/s/m2 gr-> mm/hr/m2 gr  
        self.etcmm  = etc / 2454000 * 3600     # convert J/s/m2 gr-> mm/hr/m2 gr
        self.h   = h
        self.hc  = hc
        self.hs  = hs 
        
        # return (et, ets, etc, h, hs, hc)
