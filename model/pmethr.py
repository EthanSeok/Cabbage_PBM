import math

class ETo:
    
    def __init__(self):
        self.lat   = 0.0
        self.longi = 0.0
        self.hour  = 0.0
        self.press = 0.0       # kPa
        self.doy   = 0.0
        self.Irrad = 0.0
        self.Tair  = 0.0
        self.wind  = 0.0
        self.RH    = 0.0
        self.lai   = 0.0
        
        # resulted ET
        self.pmet  = 0.0
        self.pmetc = 0.0
        self.pmets = 0.0 
        
    def setValue(self, latitude, doy, hour, Irrad, Tair, wind, RH, lai=2.0, press=100, longitude=127):
        self.lat   = math.radians(latitude)
        self.longi = math.radians(longitude)
        self.hour  = hour
        self.doy   = doy
        self.Irrad = Irrad
        self.Tair  = Tair
        self.wind  = wind
        self.RH    = RH
        self.press = press
        self.lai   = lai

    def calc(self):
        press   = self.press       # kPa
        lat     = self.lat
        longitude  = math.degrees(self.longi)
        doy     = self.doy
        hour    = self.hour
        Irrad   = self.Irrad
        wind    = self.wind
        Tair    = self.Tair
        RH      = self.RH
        lai     = self.lai
        
        #wind_c  = wind * 4.87 / math.log(67.8 * 10 - 5.42)
        slope   = 4098*(0.6108*math.exp((17.27*Tair)/(Tair+237.3)))/math.pow((Tair+237.3),2)                
        
        gamma   = 0.665 * press / 1000
        
        esTair = 0.6108 * math.exp(17.27 * Tair / (Tair + 237.3))
        ea      = esTair * RH / 100
        delsa   = esTair - ea
        
        dr      = 1 + 0.033 * math.cos(2 * math.pi / 365 * doy)
        delta   = 0.409 * math.sin((2 * math.pi / 365 * doy) - 1.39)
        delta   = max(delta, 0.0001)
        ws      = math.acos(-1 * math.tan(lat) * math.tan(delta))
        b       = 2 * math.pi * (doy - 81) / 364
        sc      = 0.1645 * math.sin(2*b) - 0.1255 * math.cos(b) - 0.025 * math.sin(b)
        
        w       = math.pi/12 *((hour + 0.5 + 0.06667 * (127-longitude) + sc) - 12)
        w1      = w - math.pi / 24
        w2      = w + math.pi / 24
        w_sin   = math.sin(w2) - math.sin(w1)
        
        v_sin   = math.sin(lat) * math.sin(delta)
        v_cos   = math.cos(lat) * math.cos(delta)
        
        Rs      = Irrad * (60*60 / 1000000)
        # Ra  MJ/hr
        Ra      = 60 * 0.082 * 12/math.pi * dr * ((w2 - w1) * v_sin + v_cos * w_sin)
        st_elev = 10
        # Rso, Rns  MJ/hr/m2,    Irrad unit W/m2 = J/s/m2
        Rso     = (0.75 + 2 * st_elev / 100000) * Ra   # station elevation ~ 10 m
        Rns     = (1 - 0.23) * Rs

        # unit correction
        B_Tair  = 0.000000004903 * math.pow((Tair + 273.16),4) / 24
        Rnl     = B_Tair * (0.34 - 0.14 * math.sqrt(ea)) * (1.35 * (Rs/Rso) - 0.35)

        # net radiation
        Rn      = Rns - Rnl
        
        # G, soil heat flux
        dayLength, sunRise, sunSet = self.solarCalc()
        if hour > sunRise and hour < sunSet:
            G   = 0.1 * Rn
        else:
            G   = 0.5 * Rn
        
        ## et calculation
        A       = slope / (slope + gamma * (1 + 0.34 * wind))  
        B       = gamma / (slope + gamma * (1 + 0.34 * wind))
        C       = (37 / (Tair + 273)) * wind

        energy  = 0.408 * A * (Rn - G)
        mass    = B * C * delsa
        eto     = energy + mass

        
        # fractionation to etc and ets (transpiration and evaporation) by simple method
        # LAIcr = 2.52   Fr = lai / LAIcr
        LAIcr   = 2.50
        Fr      = lai / LAIcr
        if lai > LAIcr:
            fr  = 0.98
        else:
            fr  = lai / LAIcr
            
        ets     = eto * (1.0 - fr)
        etc     = eto - ets
        
        self.pmet   = eto
        self.pmetc  = etc
        self.pmets  = ets
        
    # dayLength, sunRise, sunSet (hr)
    def solarCalc(self):
        lat   = self.lat
        decl  = -0.4093 * math.cos(2 * math.pi * (self.doy + 10) /365)
        dA    = math.sin(decl) * math.sin(lat)
        dB    = math.cos(decl) * math.cos(lat)
        dayLength = 24 * math.acos(-dA / dB) / math.pi
        sunSet    = 12 * math.acos(-dA / dB) / math.pi + 12
        sunRise   = 24 - sunSet
        return (dayLength, sunRise, sunSet)
        
        