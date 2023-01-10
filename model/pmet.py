import math

class ETo:
    
    def __init__(self):
        self.lat   = 0.0
        self.press = 0.0       # kPa
        self.doy   = 0.0
        self.sunhr = 0.0
        self.Tmax  = 0.0
        self.Tmin  = 0.0
        self.Tmean = 0.0
        self.wind  = 0.0
        self.RH    = 0.0
        
        # resulted ET
        self.pmet   = 0.0
        self.energy = 0.0        
        
    def setValue(self, latitude, doy, sunhr, Tmax, Tmin, wind, RH, press=100):
        self.lat   = math.radians(latitude)
        self.doy   = doy
        self.sunhr = sunhr
        self.Tmax  = Tmax
        self.Tmin  = Tmin
        self.Tmean = (Tmax + Tmin) / 2
        self.wind  = wind
        self.RH    = RH
        self.press = press
        

        
    def calc(self):
        press   = self.press       # kPa
        lat     = self.lat
        doy     = self.doy
        sunhr   = self.sunhr
        wind    = self.wind
        Tmax    = self.Tmax
        Tmin    = self.Tmin
        Tmean   = self.Tmean
        RH      = self.RH
        wind_c  = wind * 4.87 / math.log(67.8 * 10 - 5.42)
        slope   = 4098*(0.6108*math.exp((17.27*Tmean)/(Tmean+237.3)))/math.pow((Tmean+237.3),2)                
        
        gamma   = 0.665 * press / 1000
        
        A       = slope / (slope + gamma * (1 + 0.34 * wind_c))  
        B       = gamma / (slope + gamma * (1 + 0.34 * wind_c))
        C       = (900 / (Tmean + 273)) * wind_c
        
        avpTmax = 0.6108 * math.exp(17.27 * Tmax / (Tmax + 237.3))
        avpTmin = 0.6108 * math.exp(17.27 * Tmin / (Tmin + 237.3))
        es      = (avpTmax + avpTmin) / 2
        ea      = es * RH / 100
        delsa   = es - ea        
        
        dr      = 1 + 0.033 * math.cos(2 * math.pi / 365 * doy)
        delta   = 0.409 * math.sin((2 * math.pi / 365 * doy) - 1.39)
        ws      = math.acos(-1 * math.tan(lat) * math.tan(delta))
        v_sin   = math.sin(lat) * math.sin(delta)
        v_cos   = math.cos(lat) * math.cos(delta)

        
        Ra      = 24 * 60 / math.pi * 0.082 * dr * (ws * v_sin + v_cos)
        N       = 24 / math.pi * ws
        Rs      = (0.25 + 0.5 * sunhr / N) * Ra
        Rso     = (0.75 + 2 * 20 / 100000) * Ra
        Rns     = 0.77 * Rs

        
        B_Tmax  = 0.000000004903 * math.pow((Tmax + 273.16),4)
        B_Tmin  = 0.000000004903 * math.pow((Tmin + 273.16),4)
        B_mean  = (B_Tmax + B_Tmin) / 2
        Rnl     = B_mean * (0.34 - 0.14 * math.sqrt(ea)) * (1.35 * (Rs/Rso) - 0.35)
        
        Rn      = Rns - Rnl
      
        G       = 0.0
        
        energy  = 0.408 * A * (Rn - G)
        mass    = B * C * delsa
        eto     = energy + mass

        self.energy = energy
        self.pmet = eto
        