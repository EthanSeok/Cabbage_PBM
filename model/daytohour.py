import numpy as np
import pandas as pd
import datetime
import random

class DayToHour():
    
    def __init__(self, fname_in, start, end, latitude, press=1000, longitude=127):
        
        self.height    = 2.0                       # reference height (m)
        self.lat       = np.radians(latitude)      # radian
        self.longi     = np.radians(longitude)     # radian
        self.press     = press
        
        # location information
        self.dap       = 0.0                      ###""" 2021 UW 추가  for calculate DAP """
        self.doy       = 0.0
        self.Tmax      = 0.0
        self.Tmin      = 0.0
        self.Tmean     = (self.Tmax + self.Tmin) / 2
        self.sunhour   = 0.0
        self.rain      = 0.0
        self.wind      = 0.0
        self.RH        = 0.0
        
        fname  = fname_in
        # discrimiate file type (.txt, .csv, .xlsx, .xls)
        # proccess for read, set index to datetime, and missing data adjustment
        ext = fname.split('.')[1]
        if ext == 'csv':
            df = pd.read_csv(fname, parse_dates=['date']).set_index('date').sort_index()
        elif ext == 'txt':
            df = pd.read_csv(fname, sep='\t', parse_dates=['date']).set_index('date').sort_index()
        elif ext == 'xlsx' or ext == 'xls':
            df = pd.read_excel(fname).set_index('date').sort_index()
        else:
            raise Exception('This file cannot read. Change file type.')
        #df = df.resample('D').interpolate()

        start = pd.to_datetime(start, format='%Y-%m-%d')
        end = pd.to_datetime(end, format='%Y-%m-%d') #+ datetime.timedelta(days=1)
        df = df[start:end]
        idx = pd.to_datetime(df.index, format='%Y-%m-%d')                     ###  2021 calculate DAP and add column ####
        df['dap'] = (idx - start).days                                   ###  2021 calculate DAP and add culumn ####
            
        # set weather data as DataFrame type
        self.dayData = df
        
        # result of hourly data
        self.data = None
        

    def calcHour(self):
        
        res = {'timestamp':[], 'Irrad':[], 'Tair':[], 'rain':[], 'wind':[], 'RH':[], 'dap':[]}     ## 2021 add dap column  ##
        
        dayData = self.dayData
        
        for i, v in dayData.iterrows():
            doy          = i.timetuple().tm_yday
            self.doy     = doy
            self.sunhour = v.sunhr
            self.Tmax    = v.Tmax
            self.Tmin    = v.Tmin
            self.rain    = v.rain
            self.wind    = v.wind
            self.RH      = v.RH
            self.dap     = v.dap                                                      ## 2021 add dap column  ##
            dayrain      = self.calcDayRain(doy, v.rain)
            self.dayrain = dayrain
            
            for hour in range(24):
                # calculation of hourly data
                timestamp = i.replace(hour = hour, minute=0)
                hourTemp  = self.calcTemp(hour)
                hourRH    = self.calcRH(hour)
                hourWind  = self.calcWind(hour)
                hourRain  = dayrain[hour]
                hourIrrad = self.calcCloud(hour)
                
                res['timestamp'].append(timestamp)
                res['Irrad'].append(hourIrrad)
                res['Tair'].append(hourTemp)
                res['rain'].append(hourRain)
                res['wind'].append(hourWind)
                res['RH'].append(hourRH)
                res['dap'].append(self.dap)                                              ## 2021 add dap column  ##
        columns = ['timestamp', 'Irrad', 'Tair', 'rain', 'wind', 'RH', 'dap']            ## 2021 add dap column  ##
        self.data = pd.DataFrame(res, columns=columns).set_index('timestamp')

    # dayLength, sunRise, sunSet (hr)
    def solarCalc(self):
        lat   = self.lat
        decl  = -0.4093 * np.cos(2 * np.pi * (self.doy + 10) /365)
        dA    = np.sin(decl) * np.sin(lat)
        dB    = np.cos(decl) * np.cos(lat)
        dayLength = 24 * np.arccos(-dA / dB) / np.pi
        sunSet    = 12 * np.arccos(-dA / dB) / np.pi + 12
        sunRise   = 24 - sunSet
        return (dayLength, sunRise, sunSet)
    
    #### hourly variables
    # hour angle(radian) : corrected by standard meridian longitude and equation of time
    def hourAngle(self, hour):
        B   = 2 * np.pi * (self.doy - 81) / 364
        EoT = 9.87 * np.sin(2 * B) - 7.53 * np.cos(B) - 1.5 * np.sin(B)
        gammaSM = int(self.longi / (np.pi/12)) * (np.pi/12)
        corr    = (gammaSM - self.longi) / (np.pi/12)
        tcorr   = hour + corr + EoT / 60
        hourAngle = np.pi * (tcorr - 12) /12
        return hourAngle
        
    # hourly temperature from Tmax, Tmin
    def calcTemp(self, hour):
        Tmax = self.Tmax
        Tmin = self.Tmin 
        dayLength, sunRise, sunSet = self.solarCalc()
        
        # Reicosky(1989), de Wit (1978) method
        if hour < sunRise and hour >= 0:
            dTau  = np.cos(np.pi * (hour + 10) / (sunRise + 10))
            dTemp = (Tmax + Tmin) / 2 + (Tmax - Tmin) / 2 * dTau
        elif hour >= sunRise and hour <= 14:
            dTau  = np.cos(np.pi * (hour - sunRise) / (14 - sunRise))
            dTemp = (Tmax + Tmin) / 2 - (Tmax - Tmin) / 2 * dTau
        else :
            dTau  = np.cos(np.pi * (hour - 14) / (10 + sunRise))
            dTemp = (Tmax + Tmin) / 2 + (Tmax - Tmin) / 2 * dTau  
        return dTemp


    ### calculate hourly RH (modified Teh method)
    def calcRH(self, hour):
        rh    = self.RH
        rhcor = min(0.0148 * rh - 0.1669, 1)
        temp = self.calcTemp(hour)
        svp = 6.1078 * np.exp(17.269 * temp / (temp + 237.3))
        svpTmin = 6.1078 * np.exp(17.269 * self.Tmin / (self.Tmin + 237.3))
        hourRH  = svpTmin / svp * 100 * rhcor
        return hourRH
    
    ### calculate hourly wind from mean value
    def calcWind(self, hour):
        wind  = max(self.wind, 0.01)
        stdRatio = -0.175 * np.log(wind) + 0.6151 
        std   = wind * stdRatio / 2
        hourWind = wind + np.random.normal(0.0, std, 1)[0]
        hourWind = max(hourWind, 0.01)
        return hourWind    

    ### calculate hourly rain from dayrain list
    ##  average hours by doy : 6 hrs(1-59), 7 (60-90), 13 (91-120), ...
    ##  avhrs = int(-0.0004 * doy ** 2 + 0.149 * doy - 0.5987)
    def calcDayRain(self, doy, rain):
        rain  = self.rain
        doy   = self.doy
        avhrs = int(-0.0004 * doy ** 2 + 0.149 * doy - 0.5987)
        avhrs = max(avhrs, 3)
        dayrain = [0.0] *24
        time = random.sample(range(0, 23), avhrs)
        freq = [random.gammavariate(1, 1) for i in range(avhrs)]
        itsy = [rain * v/sum(freq) for v in freq]
        for i in range(avhrs):
            dayrain[time[i]] = itsy[i]
        return dayrain
        
    # calculate hourly radiation (using MRM model)
    # seasonal extra solar radiation W m-2, corrected by dist. from E and S(eccen), EoT and standard meridian
    # optical air mass
    def calcIexOam(self, hour): 
        lat   = self.lat
        press = self.press
        decl  = -0.4093 * np.cos(2 * np.pi * (self.doy + 10) /365)
        gamma = 2.0 * np.pi * (self.doy - 1) / 365
        Iex   = 1366.1 * (1.00011 + 0.034221 * np.cos(gamma) + 0.00128 * np.sin(gamma) 
                        + 0.000719 * np.cos(2 * gamma) + 0.000077 * np.sin(2 * gamma))
        tau   = self.hourAngle(hour)
        A     = np.sin(decl) * np.sin(lat)
        B     = np.cos(decl) * np.cos(lat) * np.cos(tau)
        sinbeta = max(A + B, 0)
        extraSolar = Iex * sinbeta        # W/m2
        
        theta   = np.pi / 2 - np.arcsin(sinbeta)      
        theta   = min(np.pi, max(theta, 0))       # solar zenith angle(radians)
        sza     = np.degrees(theta)               # solar zenith angle(degrees)
        C       = np.power((96.07995 - sza), -1.6364)
        optAirMass  = 1 / (sinbeta + 0.50572 * C)
        corrAirMass = optAirMass * press / 1013.25     # Po = 1013.25 mbar or 101.325 kPa
        return (extraSolar, optAirMass, corrAirMass)
    
    ### Calculation radiation inhibitors
    #  Tw (H2O vapor), Tmg(mixed gas), Toz(ozone), Tray(Rayleigh scatter), Taero (aerosol)
    def calcT(self, hour):
        lat    = self.lat 
        height = self.height
        temp   = self.Tmean
        rh     = self.RH

        Tl     = (temp + 273.15) / 100
        es     = np.exp(22.329699 - 49.140396 / Tl - 10.921853 / Tl / Tl - 0.39015156 * Tl)
        em     = es * rh / 100
        uwater = 0.493 * em / (temp + 273.15)
        Iex,m,mpr = self.calcIexOam(hour)       # Iex (extra solar), m (optical air mass)
        B      = max(1 + 119.3 * uwater * m, 0.1)
        A      = np.power(B, 0.644)
        Twater = 1 - 3.014 * m * uwater / (A + 5.814 * m * uwater)    # water coefficient
        
        Tco2 = 1 - 0.721 * mpr * 350 / ((1 + 377.89 * mpr * 350) ** (0.5855) + 3.1709 * mpr * 350)
        Tco = 1 - 0.0062 * mpr * 0.075 / ((1 + 243.67 * mpr * 0.075) ** (0.4246) + 1.7222 * mpr * 0.075)
        Tn2o = 1 - 0.0326 * mpr * 0.28 / ((1 + 107.413 * mpr * 0.28) ** (0.5501) + 0.9093 * mpr * 0.28)
        Tch4 = 1 - 0.0192 * mpr * 1.6 / ((1 + 166.095 * mpr * 1.6) ** (0.4221) + 0.7186 * mpr * 1.6)
        To2 = 1 - 0.0003 * mpr * 209500 / ((1 + 476.934 * mpr * 209500) ** (0.4892) + 0.1261 * mpr * 209500)       
        Tmix = Tco2 * Tco * Tn2o * Tch4 * To2                          # mixed gases coeff.
        
        Tozone = 1 - 0.2554 * m * 0.3 / ((1 + 6017.26 * m * 0.3) ** (0.204) + 0.471 * m * 0.3)  # from Psiloglou(2007)
        
        AA = -0.1128 * mpr ** (0.8346)
        BB = 0.9341 - mpr ** (0.9868) + 0.9391 * mpr
        Tray  = np.exp(AA * BB)
        
        beta = (0.025 + 0.1 * np.cos(lat)) * np.exp(-0.7 * height / 1000) + 0.04      #from Psiloglou(2007)
        Taero = np.exp(-1 * m * beta * (0.6777 + 0.1464 * m * beta - 0.00626 * (m * beta) ** 2) ** (-1.3))
        Taeab = 1 - 0.1 * (1 - m + m ** (1.06)) * (1 - Taero)
        Tratio = Taero / Taeab     
        return (Twater, Tmix, Tozone, Tray, Taero, Taeab, Tratio)
    
    ### calculate radiations on Clear condition (W/m2) 
    #  Itotal, Ibeam, Idiffuse
    def calcClear(self, hour):
        height = self.height
        lat    = self.lat

        Iex, m, mpr = self.calcIexOam(hour)
        Tw, Tmg, Toz, Tr, Ta, Taeab, Tratio = self.calcT(hour) 
        Ibeam = Iex * Tw * Tr * Toz * Tmg * Ta
        Ibeam = max (Ibeam, 0.0)
        
        # Idfsin, clear diffuse by single scattering of direct beam
        dayLength, sunRise, sunSet = self.solarCalc()
        Idfsin2  = Iex * Tw * Toz * Tmg * Taeab * 0.5 * (1 - Tratio * Tr)
        if hour < sunRise:
            Idfsin2 = 0.0
        elif hour > sunSet :
            Idfsin2 = 0.0
        Idfsin = Idfsin2            # diffuse by single scattering of direct beam radiation
    
        # Idfmul, clear diffuse by multiple scattering of direct beam from Psiloglou(2007)
        beta = (0.025 + 0.1 * np.cos(lat)) * np.exp(-0.7 * height / 1000) + 0.04 
        Ta166 = np.exp(-1 * 1.66 * beta * (0.6777 + 0.1464 * 1.66 * beta - 0.00626 * (1.66 * beta)**2)**(-1.3))
        alphag = 0.2
        alphar = 0.0685
        alphaa = 0.16 * (1 - Ta166)
        Idfmul2 = (Ibeam + Idfsin) * (alphag * (alphar + alphaa)) / (1 - alphag * (alphar + alphaa))
        if hour < sunRise:
            Idfmul2 = 0.0
        elif hour > sunSet :
            Idfmul2 = 0.0
        Idfmul = Idfmul2            # diffuse by multiple scattering of direct beam radiation
        Idiffuse = Idfsin + Idfmul   # diffuse total at clear condition
        Itotal = Ibeam + Idiffuse
        return (Itotal, Ibeam, Idiffuse, Idfsin, Idfmul)
    
    
    ### calculati radiations on Cloud condition (W/m2)
    #  Ictotal, Icbeam, Icdiffuse
    def calcCloud(self, hour):
        height = self.height
        lat     = self.lat
        sunhour = self.sunhour
        dayLength, sunRise, sunSet = self.solarCalc()
        
        Tcloud = 0.75 * sunhour / dayLength     # 0.9 is empirical coeff Psiloglou(2007), 0.75 from Ideriah(1981)
        Itotal, Ibeam, Idiffuse, Idfsin, Idfmul = self.calcClear(hour)
        
        Icbeam  = Ibeam * Tcloud                                               # beam radiation at cloud condition
        Icdfsin = Idfsin * Tcloud + 0.33 * (1 - Tcloud) * (Ibeam + Idfsin)     # diffuse single scattering
        
        beta = (0.025 + 0.1 * np.cos(lat)) * np.exp(-0.7 * height / 1000) + 0.04      # Psiloglou(2007)
        Ta166 = np.exp(-1 * 1.66 * beta * (0.6777 + 0.1464 * 1.66 * beta - 0.00626 * (1.66 * beta)**2)**(-1.3))
        alphag = 0.2
        alphar = 0.0685
        alphaa = 0.16 * (1 - Ta166)
        alphac = 0.4 * sunhour / dayLength   # from Psiloglou(2007)
        
        Icdfmul = (Icbeam + Icdfsin)*alphag*(alphar + alphaa + alphac)/(1-alphag*(alphar+alphaa+alphac))
        
        
        Icdiffuse = Icdfsin + Icdfmul
        Ictotal   = Icbeam + Icdiffuse
        
        return Ictotal


      