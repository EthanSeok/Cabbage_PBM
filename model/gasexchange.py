###################################################################################################
###   C3 Leaf photosynthesis module
###        input  : Irrad(total irradiance, W m-2), Ta (air temperature, C), RH (relative humidity) 
###                 wind(wind speed, m s-), CO2 (ambient CO2 concentration, ppm)
###        output : An (net assimilation umol CO2 m-2 s-1), Wc(rubisco-limited), Wj(light-limited)
###                 Ws (TPU-limited), Rd(respiration), gs(stomatal conductance H2O, mol m-2 s-1)
###                 Tleaf (leaf temp.), Ci(intercellular CO2, ubar), E(etranspiration, mol m-2 s-1)
###################################################################################################

import numpy as np

##### Constants
KELVIN = 273.15  # Kelvin temperature
R = 8.314  # ideal gas constant
SBC = 5.6697e-8  # Stefan-Boltzmann constant   W m-2 K-4
PSC = 6.66e-4  # psychrometer constant
Cp = 29.3  # specific heat of air    J mol-1 C-1
LAMBDA = 44000.  # Latent heat of vaporization at 25 C J mol-1
##### Parameters for light
scatt = 0.15  # leaf reflectance + transmittance
f = 0.15  # correction factor
conv = 4.57  # conversion from W m-2 to umol m-2 s-1
epsilon = 0.97  # leaf thermal emissivity
##### Other Common parameters
O2 = 210.  # O2 partial pressure (mbar=20.5kPa)
Kc25 = 404.9  # MM constant of Rubisco for CO2 at 25 from Bernacchi et al. (2001)
Ko25 = 278.4  # MM constant of Rubisco for O2 at 25 from Bernacchi et al. (2001)
Eac = 79430.  # Activation energy for Kc from Bernacchi et al. (2001)
Eao = 36380.  # Activation energy for Ko from Bernacchi et al. (2001)
g0 = 0.036  # residual stomatal conductance,  mol m-2 s-1
g1 = 10.0  # empirical coefficient of BWB model  이거
P = 1.013  # conversion factor from ppmv to ubar (= 1013000 / 1000000) = 101.3 kPa = 1013 mbar

## parameter for water stress
s_f = 0.5  # sensitivity
psi_f = -1.0
psi_th = -0.8  # threshold wp below which stress effect shows up

##########################################################
##### Parameters for Chinese cabbage from Experiment 2015
##########################################################
theta = 0.7  # hyperbola parameter 이거
width = 0.1  # leaf width (m) 이거
Vcm25 = 152.5  # Vcmax of CC at 25 C 이거
Jm25 = 238.6  # Jmax  of CC at 25 C 이거
TPU25 = 17.1  # 이거F
Rd25 = 1.7  # 이거
gamma25 = 42.8  #
Havcm = 65330  # Ha
# Hdvcm   = 149252.    # Hd
# Svvcm   = 486.       # Sv
Hajm = 57500  # Ha
Hdjm = 439800.  # Hd
Svjm = 1400  # Sv
Ear = 46390.  # Ha
Eagamma = 37830.  # Ha
EaTPU = 47100.  # Ha for TPU   from Kim and Leith(2003)


class GasExchange():

    def __init__(self):
        # variables
        self.Ic = 0.0
        self.Ta = 0.0
        self.RH = 0.0
        self.wind = 0.0
        self.Ca = 0.0

        ## Parameters
        self.Vcmax = 0.0
        self.Jmax = 0.0  # RuBP
        self.GammaStar = 0.0
        self.TPU = 0.0
        self.Kc = 0.0
        self.Ko = 0.0
        self.Rd = 0.0

        ## Output : assimilates
        self.An = 0.0
        self.Wc = 0.0
        self.Wj = 0.0
        self.Ws = 0.0
        self.Rd = 0.0
        ## Output : gs, Ds, E
        self.gsco = 0.0
        self.predawn = 0.0
        self.Ds = 0.0
        self.Emm = 0.0

    def setValue(self, PPFD, Ta, RH, wind, predawnPot, CO2=400):
        ## Input variables
        self.Ic = PPFD  # umol m-2 s-1, PAR (not irradiance, W m-2)
        self.Ta = Ta  # C
        self.RH = RH  # %
        self.wind = wind  # m s-1
        self.Ca = CO2 * P  # Ca(ubar),CO2 (ppm)
        self.predawn = predawnPot

    def routine(self):
        Ca = self.Ca
        Tleaf = self.Ta
        Ci = self.Ca * 0.9

        def func1(Ci, Tleaf):
            An = self.leafAssim(Tleaf, Ci)[0]
            gbc = self.gbc()
            gsco = self.gsc(An, Tleaf)[0]
            newCi = self.interCi(An, Tleaf)
            return (newCi - Ci)

        def func2(Ci, Tleaf):
            An = self.leafAssim(Tleaf, Ci)[0]
            gbc = self.gbc()
            gsco = self.gsc(An, Tleaf)[0]
            newTleaf = self.EnergyBal(An, Tleaf)[0]
            return (newTleaf - Tleaf)

        Ci, Tleaf = self.Newton_2Var(func1, func2, Ci, Tleaf)

        An, Wc, Wj, Ws, Rd = self.leafAssim(Tleaf, Ci)
        gbc = self.gbc()
        gsco, Ds = self.gsc(An, Tleaf)
        Tlnew, Emm = self.EnergyBal(An, Tleaf)
        # print("An", An, (Ca - Ci) / (1/gbc + 1/gsco) / P)
        if not np.isclose(An, (Ca - Ci) / (1/gbc + 1/gsco) / P):
            print("An", An, (Ca - Ci) / (1 / gbc + 1 / gsco) / P)
        # An = (Ca - Ci) / (1 / gbc + 1 / gsco) / P

        self.An = An
        self.Wc = Wc
        self.Wj = Wj
        self.Ws = Ws
        self.Rd = Rd
        self.gsco = gsco
        self.Emm = Emm

    #         return Tleaf, Ci, An, Wc, Wj, Ws, Rd, gs, newTleaf, E, G

    def leafAssim(self, Tleaf, Ci):
        self.calcParam()
        Vcmax = self.Vcmax
        Jmax = self.Jmax
        TPU = self.TPU
        GammaStar = self.GammaStar
        Kc = self.Kc
        Ko = self.Ko
        Rd = self.Rd

        J2 = self.Ic * (1 - scatt) * (1 - f) / 2
        J = ((J2 + Jmax) - np.sqrt((J2 + Jmax) ** 2 - 4 * J2 * Jmax * theta)) / (2 * theta)

        ## 광합성 모델식
        Wc = Vcmax * (Ci - GammaStar) / (Ci + Kc * (1 + O2 / Ko))  # Rubisco-limited
        Wj = J * (Ci - GammaStar) / (4 * (Ci + 2 * GammaStar))
        Ws = 3 * TPU
        Wp = self.minh(Wc, Wj, theta)
        W = min(Wp, Ws)  # min(Wp, Ws)
        Rd = 0.015 * Vcmax  # Dark respiration from Collatz(1991)
        An = W - Rd
        self.Ancheck = An
        # print(self.An)
        return An, Wc, Wj, Ws, Rd

    def minh(self, A, B, theta):  # hyperbolic minimum
        x = (A + B) * (A + B) - 4 * A * B * theta
        if (x < 0):
            res = min(A, B)
        else:
            res = ((A + B) - np.sqrt(x)) / (2 * theta)
        return res

    def EnergyBal(self, An, Tleaf):
        st = 1
        Ta = self.Ta
        gbw = self.gbw()
        gsw = self.gsc(An, Tleaf)[0] * 1.6  # co2 - > h2o
        gh = gbw * 0.135 / 0.147
        gv = gsw * gbw / (gsw + gbw)
        gr = 4 * epsilon * SBC * (Tleaf + KELVIN) ** 3 / Cp
        ghr = gh + gr
        thermal_air = epsilon * SBC * (Ta + KELVIN) ** 4
        es_Ta = 0.611 * np.exp(17.27 * Ta / (Ta + 237.3))  # kPa
        RH = min(100, max(self.RH, 10.0))
        D = (1 - RH / 100) * es_Ta  # Vapor pressure deficit
        Ea = RH / 100 * es_Ta  # Ambient vapor pressure
        Icabs = self.Ic / conv * (1 - scatt)
        NIR = self.Ic / 2.  # NIR = solRad - PAR
        Rabs = Icabs + 0.15 * NIR + 2 * (epsilon * SBC * (Ta + KELVIN) ** 4)
        slope = (4098 * 0.6198 * np.exp(17.269 * Ta / (Ta + 237.3))) / (Ta + 237.3) ** 2
        Tlnew = Ta + ((Rabs - thermal_air - LAMBDA * gv * D / P) / (ghr * Cp + LAMBDA * gv * slope / P))
        es_Tlnew = 0.611 * np.exp(17.27 * Tlnew / (Tlnew + 237.3))  # kPa
        E = max(0.0, 2 * gv * (es_Tlnew - Ea) / 101.3)  # mol m-2 leaf s-1 (H2O), both sides of a leaf(*2)
        Emm = E * 3600 * 18 * 10 ** (-3)  # mm m-2 leaf hr-1 (H2O)

        self.gbwcheck = gbw # 파라미터 확인용
        self.gswcheck = gsw # 파라미터 확인용
        self.gh = gh # 파라미터 확인용
        self.gv = gv # 파라미터 확인용
        self.gr = gr # 파라미터 확인용
        self.ghr = ghr # 파라미터 확인용
        self.thermal_air = thermal_air # 파라미터 확인용
        self.RH_eb = RH # 파라미터 확인용
        self.D = D # 파라미터 확인용
        self.Ea = Ea # 파라미터 확인용
        self.Icabs = Icabs # 파라미터 확인용
        self.NIR = NIR # 파라미터 확인용
        self.Rabs = Rabs # 파라미터 확인용
        self.slope = slope # 파라미터 확인용
        self.Tlnew = Tlnew # 파라미터 확인용
        self.es_Tlnew = es_Tlnew # 파라미터 확인용
        self.E = E # 파라미터 확인용
        self.Emm = Emm # 파라미터 확인용

        return Tlnew, Emm

    def interCi(self, An, Tleaf):
        gbc = self.gbc()
        gsc = self.gsc(An, Tleaf)[0]
        # print(Tleaf)
        Ca = self.Ca
        # print(An)
        Ci = Ca - An * (1 / gsc + 1 / gbc) * P
        # print(Ci)
        return Ci

    def gbw(self):  # H2O     co2 = gbw / 1.37
        wind = self.wind
        gbw = 0.147 * np.sqrt(max(0.4, wind) / (0.72 * width))  # outdoor condition C & N (1998)
        return gbw

    def gbc(self):
        gbw = self.gbw()
        # print(gbw)
        return gbw / 1.37

    def gsc(self, An, Tleaf):  # H2O  co2 = gsw / 1.6
        gbc = self.gbc()
        # print(gbc)
        Ca = self.Ca
        Ha = self.RH / 100
        Cs = Ca - An / self.gbc() * 1.013
        if Cs == 0: Cs = 0.001
        if An <= 0: An = 0.1
        stress = self.LWPeffect()
        ah = (stress * g1 * An) / Cs
        bh = g0 + gbc - (stress * g1 * An / Cs)
        # print(bh)
        ch = (-Ha * gbc) - g0
        # print(ch)
        hs = min(max((-bh + np.sqrt(bh * bh - 4 * ah * ch) / (2 * ah)), 0.1), 1.0)
        #         psi_leaf = R * Tleaf * np.log(hs) / 18     # unit : bars
        #         self.psi_leaf = psi_leaf
        # print(hs)
        gsc = min(g0 + stress * g1 * An * hs / Cs, 1.5)
        es_Tleaf = 0.611 * np.exp(17.27 * Tleaf / (Tleaf + 237.3))
        Ds = (1 - hs) * es_Tleaf
        self.cs = Cs # 파라미터 추출용
        self.ah = ah # 파라미터 추출용
        self.bh = bh # 파라미터 추출용
        self.ch = ch # 파라미터 추출용
        self.hs = hs # 파라미터 추출용
        self.gsccheck = gsc
        self.es_Tleaf = es_Tleaf # 파라미터 추출용
        self.Ds = Ds # 파라미터 추출용

        return gsc, Ds

    def LWPeffect(self):  # H2O  stress factor by evaopration
        predawn = self.predawn
        stress = min(1.0, (1 + np.exp(psi_f * s_f)) / (1 + np.exp(s_f * (psi_f - (predawn - psi_th)))))
        return stress

    #################################################
    def Newton_2Var(self, func1, func2, x0, y0, h=0.1, tolerance=1e-8, steps=1000):
        for i in range(1, steps):
            x0 = 0.01 if x0 == 0 else x0
            y0 = 0.01 if y0 == 0 else y0
            A = func1(x0, y0)
            B = (func1(x0 + h, y0) - func1(x0 - h, y0)) / (2 * h)
            C = (func1(x0, y0 + h) - func1(x0, y0 - h)) / (2 * h)
            O = func2(x0, y0)
            P = (func2(x0 + h, y0) - func2(x0 - h, y0)) / (2 * h)
            Q = (func2(x0, y0 + h) - func2(x0, y0 - h)) / (2 * h)
            aaa = C * P - B * Q
            if aaa == 0:
                print("Error: divide by 0")
                aaa = 0.001
            x = x0 + (C * O - A * Q) / (B * Q - C * P)
            y = y0 + (B * O - A * P) / (C * P - B * Q)
            error1 = abs((x - x0) / x)
            error2 = abs((y - y0) / y)
            if (error1 < tolerance) and (error2 < tolerance): break
            x0, y0 = x, y
            self.newci = x
            self.newtleaf = y
        return x, y

    #### Temperature dependence of parameters from Leuning(2002)
    def calcParam(self):
        temp = self.Ta
        Vcmax = self.tempParam2(k25=Vcm25, Ha=Havcm, temp=temp)  # , Hd=Hdvcm, Sv=Svvcm, temp=temp)
        Jmax = self.tempParam1(k25=Jm25, Ha=Hajm, Hd=Hdjm, Sv=Svjm, temp=temp)
        GammaStar = self.tempParam2(k25=gamma25, Ha=Eagamma, temp=temp)
        Kc = self.tempParam2(k25=Kc25, Ha=Eac, temp=temp)
        Ko = self.tempParam2(k25=Ko25, Ha=Eao, temp=temp)
        TPU = self.tempParam2(k25=TPU25, Ha=EaTPU, temp=temp)
        Rd = self.tempParam2(k25=Rd25, Ha=Ear, temp=temp)
        self.Vcmax = Vcmax  # rubisco
        self.Jmax = Jmax  # RuBP
        self.GammaStar = GammaStar  #
        self.TPU = TPU
        self.Kc = Kc
        self.Ko = Ko
        self.Rd = Rd

    def tempParam1(self, k25, Ha, Hd, Sv, temp):
        a = (1 + np.exp((Sv * (25 + KELVIN) - Hd) / (R * (25 + KELVIN)))) / (
                    1 + np.exp((Sv * (temp + KELVIN) - Hd) / (R * (temp + KELVIN))))
        b = np.exp((Ha / (R * (25. + KELVIN))) * (1 - (25. + KELVIN) / (temp + KELVIN)))
        return (k25 * b * a)

    def tempParam2(self, k25, Ha, temp):
        return k25 * np.exp(
            (Ha / (R * (25. + KELVIN)) * (1 - (25. + KELVIN) / (KELVIN + temp))))  # --> arrhenius function
