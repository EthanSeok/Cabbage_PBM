import numpy as np

## initial condition
# inileaf = 8            # initial leaf number, default    
# plantDensity = 4.0     # plant density per m2, LA = np.exp(1.80 * np.log(leafNumber) + 1.48) , default

## parameters 
# parameters for leaf number increasing
Rxleaf  = 1.838
Txleaf  = 34.7
Toleaf  = 20.57 

# temp. function sum for veg. and rep. stage
# tempdvs = sumtemp/satVS or 1 + (sumtemp - satVS)/(satRS - satVS) or 2
satTemp   = 30       # requirement sum of temp_rate for flowering
satRep    = 150      # requirement sum of temp_rate for harvesting
optTemp   = 15.7       # max. temp for veg. and rep, growth, (temp_rate = np.exp(-1*(np.log(temp/optTemp)**2) )   

# ver. function sum for vernalization
# verdvs = sumver / satVer or 1
satVer    = 52       # sum of ver_rate for flowering
optVer    = 6.2      # opt. temp. for vernalization , (verrate = np.exp(-1*(np.log(temp/optVer)**4) )        

# dvs = tempdvs * verdvs

#timestep 
time = 'hour'
conv = 1/24  if time == 'hour' else 1

class BD():      # balting development
 
    def __init__(self, initialLeafNumber, plantDensity, pLeafForm, daysRoot):
        self.leafNumber    = int(initialLeafNumber)
        self.plantDensity  = plantDensity
        self.pLeafForm =  pLeafForm                ###""" 2021 UW 추가 """
        self.daysRoot = daysRoot                   ###""" 2021 UW 추가 """        
        ## self.appLeafNumber = 0.0                ### 2021 UW 삭제 """
       
        self.dvs     = 0.0
        self.sumTemp = 0.0
        self.sumVer  = 0.0
        self.tempdvs = 0.0
        self.verdvs  = 0.0
        
        self.lai     = 0.0

########### 2021 UW 추가 #################################################################################
## 엽수증가식 -> 정식초기(earlyRateLN), 중기(midRateLN), 무름증 후(0.0)로 구분해서 계산
## 정식초기는 평균 정식 후 14일간, 입력받아 조정 가능함(daysRoot)
    def midRateLN(self, Ta):
        # calculation leaf number
        if (Ta > 0.0) and (Ta < Txleaf): 
            leafRate = Rxleaf *((Txleaf-Ta)/(Txleaf-Toleaf))*(Ta/Toleaf)**(Toleaf/(Txleaf-Toleaf))
            leafRate = leafRate * conv
        else:
            leafRate = 0.0
        return leafRate
    
    def earlyRateLN(self, Ta):
        return self.midRateLN(Ta)*0.05
    
########### 2021 UW 추가 #################################################################################
## plant green leaf Area calculation from green leaf area of each leaves using total leaf number
## 식물체 개체 수준의 녹엽면적만 계산(엽수 이용)

    def eachLenDistribution(self, leafnumber):  ## internal function
        eachLenDist = [] 
        if leafnumber < 8:                    ## limit of leafnumber 8장 보다 작을 경우 오류 방지
            for i in range(leafnumber):
                eachLenDist.append(1.0)
        else:
            a = 20.347*np.log(leafnumber) - 40.703
            b = 0.2086*(leafnumber) + 1.6862
            for i in range(1,leafnumber+1):
                eachLen = a * np.exp(-0.5 * ((i - b) / b)**2)
                eachLenDist.append(eachLen)
        return eachLenDist

    def eachLeafArea(self, eachLenDist):       ## internal function
        eachLeafArea = [0.3512*each**2 + 1.1328*each for each in eachLenDist]
        return eachLeafArea

    def eachBladeArea(self, eachLeafArea):     ## internal function
        ratio = [0.9217*np.exp(-0.01*order) for order in range(1, len(eachLeafArea)+1)]
        eachBladeArea = [a*b for a, b in zip(ratio, eachLeafArea)]
        return eachBladeArea

    def numGreenLeaf(self, eachBladeArea):     ## internal function
        numberGL = 0.9333*np.exp(-0.017*len(eachBladeArea))*len(eachBladeArea)
        numberGreenLeaf = int(numberGL)
        return numberGreenLeaf

    def plantGreenLeafArea(self, leafnumber):
        a = self.eachLenDistribution(leafnumber)
        b = self.eachLeafArea(a)
        c = self.eachBladeArea(b)
        d = self.numGreenLeaf(c)
        greenBladeArea = c[0:d]
        greenLeafArea = sum(greenBladeArea)
        return greenLeafArea

########### 2021 UW 추가 #################################################################################
## 추대 모형 -> 온도를 이용해서 계산  dvs = verdvs * tempdvs
## 
    def calcVerdvs(self, Ta):
        Ta  = max(Ta, 0.01)
        rate = np.exp(-1*(np.log(Ta/optVer)**4))
        self.sumVer += rate * conv
        self.verdvs = max(1, self.sumVer/satVer)
        
    def calcTempdvs(self, Ta):
        Ta  = max(Ta, 0.01)
        rate = np.exp(-1*(np.log(Ta/optTemp)**2))
        self.sumTemp += rate * conv
        sumTemp = self.sumTemp
        if sumTemp >= satRep:
            tempdvs = 2.0
        elif sumTemp > satTemp:
            tempdvs = 1 + (sumTemp - satTemp) /(satRep - satTemp)
        else:
            tempdvs = sumTemp / satTemp
        self.tempdvs = tempdvs
        

##############################################################################################################        
     
        
    def calcBD(self, Ta, dap):
        
        # calculation leaf number
        if dap < self.daysRoot:                                ###""" 2021 UW 추가 """
            leafRate = self.earlyRateLN(Ta)                    ###""" 2021 UW 추가 """
        else:                                                  ###""" 2021 UW 추가 """
            leafRate = self.midRateLN(Ta)*self.pLeafForm       ###""" 2021 UW 추가 """
        self.leafNumber += leafRate
        leafNumber = self.leafNumber

        # calculation dvs
        self.calcVerdvs(Ta)
        self.calcTempdvs(Ta)
        self.dvs = self.tempdvs * self.verdvs
        
#         # no increasing leaf number by aging
#         dvs = self.dvs
#         if dvs > 1.0:
#             self.leafNumber -= (0.5*leafRate)
            
        ###### LAI calculation
        leafNumber = self.leafNumber
        plantDensity = self.plantDensity
        plantGreenLeafArea = self.plantGreenLeafArea(int(leafNumber))
        lai  =  plantGreenLeafArea * plantDensity / 10000
        self.lai = lai
        
    

    

