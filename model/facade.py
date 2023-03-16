#     model = Facade(inifile, outfile, pltday=startdoy, iniLN=iniLN, density=density, latitude=latitude)
#    model.run(start, end, lnratio)
import pandas as pd
import numpy as np

import daytohour as dh
import disease as dl
import fraction as fr
import gasexchange as gas
import growth as gro
import stage as stg
import weather as wt


# import pmethr as pmet
# import soilwater as soil


class Facade(object):

    def __init__(self, dorh, fname_in, fname_out, start, end, iniLN, density, latitude, pLeafForm, daysRoot): ## pLeafForm, daysRoot 변경 

        # assign initial values
        self.pltday = start
        self.plantDensity = density
        self.initialLeafNumber = iniLN
        self.latitude = float(latitude)
        self.pLeafForm =  pLeafForm               ## """ 2021 UW 추가 """
        self.daysRoot = daysRoot                  ## """ 2021 UW 추가 """

        # initialize attributes:
        self.dvs = 0.0
        self.leafNumber = 0.0
        self.lai = 0.0
        iniDW = 0.0016 * iniLN ** 2.5917
        self.assimSum = iniDW * self.plantDensity
        self.assim = 0.0
        self.lossPer = 0.0
        self.hourresult = None
        self.dayresult = None
        self.fname_out = fname_out

        # create instances
        if dorh == 'd':
            dayData = dh.DayToHour(fname_in, start, end, self.latitude)
            dayData.calcHour()
            wthr = dayData
        elif dorh == 'h':
            wthr = wt.Weather(fname_in, start, end)
        else:
            raise Exception("Error: Weather data type (d or h)")
        self.wthr = wthr
        self.bd = stg.BD(initialLeafNumber=iniLN, plantDensity=density,  pLeafForm=pLeafForm, daysRoot=daysRoot)  ## 2021 변경 """
        self.frac = fr.Fractionation(latitude=self.latitude)
        self.gasEx = gas.GasExchange()
        self.growth = gro.Growth()
        self.loss = dl.Disease()

        # wet or dry regime, assume disease infects wet regime only 
        self.count = 0.0  # for count rain effect duration (7 days ?)

    def run(self, file_write=True):
        # output results as dictionary
        headers = ['datetime', 'LN', 'lai', 'dvs', 'pltDW', 'pltFW', 'lossPer', 'realFW', 'bolting']
        res = {'datetime': [], 'LN': [], 'lai': [], 'dvs': [],
               # 'netTmm':[], 'crop_st':[],'soil_st':[], 'pmet':[], 'pmetc':[], 'pmets':[],
               'pltDW': [], 'pltFW': [], 'lossPer': [], 'realFW': [], 'bolting': []}

        # run the model and after every run,
        print('\nRunning ...\n')
        
#         Web에 보일 때 필요한 내용
#         rows = [
#             headers
#         ]

        wdata = self.wthr.data
        for i, v in wdata.iterrows():
            # data retrieval from wdata dataframe
            date = i
            doy = i.timetuple().tm_yday
            hour = i.hour
            Irrad = v.Irrad
            Tair = v.Tair
            rain = v.rain
            wind = v.wind
            RH   = v.RH
            dap  = v.dap                                         ####### 2021 추가

            # calculate lai and dvs from temp
            self.bd.calcBD(Tair, dap)
            self.leafNumber = self.bd.leafNumber                 ####### 2021 워싱턴대에서 변경 ############
            LN = self.leafNumber                                                   ##    * self.lnratio  ->  삭제"""

            self.lai = self.bd.lai  # green leaf part, effective lai
            lai = self.lai
            self.dvs = self.bd.dvs
            dvs = self.dvs
            if dvs < 1.0:
                bolting = 0
            else:
                bolting = 1

            # 'datetime', 'LN', 'lai', 'dvs', 'netTmm', 'crop_st','soil_st',
            #         'pmet', 'pmetc', 'pmets','pltDW','pltFW', 'lossPer', 'realFW', 'bolting'
            # fout.write(','.join((str(x) for x in (i.timestamp(), LN, lai, dvs))))

            predawn = 0

            # total assimilate
            Ictot = Irrad * 4.57 / 2  # from Wm-2 to umol m-2 s-1 and total Irrad to PAR
            tot = self.gasEx
            tot.setValue(Ictot, Tair, RH, wind, predawn)
            tot.routine()
            aAn = tot.An

            # fractionation PAR
            # doy, hour, lai, Irrad
            self.frac.radFraction(doy=doy, hour=hour, PPFD=Ictot, LAI=self.lai)
            Icsun = self.frac.Icsun  # umol m-2 s-1
            Icsh = self.frac.Icsh  # umol m-2 s-1

            # fractionation lai
            self.frac.laiFraction(doy=doy, hour=hour, LAI=self.lai)
            laiSun = self.frac.laiSun
            laiSh = self.frac.laiSh

            # fractionation Vcmax
            # self.frac.rubFraction(doy=doy, hour=hour, LAI=self.lai, Vcmax=110)
            # VcmaxSun = self.frac.VcmaxSun
            # VcmaxSh  = self.frac.VcmaxSh

            # for sunlit leaves
            sun = self.gasEx
            sun.setValue(Icsun, Tair, RH, wind, predawn)
            sun.routine()
            aWc = sun.Wc
            aWj = sun.Wj
            aAn = sun.An

            # for shaded leaves
            shd = self.gasEx
            shd.setValue(Icsh, Tair, RH, wind, predawn)
            shd.routine()
            bWc = shd.Wc
            bWj = shd.Wj
            bAn = shd.An

            # net assimilates (umol CO2 m-2 ground s-1)
            netA = aAn * laiSun + bAn * laiSh  # from umolCO2 m-2leaf s-1 to umolCO2 m-2ground s-1
            assim = netA * 3600 * 30 * 10 ** (-6)  # convert CO2 umol m-2gr s-1 to CH2O g m-2gr hr-1
            self.assimSum += assim  # assimilate summation
            self.assim = assim  # hourly assimilate, CH2O g m-2gr hr-1

            # partitioning assimilates
            part = self.growth
            part.growCalc(Tair, assim, dvs, RDT=1.0)  # water stress
            leaf = part.wgl
            root = part.wr
            balt = part.wo
            maint = part.maint

            # DW per plant
            pltDW = leaf / self.plantDensity

            # convert DW to FW
            # FDratio = 0.1632 * leafNumber + 11.049
            # FDratio  = 6.949*np.log(leafNumber) - 6.5953
            # pltFW  = FDration * plt
            leafNumber = self.leafNumber
            FDratio = 6.949 * np.log(leafNumber) - 6.5953                                  ####### 2021 워싱턴대에서 변경 ############
            pltFW = FDratio * pltDW                                                        ##   * self.hi  -> 삭제   """

            ## loss by Heat or Disease
            # determine wet or dry regime  (for 7 days after rain and rh > 98
            if rain > 0:
                self.count = 7 * 24  # 7 days and 24 hours
            else:
                self.count -= 1
                self.count = max(self.count, 0)
            if self.count >= 1 or RH > 98:
                wet = True
            else:
                wet = False

            # calculate heat or disease loss                         ####### 2021 워싱턴대에서 변경 ######################
            if LN > 74*0.78:                                        ## 변경 : 35장 -> 57장   (2021 실험 74장 => 74*0.78 = 57.7"""
                loss = self.loss
                loss.calc(Tair, wet)
                ds = loss.ds
                lossPer = ds * 100
            else:
                lossPer = 0
            self.lossPer = lossPer
                                                                    ####### 2021 워싱턴대에서 변경 ########################
            if self.lossPer != 0:
                self.pLeafForm = 0.0                              ##  추가 : 무름증 발생 후 엽 증가속도 0.0 (2021 실험) """          
            realFW = pltFW * (1 - self.lossPer / 100)  ## FW per plant after loss

            # retrieve the model results
            lai = self.lai
            res['datetime'].append(i)
            res['LN'].append(LN)
            res['lai'].append(lai)
            res['dvs'].append(dvs)
            res['pltDW'].append(pltDW)
            res['pltFW'].append(pltFW)
            res['lossPer'].append(lossPer)
            res['realFW'].append(realFW)
            res['bolting'].append(bolting)

            self.hourresult = res  # store the model results
            

        
        ## Message for abnomal plant growth
        if self.dvs > 1.0:
            print("Be careful Balting!!")
        if self.lossPer > 10:
            print("Be careful Heat Damage!!")

        out = pd.DataFrame(self.hourresult).set_index('datetime')
        aggrigation = { 'dvs':'mean', 'LN' : 'mean', 'lai' : 'mean', 'pltDW':'mean', 'pltFW':'mean', 
                       'lossPer':'mean', 'realFW':'mean', 'bolting':'mean' }
        self.dayresult = out.groupby(out.index.date).agg(aggrigation)
        
        # write df to various file
        cols = ['dvs', 'LN', 'lai', 'pltDW', 'pltFW', 'realFW', 'lossPer', 'bolting']
        header  = ['dvs', 'leaf.number', 'LAI', 'DW(g/m2)', 'pot.FW(g)', 'real.FW(g)', 'loss(%)', 'bolting']
        
        fname = self.fname_out
        
        ext = fname.split('.')[1]
        if ext == 'xlsx':
            self.dayresult[cols].to_excel(pd.ExcelWriter(fname), float_format = '%.1f', header=header)
        elif ext == 'csv':
            self.dayresult[cols].to_csv(fname, sep=',', float_format='%.1f', header=header)
        elif ext == 'txt':
            self.dayresult[cols].to_csv(fname, sep='\t', float_format='%.1f', header=header)
        else:
            raise Exception("Only possible file name is *.xlsx, *.csv, *.txt.")
       
        print('\ndone.')  

            

