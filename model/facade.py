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
               'pltDW': [], 'pltFW': [], 'lossPer': [], 'realFW': [], 'bolting': [], 'An':[], 'Ci':[], 'Tleaf':[], 'PPFD':[]}

        check_headers = ['lai', 'Icsun', 'Icsh', 'laisun', 'laish', 'sunhgt', 'aAn', 'bAn', 'leaf', 'assim'] # 중간 점검 생성
        check = {'datetime': [],'doy':[],'dap':[],'hour':[],'LN':[] ,'lai':[], 'Icsun':[], 'Icsh':[], 'laisun':[], 'laish':[], 'sunhgt':[], 'aAn':[], 'bAn':[], 'leaf':[], 'assim':[],'assimSum':[],'dvs': [], 'fgl':[]} # 중간 점검 생성

        growth_param = {'datetime':[], 'assim_c':[], 'ddage':[], 'mgl':[], 'mo':[], 'RM':[], 'tempRM':[], 'RMpr':[], 'fgl':[], 'fr':[], 'fo':[], 'fggl':[], 'fgr':[], 'fgo':[], 'GT':[], 'available':[], 'gr_gl':[],
                        'gr_r':[], 'gr_o':[], 'gr_dl':[], 'wgl':[], 'wr':[], 'wo':[], 'maint':[], 'rootd':[]}

        gasex_param = {'datetime':[], 'An':[], 'gsc':[],'newci':[], 'newtleaf':[], 'Cs':[], 'ah':[], 'bh':[], 'ch':[], 'hs':[],'es_Tleaf':[], 'Ds':[], 'gbw':[], 'gsw':[], 'gh':[], 'gv':[], 'gr':[], 'ghr':[],
                       'thermal_air':[], 'RH_eb':[], 'D':[], 'Ea':[], 'Icabs':[], 'NIR':[], 'Rabs':[], 'slope':[], 'Tlnew':[], 'es_Tlnew':[], 'E':[], 'Emm':[], 'Ansun':[], 'Ansh':[], 'Ci':[], 'Tleaf':[], 'Vcmax':[], 'Jmax':[],
                       "Gammastar":[], 'Kc':[], 'Ko':[]}

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
            # print(doy)
            hour = i.hour
            # print(hour)
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
            # print(self.lai)
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
            # print(Ictot)`
            tot = self.gasEx
            tot.setValue(Ictot, Tair, RH, wind, predawn)
            tot.routine()
            aAn = tot.An

            # fractionation PAR
            # doy, hour, lai, Irrad
            self.frac.radFraction(doy=doy, hour=hour, PPFD=Ictot, LAI=self.lai)
            Icsun = self.frac.Icsun  # umol m-2 s-1
            # print(Icsun)
            Icsh = self.frac.Icsh  # umol m-2 s-1

            # fractionation lai
            self.frac.laiFraction(doy=doy, hour=hour, LAI=self.lai)
            laiSun = self.frac.laiSun
            laiSh = self.frac.laiSh
            # print(laiSun)
            sunhgt = self.frac.sunhgt # 중간 점검 생성

            # fractionation Vcmax
            # self.frac.rubFraction(doy=doy, hour=hour, LAI=self.lai, Vcmax=110)
            # VcmaxSun = self.frac.VcmaxSun
            # VcmaxSh  = self.frac.VcmaxSh

            #######################################################
            An = self.gasEx.Ancheck # 파라미터 추출용
            gsc = self.gasEx.gsccheck # 파라미터 추출용
            # mid_Tleaf = self.gasEx.mid_Tleaf # 파라미터 추출용
            Ci = self.gasEx.Ci
            Tleaf = self.gasEx.Tleaf
            # print(Tleaf)
            newci = self.gasEx.newci # 파라미터 추출용
            newtleaf = self.gasEx.newtleaf # 파라미터 추출용
            # gsc
            cs = self.gasEx.cs # 파라미터 추출용
            ah = self.gasEx.ah # 파라미터 추출용
            bh = self.gasEx.bh # 파라미터 추출용
            ch = self.gasEx.ch # 파라미터 추출용
            hs = self.gasEx.hs # 파라미터 추출용
            es_Tleaf = self.gasEx.es_Tleaf # 파라미터 추출용
            Ds = self.gasEx.Ds # 파라미터 추출용
            #EnergyBal
            gbw = self.gasEx.gbwcheck # 파라미터 추출용
            gsw = self.gasEx.gswcheck # 파라미터 추출용
            gh = self.gasEx.gh # 파라미터 추출용
            gv = self.gasEx.gv # 파라미터 추출용
            gr = self.gasEx.gr # 파라미터 추출용
            ghr = self.gasEx.ghr # 파라미터 추출용
            thermal_air = self.gasEx.thermal_air # 파라미터 추출용
            RH_eb = self.gasEx.RH_eb # 파라미터 추출용
            D = self.gasEx.D # 파라미터 추출용
            Ea = self.gasEx.Ea # 파라미터 추출용
            Icabs = self.gasEx.Icabs # 파라미터 추출용
            NIR = self.gasEx.NIR # 파라미터 추출용
            Rabs = self.gasEx.Rabs # 파라미터 추출용
            slope = self.gasEx.slope # 파라미터 추출용
            Tlnew = self.gasEx.Tlnew # 파라미터 추출용
            es_Tlnew = self.gasEx.es_Tlnew # 파라미터 추출용
            E = self.gasEx.E # 파라미터 추출용
            Emm = self.gasEx.Emm # 파라미터 추출용
            Vcmax = self.gasEx.Vcmax
            Jmax = self.gasEx.Jmax
            Gammastar = self.gasEx.GammaStar
            Kc = self.gasEx.Kc
            Ko = self.gasEx.Ko


            #######################################################

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
            # print(assim)
            self.assimSum += assim  # assimilate summation
            self.assim = assim  # hourly assimilate, CH2O g m-2gr hr-1
            # print(self.assimSum)

            # partitioning assimilates
            part = self.growth
            part.growCalc(Tair, assim, dvs, RDT=1.0)  # water stress
            # print(assim)
            leaf = part.wgl
            # print(leaf)
            root = part.wr # 파라미터 추출용
            balt = part.wo # 파라미터 추출용
            maint = part.maint # 파라미터 추출용

            ####################### growth.py param # 파라미터 추출용
            rootd = part.rootd
            assim_c = part.assim_c
            ddage = part.ddage
            mgl = part.mgl
            mo = part.mo
            RM = part.RM
            tempRM = part.tempRM
            RMpr = part.RMpr
            fgl = part.fgl
            fr = part.fr
            fo = part.fo
            fggl = part.fggl
            fgr = part.fgr
            fgo = part.fgo
            GT = part.GT
            available = part.available
            gr_gl = part.gr_gl
            gr_r = part.gr_r
            gr_o = part.gr_o
            gr_dl = part.gr_dl
            #######################


            # DW per plant
            pltDW = leaf / self.plantDensity

            # print(leaf)
            # convert DW to FW
            # FDratio = 0.1632 * leafNumber + 11.049
            # FDratio  = 6.949*np.log(leafNumber) - 6.5953
            # pltFW  = FDration * plt
            leafNumber = self.leafNumber
            FDratio = 6.949 * np.log(leafNumber) - 6.5953                                  ####### 2021 워싱턴대에서 변경 ############
            pltFW = FDratio * pltDW                                                        ##   * self.hi  -> 삭제   """
            # print(pltDW)

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
            res['An'].append(An)
            res['Ci'].append(Ci)
            res['Tleaf'].append(Tleaf)
            res['PPFD'].append(Ictot)

            # 중간 점검 생성
            check['datetime'].append(i)
            check['doy'].append(doy)
            check['dap'].append(dap)
            check['hour'].append(hour)
            check['LN'].append(LN)
            check['lai'].append(lai)
            check['Icsun'].append(Icsun)
            check['Icsh'].append(Icsh)
            check['laisun'].append(laiSun)
            check['laish'].append(laiSh)
            check['sunhgt'].append(sunhgt)
            check['aAn'].append(aAn)
            check['bAn'].append(bAn)
            check['leaf'].append(leaf)
            check['assim'].append(assim)
            check['dvs'].append(dvs)
            check['assimSum'].append(self.assimSum)
            check['fgl'].append(fgl)


            growth_param['datetime'].append(i)
            growth_param['rootd'].append(rootd)
            growth_param['assim_c'].append(assim_c)
            growth_param['ddage'].append(ddage)
            growth_param['mgl'].append(mgl)
            growth_param['mo'].append(mo)
            growth_param['RM'].append(RM)
            growth_param['tempRM'].append(tempRM)
            growth_param['RMpr'].append(RMpr)
            growth_param['fgl'].append(fgl)
            growth_param['fr'].append(fr)
            growth_param['fo'].append(fo)
            growth_param['fggl'].append(fggl)
            growth_param['fgr'].append(fgr)
            growth_param['fgo'].append(fgo)
            growth_param['GT'].append(GT)
            growth_param['available'].append(available)
            growth_param['gr_gl'].append(gr_gl)
            growth_param['gr_r'].append(gr_r)
            growth_param['gr_o'].append(gr_o)
            growth_param['gr_dl'].append(gr_dl)
            growth_param['wgl'].append(leaf)
            growth_param['wr'].append(root)
            growth_param['wo'].append(balt)
            growth_param['maint'].append(maint)

            gasex_param['datetime'].append(i)
            gasex_param['An'].append(An)
            gasex_param['gsc'].append(gsc)
            gasex_param['newci'].append(newci)
            gasex_param['newtleaf'].append(newtleaf)
            gasex_param['Cs'].append(cs)
            gasex_param['ah'].append(ah)
            gasex_param['bh'].append(bh)
            gasex_param['ch'].append(ch)
            gasex_param['hs'].append(hs)
            gasex_param['es_Tleaf'].append(es_Tleaf)
            gasex_param['Ds'].append(Ds)
            gasex_param['gbw'].append(gbw)
            gasex_param['gsw'].append(gsw)
            gasex_param['gh'].append(gh)
            gasex_param['gv'].append(gv)
            gasex_param['gr'].append(gr)
            gasex_param['ghr'].append(ghr)
            gasex_param['thermal_air'].append(thermal_air)
            gasex_param['RH_eb'].append(RH_eb)
            gasex_param['D'].append(D)
            gasex_param['Ea'].append(Ea)
            gasex_param['Icabs'].append(Icabs)
            gasex_param['NIR'].append(NIR)
            gasex_param['Rabs'].append(Rabs)
            gasex_param['slope'].append(slope)
            gasex_param['Tlnew'].append(Tlnew)
            gasex_param['es_Tlnew'].append(es_Tlnew)
            gasex_param['E'].append(E)
            gasex_param['Emm'].append(Emm)
            gasex_param['Ansun'].append(aAn)
            gasex_param['Ansh'].append(bAn)
            gasex_param['Ci'].append(Ci)
            gasex_param['Tleaf'].append(Tleaf)
            gasex_param['Vcmax'].append(Vcmax)
            gasex_param['Jmax'].append(Jmax)
            gasex_param['Gammastar'].append(Gammastar)
            gasex_param['Kc'].append(Ko)
            gasex_param['Ko'].append(Ko)


            self.hourresult = res  # store the model results
            self.checkresult = check # 중간 점검 생성
            self.growth_parameter = growth_param
            self.gasex_parameter = gasex_param


        
        ## Message for abnomal plant growth
        if self.dvs > 1.0:
            print("Be careful Balting!!")
        if self.lossPer > 10:
            print("Be careful Heat Damage!!")

        out = pd.DataFrame(self.hourresult).set_index('datetime')
        out_test = pd.DataFrame(self.checkresult).set_index('datetime') # 중간 점검 생성
        out_growth_param = pd.DataFrame(self.growth_parameter).set_index('datetime') # 파라미터 추출용
        out_gasex_param = pd.DataFrame(self.gasex_parameter).set_index('datetime') # 파라미터 추출용

        out_test.to_csv('output/check_hour.csv') # 파라미터 추출용
        out_growth_param.to_csv('../ipynb/growth/growth_param.csv') # 파라미터 추출용
        out_gasex_param.to_csv('../ipynb/gasexchange/gasex_param.csv') # 파라미터 추출용

        aggrigation = { 'dvs':'mean', 'LN' : 'mean', 'lai' : 'mean', 'pltDW':'mean', 'pltFW':'mean', 
                       'lossPer':'mean', 'realFW':'mean', 'bolting':'mean' ,'An':'mean', 'Ci':'mean', 'Tleaf':'mean', 'PPFD':'mean'}
        self.dayresult = out.groupby(out.index.date).agg(aggrigation)
        # self.timeresult = out_test # 중간 점검 생성


        # write df to various file
        cols = ['dvs', 'LN', 'lai', 'pltDW', 'pltFW', 'realFW', 'lossPer', 'bolting', 'An', "Ci", "Tleaf", "PPFD"]
        header  = ['dvs', 'leaf.number', 'LAI', 'DW(g/m2)', 'pot.FW(g)', 'real.FW(g)', 'loss(%)', 'bolting','An', "Ci", "Tleaf", 'PPFD']

        # check_cols = ['lai', 'Icsun', 'Icsh', 'laisun', 'laish', 'sunhgt', 'aAn', 'bAn', 'leaf', 'assim'] # 중간 점검 생성
        # check_headers = ['lai', 'Icsun', 'Icsh', 'laisun', 'laish', 'sunhgt', 'aAn', 'bAn', 'leaf', 'assim'] # 중간 점검 생성
        fname = self.fname_out

        ext = fname.split('.')[1]
        if ext == 'xlsx':
            self.dayresult[cols].to_excel(pd.ExcelWriter(fname), float_format = '%.1f', header=header)
            # self.timeresult[check_cols].to_excel(pd.ExcelWriter(fname), float_format = '%.1f', header=check_headers) # 중간 점검 생성
        elif ext == 'csv':
            self.dayresult[cols].to_csv(fname, sep=',', float_format='%.1f', header=header)
            # self.timeresult[check_cols].to_csv(fname, sep=',', float_format='%.1f', header=check_headers) # 중간 점검 생성
        elif ext == 'txt':
            self.dayresult[cols].to_csv(fname, sep='\t', float_format='%.1f', header=header)
            # self.timeresult[check_cols].to_csv(fname, sep='\t', float_format='%.1f', header=check_headers) # 중간 점검 생성
        else:
            raise Exception("Only possible file name is *.xlsx, *.csv, *.txt.")
       
        print('\ndone.')  

            

