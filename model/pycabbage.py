a = """
===============================================================================
===============================================================================  
=============           Chinese cabbage model (v.2.0)          ================
===============================================================================
                           - for prediction of Chinese cabbage growth and yield.

How to use:
    1. to double click "pycabbage.exe"
    2. to type requested data
    3. Weatherdata file should be same directory of "pacabbage.exe".
    4. Output file will be made at same directory.
    
                                                                    by K H Moon                                      
===============================================================================
===============================================================================
"""
b = """
===============================================================================
How to input:
  == Set simulation coditions ==
    ? latitude(default = 37.0) : -> type latitude or just enter
    ? weather file name : -> type weather file name
      * (.csv, .txt, .xlsx, .xls) file types are possible 
      * Data file shoud have specific name and format
        - daily name(format)  : date(YYYY-mm-dd), sunhr, Tmax, Tmin, rain,
                                wind, RH (all float number)
        - hourly name(format) : timestamp(YYYY-mm-dd HH:MM), Irrad, Tair,
                                rain, wind, RH (all float number)             
        - unit: sunhr(hours), Tmax(C), rain(mm), wind(m/s), RH(%), Irrad(W/m2)
    ? is weather data daily(d) or hourly(h) : -> type 'd' or 'h'
    ? output file name : -> type output file name
      * (*.csv, *.txt, *.xlsx) file names are possible    
    ? planting date(YYYY-mm-dd) : -> type YYYY-mm-dd 
    ? end date(YYYY-mm-dd) : -> type YYYY-mm-dd
    
  == Set initial conditions ==
    ? planting density : -> type greater than 3.6 or just enter
    ? leaf number at planting : -> type leaf number(default=8) or just enter
===============================================================================
"""
"""
Example: 
    pycabbage.exe (double click)
    
       === Set simulation conditions. ===
       ? latitude : 33.2 (enter)
       ? weather file name : jeju2018.csv (enter)
       ? is weather data daily(d) or hourly(h) : d (enter)
       ? output file name(including extension) : jejuoutput.xlsx (enter)
       ? planting date(YYYY-mm-dd) : 2018-4-3 (enter)
       ? end date(YYYY-mm-dd) : 2018-7-12 (enter)
       
       === Set initial condition.===
       ? planting density(plants/m2, default = 3.8) : (enter)
       ? leaf number at planting(default = 6) : 6.0 (enter)
       ? days for rooting after transplanting(default = 10) : 10 (enter)
       ? parameter for daily leaf formation rate(0 ~ 1, default = 0.69) : 0.69 (enter)
       
    then model will be working and make a output file (named "jejuoutput.xlsx")
    
    It means the model, located at 33.2 latitude, wether file is 'jeju2018.csv', and output file is 'jejuoutput.xlsx',
    and both these files are stored in same folder. The model is to be run from 2018-4-3 to 2018-7-12.

@author Kyung Hwan Moon (RICCA, RDA, Korea)

"""
import sys
from facade import Facade   
import pandas as pd
from matplotlib import dates
import matplotlib.pyplot as plt

def run():
    """Main entry point for the program.

    Run the model simulation (hourly or daily)

    """

    # set the accepted flags, and parse the options and arguments:
    inifile = outfile = None
    # print(a)
    # print(b)
    # print("\nLet's try model working!!")
    # print("\n=== Set simulation conditions. === ")
    # latitude  = float(input("? latitude(default = 37.0) : ") or "37.0")
    # inifile   = input("? weather file name(including extension) : ")
    # dorh      = input("? is weather data daily(d) or hourly(h) : ")
    # outfile   = input("? output file name(including extension) : ")
    # start     = input("? planting date(YYYY-mm-dd) : ")
    # end       = input("? end date(YYYY-mm-dd) : ")
    #
    # sp_start = start.split('-')
    # sp_end   = end.split('-')
    # if len(sp_start) < 3 and len(sp_end) < 3 :
    #     raise Exception("Planting date format is YYYY-mm-dd ")
    #
    # # set growing condition
    # print("\n=== Set initial condition.===")                                                               ####### 2021 워싱턴대에서 변경 ############
    # density  = float(input('? planting density(plants/m2, default = 3.8) : ') or "3.8")             ## 4.0 -> 3.8 로 변경 (고랭지 평균)
    # iniLN    = float(input('? leaf number at planting(default = 6) : ') or "6.0")                          ## 8.0 -> 6.0 으로 변경
    # daysRoot = float(input('? days for root adaptation after transplanting(default = 10) : ') or "10.0")           ## 정식 후 활착 기간 (평균 14일)
    # pLeafForm = float(input('? suppression parameter for leaf formation rate(default = 0.69): ') or "0.69")     ## 일엽생성속도 장해 계수(0~1), 춘광 0.69

    # calculation condition
#    print("\n=== Set ratio of apparent leaf number.===")
#    lnratio = float(input('? ratio of leaf number(default = 0.78) : '))

    latitude = 37.0
    inifile = 'input/ricca17.csv'
    dorh = 'h'
    outfile = 'output/ricca17_test_result.csv'
    start = '2017-9-20'
    end = '2017-12-30'

    sp_start = start.split('-')
    sp_end = end.split('-')
    if len(sp_start) < 3 and len(sp_end) < 3:
        raise Exception("Planting date format is YYYY-mm-dd ")

    density = 3.8
    iniLN = 6
    daysRoot = 10
    pLeafForm = 0.69


    if None in [inifile, outfile]:
        print('One or more parameters are missing.')
        print(__doc__)
        sys.exit()
                                                                                                         ####### 2021 워싱턴대에서 변경 ######
                                                                                                         ## lnratio=0.78, hi=0.6) - 삭제
    model = Facade(dorh, inifile, outfile, start, end, iniLN=iniLN, density=density, latitude=latitude, \
                   pLeafForm=pLeafForm, daysRoot=daysRoot) 
    model.run()

    # inputfile = pd.read_csv(inifile)
    origin = pd.read_csv('output/test_normal.csv')
    diff = pd.read_csv(outfile)

    # print(origin)
    out = origin.merge(diff, how='left', on='Unnamed: 0').set_index('Unnamed: 0')

    fig, ax = plt.subplots(3, 2, figsize=(10, 10),gridspec_kw={"wspace":0.15, "hspace":0.2}, sharex=True)
    # plt.subplot(3, 2, 1).title.set_text('dvs')
    ax[0, 0].title.set_text('dvs')
    ax[0, 0].plot(out[['dvs_x', 'dvs_y']])
    ax[0, 0].set_xticks(ticks=out.index)
    # ax[2, 0].tick_params(labelrotation=45, labelsize=7)
    for tick in ax[2,0].get_xticklabels():
        tick.set_rotation(45)
    for tick in ax[2,1].get_xticklabels():
        tick.set_rotation(45)
    ax[2, 0].xaxis.set_major_locator(dates.WeekdayLocator(interval=3))
    ax[2, 1].xaxis.set_major_locator(dates.WeekdayLocator(interval=3))

    plt.subplot(3, 2, 2).title.set_text('leaf_number')
    plt.plot(out[['leaf.number_x', 'leaf.number_y']])

    plt.subplot(3, 2, 3).title.set_text('LAI')
    plt.plot(out[['LAI_x', 'LAI_y']])

    plt.subplot(3, 2, 4).title.set_text('DW')
    plt.plot(out[['DW(g/m2)_x', 'DW(g/m2)_y']])

    plt.subplot(3, 2, 5).title.set_text('FW')
    plt.plot(out[['pot.FW(g)_x', 'pot.FW(g)_y']])

    plt.subplot(3, 2, 6).title.set_text('bolting')
    plt.plot(out[['bolting_x', 'bolting_y']])

    # print()

    plt.savefig('graph/결과비교.png')

    out.to_csv('output/difference.csv')


if __name__ == '__main__':
    run()

