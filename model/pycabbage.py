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
import datetime
from facade import Facade   

def run():
    """Main entry point for the program.

    Run the model simulation (hourly or daily)

    """

    # set the accepted flags, and parse the options and arguments:
    inifile = outfile = None
    print(a)
    print(b)
    print("\nLet's try model working!!")
    print("\n=== Set simulation conditions. === ")
    latitude  = float(input("? latitude(default = 37.0) : ") or "37.0")
    inifile   = input("? weather file name(including extension) : ")
    dorh      = input("? is weather data daily(d) or hourly(h) : ")
    outfile   = input("? output file name(including extension) : ")
    start     = input("? planting date(YYYY-mm-dd) : ")
    end       = input("? end date(YYYY-mm-dd) : ")
    
    sp_start = start.split('-')
    sp_end   = end.split('-')
    if len(sp_start) < 3 and len(sp_end) < 3 : 
        raise Exception("Planting date format is YYYY-mm-dd ")
    
    # set growing condition
    print("\n=== Set initial condition.===")                                                               ####### 2021 워싱턴대에서 변경 ############
    density  = float(input('? planting density(plants/m2, default = 3.8) : ') or "3.8")             ## 4.0 -> 3.8 로 변경 (고랭지 평균)
    iniLN    = float(input('? leaf number at planting(default = 6) : ') or "6.0")                          ## 8.0 -> 6.0 으로 변경
    daysRoot = float(input('? days for root adaptation after transplanting(default = 10) : ') or "10.0")           ## 정식 후 활착 기간 (평균 14일)
    pLeafForm = float(input('? suppression parameter for leaf formation rate(default = 0.69): ') or "0.69")     ## 일엽생성속도 장해 계수(0~1), 춘광 0.69

    # calculation condition
#    print("\n=== Set ratio of apparent leaf number.===")
#    lnratio = float(input('? ratio of leaf number(default = 0.78) : '))

    
    if None in [inifile, outfile]:
        print('One or more parameters are missing.')
        print(__doc__)
        sys.exit()
                                                                                                         ####### 2021 워싱턴대에서 변경 ######
                                                                                                         ## lnratio=0.78, hi=0.6) - 삭제
    model = Facade(dorh, inifile, outfile, start, end, iniLN=iniLN, density=density, latitude=latitude, \
                   pLeafForm=pLeafForm, daysRoot=daysRoot) 
    model.run()

if __name__ == '__main__':
    run()
