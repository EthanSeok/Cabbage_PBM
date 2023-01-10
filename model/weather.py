# weather.py module
#
# function : read weather data file (txt, csv, xlsx, xls files)
#               input files should have 'timestamp' column  (%Y-%m-%d %H:%m
# input   : filename
# output : dataframe object
#
# auther : K. H. Moon (2018)

import datetime
import pandas as pd

class Weather():
    
    def __init__(self, fname, start, end):
        
        # discrimiate file type (.txt, .csv, .xlsx, .xls)
        # proccess for read, set index to datetime, and missing data adjustment
        ext = fname.split('.')[1]
        if ext == 'csv':
            df = pd.read_csv(fname, parse_dates=['timestamp']).set_index('timestamp').sort_index()
        elif ext == 'txt':
            df = pd.read_csv(fname, sep='\t', parse_dates=['timestamp']).set_index('timestamp').sort_index()
        elif ext == 'xlsx' or ext == 'xls':
            df = pd.read_excel(fname).set_index('timestamp').sort_index()
        else:
            raise Exception('This file cannot be readed. Change file type.')
        
        
        df = df[~df.index.duplicated()]
        df = df.resample('H').interpolate()

        start = pd.to_datetime(start, format='%Y-%m-%d %H:%M')
        end = pd.to_datetime(end, format='%Y-%m-%d %H:%M') + datetime.timedelta(hours=23)
        df = df[start:end]
        start_d = pd.to_datetime(start, format='%Y-%m-%d')                    ###  2021 calculate DAP and add column ####
        idx = pd.to_datetime(df.index, format='%Y-%m-%d')                     ###  2021 calculate DAP and add column ####
        df['dap'] = (idx - start_d).days                                      ###  2021 calculate DAP and add culumn ####        
            
        # set weather data as DataFrame type
        self.data = df