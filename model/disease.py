"""
   Disease class : calculation of disease loss by rate
"""

import numpy as np

# parameters  
# rate = beta function with, Tx, Tn, To
# disease severity(ds) = 1 - exp(-((rate_sum/33.0)**5.0)), 0 - 1
#Tx = 40     # maximum temp for beta function
#To = 24     # optimum temp for beta function
#Tn = 18     # minimum temp for beta function
#A  = 33
#B  = 5.0
A = 2.0
B = 4.0

# time calibration
time   = 'hour'
conv   = 1 / 24  if time == 'hour' else 1

class Disease(object):
    
    def __init__(self):
        self.ds = 0.0                    # Severity (0-1)
        self.rsum = 0.0                  # rate sum

    def calc(self, Tair, wet):      
        # calculation
        if wet and Tair > 0:
            # use modified lognormal function
            rate  = np.exp(-1*(np.log(Tair/100))**4) * conv
            # print(rate)
        else :
            rate  = 0

        # use beta function
#        if (Tair > Tn and Tair < Tx) and wet :
#            ax   = (To-Tn)/(Tx-To)
#            bx   = (Tx-Tair)/(Tx-To)
#            cx   = (Tair-Tn)/(To-Tn)
#            rate = 1 * bx * (cx ** ax) * conv
#        else:
#            rate = 0

        self.rsum += rate
        dx   = -1 * (self.rsum / A)**B      
        ds = 1 - np.exp(dx)
        self.ds = ds
