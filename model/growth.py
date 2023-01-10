'''
  Growth class for partitioning assimilates to shoot, root
'''

import numpy as np

### Parameters for Chinese cabbage
Tbase  = 0.0           # base temperature
dm     = 0.6           # max root depth
dg     = 0.012         # root elongation rate (m/day)
hgt    = 0.4           # mean plant height (m)
width  = 0.20          # mean leaf width (m)
laicr  = 3             # lai critical for leaf death

### Constants for maintenace respiration at 25C for various plant parts
kgl  = 0.03            # for green leaves (gCH2O/gDM/day)
kr   = 0.015           # for root (gCH2O/gDM/day)
ko   = 0.020           # for rep. org

### Glucose requirement for synthesis of various plant parts 
ggl  = 1.463           # for green leaves (gCH2O/gDM)
gr   = 1.444           # for root (gCH2O/gDM)
go   = 1.463

# ## Table for partitioning
istage = [0.00, 0.50, 0.95, 1.00, 1.10, 2.00]         # stage
ifgl   = [0.92, 0.93, 0.94, 0.70, 0.00, 0.00]         # green leaf ratio
ifr    = [0.08, 0.07, 0.05, 0.00, 0.00, 0.00]         # root ratio
ifo    = [0.00, 0.00, 0.01, 0.30, 1.00, 1.00]

## ratio to DM / CH2O
cd     = 1.125         # factor convert mass to CH2O weight, = mass * 45%C / C_MW * CH2O_MW

## ratio of fresh weight to dry weight of leaf
FDR    = 16.76

## time conversion
time  = 'hour'
conv  = 1/24 if time == 'hour' else 1

class Growth(object):
    
    def __init__(self):
        self.wgl = 1.0            # green leaves weight (DW g m-2)
        self.wr = 0.2             # root weight (DW g m-2)
        self.wo = 0.0             # reproductive organ (stem + flower + seed)
        self.maint = 0.0          # maintenace respiration
        self.rootd  = 0.0         # root depth (m)
        
        
    def growCalc(self, Ta, assim, dvs, RDT=1.0):
#        import pdb; pdb.set_trace()
        corr = 1 / 24
        assim_c = assim * RDT     # reduced by water stress

        # print(assim_c)
    
        ## root depth by elongation (m)
        rgrow = min(dm, dg * RDT ) * corr
        self.rootd += rgrow
        self.assim_c = assim_c  # 파라미터 추출용
        
        ## leaf death rate by old age(dvs > 1.0)
        ddage = 0.0
        if (dvs > 0.9 and dvs < 1.9):
            ddage = 0.05 / (2.0 - dvs)
        elif dvs >= 1.9:
            ddage = 0.05 / 0.1
        else:
            ddage = 0.0
        # print(ddage)
        ddage = ddage * conv
        self.ddage = ddage # 파라미터 추출용
        # print(ddage)


        ## calculaltion green leaves, death laeves, stem, root, storage (g CH2O)
        wgl = self.wgl
        # print(wgl)
        wr = self.wr
        wo = self.wo
        
        #### calculation of maintenances respiration (g CH2O)
        mgl = wgl * kgl
        # print(wgl)
        mr  = wr * kr
        mo  = wo * ko
        RM = mgl + mo
        # print(RM)
        tempRM = RM * 2 **((Ta - 20)/10)  # temperature   from Teh
        # print(tempRM)
        tempRM = tempRM * corr            # from per day to per hour 
        RMpr = min(tempRM, assim_c)
        # print(RMpr)
        self.mgl = mgl # 파라미터 추출용
        self.mo = mo # 파라미터 추출용
        self.RM = RM # 파라미터 추출용
        self.tempRM = tempRM # 파라미터 추출용
        self.RMpr = RMpr # 파라미터 추출용
       
        
        ### calculation of growth respiration (g CH2O)
        fgl = np.interp(dvs, istage, ifgl)
        # print(fgl)
        fr  = np.interp(dvs, istage, ifr)
        # print(fr)
        fo  = np.interp(dvs, istage, ifo)
        # print(fo)

        self.fgl = fgl # 파라미터 추출용
        self.fr = fr # 파라미터 추출용
        self.fo = fo # 파라미터 추출용

        # nomalization of all partitioning
        fggl = fgl * ggl
        fgr  = fr  * gr
        fgo  = fo  * go
        GT = fggl + fgr + fgo
        # print(GT)

        self.fggl = fggl # 파라미터 추출용
        self.fgr = fgr # 파라미터 추출용
        self.fgo = fgo # 파라미터 추출용
        self.GT = GT # 파라미터 추출용

        available = (assim_c - RMpr) / GT      # available assim per hour (g DW)
        # print(available)

        available *= 0.84

        self.available = available # 파라미터 추출용
        # print(self.available)

        gr_gl = fgl * available     # for green leaves DW
        # print(gr_gl)
        gr_r  = fr  * available     # for root DW
        gr_o  = fo  * available     # for rep. organ
        gr_dl = ggl * ddage         # death leaf weight by aging

        self.gr_gl = gr_gl # 파라미터 추출용
        self.gr_r = gr_r # 파라미터 추출용
        self.gr_o = gr_o # 파라미터 추출용
        self.gr_dl = gr_dl # 파라미터 추출용


        self.wgl += (gr_gl - gr_dl)
        self.wr  += gr_r
        self.wo  += gr_o
        self.maint = assim_c - available 
        # print(self.wgl)