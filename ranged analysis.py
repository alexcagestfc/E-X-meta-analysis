from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
pd.set_option('display.max_rows', 300)
pd.set_option('display.max_columns', 15)

from scipy.optimize import curve_fit 

sys.path.append('C:/Users/yth69564/OneDrive - Science and Technology Facilities Council/Documents/Alexinfo/working_scripts')
from EnginXMetaAnalysisBlocks import *
"""
SECTION1: identification

This script is very similar to main analysis script, This is for preforming analyses over particular bounds for I/C

"""

VNbDf = identifyexperiments()
VNbDf = getprefix(VNbDf)
criterion = '4 x|4x|x4|4_|4h|4c'

#VNbDf = ObjectSelections(criterion,VNbDf)

newindexing(VNbDf) #I'm about to sort the dataframe inorder of charge, newindex will note chronology 
VNbDf =AmpHourSort(VNbDf)
minAhr=300
sufficientamphour(VNbDf,minAhr)




"""
SECTION2: load experiments

Mantid has a maxwell for loading an experiment, my library uses it to sucessivly to load all experiments in a data frame made by identifyexperiments()
to load only useable experiments, VNbDf2= VNbDf.loc[VNbDf['use']==True]
"""

spaces=0
detector=2402
VNbDf= VNbDf.loc[VNbDf['use']==True]
spaces=loadN(VNbDf,spaces,detector,detector)



"""
following hasn't been tidied yet




SECTION3: reading data sets
extracts
"""
rebinall(VNbDf,'1')
experimentallist=ReadAll(VNbDf)
"""
analyse:
"""

S=0
E=10000
index=1









dates, arealist=BoundIntegralAllSpectra(VNbDf, experimentallist,S,E)

EnginXBDay='01-JUN-2003'
Day1=dayssince_1_1_2000(EnginXBDay)
Days_2 = list(np.asarray(dates) - Day1)


InfoFrame = Framed(Days_2,arealist,experimentallist,VNbDf,[])
I = InfoFrame['I'].to_numpy()
C = InfoFrame['Charge'].to_numpy()
Xd = InfoFrame['Days'].to_numpy()

Y2=chargediv(I,C)


InfoFrame2 = Framed2(Xd,Y2,experimentallist,VNbDf,[])

IoC = InfoFrame2['I/C'].to_numpy()
S3Z = InfoFrame2['Charge'].to_numpy()
S3X = InfoFrame2['Days'].to_numpy()
 

rslt_df = InfoFrame2[InfoFrame2['Days'] > 1900]
 

NoP = rslt_df['I/C'].to_numpy()
uAh = rslt_df['Charge'].to_numpy()
date = rslt_df['Days'].to_numpy()

plt.figure()
date=d2y(date)
NoP = Proportionise(NoP,date)
NoP=NoP*100
date,NoP = RangeLimit(date,NoP,0,200)
m, c = np.polyfit(date, NoP, 1)
P2Yfit = vars_simplefit(date,NoP)

plt.plot(date,P2Yfit,'r')
plt.scatter(date,NoP)
print(m)
m=round(m, 4)
c=round(c)
plt.title('neutrons per uAh vs date, TOF range: '+str(S)+' to '+str(E)+' ,fit: y='+str(m)+'x +'+str(c))
plt.xlabel('Year')
plt.ylabel(' monitor count per unit charge (as percentage of earliest)')


