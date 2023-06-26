"""
This program will use scripts from EnginXMetaAnalysisBlocks  for the purpose of beamline performance analysis

the results output will be six plots as well as statistical info which will be printed
statistical info will include covarance info, correlation info and standard deviations

of the plots, the sixth plot will present data that has been manipulated for the purpose of legibility

plot explanations:
figure 1: un edited counts vs number of days since debut of instrument
figure 2: 3D plot of counts vs number of days vs associated charge

figure 3: 3D plot of counts/charge vs number of days vs charge
(atempt to control for the relation between counts and charge)

figure 4: same as figure 3, now with line of best fit as made by principal component analysis
figure 5: counts/chage vs time since debut with line of best fitting_algorithm

figure 6: numbers edited for legibility counts per unit charge over years
"""


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

I'm provided summary.txt a list of experiment conducted by engin-x

this phase identifies from summary which experiments should be used

VNbDf stands for Vanadium-NioBium DataFrame, this process uses the callibration experiments of a VNb sample

identifyexperiments() reads through the actual CSV 
getprefix(VNbDf) appends information pertainant to loading the data (which is stored separatly to the info

newindexing(VNbDf) applies a new column that numbers within Df and not their number within the summary
AmpHourSort(VNbDf) organises by largest detection numbers
sufficientamphour labels those that have sufficient detections

"""

VNbDf = identifyexperiments() #read through summary.txt
VNbDf = getprefix(VNbDf) #assignsnumbers for locating data file
criterion = '4 x|4x|x4|4_|4h|4c' #only use experiments with certain perameters, I don't use this atm
#"|" is basically an "or" statement 

#VNbDf = ObjectSelections(criterion,VNbDf)

newindexing(VNbDf) #I'm about to sort the dataframe inorder of charge, newindex will note chronology 
VNbDf =AmpHourSort(VNbDf) # sort max to min charge values
minAhr=300 # min charge value used in analysis
sufficientamphour(VNbDf,minAhr) # marks those less than min charge value for removal




"""
SECTION2: load experiments

Mantid has a maxwell for loading an experiment, my library uses it to sucessivly to load all experiments in a data frame made by identifyexperiments()
to load only useable experiments, VNbDf2= VNbDf.loc[VNbDf['use']==True]
"""

spaces=0
detector=2401
VNbDf= VNbDf.loc[VNbDf['use']==True] #remove those that don't meet criteria
spaces=loadN(VNbDf,spaces,detector,detector) # loads data as workspace



"""
SECTION3: reading data sets
extracts data

readall works by saving data as a csv and then reading as csv. 
you'll need to edit the script's decorator: "readexperiment"
so that it saves to a file on your computer
it is the first script in section 3. 
Just create an empty note pad doc, copy the file path, put .csv instead of .txt


"""
rebinall(VNbDf,'1')#set bin width. the width itself doesn't matter too much, only that they're equal
experimentallist=ReadAll(VNbDf)#acquire data as pandas data frame. 

"""
Section4 analysis:
 
 
pt1:
integrates all spectra, notes their dates so they can be plottted wrt time
"""
points=401
days, arealist = SelectiveFullAll(points,experimentallist,VNbDf)

Days_2,rangedsums=RangeLimit(days,arealist,0*10**3,1*10**12)#remove anomalously large or small results


EnginXBDay='01-JUN-2003'
Day1=dayssince_1_1_2000(EnginXBDay)
Days_2 = list(np.asarray(Days_2) - Day1)# reformat date interms of time since Engin-X debut

"""
plotting Integral against T (days)
"""
plt.figure()
plt.xlabel('time (days) since debut')
plt.ylabel('integrated monitor trace (normalised)')

m, c = np.polyfit(Days_2, rangedsums, 1)
P2Yfit = vars_simplefit(Days_2,rangedsums)

plt.plot(Days_2,P2Yfit,'r')
plt.scatter(Days_2,rangedsums)

m=round(m, 2)
c=round(c)
plt.title('integral vs time ,fit: y='+str(m)+'x +'+str(c))

"""
pt2:
establish data frame for storage of processed data
"""

InfoFrame = Framed(Days_2,rangedsums,experimentallist,VNbDf,[])



IntWOcharge = InfoFrame['I'].to_numpy()
S3Y = InfoFrame['Charge'].to_numpy()
S2X = InfoFrame['Days'].to_numpy()

S3X=d2y(S2X)
"""
3D plot I vs T vs Charge
"""
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plt.title('Integral vs Charge vs time ')

ax.scatter(np.round(S3X,0), S3Y, IntWOcharge)

ax.set_xlabel('date (year)')
ax.set_zlabel('Neutron count')
ax.set_ylabel('Charge')

"""
pt3: charge adjust

the process will now divide integral by respectivr ccharge to assertain intensity ber unit charge
"""
IntByCharge=rangedsums

IntByCharge=chargediv(rangedsums,S3Y)
"""
pt4: intesity vs charge analysis
info is stored for ease of access
"""

Days,IntByCharge = RangeLimit(Days_2,IntByCharge,0*10*3,1.0*10**12)

InfoFrame2 = Framed2(Days,IntByCharge,experimentallist,VNbDf,[])

IoC = InfoFrame2['I/C'].to_numpy()
S3Z = InfoFrame2['Charge'].to_numpy()
S3X = InfoFrame2['Days'].to_numpy()


rslt_df = InfoFrame2[InfoFrame2['Days'] > 2100]


S3Y = rslt_df['I/C'].to_numpy()
S3Z = rslt_df['Charge'].to_numpy()
S2X = rslt_df['Days'].to_numpy()
S3X=d2y(S2X)


"""
3D plotting:
"""
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plt.title(str(detector)+': I/C vs Charge vs time, after change, pc1')
ax.scatter(S3X, S3Y, S3Z)

ax.set_xlabel('time (days)')
ax.set_ylabel('Integral per unit charge')
ax.set_zlabel('Charge')


"""
3D line fit:

easiest method of fitting for points with three perameters is pca
"""


datalist=np.concatenate((S3X[:,np.newaxis],S3Y[:,np.newaxis],S3Z[:,np.newaxis]),axis=1)


datamean = datalist.mean(axis=0)
uu, dd, vv = np.linalg.svd(datalist - datamean)

p1=+1000000/2
p2=-1*1000000/2
#these two points work well when they have the same magnitude the appropriate mag varies
linepts = vv[0] * np.mgrid[p1:p2:2j][:, np.newaxis]


linepts += datamean



import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d

ax = m3d.Axes3D(plt.figure())
ax.scatter3D(*datalist.T)
ax.plot3D(*linepts.T)


ax.set_xlabel('time (years)')
ax.set_ylabel('Integral per unit charge')
ax.set_zlabel('Charge')
ax.set_title('1srtprincipal compent')

covarxy=np.cov(S3X,S3Y)
covaryz=np.cov(S3Y,S3Z)
covarxz=np.cov(S3X,S3Z)
print('covrance matrix date and integral/charge:')
print(np.round(covarxy,1))
print('covrance matrix date and charge:')
print(np.round(covarxz,1))
print('covrance matrix integral/charge and charge:')
print(np.round(covaryz,1))

svx=np.std(S3X)
svy=np.std(S3Y)
svz=np.std(S3Z)



XPC=linepts.T[0]
YPC=linepts.T[1]
ZPC=linepts.T[2]

Xavr=np.mean(XPC)
Yavr=np.mean(YPC)
Zavr=np.mean(ZPC)

XRS=svx/Xavr
YRS=svy/Yavr
ZRS=svz/Zavr
XRSpc=XRS*100
YRSpc=YRS*100
ZRSpc=ZRS*100

#print('relative standard deviation in date, calculated as (sigma)/mean '+str(round(XRSpc))+'%')
print('relative standard deviation in integral/charge '+str(round(YRSpc))+'%')
#print('relative standard deviation in charge '+str(round(ZRSpc))+'%')




xycor=np.corrcoef(S3X,S3Y)
xzcor=np.corrcoef(S3X,S3Z)
yzcor=np.corrcoef(S3Y,S3Z)
xycor=xycor[1]
xycor=xycor[0]

xzcor=xzcor[1]
xzcor=xzcor[0]

yzcor=yzcor[1]
yzcor=yzcor[0]
print('correlation between date and integral/charge '+str(round(xycor,2)))
print('correlation between date and charge '+str(round(xzcor,2)))
print('correlation between integral/charge and charge is '+str(round(yzcor,2)))




plt.figure()
plt.xlabel('year')
plt.ylabel('integrated monitor count per unit charge')
S3X2=d2y(S3X)

#S3X,S3Y=RangeLimit(S3X,S3Y,5*10**2,1.25*10**3)
plt.scatter(S3X2,S3Y)

m, c = np.polyfit(S3X2, S3Y, 1)

S3Yfit = vars_simplefit(S3X2,S3Y)
plt.plot(S3X2,S3Yfit,'r')
m=round(m, 7)
c=round(c)
plt.title('monitor 2402: Neutrons per microAmpHr vs date,'+',fit: y='+str(m)+'x +'+str(c))

"""
graph proportially

it doesn't change anythig but the data is more intuitive this way
"""
S2X = rslt_df['Days'].to_numpy()
S3X2=d2y(S2X)

Ydata2 = Proportionise(S3Y,S3X2)
Ydata2=Ydata2*100

S3X2,Ydata2=RangeLimit(S3X2,Ydata2,50,150)
plt.figure()

plt.scatter(S3X2,Ydata2)

m, c = np.polyfit(S3X2, Ydata2, 1)

S3Yfit = vars_simplefit(S3X2,Ydata2)
plt.plot(S3X2,S3Yfit,'r')
m=round(m,1)
c=round(c)
plt.xlabel('Year')
plt.ylabel('monitor count per pulse (as percentage of earliest)')
plt.title('Neutrons per pulse vs date'+', gradient: '+str(m)+' percent per year')

"""
"""

