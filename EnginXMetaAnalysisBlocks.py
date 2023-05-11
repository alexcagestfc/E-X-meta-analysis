# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

from mantid.simpleapi import *
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 200)
from scipy.optimize import curve_fit 
import math


"""
introduction statement:

These scripts were written for a "historical analysis project"
We suspected that Engin-X had been depreciating over time and wanted to consider the possibility that it had been producing less intense neutron beams

-Purpose and goals
My objective was to analyse various experiments that occured since Engin-X's debut. considering the spectra from the two monitors along the beam line.
control for all extraneous varables.
compare spectra to one another and conclude as to if they're registering less flux as time goes on.


-methods:
I've formatted this library into 4 sections:

    -1- acquiring and selecting experiments:
    
My first step is to select the experiments from which I'll load data
I've been provided a doccument (SUMMARY.txt) listing all of engin-x's experiments with info about the experiment, the scripts in this section do the following:
open the .txt file,
identify experiments that will be useful,
read the information,
make that information accessible for later processing

This was acheved using the pandas data analysis library, The info would be tabulated in a data frame
The experiments I would choose would be calibration experiments using V-Nb samples, 
these data frame I would use would be named VNbDf as in "V-Nb Dataframe"

    -2- sequenctial loading experiment data:
    
from the info I've now aquired from SUMMARY.txt, I will need to load the relevant data from each experiment.
I'm using Mantid IDE and api which provides the script for loading experiments. 
mantid.simpleapi.load() does require a file path. This must be done by File ->manage user directory-> navigate to folder


I trialed multiple methods for efficiently loading multiple experiments, used together they can:
Sequentially access experiments listed in VNbDF
derive each experiment's file name from VNbDf info
succesivly use each experiment's file name in mantid.simpleapi.load()

    -3- read loaded data:
    
The experiments would be loaded as workspaces in Mantid's IDE. these were not ideal as I need the info in arrays that can subjected to mathematical processes
This section accomplished created these arrays like so:
succesivly saves the data of each experiment as a CSV.
then uses the pandas function to access the very same CSV and forms a dataframe for each experiments data
each experiments data frame is appended in a list named experimental list



    -4- processing
only now do I have all the data with which to analyse.
the processing section is a list of scripts I use to analyse the data the majority of the lines of code are in section 4. and the individual scripts vary a lot in function
"""
"""
SECTION 1: acquire information
"""
def identifyexperiments():
    """
    e.g. variable = identifyexperiments()
    
    variable is set as a pandas dataframe:
    columns: ['Number','Object','date','time','secs','experiment']
    these will be read directly from SUMMARY.TXT. Many aspects will require certain processes before being usable.
        one example of this sort of processing is 'Number'. which takes the format ('ENG' +str(opqrs)) where lmnop are integers different integers
        'Number' is queried to find the relevant filename. the file name will have the format ('ENGINX'+str(lmn)+str(opqrs)) qrs is identified in later scripts
        
    It is a requirement that the program is able to access SUMMARY.TXT. 
       
       

    """
    col_names = ['Number','Object','date','time','secs','experiment']
    col_width = [(0,28),(28,52),(52,63),(63,72),(72,80),(80,88)]
    data = pd.read_fwf('C:\SUMMARY.TXT',names = col_names, colspecs = col_width, na_filter = False, encoding = 'cp1252')
    df = pd.DataFrame(data)
    VNbDf = df[df['Object'].str.contains('vnb|VNb|v-nb|V-Nb|vanadium|Vanadium|VdNb|V_Nb|Nb_V')]
    return VNbDf

def getprefix(VNbDf):
    """
    e.g. variable2  = getprefix(variable1)
    variable1 will be required to have a column called 'date'. this condition will be met by dataframes set by identifyexperiments()
    variable2 will be identical to variable1 except for having an extra column called 'prefixes'. 
    each entry is a string formed of three integers, which I refer to in docstrings as lmn
    the afformentioned prefix so that the file name can be formed from numbers in the dataframe
    
    
    numbers lmn will be one of: 001,002,003,999
    the former three are genuine prefiexes that are used file accessing
    999 marks it as using the pre 2008 file name format. My loading script will load it without issue.
    """
    date_aframe = VNbDf['date'] #essentially makes a pandas series 

    rows = len(date_aframe) #find number of experiments
    i=1
    date1 = datetime(2008,1,1)
    date2 = datetime(2013,2,27)
    date3 = datetime(2018,11,24)
    n=[]

    for i in range(rows):
        date = datetime.strptime(date_aframe.iloc[i],'%d-%b-%Y')
        if date < date1:
            prefix='999'
            n.append(prefix)
            #marks out  pre 2008 experiments
        if date1<date and date<date2:
            prefix='001'
            n.append(prefix)
            #2008-2012 marked 001 etc...
        if date2<date and date<date3:
            prefix='002'
            n.append(prefix)
        if date3<date:
            prefix='003'
            n.append(prefix)
        i=i+1
    VNbDf['prefixes'] = n
    return VNbDf
    
    
def ObjectSelections(criterion,VNbDf):
    """
    removes all lines that don't meet requirements
    
    MakeFurtherSelections==false actually causes an error.
        I will fix this soon aswell. until then simply comment the line out when not in use.
        
    """
    
    VNbDf = VNbDf[VNbDf['Object'].str.contains(criterion)]
    return(VNbDf)
    
def AmpHourSort(VNbDf):
    """
    data frame of experiment info is changed to be in order of most counts. 
    it would otherwise be chronological
    
    """
    pd.to_numeric(VNbDf['secs'])
    VNbDf = VNbDf.astype({"secs": float})
    VNbDf = VNbDf.sort_values(['secs'],ascending=False)
    return(VNbDf)
    
def newindexing(VNbDf):
    """
    creates a column that labels the order of rows in current state. 
    I used this to number the rows in choronological order and keep track of this order after sorting for amphour
    
    """
    b = (len(VNbDf))
    index2 = np.arange(0,b,1)
    VNbDf['index2'] = index2
    return(VNbDf)   
    
    
    
    

    
    
    
"""
SECTION 2 load experiments
"""
def loadfromexn(expnumber,min,max,name,explist):
    """
    Loader function 1: loads experiment from 5 digit number 
    
    inputs are  expnumber: a 5 digit code associated with the experiment in question
                min: the first spectrum from the experiment that'll be used
                max: the final spectrum
                name: a string for which will be used for the workspace name      
           
    this number will be the very same 5 digit number the dataframe,'Number' characters[3:8] (in python's indexing) 
    
    
    function is very simple given that mantid api was already able to take the correct file name and load based on that information
    loadfromexn constructs the file name based on info provided by summary.txt   
    """
    
    exp = ('ENG'+str(expnumber))
    #print(exp)
    row1 = explist[explist['Number'].str.contains(exp)] #finds experiments information from the data drame
    pref=row1['prefixes']# experiment file will have 8 digits, pref is the first three that aren't included in summary.txt   
    filename = ('ENGINX'+pref+str(expnumber)+'.RAW')#file name has this structure e.g. 'ENG00113579.RAW' where expnumber is '13579'
    filename = (str(filename))#file name converted from object to string
    filename= filename[10:28] #extra information removed
    pref = (str(pref)) #convert to string
    pref= pref[10:13]# extra information removed
    #if experiment is pre2008, it won't load
    if pref !=999:
        
        filename = ('ENGINX'+pref+str(expnumber)+'.RAW')
        filename = (str(filename))
        filename= filename#[0:28]
        #reasserts filename
        Load(filename, OutputWorkspace=name,SpectrumMin=min, SpectrumMax=max)



    
def loadN(VNbDf,start,N1,N2):
    """
    This function will analyse rows of the dataframe of experimental info

	for each row it will figure out the coresponding filename and use Mantid API's load function to load it.
	It also returns the number of experiments loaded. This allows for the script to easily be used to load multiple data frames of info without over 	writing.  

    """
    a=start
    b=(len(VNbDf))
    b=a+b
    min=N1
    max=N2
    c=0
    while a !=b:
        
        print(a)
        try:
            row1=VNbDf.iloc[c]
            expname=row1['Number']
            expnumber = (expname[3:8])
            pref=row1['prefixes']
            filename = ('ENGINX'+pref+str(expnumber)+'.RAW')
            filename = (str(filename))
            filename= filename[10:28]
            try:
                    
                filename = ('ENGINX'+pref+str(expnumber)+'.RAW')
                filename = (str(filename))
                filename= filename#[0:28]
                loadfromexn(expnumber,min,max,str(a),VNbDf)
            except:
                try:
                        
                    filename = ('ENG'+str(expnumber)+'.RAW')
                    filename = (str(filename))
                    filename= filename#[0:28]
                    Load(filename, OutputWorkspace=str(a),SpectrumMin=min, SpectrumMax=max)
                    
                except:
                    print('pre 2008 load failed')
                    pass
        except:
            print('load failed')
            pass
        c=c+1
        a=a+1
    return(a)
    
"""
SECTION 3: read experimental data


	having load the experimental data. It is not instantly in a state to be processed and analysed.

	these scripts will save data as an CSV and then read the contents of the file to produce a pandas data frame.

	this data frame is much easier to analyse
"""

def readexperiment(wkspce,exprimentaldict):
    """
    reads through data in workspace , acquires the data for processing and graphing. 
    each workspace will produce a data frame of information
    
    reading the information makes use of and edits info in a CSV file called 'atmpt1'
    
    others who may use this this library will have to edit the file paths used by SaveAscii and read_csv
        
        
    
    """
    o=str(wkspce)
    SaveAscii(InputWorkspace=o, Filename='C:/Users/yth69564/atmpt1.CSV',Separator="CSV")
    column_names='X' ,'Y', 'E'
    col_width = [(0,5),(6,16),(17,215)]
    currentdat = pd.read_csv('C:/Users/yth69564/atmpt1.CSV',names = column_names, na_filter = False, encoding = 'cp1252')
    return(currentdat)
    currentdat=[]
    
def ReadAll(VNbDf):
    """
    Sequentially use readexperiment to aquire data from each experiment described in the info dataframe, its use is contingent on priored data loading using section 2 scripts
    """
    experimentallist=[]
    a=0
    b=len(VNbDf)-1
    while a!=b:
        
        dataset=readexperiment(a,experimentallist)
        experimentallist.append(dataset)
        a=a+1
    return(experimentallist)



"""
SECTION 4: processing and plotting

The purpose of these scripts are to find the areas under each spectrum and then plot the areas(the integral) against the dates of experiment.

Note that many Section 4 functions make referance to X and Y values 
these X Y values may be the respective: TOF vals and counts per second of individual spectra 

X and Y at times also might refer to dates of experiements and the integrals of each experimental spectrum. 

I will endevour to make it clear which the script refers to.

By the way I use area and integral interchangeably  
"""
def rebinall(VNbDf,width):
    """
    Rebins all spectra 
    width should be a string that contains a number
	I invaribly use a bin width of 1.
    """
    a=0
    b=len(VNbDf)-1
    while a!=b:
        stra=str(a)
        Rebin(InputWorkspace=stra, OutputWorkspace=stra, Params=width)
        a=a+1


    
    

    
    
def dayssince_1_1_2000(x):
    """
    summarys date formatting is incompatiable with typical datetime library, I made this to operate on it instead
    """
    try:
        x1=x[0:2]
        count1= int(x1)
        #print(count1)
        X2=x[3:6]

        if X2=='JAN':
            count2=0
        if X2=='FEB':
            count2=31
        if X2=='MAR':
            count2=59
        if X2=='APR':
            count2=90
        if X2=='MAY':
            count2=120
        if X2=='JUN':
            count2=151
        if X2=='JUL':
            count2=181
        if X2=='AUG':
            count2=212
        if X2=='SEP':
            count2=243
        if X2=='OCT':
            count2=273
        if X2=='NOV':
            count2=304
        if X2=='DEC':
            count2=334
        
        
        if X2=='Jan':
            count2=0
        if X2=='Feb':
            count2=31
        if X2=='Mar':
            count2=59
        if X2=='Apr':
            count2=90
        if X2=='May':
            count2=120
        if X2=='Jun':
            count2=151
        if X2=='Jul':
            count2=181
        if X2=='Aug':
            count2=212
        if X2=='Sep':
            count2=243
        if X2=='Oct':
            count2=273
        if X2=='Nov':
            count2=304
        if X2=='Dec':
            count2=334
        
        x3=x[9:11]
        year=int(x3)
        count3=year*365
        xnew=(count1+count2+count3)
        return(xnew)
    except:
        print('analyse date string failed:'+ x)


    
def vars_Time_Sum(VNbDf,experimentallist):
    """
	

    produces list X and list Y

	Each element in X should be a date as formatted in the SUMMARY.txt
    X relates to how long the instrument has been running, though it simply uses and arbitraty datetime64
    
	Each element in Y is an integral of a spectrum
    Y is the sum of the Y norm column for each experiment
    
    typical method would use the datetime modules strptime command 
    however SUMMARY.txt's conventions do not allow this. I've made a module called dayssince_1_1_2000 That figures it out.
	 
    
    """
    X=[]
    Y=[]
    VNbDf['X']=0
    VNbDf['Y']=0
    a=0
    b=len(VNbDf)-1
    c=0
    while a!=b:
        
        row1=VNbDf.iloc[a]
        if row1['use']==True:
            
            datestring = row1['date']
            n=dayssince_1_1_2000(datestring)
            
            subject = experimentallist[a]
            subject = subject.iloc[1: , :]
             
            Yvalue =subject['Y']
            Yvalue = Yv.iloc[1:]
            Yvalue = pd.to_numeric(Yvalue, downcast="float")
            m=Yvalue.sum()
            
            X.append(n)
            Y.append(m)
            c=c+1
        a=a+1
    return(X,Y)
    
    
    

    
    
    
def vars_simplefit(X,Y):
    """
	X,Y are lists of Dates and integrals
    numpy polyfit used for line of best fit
    
    note that plotting does not occur here, it provides Yfit that may be fitted in a separate line_buffering
    """
    
    fitperam1, fitperam2 = np.polyfit(X, Y, 1)
    Yfit =[]
    a=0
    b=len(X)
    while a != b:
        fitpoint_a= fitperam1*X[a]+fitperam2
        Yfit.append(fitpoint_a)
        a=a+1

    return(Yfit)
    
    
    
def sufficientamphour(VNbDf,n):
    """
	ignores experiments with low pulse number (note to self: might be redundent)
	"""
    a=0
    use=[]
    b=len(VNbDf)


    while a !=b:
        row  = VNbDf.iloc[a]
        s = float(row['secs'])
        
        if s>n:
            viable=True
            use.append(viable)
        if s<n:
            viable=False
            use.append(viable)
        a=a+1
    VNbDf['use']=use
    return(VNbDf)
    
    
def conditonappend(X,Y,VNbDf):
    """
	new col in VNbDf with sum of each experiments spectrum 
	possibly redundant
	"""
    a=0
    c=0
    b=len(VNbDf)-1
    VNbDf['X']=0
    VNbDf['Y']=0
    while a!=b:
        row=VNbDf.iloc[a]
        if row['use']==True:
            
            row['X']=X[c]
            row['Y']=Y[c]
            VNbDf.iloc[a]=row
            c=c+1
        a=a+1
    return(VNbDf)
    
  
def RangeLimit(X,Y,g,h):
    """
    X,Y are Dates,Integrals
    designed to take a list of X and Y values that would be plotted and excluded the data points that are outside a certain min and max range e.g:

    RangeMin=10
    RangeMax=100
    X2,Y2=RangeLimit(Y,X,RangeMin,RangeMax)

     Y2 will now only contain only the elements of Y that above 10 and below 100
     X2 will be the elements of X that had the same index of the Y values that met the criteria 


    This function can be used to controll domain instead of range by simily running Y2,X2=RangeLimit(Y,X,domainMin,DomainMax)



    X and Y are equal length lists 
    g and h are floating points

    selects max (h) and min (g) values in list Y. Then returns X2 and Y2.
    Y2 is a list of elements in Y that are between the min and max vals
    X2 is the list of values appended from X that corespond to the values in Y2

    """
    a=0
    X2=[]
    Y2=[]
    b=len(Y)
    while a!=b:
        c=Y[a]
        if g<c and c<h:
            Y2.append(c)
            X2.append(X[a])
        
        a=a+1

    return(X2,Y2)
        
        
"""
one of my prototypes tried fitting a maxwell boltzmann (MWBM) dist to the spectrum and analysing that. 
It worked only for 2402 and didn't help me much. I might leave it in the final version if it seems useful to others.
"""
        
##processing mk2

def maxwell(X,a,b,c):
    """
    the MWBM function fitted
    X is the TOF of some particle 
    Y is a function of X and is the coresponding count rate
    """
    return a*(1*b*X**2)*np.exp(-1*c*X**2)

def fitmaxwell(points,index,experimentallist):
    """
    fits a designated spectrum to MWBM with a predetermined number of points
    X,Y are TOF values, Count rates of each TOF
    """
    plt.figure()
    subject = experimentallist[index]
    subject = subject.iloc[2: , :]

    Xdata = subject['X']
    Xd=np.linspace(0,100,points)


    Ydata =subject['Y']
    n=len(Ydata)/points
    n=round(n)
    Yd=Ydata[Ydata.index % n == 0]

    while len(Yd) >  len(Xd):
        Yd.drop(index=Yd.index[-1],axis=0,inplace=True)

    while len(Xd) > len(Yd):
        Xd.drop(index=Xd.index[-1],axis=0,inplace=True)

    Yd = pd.to_numeric(Yd, downcast="float")
    popt, pcov = curve_fit(maxwell, Xd, Yd, bounds=(0, [50., 10., 0.5]))
    # popt, pcov = curve_fit(maxwell, Xdata, Ydata)
    Yf= maxwell(Xd, *popt)
    print(popt)
    print(pcov)
    plt.plot(Xd,Yf , 'r-')
    plt.plot(Xd,Yd)
    return(Xd,Yd,Yf,popt,pcov)
    
def NVfitmaxwell(points,index,experimentallist):
    """
    MWBM fit but without plotting
    """
    subject = experimentallist[index]
    subject = subject.iloc[2: , :]

    Xdata = subject['X']
    Xd=np.linspace(0,100,points)


    Ydata =subject['Y']
    n=len(Ydata)/points
    n=round(n)
    Yd=Ydata[Ydata.index % n == 0]

    while len(Yd) >  len(Xd):
        Yd.drop(index=Yd.index[-1],axis=0,inplace=True)

    while len(Xd) > len(Yd):
        Xd=Xd[:-1]

    Yd = pd.to_numeric(Yd, downcast="float")
    popt, pcov = curve_fit(maxwell, Xd, Yd, bounds=(0, [50., 10., 0.5]))
    # popt, pcov = curve_fit(maxwell, Xdata, Ydata)
    Yf= maxwell(Xd, *popt)
    return(Xd,Yd,Yf,popt,pcov)

    
def XAxisCorrection(Xdata,experimentallist,index):
    """
    I use a number of functions that will sample integral data points equidistant on the X axis (TOF) from a given spectrum
    This function will asign the sampled data the correct X value. (instead of simply designating e.g. 1-100 step: 1)


    e.g.
    Xdata= np array(1,2,3...100)
    Xdata2= (5000,5001 ...79999)
    """
        
    subject = experimentallist[index]
    subject = subject.iloc[2: , :]
    Xinfo = subject['X']
    factor=(len(Xinfo))/100
    Xdata2=Xdata*factor
    return(Xdata2)

    
def Framed(Xinfo,Yinfo,experimentallist,VNbDf,failiures):
    """
    establishes a data frame of results from each spectra's analysis
    note that X and Y are dates and integrals 



    using inputs:
    Xinfo: which will be designated as the dates of each experiment
    Yinfo: the Integral of each experiments produced spectrum
    VNbDf: which contains the charge info 

    produces data frame listing locations in 3D space to plot points:
    the info in each axis being
    Integral of experiment
    date of experiment
    associated charge of experiment

    Note to self: I dont recall why failiures and experimentalist are inputs remove in next version

    """
    allcharge = VNbDf['secs']   
    a=0
    b=len(allcharge)
    charge=[]

    while a!=b:
        if a not in failiures:
            target = allcharge.iloc[a]
            
            charge.append(target)
        a=a+1
        
    charge=charge[0:(len(Xinfo))]
    Info3D = list(zip(Xinfo,charge,Yinfo))
    InfoFrame = pd.DataFrame(Info3D, columns=['Days','Charge','I'])
    return InfoFrame

    
    
def Framed2(Xinfo,Yinfo,experimentallist,VNbDf,failiures):
    """
    Identical to framed1 except it will store Integral/Charge instead of Integral
    """
    allcharge = VNbDf['secs']   
    a=0
    b=len(allcharge)
    charge=[]

    while a!=b:
        if a not in failiures:
            target = allcharge.iloc[a]
            
            charge.append(target)
        a=a+1
        
    charge=charge[0:(len(Xinfo))]
    Info3D = list(zip(Xinfo,charge,Yinfo))
    InfoFrame = pd.DataFrame(Info3D, columns=['Days','Charge','I/C'])
    return InfoFrame
    
def Selections(points,index,experimentallist):
    """
    fascilitates faster processing by selecting N equidistant data points from each spectrum and analysing those

    The X and Y data are the data in one spectrum: TOF vals and count rates respective
    The spectrum is slected based on its index number

    """

    plt.figure()
    subject = experimentallist[index]
    subject = subject.iloc[2: , :]

    Xdata = subject['X']
    Xd=np.linspace(0,100,points)
    Xd = pd.Series(Xd)

    Ydata =subject['Y']
    n=len(Ydata)/points
    n=round(n)
    Yd=Ydata[Ydata.index % n == 0]

    while len(Yd) >  len(Xd):
        Ydlen=len(Yd)
        TargetIndex=math.floor(Ydlen/2)
        Yd.drop(index=Yd.index[TargetIndex],axis=0,inplace=True)

    while len(Xd) > len(Yd):
        Xd.drop(index=Xd.index[-1],axis=0,inplace=True)

    Yd = pd.to_numeric(Yd, downcast="float")
    plt.plot(Xd,Yd)
    return(Xd,Yd)    
    
def SelectionsNV(points,index,experimentallist):
    """
    Selections but no plotting
    """
    subject = experimentallist[index]
    subject = subject.iloc[2: , :]

    Xdata = subject['X']
    Xd=np.linspace(0,100,points)
    Xd = pd.Series(Xd)

    Ydata =subject['Y']
    n=len(Ydata)/points
    n=round(n)
    Yd=Ydata[Ydata.index % n == 0]

    while len(Yd) >  len(Xd):
        Ydlen=len(Yd)
        TargetIndex=math.floor(Ydlen/2)
        Yd.drop(index=Yd.index[TargetIndex],axis=0,inplace=True)

    while len(Xd) > len(Yd):
        Xd.drop(index=Xd.index[-1],axis=0,inplace=True)

    Yd = pd.to_numeric(Yd, downcast="float")
    return(Xd,Yd) 

    
def SelectiveFullIntegralNV(Xdata,Yd):
    """
    for some spectrums data: X and Y referring to TOF and count rate

    area is calculated and returned
    """

    area=0
    a=0
    fin=len(Xdata)
    b=fin
    c=0


    Yranged=[]
    while a<b:
        yi=Yd[a]
        area=area+yi
        a=a+1
        c=c+1
    return(area)


    

def SelectiveFullIntegralNV(Xdata,Yd):
    """
    calcualtes area without plotting
    """
    area=0
    a=0
    fin=len(Xdata)
    b=fin
    c=0
    dx=Xdata[1]-Xdata[0]

    Yranged=[]
    while a<b:
        yi=Yd[a]
        
        
        area=area+(yi*dx)
        a=a+1
        c=c+1
    return(area) 

def SelectiveFullAll(res,experimentallist,VNbDf):
    """
    possibly identical to SelectiveAll


    """
    a=0
    b=len(experimentallist)
    points=res
    arealist=[]
    dates=[]
    failiures=[]
    while a!=b:
        Xdata,Ydata =SelectionsNV(points,a,experimentallist)
        Xdata = pd.Series(Xdata)
        Xdata2=XAxisCorrection(Xdata,experimentallist,a)
        
        EnginXBDay='01-JUN-2003'
        Day1=dayssince_1_1_2000(EnginXBDay)
        Xdata2 = list(np.asarray(Xdata2) - Day1)
        Ydata=Ydata.to_numpy()
        area=SelectiveFullIntegralNV(Xdata2,Ydata)
        
        arealist.append(area)
        row1=VNbDf.iloc[a]
        datestring = row1['date']
        n=dayssince_1_1_2000(datestring)
        dates.append(n)
        a=a+1
    return dates, arealist
 
def Proportionise(Ydata,Xdata):
    """
    X,Y are dates, areas

    sets earliest experiment's area as one and all after as its proportional size in relation to first.

    used for making data more intuitive
    """

    StArg=np.argmin(Xdata)
    standard=Ydata[StArg]

    a=0
    b=len(Ydata)
    Ydata2=[]
    while a!=b:
        appendix=Ydata[a]/standard
        Ydata2.append(appendix)
        a=a+1
    Ydata2=np.array(Ydata2)
    return Ydata2
    
def d2y(D):
    """
    dates will be a number of days.
    This will calculate the year it occured for ease of reading

    """
    a=0
    b=len(D)
    Year=[]
    while a!=b:
        YearsSinceDebut=D[a]/365
        Y=YearsSinceDebut+2003.5
        Year.append(Y)
        a=a+1
    Year = np.array(Year)
    return(Year)

    
    
def chargediv(I,C):
    """
    I is a list of integrals
    C is a list of charges

    returns list of same size containging I/C
    """
    a=0
    b=len(I)
    IoC=[]
    while a!=b:
        entrant1=I[a]
        charge=C[a]
        result=entrant1/charge
        IoC.append(result)
        a=a+1
    return(IoC)
    
def BoundIntegral(experimentallist,index,S,E):
    """
    integrate the experimental data of one experiment between [lowest TOF bound]S and [hightest TOF bound]E

    """
    subject=experimentallist[index]
    data = subject.iloc[2: , :]

    Xd=data['X']
    Yd=data['Y']

    Xd = pd.to_numeric(Xd, downcast="float")
    Yd = pd.to_numeric(Yd, downcast="float")

    Xd=Xd.to_numpy()
    Yd=Yd.to_numpy()

    dx=Xd[42]-Xd[41]
    a=0
    b=len(Xd)
    area=0
    Yb=[]
    while a != b:
        if Xd[a]>S and Xd[a]<E:
            yi=Yd[a]
            Yb.append(yi)
            ai=yi*dx
            area=area+ai
        else:
            Yb.append(0)

        a=a+1
        

    return(area,Yb,Xd,Yd)
    
def boundplot(Xd,Yd,Yb):
    """
    code generates a figure showing the spectrum and a filled in area that boundIntegral analysed
    """
    fig, ax = plt.subplots()
    ax.plot(Xd, Yd, color='blue', alpha=1.0)
    ax.fill_between(Xd, 0, Yb, color='m', alpha=0.3)
    ax.set_xlabel('Time-of-flight ($\\mu s$)')
    ax.set_ylabel('Counts ($\\mu s$)$^{-1}$')
    ax.set_title('TOF energy spectra')
    
def BoundIntegralAllSpectra(VNbDf,experimentallist,S,E):
    """
    applied bound integral to each entry in experimentalist using the same bounds
    """
    index=0
    final =len(experimentallist)
    arealist=[]
    dates=[]
    while index <final:
        Area,Yb,Xd,Yb = BoundIntegral(experimentallist,index,S,E)
        print('index: '+str(index)+' integration successful')
        arealist.append(Area)
        
        row1=VNbDf.iloc[index]
        datestring = row1['date']
        n=dayssince_1_1_2000(datestring)
        dates.append(n)
        index=index+1
    return dates, arealist
    
    
    
    
    def chargediv(I,C):
    a=0
    b=len(I)
    IoC=[]
    while a!=b:
        entrant1=I[a]
        charge=C[a]
        result=entrant1/charge
        IoC.append(result)
        a=a+1
    return(IoC)