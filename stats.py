from numpy.linalg import inv
import numpy as np
import pykov
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
from itertools import cycle
import matplotlib
from matplotlib.pyplot import *
import brewer2mpl
import seaborn as sns
from scipy.stats import ks_2samp
from scipy.stats import mode
import statistics as st
import pandas as pd
import seaborn as sns



df = pd.read_csv('MCresults.csv')
data=[]
for i in range(3,10):
    nsfData=df.loc[(df['N'] == i) & (df['PType'] =='nsf')]
    nonsfData=df.loc[(df['N'] == i) & (df['PType'] =='no_nsf')]
    pNSF=nsfData['Probability']
    pNoNSF=nonsfData['Probability']
    
    ks=ks_2samp(pNSF,pNoNSF)
    
    modeNSF=mode(pNSF)
    maxNSF=max(pNSF)
    minNSF=min(pNSF)
    meanNSF=st.mean(pNSF)
    stdNSF=st.pstdev(pNSF)
    
    modeNoNSF=mode(pNoNSF)
    maxNoNSF=max(pNoNSF)
    minNoNSF=min(pNoNSF)
    meanNoNSF=st.mean(pNoNSF)
    stdNoNSF=st.pstdev(pNoNSF)
    
    row=[]
    row.append(i)
    row.append(minNSF)
    row.append(minNoNSF)
    #row.append(modeNSF)
    #row.append(modeNoNSF)
    row.append(maxNSF)
    row.append(maxNoNSF)
    row.append(meanNSF)
    row.append(meanNoNSF)
    row.append(ks[1])
    data.append(row)
df = pd.DataFrame(data, columns=["N","minNSF","minNoNSF","maxNSF","maxNoNSF","meanNSF","meanNoNSF","p-value"])
print(df)
df.to_csv("stats.csv")

    #print('N:'+str(i)+':minNSF:'+str(minNSF)+':minNoNSF:'+str(minNoNSF)+':modeNSF:'+str(modeNSF[0][0])+':modeNoNSF:'+str(modeNoNSF[0][0])+':maxNSF:'+str(maxNSF)+':maxNoNSF:'+str(maxNoNSF)+':meanNSF:'+str(meanNSF)+':meanNoNSF:'+ str(meanNoNSF)+':ks:'+str(ks))
