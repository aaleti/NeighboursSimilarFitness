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
import pandas as pd

files = ["3","4","5","6","7","8","9"]
Ptype=["nsf","no_nsf"]
data=[]
datafitnesses=[]
for fNumber in files:
    iterations=0
    if(fNumber=="9"):
        iterations=100
    else:
        iterations=1000
    for ptype in Ptype:
        for i in range(iterations):
            maxF=0
            with open("local-search-july-2017/"+ptype+fNumber) as f:
                lines = f.readlines()
                #reading the files
                for k in range(0, len(lines)):
                    line = lines[k]
                    if(str(i)+") - gen" in line):
                        k=k+3
                        line = lines[k]
                        linex = line.split(",")
                        fitnesses = []
                        #reading the search space
                        #1 - 0, 2 - 3, 3 - 0, 4 - 3, 5 - 3, 6 - 1, 7 - 0, 8 - 1
                        for item in linex:
                            itemx = item.split("-")
                            fitnesses.append(float(itemx[1]))
                            fdata=[]
                            fdata.append(fNumber)
                            fdata.append(ptype)
                            fdata.append(float(itemx[1]))
                            datafitnesses.append(fdata)
                        #calculation of good enough fitness
                        modeF=mode(fitnesses)
                        maxF=max(fitnesses)
                        minF=min(fitnesses)
                        vge=modeF[0]+(maxF-modeF[0])/2
                
                    #reading the transition probabilities
                    if("it("+str(i)+");" in line):
                        s1=line.split(" ")
                        mSize=int(s1[1])
                        P= np.array([]).reshape(0,mSize)
                        for j in range(mSize):
                            line = lines[k+j+1]
                            line=line.rstrip()
                            row = line.split(" ")
                            a = np.array([])
                            for item in row:
                                itt = float(item)
                                a = np.append(a, itt)
                            P = np.vstack([P,a])

                        lenP=len(P)
                        rm= []
                        nvge=[]
                        allRm=[]
                        listS=[]
                        #Find absorbing states and optima
                        for j in range(lenP):
                            flag=0
                            ff = 0
                            for s in range(lenP):
                                # if there are no outgoing probabilities, then this is a local/global optimum.
                                if(P[j,s]>0):
                                    ff = 1
                                
                                if(j not in listS and s not in listS):
                                    # plateoux of two solutions
                                    if(P[j,s]==1.0 and P[s,j]==1.0):
                                        flag=1
                                        listS.append(j)
                                    
                                    # absorbing state
                                    if(P[j,s]==1.0 and j==s):
                                        flag=1
                                        listS.append(j)
                                    
                                    for k in range(lenP):
                                        if(k not in listS):
                                            # plateoux of three solutions
                                            if(P[j,s]==1.0 and P[s,k]==1.0 and P[k,j]==1.0):
                                                flag=1
                                                listS.append(j)
                                            
                                            # plateoux of four solutions
                                            if(P[j,s]==1.0 and P[s,j]>0 and P[s,k]>0 and (P[s,j]+P[s,k])==1.0 and P[k,s]==1.0):
                                                flag=1
                                                listS.append(j)
                                    
                                            if(P[j,s]==1.0 and P[s,j]>0 and P[s,k]>0 and (P[s,j]+P[s,k])==1.0 and P[k,j]==1.0):
                                                flag=1
                                                listS.append(j)
                        
                            # list that keep track of absorbing states and local/global optima
                            if(flag==1 or ff==0):
                                rm.append(j)
                                allRm.append(j)
                            if(fitnesses[j]<vge):
                                nvge.append(j)
                                allRm.append(j)
                            
                        
                        keptFitnesses = []
                        removedFitnesses = []
                        nvgeFitnesses = []
                        keep=[]
                        
                        for j in range(lenP):
                            if(j in nvge):
                                nvgeFitnesses.append(fitnesses[j])
                            if(j not in rm and j not in nvge):
                                keptFitnesses.append(fitnesses[j])
                                keep.append(j)
                            if(j in rm):
                                removedFitnesses.append(fitnesses[j])
                        
                        R=np.zeros((len(keep),len(rm)), dtype='float')
                        #create a vector of 1s for calculating number of visits
                        mat1=[]
                        # canonical representation by removing absorbing states and local
                        for j in range(len(keep)):
                            mat1.append(1)
                            for s in range(len(rm)):
                                R[j,s]=P[keep[j],rm[s]]
                                    #removing
                        P=np.delete(P, allRm, axis=1)
                        P=np.delete(P, allRm, axis=0)
        
                        sm=0.0
                        sb=0.0
                        try:
                            if(len(P)>0):
                                iM=np.identity(len(P))
                                mM=iM-P
                                # Fundamental matrix
                                N = inv(mM)
                                # probability of reaching an absorbing state from any point
                                M=np.dot(N,R)
                                # expected number of steps to absorbion from any state
                                B=np.dot(N,mat1)
                                colsM = M.shape[1]
                                nrows=N.shape[0]
                                
                                # calculating the probability of reaching a global optima
                                globalC=0
                                for j in range(colsM):
                                    # if the absorbing state or optimum is a global optimum
                                    if(removedFitnesses[j]==maxF):
                                        globalC=globalC+1
                                        sumTemp=sum(row[j] for row in M)
                                        avgTemp=sumTemp/nrows
                                        sm=sm+avgTemp
                                sm=sm/globalC
                                '''
                                colsN = N.shape[1]
                                for j in range(colsN):
                                    if(keptFitnesses[j]==max):
                                        tempf=0
                                        for s in range(colsM):
                                            if(M[j,s]>0.0):
                                                tempf=1
                                        if(tempf==0):
                                            sumTemp=sum(row[j] for row in N)
                                            avgTemp=sumTemp/nrows
                                            if(avgTemp>=1.0):
                                                avgTemp=1.0
                                            sm=sm+avgTemp
                                '''

                            else:
                                countO=0
                                colsR = R.shape[1]
                                
                                for j in range(colsR):
                                    # if the absorbing state or optimum is a global optimum
                                    if(removedFitnesses[j]==maxF):
                                        countO=countO+1
                                sm=countO/colsR
                            nrows=B.shape[0]
                            globalC=0
                            for j in range(nrows):
                                if(removedFitnesses[j]==maxF):
                                    globalC=globalC+1
                                    sb=sb+B[j]
                            sb=sb/globalC
                            recD=[]
                            recD.append(fNumber)
                            recD.append(ptype)
                            #probability reaching global optimum
                            recD.append(sm)
                            #number of steps
                            recD.append(sb)
                            recD.append(globalC)
                            data.append(recD)
                        except:
                            print("error"+fNumber)

    # drawing the boxplots
    
df = pd.DataFrame(data, columns=["N","PType","Probability","Steps","NGlobal"])
df.to_csv("MCresults.csv")
df2 = pd.DataFrame(datafitnesses, columns=["N","PType","Fitness"])
df2.to_csv("MCfitnesses.csv")
