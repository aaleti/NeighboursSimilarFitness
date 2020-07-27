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

files = ["3","4","5","6","7","8","9"]
Ptype=["nsf","no_nsf"]
data = []
for fNumber in files:
    iterations=0
    print(fNumber)
    if(fNumber=="9"):
        iterations=100
    else:
        iterations=1000
    for ptype in Ptype:
        for i in range(iterations):
            maxF=0
            with open("local-search-july-2017/"+ptype+fNumber) as f:
                lines = f.readlines()
                errorLine = ""
                #reading the files
                for k in range(0, len(lines)):
                    line = lines[k]
                    if(str(i)+") - gen" in line):
                        errorLine = line
                        k=k+3
                        line = lines[k]
                        linex = line.split(",")
                        fitnesses = []
                        #reading the search space
                        #1 - 0, 2 - 3, 3 - 0, 4 - 3, 5 - 3, 6 - 1, 7 - 0, 8 - 1
                        for item in linex:
                            itemx = item.split("-")
                            fitnesses.append(float(itemx[1]))
                        #calculation of good enough fitness
                        modeF=mode(fitnesses)
                        maxF=max(fitnesses)
                        minF=min(fitnesses)
                        stdF=st.pstdev(fitnesses)
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
                
                        removedFitnesses = []
                        globalC=0
                        for j in range(lenP):
                            if(j in rm):
                                removedFitnesses.append(fitnesses[j])
                                if(fitnesses[j]==maxF):
                                    globalC=globalC+1
                            
                        dat=[fNumber,ptype,maxF,minF,modeF[0][0],modeF[1][0],stdF,globalC,len(removedFitnesses)]
                        data.append(dat)

                                
df = pd.DataFrame(data, columns=["N","PType","MaxF","MinF","ModeF","NMode","STD","NglobalOptima","NAbsorbing"])
df.to_csv("resultsAll2.csv")
