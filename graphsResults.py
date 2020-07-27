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

#df = pd.read_csv('MCfitnesses.csv')

#plt.figure(figsize=(4, 2.5))
#g=sns.boxplot(y='Fitness', x='N',data=df, palette="Pastel1", hue='PType').set(xlabel=' ', ylabel=' ')
#plt.legend(title='',loc='upper left', ncol=3, fancybox=False)
#plt.savefig('AllFitnesses.pdf')

df = pd.read_csv('MCNewresults.csv')

plt.figure(figsize=(4, 2.5))
g=sns.boxplot(y='Probability', x='N',data=df, palette="Pastel1", hue='PType').set(xlabel=' ', ylabel=' ')
plt.legend(title='',loc='upper left', ncol=3, fancybox=False)
plt.savefig('Probability.pdf')
