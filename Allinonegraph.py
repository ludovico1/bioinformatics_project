import pandas as pd
#import numpy as np
import os
import matplotlib.pyplot as plt
def Allinonegraph(kmchisto,sakeimahisto,sprisshisto,kmergeniehisto, title,savedir):

    abkmc = pd.read_csv(kmchisto, delimiter=' ', names=['Abundance','Numberofkmers'])
    absakeima =pd.read_csv(sakeimahisto, delimiter=' ', names=['Abundance','Numberofkmers'])
    abspriss = pd.read_csv(sprisshisto, delimiter=' ', names=['Abundance','Numberofkmers'])
    abkmergenie= pd.read_csv(kmergeniehisto, delimiter='\t',skiprows=2, names=['Abundance','Numberofkmers'])
    #da riattivare quando a posto
    fig, ax =   plt.subplots()
    ax.loglog('Abundance', 'Numberofkmers', data=abkmergenie)

    ax.loglog('Abundance','Numberofkmers', data=absakeima)
    ax.loglog('Abundance','Numberofkmers', data=abspriss)
    ax.loglog('Abundance', 'Numberofkmers', data=abkmc)

    #ax.set(title=title, ylabel='Number of kmers', xlabel='Abundance')
    ax.legend(['KMERGENIE-ntcard', 'SAKEIMA','SPRISS','KMC'])
   # ax.legend(['KMERGENIE', 'KMC', 'SPRISS'])
    ax.set(title=title, xlabel='Abundance', ylabel='Number of kmers')
    #ax.legend(['KMC', 'SAKEIMA', 'SPRISS'])
   # fig.tight_layout()
    fig.set_size_inches(18, 10)
    # plt.plot()

    pngfile = os.path.join(savedir, title+".png")
    plt.savefig(pngfile)
    plt.show()