import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
def abundancehistograph(filehisto, title):
    abundance = pd.read_csv(filehisto,delimiter=' ', names=['Abundance','Numbersofkmers'])
    sum_abundance=abundance['Abundance'].sum()
    print(sum_abundance)
    fig, ((ax1, ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2)
    ax1.semilogx(abundance['Abundance'],abundance['Numbersofkmers'] )
    #ax1.set(title='semilogx '+title, ylabel='Abundance', xlabel='Frequence')
    ax1.set(title='semilogx ', xlabel='Abundance', ylabel='Number of kmers')
    ax1.grid()
    ax2.semilogy(abundance['Abundance'],abundance['Numbersofkmers'] )
    #ax2.set(title='semilogy '+title , ylabel='Abundance', xlabel='Frequence')
    ax2.set(title='semilogy ' , xlabel='Abundance', ylabel='Number of kmers')
    ax2.grid()
    ax3.loglog(abundance['Abundance'],abundance['Numbersofkmers'] )
    #ax3.set(title='loglog '+title, ylabel='Abundance', xlabel='Frequence')
    ax3.set(title='loglog ' , xlabel='Abundance', ylabel='Number of kmers')
    ax3.grid()
    #ax4.loglog(abundance['FreqNorm'], abundance['Abundance'] )
    #ax4.set(title='Norm. Freq. loglog', ylabel='frequenze', xlabel='abundance')
    #ax4.grid()
    #ax5.loglog(abundance['FreqNorm'], (abundance['Abundance']/sum_abundance) )
    #ax5.set(title='Norm. Freq. Norm. Abundance loglog', ylabel='frequenze', xlabel='abundance')
    #ax5.grid()
    #plt.semilogx(abundance['FreqKmer'], abundance['Abundance'],color='blue')
    plt.title=title
    fig.tight_layout()

    #plt.plot()
    plt.show()