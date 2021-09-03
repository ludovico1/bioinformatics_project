import pandas as pd

import os
import math
import matplotlib.pyplot as plt
import statistics

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)

def dist_error_2graph(refhisto, histo, histo2,tool1,tool2, title1,savedir):
    colref = "red"

    if (tool1 == "SAKEIMA"):
        colortool1 = "orange"
    elif (tool1 == "SPRISS"):
        colortool1 = "green"
    elif (tool1 == "KMERGENIE"):
        colortool1 = "blue"

    if (tool2 == "SAKEIMA"):
        colortool2 = "orange"
    elif (tool2 == "SPRISS"):
        colortool2 = "green"
    elif (tool2 == "KMERGENIE"):
        colortool2 = "blue"


    ref = pd.read_csv(refhisto, delimiter=' ', names=['Abundance', 'Numberofkmers'])
    #refhisto è sempre kmc

    if (tool1 == "SAKEIMA" or tool1 == "SPRISS"):
        conf = pd.read_csv(histo, delimiter=' ', skiprows=1, names=['Abundance', 'Numberofkmers'])
    elif (tool1 == "KMERGENIE"):
        conf = pd.read_csv(histo, delimiter='\t', skiprows=2, names=['Abundance', 'Numberofkmers'])

    if (tool2 == "SAKEIMA" or tool2 == "SPRISS"):
        conf2 = pd.read_csv(histo2, delimiter=' ', skiprows=1, names=['Abundance', 'Numberofkmers'])
    elif (tool2 == "KMERGENIE"):
        conf2 = pd.read_csv(histo2, delimiter='\t', skiprows=2, names=['Abundance', 'Numberofkmers'])



    merge1 = ref.merge(conf, left_on='Abundance',right_on='Abundance', suffixes=('_ref', '_conf1'))
    merge2 = ref.merge(conf2, left_on='Abundance', right_on='Abundance', suffixes=('_ref', '_conf2'))


    merge1['diff1'] = ((merge1['Numberofkmers_ref'] - merge1['Numberofkmers_conf1']))/merge1['Numberofkmers_ref']
    merge2['diff2'] = ((merge2['Numberofkmers_ref'] - merge2['Numberofkmers_conf2']))/merge2['Numberofkmers_ref']

    lista1 = []
    lista2 = []
    pearson_merge1=[]
    mean_merge1=[]
    ablist=[]
    for index, row in merge1.iterrows():
        print(row['Abundance'],row['Numberofkmers_ref'], row['Numberofkmers_conf1'])
        lista1.append(row['Numberofkmers_ref'])
        lista2.append(row['Numberofkmers_conf1'])
        ra=row['Abundance']
        if (index % 10 == 0  and index!=0) or index == len(merge1):
            p=pearson_def(lista1,lista2)
            print(p)
            pearson_merge1.append(p)
            mean_merge1.append(abs(statistics.mean(lista1)- statistics.mean(lista2)))
            ablist.append(ra)
            lista1 = []
            lista2 = []

    lista1 = []
    lista2 = []
    pearson_merge2 = []
    ablist2 = []
    mean_merge2 = []

    for index, row in merge2.iterrows():
        print(row['Abundance'], row['Numberofkmers_ref'], row['Numberofkmers_conf2'])
        lista1.append(row['Numberofkmers_ref'])
        lista2.append(row['Numberofkmers_conf2'])
        ra = row['Abundance']
        if (index % 10 == 0 and index != 0) or index == len(merge2):
            p = pearson_def(lista1, lista2)
            print(p)
            pearson_merge2.append(p)
            mean_merge2.append(abs(statistics.mean(lista1) - statistics.mean(lista2)))
            ablist2.append(ra)
            lista1 = []
            lista2 = []

    ##### Inizio Grafico Pearson
    fig, ax = plt.subplots()
    ax.grid(b=True, which='major', color='0.9', axis='y', linestyle='-')
    ax.grid(b=True, which='minor', color='0.9', axis='x', linestyle='-')
    ax.set(title=title1+" Pearson "+tool1+","+tool2, xlabel='Abundance', ylabel='Diff')
    shape = '.'
    fig.set_size_inches(15, 18)
    ax.semilogx(ablist, pearson_merge1,color=colortool1)
    ax.semilogx(ablist2, pearson_merge2, color=colortool2)
    ax.legend([tool1 + ' Pearson vs kmc', tool2 + ' Pearson vs kmc'])
    pngfile = os.path.join(savedir, title1+"_"+tool1+"_"+tool2  + "_Pearson.png")
    plt.savefig(pngfile)
    plt.show()
    ##### Inizio Grafico Mean
    fig2, ax2 = plt.subplots()
    ax2.set(title=title1 + " Mean 10 vs KMC " + tool1 + "," + tool2, xlabel='Abundance', ylabel='Mean')
    fig2.set_size_inches(15, 18)
    ax2.loglog(ablist, mean_merge1, color=colortool1)
    ax2.loglog(ablist2, mean_merge2, color=colortool2)
    ax2.legend([tool1 + ' Mean vs kmc', tool2 + 'Mean vs kmc'])
    ax2.grid(b=True, which='minor', color='0.9', axis='y', linestyle='-')
    ax2.grid(b=True, which='minor', color='0.9', axis='x', linestyle='-')
    ax2.set(title=title1 + " Mean 10 vs KMC " + tool1 + "," + tool2, xlabel='Abundance', ylabel='Mean')
    pngfile = os.path.join(savedir, title1+"_"+tool1+"_"+tool2  + "_mean_10.png")
    plt.savefig(pngfile)
    plt.show()

    ##### Inizio Grafico differenza
    fig3, ax3 = plt.subplots()
    ax3.set(title=title1 + " Errore Percentuale " + tool1 + " vs " + tool2, xlabel='Abundance', ylabel='Errore Percentuale')
    fig3.set_size_inches(15, 18)
    ax3.plot('Abundance', 'diff1', data=merge1, ms=5, c=colortool1,  marker=shape,  zorder=9, lw=0.1, mew=0.0)
    ax3.plot('Abundance', 'diff2', data=merge2, ms=5, c=colortool2, marker=shape, zorder=9, lw=0.1, mew=0.0)
    ax3.legend([tool1 + ' Mean vs kmc', tool2 + 'Mean vs kmc'])
    ax3.grid(b=True, which='major', color='0.9', axis='y', linestyle='-')
    ax3.grid(b=True, which='major', color='0.9', axis='x', linestyle='-')
    ax3.set(title=title1 + " Errore Percentuale " + tool1 + " vs " + tool2, xlabel='Abundance', ylabel='Errore Percentuale')
    pngfile = os.path.join(savedir, title1+"_"+tool1+"_"+tool2 + "_Errore_Percentuale.png")
    plt.savefig(pngfile)
    plt.show()
