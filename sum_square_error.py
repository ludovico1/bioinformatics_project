import pandas as pd
# import numpy as np
import os
import matplotlib.pyplot as plt


def calc_error_cresc(refhisto, histo, tool):
    ref = pd.read_csv(refhisto, delimiter=' ', names=['Abundance', 'Numberofkmers'])
    if (tool == "SAKEIMA" or tool == "SPRISS"):
        data = pd.read_csv(histo, delimiter=' ', skiprows=1, names=['Abundance', 'Numberofkmers'])
    elif (tool == "KMERGENIE"):
        #nota salto il prio valore di kmergenie perche sballa il grafico
        data = pd.read_csv(histo, delimiter='\t', skiprows=3, names=['Abundance', 'Numberofkmers'])

    ref_tuples = list(ref.itertuples(index=False, name=None))
    data_tuples = list(data.itertuples(index=False, name=None))
    #errore crescente
    error_1 = []
    error_1.append((0, 0))

    i = 0  # counter per ref_tuples
    j = 0  # counter per data_tuples

    while i < len(ref_tuples) and j < len(data_tuples):
        if int(ref_tuples[i][0]) == int(data_tuples[j][0]):
            #tmp = error_1[len(error_1) - 1][1] + (((ref_tuples[i][1] - data_tuples[j][1])/(ref_tuples[i][1] + data_tuples[j][1])) ** 2)
            #tmp = error_1[len(error_1) - 1][1] + (abs((ref_tuples[i][1] - data_tuples[j][1]) / (ref_tuples[i][1] + data_tuples[j][1])))
            tmp = error_1[len(error_1) - 1][1] + ((ref_tuples[i][1] - data_tuples[j][1]) ** 2)
            error_1.append((int(ref_tuples[i][0]), tmp))
            i += 1
            j += 1
        else:
            if int(ref_tuples[i][0]) < int(data_tuples[j][0]):
                i += 1
            else:
                j += 1

    error_1_file = histo + '_RSS_crescente.HISTO'
    with open(error_1_file, 'w') as hf:
        for j in error_1[1:]:
            hf.write("{} {}\n".format(j[0], j[1]))
        hf.close()

    return error_1[1:]


def calc_error_decr(refhisto, histo, tool):
    ref = pd.read_csv(refhisto, delimiter=' ', names=['Abundance', 'Numberofkmers'])
    if (tool == "SAKEIMA" or tool == "SPRISS"):
        data = pd.read_csv(histo, delimiter=' ', skiprows=1, names=['Abundance', 'Numberofkmers'])
    elif (tool == "KMERGENIE"):
        # nota salto il prio valore di kmergenie perche sballa il grafico
        data = pd.read_csv(histo, delimiter='\t', skiprows=3, names=['Abundance', 'Numberofkmers'])

    ref_tuples = list(ref.itertuples(index=False, name=None))
    data_tuples = list(data.itertuples(index=False, name=None))
    # errore crescente
    error_1 = []
    error_1.append((0, 0))

    i = len(ref_tuples) -1 # counter per ref_tuples
    j = len(data_tuples) -1 # counter per data_tuples

    while i >=0 and j >=0 :
        if int(ref_tuples[i][0]) == int(data_tuples[j][0]):
            # tmp = error_1[len(error_1) - 1][1] + (((ref_tuples[i][1] - data_tuples[j][1])/(ref_tuples[i][1] + data_tuples[j][1])) ** 2)
            # tmp = error_1[len(error_1) - 1][1] + (abs((ref_tuples[i][1] - data_tuples[j][1]) / (ref_tuples[i][1] + data_tuples[j][1])))
            tmp = error_1[len(error_1) - 1][1] + ((ref_tuples[i][1] - data_tuples[j][1]) ** 2)
            error_1.append((int(ref_tuples[i][0]), tmp))
            i -= 1
            j -= 1
        else:
            if int(ref_tuples[i][0]) < int(data_tuples[j][0]):
                j -= 1
            else:
                i -= 1

    error_1_file = histo + '_RSS_decrescente.HISTO'
    with open(error_1_file, 'w') as hf:
        for j in error_1[1:]:
            hf.write("{} {}\n".format(j[0], j[1]))
        hf.close()

    return error_1[1:]

def diff_error(refpah,crespath,decrpath,typecres,typedecr,graphdir,datasetType):
    if (typecres == "SAKEIMA"):
        colorcres="orange"
    elif(typecres == "SPRISS"):
        colorcres = "green"
    elif (typecres == "KMERGENIE"):
        colorcres = "blue"

    if (typedecr == "SAKEIMA"):
        colordecr="orange"
    elif(typedecr == "SPRISS"):
        colordecr = "green"
    elif (typedecr == "KMERGENIE"):
        colordecr = "blue"

    title = 'RSS error '+typecres+' vs '+typedecr+' ' + datasetType
    histo1 = calc_error_cresc(refpah, crespath, typecres)
    histo2 = calc_error_decr(refpah, decrpath, typedecr)

    fig, ax = plt.subplots()
    ax.loglog([x[0] for x in histo1], [y[1] for y in histo1],color=colorcres)
    ax.loglog([x[0] for x in histo2], [y[1] for y in histo2], color=colordecr)
    ax.legend([typecres+' RSS error crescente', typedecr+' RSS error decrescente'])
    ax.set(title=title, xlabel='Abundance' ,ylabel='RSS Error')

    pngfile = os.path.join(graphdir, title + ".png")
    plt.savefig(pngfile)
    plt.show()