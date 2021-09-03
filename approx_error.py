import pandas as pd
# import numpy as np
import os
import matplotlib.pyplot as plt


def calc_error(refhisto, histo, tool):
    ref = pd.read_csv(refhisto, delimiter=' ', names=['Abundance', 'Numberofkmers'])
    if (tool == "SAKEIMA" or tool == "SPRISS"):
        data = pd.read_csv(histo, delimiter=' ', skiprows=1, names=['Abundance', 'Numberofkmers'])
    elif (tool == "KMERGENIE"):
        # nota salto il prio valore di kmergenie perche sballa il grafico
        data = pd.read_csv(histo, delimiter='\t', skiprows=2, names=['Abundance', 'Numberofkmers'])

    ref_tuples = list(ref.itertuples(index=False, name=None))
    data_tuples = list(data.itertuples(index=False, name=None))
    # errore
    error_1 = []

    i = 0  # counter per ref_tuples
    j = 0  # counter per data_tuples

    while i < len(ref_tuples) and j < len(data_tuples):
        if int(ref_tuples[i][0]) == int(data_tuples[j][0]):
            tmp = abs((ref_tuples[i][1] - data_tuples[j][1]) / ref_tuples[i][1])
            error_1.append((int(ref_tuples[i][0]), tmp))
            i += 1
            j += 1
        else:
            if int(ref_tuples[i][0]) < int(data_tuples[j][0]):
                i += 1
            else:
                j += 1

    error_1_file = histo + '_AE_.HISTO'
    with open(error_1_file, 'w') as hf:
        for j in error_1:
            hf.write("{} {}\n".format(j[0], j[1]))
        hf.close()

    return error_1


def appr_error(refpath, crespath, decrpath, typecres, typedecr, graphdir, datasetType):
    if (typecres == "SAKEIMA"):
        colorcres = "orange"
    elif (typecres == "SPRISS"):
        colorcres = "green"
    elif (typecres == "KMERGENIE"):
        colorcres = "blue"

    if (typedecr == "SAKEIMA"):
        colordecr = "orange"
    elif (typedecr == "SPRISS"):
        colordecr = "green"
    elif (typedecr == "KMERGENIE"):
        colordecr = "blue"

    title = 'Approximation  error ' + typecres + ' vs ' + typedecr + ' ' + datasetType

    histo1 = calc_error(refpath, crespath, typecres)  # non devono essere fissi?
    histo2 = calc_error(refpath, decrpath, typedecr)

    #histo1 = calc_error(refpath, crespath, "SPRISS")  #non devono essere fissi?
    #histo2 = calc_error(refpath, decrpath, "KMERGENIE")

    fig, ax = plt.subplots()
    ax.loglog([x[0] for x in histo1], [y[1] for y in histo1], color=colorcres)
    ax.loglog([x[0] for x in histo2], [y[1] for y in histo2], color=colordecr)
    ax.legend([typecres + ' approx error ', typedecr + ' approx error'])
    ax.set(title=title, xlabel='Abundance', ylabel='Approximation  Error')

    pngfile = os.path.join(graphdir, title + ".png")
    plt.savefig(pngfile)
    plt.show()
