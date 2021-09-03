import pandas as pd
import os
import time

import gc

def run_makeabundance(filekmercount, dataset, tool, savedir):
    # type 0 KMC kmer e count
    # type 1 Spriss sakeima kmer count e freqnorm
   # auxpath =os.path.split(filekmercount)
    gc.collect()
    conc=1
    start_pd = time.time()
    print("Inizio caricamento PD"+dataset+" "+tool+" \n" + str(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))) + "\n")
    abundance_histo = {1: 0}

    cz=1000000
    skipr=0
    delimiterchar=' '
    #if tool=='KMERGENIE': #non serve per kmergenie, fa già abundance
     #   skipr=2 #salto righe F0 e F1
      #  delimiterchar='\t'

    for ichunk in pd.read_csv(filekmercount, delimiter=delimiterchar, header=None, skiprows=skipr, chunksize=cz):

        print(ichunk.shape," ",str(conc)," ",str(cz*conc))
        for i in ichunk[1] :
            if i in abundance_histo:
                abundance_histo[i] += 1
            else:
                abundance_histo[i] = 1
        conc = conc + 1
    gc.collect()

    end_pd = time.time()
    print("Finito caricamento PD "+dataset+" "+tool+" \n" + str(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))) + "\n")
    print("Tempo di caricamento PD: " + str(end_pd - start_pd) + "\n")


    ol = sorted(abundance_histo.keys())
    histofile = filekmercount+'.HISTO'
    histofile = os.path.join(savedir, histofile)
    totkmershisto = 0
    with open(histofile, 'w') as hf:
        for j in sorted(abundance_histo.keys()):
            # faccio già ordinamento per numero di kemers
            ab = abundance_histo[j]
            totkmershisto = totkmershisto + ab
            hf.write("{} {}\n".format(j, ab))
        hf.close()
    print(">>>finita creazione  " + histofile + " " + "\n kmer totali : " + str(totkmershisto) + "\n")
    return histofile