import pandas as pd
import os
import time

import gc

def run_makeabundance(filekmercount, dataset, tool, savedir,sample_size):
    # type 0 KMC kmer e count
    # type 1 Spriss sakeima kmer count e freqnorm
   # auxpath =os.path.split(filekmercount)
    gc.collect()
    conc=1
    start_pd = time.time()
    print("Inizio caricamento PD MOD "+dataset+" "+tool+" \n" + str(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))) + "\n")
    abundance_histo = {1: 0}

    cz=1000000
    for ichunk in pd.read_csv(filekmercount, delimiter=' ', header=None, chunksize=cz):
        print(ichunk.shape," ",str(conc)," ",str(cz*conc))
        for i in ichunk[1] :
            j=round(i*(1/sample_size)) #agg. 04/06 per adattare abundance
            if j in abundance_histo:
                abundance_histo[j] += 1
            else:
                abundance_histo[j] = 1
        conc = conc + 1
    gc.collect()

    end_pd = time.time()
    print("Finito caricamento PD MOD "+dataset+" "+tool+" \n" + str(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))) + "\n")
    print("Tempo di caricamento PD MOD: " + str(end_pd - start_pd) + "\n")


    ol = sorted(abundance_histo.keys())
    histofile = filekmercount+'.HISTO2'
    histofile = os.path.join(savedir, histofile)
    totkmershisto =0
    with open(histofile, 'w') as hf:
        for j in sorted(abundance_histo.keys()):
            # faccio già ordinamento per numero di kemers
            ab=abundance_histo[j]
            totkmershisto = totkmershisto + ab
            #ab=ab*(1/sample_size)
            #ab=round(ab)
            hf.write("{} {}\n".format(j, ab))
        hf.close()
    print(">>>finita creazione  " + histofile+ " "+ "\n kmer totali : "+str(totkmershisto)+ "\n")
    return histofile