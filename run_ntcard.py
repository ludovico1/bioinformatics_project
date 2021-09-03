#deve essere presente ntcard nella stessa dir passata a run_all con parametro -p
#dentro alla directory ntCard comppilato e funzonante
# viene utilizzato al posto di kmergenie perche' la parte che
# ci interessa per il conteggio è fatto esclusivamente da ntCard
##todo: togliere le prime righe di F0 e F1
import os
import time

def run_ntcard(fastqfile, dataset, k,working_dir,  exec_fold):
    outlogfile = os.path.join(working_dir, "log_ntcard" + dataset + "_" + str(k) + ".txt")
    of = open(outlogfile, 'w')
    print("ntCard : "+dataset)
    print("K: "+ str(k))
    start_mining = time.time()


    kmcexec = os.path.join(exec_fold, "ntCard/ntcard") #limito ai primi 10.000
    kmer_counts_file = working_dir + "/KmersCount_" + dataset + "ntcard_k" + str(k) +  ".txt"
    cmd = kmcexec + " -k" + str(k) +  " -c 10000 -p" + kmer_counts_file + " "+fastqfile
    print(">>COMANDO ntcard LANCIATO: " + cmd + "\n")
    os.system(cmd)
    end_mining = time.time()
    counting_time = end_mining - start_mining

    print("Total_time ntcard = " + str(counting_time))
    of.write("Total_time ntcard= " + str(counting_time)+" \n")

    of.write(dataset + " \n")
    of.close()
    #lo rinomina
    kmer_counts_file =kmer_counts_file+ "_k"+str(k) +".hist"
    return kmer_counts_file