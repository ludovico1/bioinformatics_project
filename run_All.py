import pyfastx
import time
import pandas as pd
from time import gmtime, strftime
import argparse
import sys
import os
import psutil
import run_kmc as kmc
import run_spriss as spriss
import run_sakeima as sakeima
import make_abundance3 as ma3
import make_abundance3_new as manew
import AbundanceHistoGraph as ab
import Allinonegraph as ag
import run_ntcard as nt
import sum_square_error as sse
import approx_error as approx
import mean_squared_error as ms
import distgraph_2graphv2 as d2gv2



# parameter parsing
parser = argparse.ArgumentParser()
parser.add_argument("-k", type=int, help="length of k-mers (>0)")
# capire se forzare il fatto che sia inserito k
parser.add_argument("-t", "--theta", type=float, help="frequency threshold (in (0,1))")
# theta parametri per sakeima
parser.add_argument("-o", "--output", help="path to output file (counts of frequent k-mers)")
# cartellina dove vengono creati output
parser.add_argument("-D", "--dataset", help="Mnemonic name for dataset")
# nome del dataset

parser.add_argument("-f", help="path to input fastq file (dataset of reads)")

parser.add_argument("-thr", type=int, help="Number of threads to use for counting (>0, def. 1)")
parser.add_argument("-v", "--verbose", action="store_true", help="Verbose/debug Mode ")
parser.add_argument("-p", "--path", help="tools Path dir ")
# aggiunto la pth dir per i tools
# mettiamo anche il parametro per decidere autonomamente quanti thread usare, default verrà scelto con calcolo psutil

# parametri che seguono sono opzionali vedere se ha senso lasciarli
parser.add_argument("-l", "--lambda", type=float,
                    help="desired fraction between sample size and dataset size (in (0,2))")
parser.add_argument("-e", "--epsilon", type=float,
                    help="approximation accuracy parameter (in (0,1), def. theta + 2/dbtot)")
parser.add_argument("-d", "--delta", type=float, help="approximation confidence parameter (in (0,1), def. 0.1)",
                    default=0.1)
args = parser.parse_args()
if not args.path:
    print(" Path dir for tools is needed")
    print("-----------------------------------")
    parser.print_help(sys.stderr)
    exit()
install = args.path

if not args.k:
    print(" Error: Argument k is needed")
    print("-----------------------------------")
    parser.print_help(sys.stderr)
    exit()
else:
    if args.k <= 0:
        print("Error: Argument k needs to be >= 0")
        print("-----------------------------------")
        parser.print_help(sys.stderr)
        exit()
if not os.path.isfile(str(args.f)):
    print("Error: path to dataset dirbase is not correct of file is not found !")
    print("----------------------------------------------")
    parser.print_help(sys.stderr)
    exit()
if not args.dataset:
    print(" Error: DATASET name is needed")
    print("-----------------------------------")
    parser.print_help(sys.stderr)
    exit()


def print_parameters():
    # prova per vedere se tutti i parametri sono stati messi correttamente

    print("Output   :" + str(args.output))
    print("K        :" + str(args.k))


def print_fastq_statistics():
    print("Stampa delle statistiche")
    print("argf " + str(args.f))
    fq = pyfastx.Fastq(str(args.f))

    print("File name      : " + str(fq.file_name))
    print("Dataset name   : " + str(args.dataset))
    print("\n")
    print("================================= DATASET STATISTICS ====================================")
    print("Dataset Len Reads: " + str(len(fq)))
    totreads = len(fq)
    print("Dateset Size     : " + str(fq.size))
    datasetssize = fq.size
    print("Avg Lenght       : " + str(fq.avglen))

    avgreadlength = fq.avglen
    print("Max Lenght       : " + str(fq.maxlen))

    maxreadlength = fq.maxlen
    print("Min Lenght       : " + str(fq.minlen))
    print("Encoding         : " + str(fq.encoding_type))
    print("Max quality      : " + str(fq.maxqual))
    print("Min quality      : " + str(fq.minqual))
    print("Composition      : " + str(fq.composition))
    print("GC content       : " + str(fq.gc_content))
    print("=========================================================================================")
    return totreads, datasetssize, avgreadlength, maxreadlength


def print_cpu_mem_stats():
    # questa parte per il calcolo dei thread disponibili e della memoria da utilizzare
    print("=============================== PC resource statistics ==================================")

    print("Virtual memory: " + str(psutil.virtual_memory()))
    memavailable = psutil.virtual_memory().available
    memfree = ((
                   psutil.virtual_memory().free) * 0.8)  # da capire come impostarlo correttamente 80% completamento della memoria rimanente
    memfree = int(memfree / 1024 / 1024 / 1024)  # tasformo in giga e arrotondo 2 cifre
    if memfree == 0:
        memfree = 1  # caso peggiore che non si abbia giga libero maKMC per esempio richiede almeno 1 GB

    print("memoria available: " + str(memavailable / 1024 / 1024 / 1024) + "free: " + str(memfree))
    print("Numero Core Totali disponibili :" + str(psutil.cpu_count()))
    numthread = psutil.cpu_count() - 1  # fissa come thread possibili i numeri di core a disposizione meno 1
    print("===========================================================================================")
    return numthread, memfree


def create_dirs():
    # verifica se ci sono le directory per i vari script altrimenti le crea
    tools = ["KMC", "SAKEIMA", "SPRISS", "KMERGENIE","Graph"]

    parent_dir = str(args.output)
    dataset = str(args.dataset)
    working_dir = []

    if os.path.isdir(parent_dir):
        print("ok directory esistente:" + parent_dir)
        parent_dir = os.path.join(parent_dir, dataset)
        if not os.path.isdir(parent_dir):
            try:
                print("Attempting create dataset: " + str(dataset))
                os.mkdir(parent_dir)
            except OSError as error:
                print("Creation problem: " + str(error))
                return False, working_dir
        else:
            print("ok directory esistente:" + parent_dir)
        for toolname in tools:
            path = os.path.join(parent_dir, toolname)
            if not os.path.isdir(path):
                working_dir.append(str(path))
                try:
                    print("Attempting create " + str(path))
                    os.mkdir(path)
                except OSError as error:
                    print("Creation problem: " + str(error))
                    return False, working_dir
            else:
                working_dir.append(str(path))

        return True, working_dir
    else:
        print(" Error: " + str(args.output) + " Output parameters must be a directory")
        return False, working_dir


def set_progress_status():
    parent_dir = str(args.output)
    dataset = str(args.dataset)
    progress_file = os.path.join(parent_dir, "progress_" + dataset + ".txt")
    pv = -1
    try:
        with open(progress_file, "r") as pf:
            print("Progress file opened")
            pv = int(pf.read())
            pv = pv + 1  # incremento
            print(" Progress found,  Next Step :" + str(pv))
            pf.close()

    except IOError:
        print("Progress File not accessible")
        pv = -1
        exit(-1)  # faccio uscire per far gestire manualmente

    try:
        with open(progress_file, "w") as pf:
            print("Progress file opened")
            pf.write(str(pv))
            print(" Progress Status, Step:" + str(pv))
            pf.close()

    except IOError:
        print("Progress File not accessible")
        pv = -1

    # setta stato
    return pv


def read_progress_status():
    # utile per la gestione dello stato avanzamento per eventuali ripartenze
    # è un file progres_status.txt dentro alla directory padre del dataset
    # se scritto 1 vuol dire che la fase 1 è finita e conclusa, aggiornamento va alla fine
    # se il primo giro non trova il file lo carica con valore 0

    parent_dir = str(args.output)
    dataset = str(args.dataset)

    progress_file = os.path.join(parent_dir, "progress_" + dataset + ".txt")
    pv = -1
    if os.path.isfile(progress_file):
        try:
            with open(progress_file, "r") as pf:

                pv = int(pf.read())
                pf.close()

        except IOError:
            print("Progress File not accessible")
            pv = -1
            exit(-1)
    else:

        try:
            with open(progress_file, "w") as pf:
                print("Progress file creation")
                pv = 0
                pf.write(str(pv))
                print(" Progress Step found:" + str(pv))
                pf.close()
                # il file non era presente è stato creato e messo a 0
        except IOError:
            print("Creation File error")
            pv = -1
            exit(-1)
    return pv


def write_log_exec_result(line,printconsole):
    # scrive sul log dei risultati per mantenere tempistiche e altre info da fuori dei tools
    # crea un file assieme a progress
    #se printconsole =0 non scrive a video se è 1 sì
    parent_dir = str(args.output)
    dataset = str(args.dataset)

    logfile = os.path.join(parent_dir, "Log_exec_result_" + dataset + ".txt")

    try:
        with open(logfile, "a") as lf:

            lf.write(str(line) + "\n")
            if printconsole == 1:
                print(str(line))
            lf.close()
            # il file non era presente è stato creato e messo a 0
    except IOError:
        print("Creation File error")

        exit(-1)
    return 0

def write_histo_path(pathhisto,printconsole):
    # scrive sul log le path histo per passarli all'ultima fase anche con ripartenza
    #se printconsole =0 non scrive a video se è 1 sì
    parent_dir = str(args.output)
    dataset = str(args.dataset)

    logfile = os.path.join(parent_dir, "Path_histo" + dataset + ".txt")

    try:
        with open(logfile, "a") as lf:
            #print("LOG FILE writing")

            lf.write(str(pathhisto) + "\n")
            if printconsole == 1:
                print("Scritto nel file Histo:" + str(pathhisto))
            lf.close()
            # il file non era presente è stato creato e messo a 0
    except IOError:
        print("Creation File error")

        exit(-1)
    return 0

def read_histo_path(numtool):
    # Legge le histo path dai tool, viene passato solo il numero dell'ultimo tool
    parent_dir = str(args.output)
    dataset = str(args.dataset)

    logfile = os.path.join(parent_dir, "Path_histo" + dataset + ".txt")
    listpath=[]
    try:
        ref = pd.read_csv(logfile, delimiter=',', names=['numtool', 'histopath'])
        for index, row in ref.iterrows():
            print("tool index"+str(index))
            if (index+1) == row['numtool']:
                print("path corretto: "+str(row['histopath']))
                listpath.append(str(row['histopath']))
                if index+1 >numtool:
                    print("Errore file di path histo!")
                    break


    except IOError:
        print("Creation File error")

        exit(-1)

    return listpath


def main_run_all():
    working_dir = []
    dataset = str(args.dataset)
    datafile = str(args.f)
    K = args.k
    if not args.delta:
        delta = 0.1  # dovrebbe essere già di default a 0
    else:
        delta = args.delta  # se non esiste delta imposto a 0
    theta = args.theta

    if args.verbose:
        print_parameters()
    # print_fastq_statistics()

    tot_reads, datasets_size, avg_read_length, max_read_length = print_fastq_statistics()

    num_thread, mem_free = print_cpu_mem_stats()
    print("numero di thread UTILIZZABILI per operazioni:" + str(num_thread))
    print("Memoria UTILIZABILE per operazioni:" + str(mem_free))

    laststep = 6  # INDICARE ULTIMO PASSO

    if read_progress_status() < laststep:  # nel caso iniziale parte da stato 0 e set_prog lo aumenta, nel caso di ripartenza
        # lo legge , controlla le directory in ogni caso
        # e parte da quello successivo
        check, working_dir = create_dirs()
        if not check:
            print(" Directory creations error aborting...")
            exit(1)
        set_progress_status()

    if read_progress_status() == 1:
        print("starting form step 1 - KMC")
        write_log_exec_result("Inizio tool KMC per Dataset: " + str(args.dataset),1)
        start_time = time.time()
        write_log_exec_result("Start time KMC:" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)

        print("lancio KMC")
        kfile = kmc.run_kmc(datafile, dataset, datasets_size, K, working_dir[0], install, num_thread, mem_free)
        end_time = time.time()

        write_log_exec_result("End time KMC :" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)
        duration = (end_time - start_time)
        write_log_exec_result("Total duration KMC " + dataset + ": " + str(duration),1)

        write_log_exec_result("Inizio abundance file per tool KMC e Dataset: " + str(args.dataset),1)
        write_log_exec_result("Start time abundance :" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)
        start_ab = time.time()
        kmchisto = ma3.run_makeabundance(kfile, dataset, 'KMC', working_dir[0])
        #print("Fine abundance file per tool KMC e Dataset: " + str(args.dataset))

        write_log_exec_result("End time abundance file KMC,Dataset: " + str(args.dataset) + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)
        end_ab = time.time()
        durationab =end_ab-start_ab
        write_log_exec_result("Total duration abundance KMC, Dataset " + str(args.dataset) + ": " + str(durationab),1)
        ab.abundancehistograph(kmchisto,"KMC  dataset: "+dataset+" K:"+str(K)+"theta:"+str(theta))

        write_histo_path("1," + kmchisto, 1)
        set_progress_status()

    if read_progress_status() == 2:
        print("starting form step 2 - Sakeima")
        write_log_exec_result("Starting tool Sakeima for Dataset: " + str(args.dataset),1)
        start_time = time.time()
        write_log_exec_result("Start time:" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)

        print("lancio Sakeima")

        safile, sample_sakeima = sakeima.run_sakeima(datafile, dataset, datasets_size, tot_reads, avg_read_length, max_read_length, K,delta, theta, working_dir[1], install, num_thread, mem_free)

        end_time = time.time()
        write_log_exec_result("End time:" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)
        duration = (end_time - start_time)
        write_log_exec_result("Total duration SAKEIMA " + dataset + ": " + str(duration),1)
        print("Crea HistoTab per Sakeima")
        print("Inizio abundance file per tool SAKEIMA e Dataset: " + str(args.dataset))
        write_log_exec_result("Inizio abundance file per tool SAKEIMA e Dataset: " + str(args.dataset), 1)
        write_log_exec_result("Start time abundance :" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()), 1)
        start_ab = time.time()
        sakeimahisto =ma3.run_makeabundance(safile, dataset, 'SAKEIMA', working_dir[1])
        # print("Fine abundance file per tool KMC e Dataset: " + str(args.dataset))

        write_log_exec_result("End time abundance file SAKEIMA,Dataset: " + str(args.dataset) + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)
        end_ab = time.time()
        durationab = end_ab - start_ab
        write_log_exec_result("Total duration abundance SAKEIMA, Dataset " + str(args.dataset) + ": " + str(durationab), 1)
        ab.abundancehistograph(sakeimahisto, "SAKEIMA  dataset: "+dataset+" K:"+str(K)+"theta:"+str(theta))

        write_log_exec_result("Inizio abundance MODIFICATO file per tool SAKEIMA e Dataset: " + str(args.dataset), 1)
        write_log_exec_result("Start time abundance MODIFICATO :" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()), 1)
        start_ab = time.time()

        sakeimahisto2 = manew.run_makeabundance(safile, dataset, 'SAKEIMA_MOD', working_dir[1],
                                               sample_sakeima)
        # print("Fine abundance file per tool KMC e Dataset: " + str(args.dataset))
        write_log_exec_result("End time abundance file SAKEIMA MODIFICATO,Dataset: " + str(args.dataset) + strftime(
            "%a, %d %b %Y %H:%M:%S", gmtime()), 1)
        end_ab = time.time()
        durationab = end_ab - start_ab
        write_log_exec_result( "Total duration abundance SAKEIMA MODIFICATO, Dataset " + str(args.dataset) + ": " + str(durationab), 1)
        # aggiungiamo il finale per avere sia histo che histo2 shiftato
        write_histo_path("2," + sakeimahisto, 1)

        set_progress_status()

    if read_progress_status() == 3:
        print("starting form step 3 - SPRISS")
        write_log_exec_result("Inizio tool SPRISS per Dataset: " + str(args.dataset),1)
        start_time = time.time()
        write_log_exec_result("Start time:" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)

        print("nella caltella di lavoro: " + working_dir[2])
        spfile,sample_spriss=spriss.run_spriss(datafile, dataset, datasets_size, tot_reads, avg_read_length, max_read_length, K, delta, theta, working_dir[2], install, num_thread, mem_free)
        end_time = time.time()
        write_log_exec_result("End time:" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)
        duration = (end_time - start_time)
        write_log_exec_result("Total duration SPRISS " + dataset + ": " + str(duration),1)
        print("Crea HistoTab per SPRISS")

        print("Inizio abundance file per tool SPRISS e Dataset: " + str(args.dataset))
        write_log_exec_result("Inizio abundance file per tool SPRISS e Dataset: " + str(args.dataset), 1)
        write_log_exec_result("Start time abundance :" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()), 1)
        start_ab = time.time()
        sprisshisto= ma3.run_makeabundance(spfile, dataset, 'SPRISS', working_dir[2])
        # print("Fine abundance file per tool KMC e Dataset: " + str(args.dataset))
        write_log_exec_result("End time abundance file SPRISS,Dataset: " + str(args.dataset) + strftime("%a, %d %b %Y %H:%M:%S", gmtime()), 1)
        end_ab = time.time()
        durationab = end_ab - start_ab
        write_log_exec_result("Total duration abundance SPRISS, Dataset " + str(args.dataset) + ": " + str(durationab),1)
        #*************** variazione 5/06
        write_log_exec_result("Inizio abundance MODIFICATO file per tool SPRISS e Dataset: " + str(args.dataset), 1)
        write_log_exec_result("Start time abundance MODIFICATO :" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()), 1)
        start_ab = time.time()

        sprisshisto2 = manew.run_makeabundance(spfile, dataset, 'SPRISS_MOD', working_dir[2],sample_spriss)
        # print("Fine abundance file per tool KMC e Dataset: " + str(args.dataset))
        write_log_exec_result("End time abundance file SPRISS MODIFICATO,Dataset: " + str(args.dataset) + strftime("%a, %d %b %Y %H:%M:%S",gmtime()), 1)
        end_ab = time.time()
        durationab = end_ab - start_ab
        write_log_exec_result("Total duration abundance SPRISS MODIFICATO, Dataset " + str(args.dataset) + ": " + str(durationab),1)
        #*********


        ab.abundancehistograph(sprisshisto, "SPRISS  dataset: "+dataset+" K:"+str(K)+"theta:"+str(theta))
        ab.abundancehistograph(sprisshisto2, "SPRISS MOD  dataset: " + dataset + " K:" + str(K) + "theta:" + str(theta))
        write_histo_path("3," + sprisshisto, 1)
        #aggiungiamo il finale per avere sia histo che histo2 shiftato
        set_progress_status()

    if read_progress_status() == 4:


        print("starting form step 4 - KMERGENIE-ntcard")
        write_log_exec_result("Inizio tool KMERGENIE-ntcard per Dataset: " + str(args.dataset),1)
        start_time = time.time()
        print("lancio ntcard")
        write_log_exec_result("Start time ntcard:" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)

        ntcardhisto = nt.run_ntcard(datafile, dataset, K, working_dir[3], install)


        end_time = time.time()
        write_log_exec_result("End time:" + strftime("%a, %d %b %Y %H:%M:%S", gmtime()),1)
        duration = (end_time - start_time)
        write_log_exec_result("Total duration KMERGENIE-ntcard " + dataset + ": " + str(duration),1)
        print("Crea HistoTab per KMERGENIE-ntcard")
        write_histo_path("4," + ntcardhisto, 1)

        set_progress_status()

    if read_progress_status() == 5:

        #questo funziona solo se eseguito di seguito non cerca le .histo
        print("Elaborazione e creazione Grafici")
        #ag.Allinonegraph(kmchisto , sakeimahisto2, sprisshisto2, ntcardhisto, 'histo2 Dataset: '+dataset+" K "+ str(K)+" Theta: "+str(theta),working_dir[4])


        kmchisto, sakeimahisto, sprisshisto, ntcardhisto = read_histo_path(4)
        #forzo ad affiungere solo il 2 finale per non dover usare 2 file di caricamento
        #nel caso di faccia ripartenza solo da grafici
        sakeimahisto2 =sakeimahisto+'2'
        sprisshisto2=sprisshisto+'2'


        ag.Allinonegraph(kmchisto, sakeimahisto2, sprisshisto2, ntcardhisto,'histo2_' + dataset , working_dir[4])

        #ag.Allinonegraph(kmchisto, sakeimahisto, sprisshisto, ntcardhisto, 'histo1 Dataset: ' + dataset + " K " + str(K) + " Theta: " + str(theta),working_dir[4])
        ag.Allinonegraph(kmchisto, sakeimahisto, sprisshisto, ntcardhisto, 'histo1_' + dataset , working_dir[4])
        print("Grafici differenze:TODO")
        #d2g.dist_error_2graph(kmchisto,sakeimahisto2,ntcardhisto,"SAKEIMA","KMERGENIE",'prova1',working_dir[4])
        d2gv2.dist_error_2graph(kmchisto,sakeimahisto2,ntcardhisto,"SAKEIMA","KMERGENIE",dataset,working_dir[4])
        d2gv2.dist_error_2graph(kmchisto, sprisshisto2, ntcardhisto, "SPRISS", "KMERGENIE", dataset, working_dir[4])
        #vedere se mettere dataset anche qui
        #verificare se salva grafici

        print("Grafico SSE")
        #confronto sempre sakeima o spriss vs kmergenie
        #sse.diff_error(kmchisto,sakeimahisto2,ntcardhisto,"SAKEIMA","KMERGENIE",working_dir[4],dataset)
        sse.diff_error(kmchisto, ntcardhisto,sakeimahisto2, "KMERGENIE", "SAKEIMA", working_dir[4], dataset)
        #sse.diff_error(kmchisto, sprisshisto2, ntcardhisto, "SPRISS", "KMERGENIE", working_dir[4], dataset)
        sse.diff_error(kmchisto, ntcardhisto, sprisshisto2, "KMERGENIE", "SPRISS", working_dir[4], dataset)
        print("Grafici Approx")
        approx.appr_error(kmchisto,sakeimahisto2,ntcardhisto,"SAKEIMA","KMERGENIE",working_dir[4],dataset)
        approx.appr_error(kmchisto,sakeimahisto2,sprisshisto2,"SAKEIMA","SPRISS",working_dir[4],dataset)
        approx.appr_error(kmchisto, sprisshisto2, ntcardhisto, "SPRISS", "KMERGENIE", working_dir[4], dataset)
        print(" Grafici Sum Square Error")

        #ms.sd_error(kmchisto, sakeimahisto2, ntcardhisto, "SAKEIMA", "KMERGENIE", working_dir[4], dataset)
        ms.sd_error(kmchisto, ntcardhisto, sakeimahisto2, "KMERGENIE", "SAKEIMA", working_dir[4], dataset)
        #ms.sd_error(kmchisto, sprisshisto2, ntcardhisto, "SPRISS", "KMERGENIE", working_dir[4], dataset)
        ms.sd_error(kmchisto, ntcardhisto, sprisshisto2, "KMERGENIE", "SPRISS", working_dir[4], dataset)
        set_progress_status()
    if read_progress_status() >= laststep:
        print("Raggiunta la fine")

main_run_all()
