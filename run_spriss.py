import os
import time
import math


def run_spriss(fastqfile, dataset, datasets_size, tot_reads,avg_read_length,max_read_length, k, delta, thetastr, working_dir,  exec_fold, thread, memory):

    sample_folder = os.path.join(working_dir, "sample/")
    if not os.path.isdir(sample_folder):
        try:
            print("Attempting create " + str(sample_folder))
            os.mkdir(sample_folder)
        except OSError as error:
            print("Creation problem: " + str(error))

    outlogfile= os.path.join(working_dir, "log_spriss_approx_mining_time_"+dataset+"_"+str(thetastr)+".txt")
    of = open(outlogfile, 'w')

    # aggiunta theta in formato scientifico non puo' essere semore con la base 10e-8
    theta = float(thetastr)

    print(dataset)
    print(str(theta))

    of.write(dataset + " \n")


    print(">>PRENDO file:" + str(fastqfile) + "\n")

    of.write(str(theta) + " \n")
    start_sample = time.time()

    output_file =os.path.join(working_dir, dataset+ "_kmc_" + str(k) + "-mers_db")
    epsilon = theta - 2.0 / datasets_size
    l = math.floor((0.9 / theta) / (avg_read_length - k + 1))
    m = math.ceil((2 / ((epsilon * l * (avg_read_length - k + 1)) ** 2)) * (math.floor(math.log2(min(2 * l * (max_read_length - k + 1), 4 ** 31))) + math.log(2.0 / delta)))
    ml = int(m * l)
    sample_size = float(ml) / float(tot_reads)
    print("Sample_size= " + str(sample_size))
    of.write("Sample_size= " + str(sample_size) + " \n")
    of.write("Dataset Size:"+ str(datasets_size) + " \n")
    of.write("Total Reads:"+ str(tot_reads) + " \n")
    of.write("l: "+str(l)+ " \n")
    of.write("m: " + str(m) + " \n")
    of.write("ml: " + str(ml) + " \n")
    of.write("epsilon: " + str(epsilon) + " \n")
    of.write("Theta: " + str(theta) + " \n")
    sample_file =os.path.join(sample_folder,dataset+ "_kmc_sample.fastq" )

    print(">>CREO FILE :" + str(fastqfile) + "\n")

    cmd_sample_dir = os.path.join(exec_fold , "SPRISS/scripts/")
    cmd_sample_exec = os.path.join(cmd_sample_dir , "create_sample ")

    if not os.path.isfile(cmd_sample_exec):
        try :
            os.system("g++ -o"+ cmd_sample_dir+"create_sample"+" "+ cmd_sample_dir+"create_sample.cpp")
        except OSError as error:
            print("Create sample compilation problem: " + str(error))

    cmd = cmd_sample_exec+" " + fastqfile + " " + sample_file + " " + str(int(tot_reads)) + " " + str(ml)
    print(cmd)
    os.system(cmd)
    end_sample = time.time()
    print("Time_sample_creation= " + str(end_sample - start_sample))
    of.write("Time_sample_creation= " + str(end_sample - start_sample) + " \n")
    start_mining = time.time()
    MemMax = memory  # Imposto Massimo 4 GB

    kmcexec = os.path.join(exec_fold, "SPRISS/bin/kmc")
    cmd = kmcexec +" -v -k"+ str(k) + " -cs" + str(datasets_size) + " -m" + str(MemMax) + " -ci1 -t"+str(thread) +" "+ sample_file + " " + output_file + " "+ sample_folder
    print(">>COMANDO KMC LANCIATO: " + cmd + "\n")
    os.system(cmd)
    end_mining = time.time()
    counting_time = end_mining - start_mining
    start_dump1 = time.time()
    denominator = m * l * (avg_read_length - k + 1)

    kmer_counts_file = working_dir + "/KmersCount_"+dataset+"SPRISS_k"+str(k)+"_theta"+str(thetastr)+".txt"


    kmcdumpexec = os.path.join(exec_fold, "SPRISS/bin/kmc_dump")
    cmd = kmcdumpexec +" -ci1 -theta" + str(theta) + " -epsilon" + str(epsilon) + " -n_bags" + str(
        int(m)) + " -denominator" + str(denominator) + " " + output_file + " " + kmer_counts_file
    print(">>COMANDO KMC_DUMP LANCIATO: " + cmd + "\n")
    print(cmd)
    os.system(cmd)
    end_dump1 = time.time()
    dump1_time = end_dump1 - start_dump1
    print("Total_time= " + str(counting_time + dump1_time + (end_sample - start_sample)))
    of.write("Total_time= " + str(counting_time + dump1_time + (end_sample - start_sample)) + " \n")
    of.close()
    return kmer_counts_file, sample_size #agg. sample_size in uscita per correggere abundance
