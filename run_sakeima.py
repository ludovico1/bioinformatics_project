import os
import time


def run_sakeima(fastqfile, dataset, datasets_size, tot_reads, avg_read_length, max_read_length, k, delta, thetastr,
                working_dir, exec_fold, thread, memory):
    outlogfile = os.path.join(working_dir, "log_sakeima_approx_time_" + dataset + "_" + str(thetastr) + ".txt")

    of = open(outlogfile, 'w')

    theta = float(thetastr)

    print(dataset)
    print(str(theta))
    of.write("Theta: "+str(theta) + " \n")
    of.write("dataset size: " + str(datasets_size) + " \n")

    start = time.time()

    output_file = os.path.join(working_dir, "KmersCount_tmp_" + dataset + "SAKEIMA_k" + str(k) + "_theta" + str(thetastr) + ".txt")
  #  output_file2 = working_dir + "/KmersCount_tmp_" + dataset + "SAKEIMA_k" + str(k) + "_theta" + str(thetastr) + ".txt"

    exec_jelly = os.path.join(exec_fold, "SAKEIMA/")
    cmd = "python3 run_SAKEIMA_2.py -v true -k " + str(k) + " -db " + fastqfile + " -o " + output_file + " -thr 1 -t " + str(
        theta) + " -dt " + str(datasets_size) + " -w " + working_dir + " -bin " + exec_jelly +" -d " + str(delta)
    #esecutione senza Theta di sakeima
    #cmd = "python3 run_SAKEIMA_2.py -v true -k " + str(k) + " -db " + fastqfile + " -o " + output_file + " -thr 1 -t  -dt " + str(datasets_size) + " -w " + working_dir + " -bin " + exec_jelly + " -d " + str(delta)
    print(cmd)
    os.system(cmd)
    end = time.time()
    of = open(outlogfile, 'a')
    print("Count_time: " + str(end - start))
    of.write("Count_time: " + str(end - start) + " \n")
    # write frequent kmers
    fallkamer = os.path.join(working_dir,"KmersCount_tmp_" + dataset + "SAKEIMA_k" + str(k) + "_theta" + str(thetastr) + ".txt")
    allkmers = open(fallkamer, 'r')
    ffreqkmers = os.path.join(working_dir,"KmersCount_" + dataset + "SAKEIMA_k" + str(k) + "_theta" + str(thetastr) + ".txt")
    freqkmers = open(ffreqkmers, 'w')
    line = allkmers.readline()
    while (line):
        splitted = line.split(';')
        kmer = splitted[1]
        support = float(splitted[0])
        unbiased_frequency = float(splitted[2])
        freqkmers.write(kmer + " " + str(support) + " " + str(unbiased_frequency) + " \n")
        line = allkmers.readline()
    allkmers.close()
    freqkmers.close()
    kmer_counts_file = os.path.join(working_dir, "KmersCount_" + dataset + "SAKEIMA_k" + str(k) + "_theta" + str(thetastr) + ".txt")
    #return os.path.join(working_dir+"/", "/KmersCount_" + dataset + "SAKEIMA_k" + str(k) + "_theta" + str(thetastr) + ".txt")


    ## aggiunta per il sampling size
    #leggo dal file:
    ssfile = os.path.join(working_dir, "work_dir","samplesize.txt")

    try:
        with open(ssfile, "r") as pf:
            # print("Progress file opened")
            ratiosize = float(pf.readline())
            sample_size= float(pf.readline())
            # da cancellare ------------------------------------------------
            # print(" Progress Step found:" + str(pv))
            pf.close()

    except IOError:
        print("Progress File not accessible")
        pv = -1
        exit(-1)
    of.write("sample size: " + str(sample_size) + " \n")
    of.close()
    #ratiosize = sample_size
    return kmer_counts_file, ratiosize
