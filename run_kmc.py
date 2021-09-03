import os


def run_kmc(fastqfile, dataset, datasets_size, k, working_dir,  exec_fold, thread, memory):
    temp_folder=os.path.join(working_dir,"work_dir_exact/")
    if not os.path.isdir(temp_folder):
        try:
            print("Attempting create " + str(temp_folder))
            os.mkdir(temp_folder)
        except OSError as error:
            print("Creation problem: " + str(error))
            #exit(1)
    kmcexec = os.path.join(exec_fold, "SPRISS/bin/kmc")
    outputfilekmc =os.path.join(temp_folder,dataset)
    cmd = kmcexec + " -v -k" + str(k) + " -cs" + str(datasets_size) + " -m" + str(memory) + " -ci1 -t" + str(thread) + " " + fastqfile + " "  + outputfilekmc +" " + temp_folder# " work_dir_exact/"
    #cmd = exec_fold+"SPRISS/bin/kmc -v -k" + str(k) + " -cs" + str(datasets_size) + " -m"+str(memory) +" -ci1 -t"+str(thread) +" " + dataset + " " + working_dir + " "+temp_folder
    print(cmd)
    print("Starting step 1 kmc")
    os.system(cmd)
    kmer_counts_file = working_dir + "/KmersCount_"+dataset+"_kmc_k"+str(k)+"_exact_ci1.txt"
    kmcdumpexec = os.path.join(exec_fold, "SPRISS/bin/kmc_dump")
    cmd = kmcdumpexec +" -ci1 " + outputfilekmc + " " + kmer_counts_file
    print(cmd)
    print("Starting step 2 kmcdump")
    os.system(cmd)
    return(kmer_counts_file)
    #<carica file fasta considerato pescando da directory origine
    #imposta file in directory uscita(verrà creata dal run_all in create_dirs)
    #<rilevazione tempo di inizio tempo KMC
    #esegui kmc   e kmcdump
    #metti i file con kmers e frequenze
    #rilevazione tempo di fine KMC