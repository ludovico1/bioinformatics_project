# Comparison of sampling techniques for k-mer count approximation 

## Prerequisiti
- Clonare la repository di Spriss e Sakeima reperibile al link: https://github.com/VandinLab/SPRISS e installare secondo le istruzioni presenti in essa
- Installare Kmergenie che si trova al link http://kmergenie.bx.psu.edu/ e copiare la cartella ntCard nella directory in cui si trovano Sakeima e Spriss
- isntallare pyfastx usando il comando  ` pip install pyfastx`

## Utilizzo

Per avviarelo usare il comando `python3 run_All.py` seguito dai segueti parametri obbligatori:
-  `--output` percoso dove salvare i risultati della run
-  `--verbose` per avere maggior dettaglio nei log
-  `-f` posizione del file fastq 
-  `-k` valore di K
-  `-D` nome dataset (viene utilizzato per nominare la directory)
-  `-p` directory di installazione contenete Spriss, Sakeima, KMC e ntCard
-  `-t` parametro theta

`--output /home/Software/Algobio/Gruppo1/Datasets/ --verbose -f /home/Software/Algobio/Gruppo1/Datasets/Rhodobacter_frag.fastq 
-k 31 -D Rhodobacter_frag_k31_t5e-6 -p /home/Software/Algobio/SPRISS -t 5e-6`
## Output

Verrà creata una directory principale con il nome del dataset. All'interno sarà presente la 
directory Graph con tutti i risultati prodotti. Sempre nella stessa directory per ogni tool
viene creata una directory contenenti file prodotti dal tools e elaborazioni del nostro tools
Il file log_exec_result_[nomedataset].txt conterrà tutti i tempi rilevati.

Sono previste 6 step:
- 1 Conteggio KMC
- 2 Conteggio Sakeima
- 3 Conteggio SPRISS
- 4 Conteggio Kmergenie-ntCard
- 5 Elaborazione e creazione Grafici

Nella directory dove viene creata la directory principale del dataset viene aggiornato 
progress_[datasetname].txt file per gestire le ripartenze. Deve contenere l'ultima fase 
da considerare completata, riparte da quella successiva.

Per ulteriori parametri eseguire `python3 run_All.py -h` per ottenere l'help
