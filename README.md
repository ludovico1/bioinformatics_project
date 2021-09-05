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

Per ulteriori parametri eseguire `python3 run_All.py -h` per ottenere l'help
