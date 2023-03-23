
# biodb

Quick CLI tool to help save disk usage, and easily retrieve the necessary sequences.

## RocksDB Keys

There are two different key formats within the RocksDB:

hmmsearch:XXX - Where xxx is the id# of the hmm search, starting at 0.  Contains JSON objects including blast results.
HEADER_TRANSLATE - The base header including any frame info.  All spces are replaced with an underscore (_).  For example:
    NODE_100000_length_301_[revcomp]:[translate(2)]

## Usage

Run "./biodb --help" for list of all available options.  Global flags that should be included with every command are:

    --action [-a] - The action to perform, one of the below.
    --input [-i] - The input directory (ie. orthograph_results/species_name).

#### upgrade-db

First, upgrade existing SQLite database.  Please note, this will not delete the SQLite database from the 
machine, but upon successful run it may be manually deleted if desired.
    ./biodb -a upgrade-db -i Syrphidae/orthograph_results/Acroceridae/SRR6453524.fa

This will create a ~/rocksdb/ sub-directory within the input directory as well.  If you need to reset the RocksDB database, simply delete that directory.

#### get-sequence

Next, get sequences.  For example, full est sequence of base header:
    ./biodb -a get-sequence -i -i Syrphidae/orthograph_results/Acroceridae/SRR6453524.fa -h "NODE_2347144_length_252"

Get nt sequence, coords 2 - 63:
    ./biodb -a get-sequence -i -i Syrphidae/orthograph_results/Acroceridae/SRR6453524.fa -h "NODE_2347144_length_252" -c 2-63

Get aa sequence, coords 2 - 41.  This is simply the nt sequence translated via the fasatranslate algorithm:
    ./biodb -a get-sequence -i -i Syrphidae/orthograph_results/Acroceridae/SRR6453524.fa -h "NODE_2347144_length_252" -c 2-41 -t aa

#### get-hmmsearches

Retrieve a list of hmm searches:
    ./biodb -a get-hmmsearches 

Retrive the next 500 hmm searches starting at the 1000 result:
    ./biodb -a get-hmmsearches -s 1000 -l 500

#### get-hmmsearch

Only applicable if you know the exact hmm search id# you wish to obtain.  For example, id# 9116:
    ./biodb -a get-hmmsearch -h 9116


## Todo

Internalize exonerate, so the sequences provided are already processed through the exonerate algorithm.

