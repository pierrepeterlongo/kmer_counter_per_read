# Kmer Counter Per Read

## Description
Kmer Counter Per Read is a tool designed to count k-mers in sequencing reads. 

Input: 
- a read set (fasta, fastq) [.gz|.zstd]
- k (kmer size)
- t threshold 

Output: 
For each input read, prints the header + kmers plus their counts for  kmers occurring at least $t$ times.

## Install
To install the tool, run:

```bash
cargo install --path .
```

## Test
```bash
kmer_counter_per_read --input example/head_ERR13885951.fa.gz -k 9 -t 2
```

