use std::collections::HashMap;
use std::error::Error;
use clap::Parser;
use fxread::initialize_reader;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file path (FASTA or FASTQ, gzipped or zstd supported)
    #[arg(short, long)]
    input: String,

    /// K-mer size
    #[arg(short, long, default_value_t = 5)]
    k: usize,

    /// Minimum occurrence threshold for reporting
    #[arg(short, long, default_value_t = 2)]
    threshold: usize,
}

fn count_kmers(sequence: &str, k: usize) -> HashMap<String, usize> {
    let mut kmer_counts = HashMap::new();

    if sequence.len() < k {
        return kmer_counts;
    }

    for window in sequence.as_bytes().windows(k) {
        if window.iter().all(|&b| b.is_ascii()) {
            let kmer = std::str::from_utf8(window).unwrap().to_string();
            *kmer_counts.entry(kmer).or_insert(0) += 1;
        }
    }

    kmer_counts
}


fn process_file(filename: &str, kmer_size: usize, threshold: usize) -> Result<(), Box<dyn Error>> {
    let mut reader = initialize_reader(&filename)?;
    loop {
        let Some(mut record) = reader.next_record()? else {
            break;
        };
        record.upper();
        let acgt_sequence = record.seq();

        // for each kmer of the sequence, insert it in the kmer_set
        if acgt_sequence.len() < kmer_size {
            continue;
        }
        let kmers = count_kmers(std::str::from_utf8(acgt_sequence)?, kmer_size);
        println!(">{}", std::str::from_utf8(record.id())?);
        for (kmer, count) in kmers.iter().filter(|&(_, &count)| count >= threshold) {
            println!("\t {} {}", kmer, count);
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    process_file(&args.input, args.k, args.threshold)
}