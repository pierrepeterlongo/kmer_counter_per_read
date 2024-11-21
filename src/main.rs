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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_kmers() {
        let sequence = "ACGTACGT";
        let k = 3;
        let result = count_kmers(sequence, k);
        let mut expected = HashMap::new();
        expected.insert("ACG".to_string(), 2);
        expected.insert("CGT".to_string(), 2);
        expected.insert("GTA".to_string(), 1);
        expected.insert("TAC".to_string(), 1);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_count_kmers_with_non_ascii() {
        let sequence = "ACGTACGT😊";
        let k = 3;
        let result = count_kmers(sequence, k);
        let mut expected = HashMap::new();
        expected.insert("ACG".to_string(), 2);
        expected.insert("CGT".to_string(), 2);
        expected.insert("GTA".to_string(), 1);
        expected.insert("TAC".to_string(), 1);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_count_kmers_with_short_sequence() {
        let sequence = "AC";
        let k = 3;
        let result = count_kmers(sequence, k);
        let expected: HashMap<String, usize> = HashMap::new();
        assert_eq!(result, expected);
    }

//     #[test]
//     fn test_process_file() -> std::result::Result<(), anyhow::Error> {
//         let filename = "example/head_ERR13885951.fa.gz";
//         let kmer_size = 9;
//         let threshold = 2;

//         let mut cmd = assert_cmd::Command::cargo_bin("kmer_counter_per_read")?;

//         cmd.args([
//             "--input",
//             &format!("{}", filename),
//             "-k",
//             &format!("{}", kmer_size),
//             "-t",
//             &format!("{}", threshold),
//         ]);
//         // check that the functions prints the following output
//         let truth: &'static str = ">gnl|SRA|ERR13885951.1.100055195-06b3-41c2-b713-78618bd15975/a9f35072049269c4116497fa804f5096351ed9d6 Biological (Biological)
// \t CAGATGTGT 2
// \t AAAAAAAAA 7
// >gnl|SRA|ERR13885951.2.10010f4b6-ae0b-4fea-a588-ec21325ceeae/a9f35072049269c4116497fa804f5096351ed9d6_1 Biological (Biological)
// \t CAAGGTACT 2
// \t CCCTCCCCC 2
// \t GCTCACAGC 2
// \t CCCCAAGTC 2
// \t CCCCTCCCC 2
// \t CAGCTGGCC 2
// \t TCTAAACGT 2
// \t AAACCTCTG 2
// \t CTAAACGTC 2
// \t AAGGTACTT 2
// \t ACAGCTGGC 2
// \t CAAACCTCT 2
// \t TAAACGTCT 2
// \t CCCAAGTCC 2
// \t AGTGGGCAA 2
// >gnl|SRA|ERR13885951.3.10010f4b6-ae0b-4fea-a588-ec21325ceeae/a9f35072049269c4116497fa804f5096351ed9d6_2 Biological (Biological)
// \t AGGTCCTGG 2
// \t AAAAAAAAA 51";
//         let assert = cmd.assert();
//         assert.success().stdout(truth);
//         Ok(())
//     }
}