use bio::io::{fasta::{Reader as FastaReader, Record, Writer as FastaWriter}};
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use regex::{Regex, Captures};
use std::{path::PathBuf, fs::{self, File}, collections::{HashSet, HashMap}, borrow::Cow};
use std::time::{Duration, Instant};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::mpsc::channel;
use clap::Parser;
use std::io;
use serde_json::json;


#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Input FASTA file
    #[arg(short = 'i', long = "input")]
    input_path: std::path::PathBuf,
    /// Input regex TSV file
    #[arg(short = 'r', long = "regex-file")]
    regex_path: std::path::PathBuf,
    /// Line number to select in the regex file, starts at 1
    #[arg(short = 'l', long = "select-line")]
    regex_line: usize,
    /// Optional output stats/log json file with a bunch of metrics
    #[arg(short = 's', long = "stats-file")]
    stats_path: Option<std::path::PathBuf>,
    /// Flag to compress output to gzip
    #[arg(short = 'g', long = "gzip-output")]
    is_gzip_output: bool,
    /// Flag to avoid doing header string replace
    #[arg(short = 'n', long = "no-replace")]
    avoid_replace: bool,
    /// Flag to only process the 1st record and outputs explanation on STDERR
    #[arg(short = 'e', long = "explain")]
    explain: bool,
    /// Flag to ignore and discard records that do not match the given regex
    #[arg(short = 'y', long = "ignore-no-match")]
    ignore_no_match: bool,
    /// Flag to ignore and discard records with empty sequence
    #[arg(short = 'x', long = "ignore-empty")]
    ignore_empty: bool,
    /// Flag to only use the first sequence for each taxon
    #[arg(short = 'f', long = "first-seq-only")]
    first_seq_only: bool
}


// TODO:
// gzip input (auto detect)         --      --      --      --      --      --      --      --      --  DONE
// + output gzip (command option),  --      --      --      --      --      --      --      --      --  DONE
// + replace is option      --      --      --      --      --      --      --      --      --      --  DONE
// input from STDIN ? NOPE
// precise > header stuff in readme --      --      --      --      --      --      --      --      --  DONE
// allow 2 regex even with no-replace       --      --      --      --      --      --      --      --  DONE
// add --explain thing      --      --      --      --      --      --      --      --      --      --  DONE
// count comment and empty lines in regex file      --      --      --      --      --      --      --  DONE
// ignore records that don't match regex, no written to output, add to stats, option        --      --  DONE
// count records with empty sequence to stats       --      --      --      --      --      --      --  DONE
// option to ignore empty sequence                  --      --      --      --      --      --      --  DONE
// ah en fait dans les stats tu peux meme mettre le path du fichier input, le path du ficher output,
// la version du software, les paranetres et leur values            --      --      --      --      --  DONE
// stats for sequences per taxons for duplicate only        --      --      --      --      --      --  DONE
// pretty print json        --      --      --      --      --      --      --      --      --      --  DONE
// option to only consider the first sequence for each taxon        --      --      --      --      --  DONE
// git the damn thing       --      --      --      --      --      --      --      --      --      --  DONE
// update readme


fn main() {
    let args = Cli::parse();

    let selected_regex_line = parse_regex_file(
        &args.regex_path, 
        args.regex_line, 
        args.avoid_replace);

    // manage input
    let path_string = args.input_path.to_str().expect("Could not parse input file path");
    let file = File::open(path_string).expect("Problem opening input file");
    let bufreader: Box<dyn BufRead> = if path_string.ends_with("gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    }
    else {
        Box::new(BufReader::new(file))
    };
    let fastareader = FastaReader::from_bufread(bufreader);
    let mut records = fastareader.records();
    
    //let mut seq_map: AHashMap<String, AHashSet<Vec<u8>>> = AHashMap::default(); 
    let mut aggreg_hashset: HashSet<String> = HashSet::new();

    // Create a channel to communicate with the writer thread
    let (tx, rx) = channel::<Record>();

    // Spawn a new thread to handle the writing
    let handle = std::thread::spawn(move || {
        // manage output
        let writer: Box<dyn Write> = if args.is_gzip_output {
            Box::new(GzEncoder::new(io::stdout(), Compression::default()))
        }
        else {
            Box::new(io::stdout())
        };
        let mut fastawriter = FastaWriter::from_bufwriter(BufWriter::new(writer));

        while let Ok(record) = rx.recv() {
            //println!("inside thread: {}", &record);
            fastawriter.write_record(&record)
                .expect(format!("Error while writing record with id: {}", record).as_str());
        }
        fastawriter.flush().expect("Error while flushing fasta writer at the end");
    });

    // collect duplicate stats for each taxons
    let mut duplicate_seq_map: HashMap<String, u32> = HashMap::new();
    let mut global_taxo_set: HashSet<String> = HashSet::new();

    //let mut record = fasta::Record::new();
    let mut total_count = 0;
    let mut identical_count = 0;
    let mut empty_seq_count = 0;
    let mut unmatched_record_count = 0;
    let mut output_entries_count = 0;
    let start = Instant::now();
    /*let mut sum_duration2 = Duration::default();
    let mut sum_duration3 = Duration::default();
    let mut sum_duration4 = Duration::default();
    let mut sum_duration5 = Duration::default();*/
    while let Some(Ok(record)) = records.next() {
        //let loop_start = Instant::now();

        // to check or not to check ? seemingly no performance impact
        record.check().expect(format!("Error while checking fasta record with rust-bio lib:\n{}", record).as_str());
        
        // DEBUG stuff
        /*if total_count % 100_000 == 0 {
            match args.stats_path {
                Some(ref path) => {
                    let content = format!(
                    "{} {}\ndur2 {:?}\ndur3 {:?}\ndur4 {:?}\ndur5 {:?}\n", 
                    total_count, identical_count, sum_duration2, sum_duration3, sum_duration4, sum_duration5);
                    fs::write(path, content).expect(format!("Error while writing stats file").as_str());
                },
                None => {}
            }
        }*/
        /*if total_count == 500_000 {
            break;
        }*/

        // ensure the header we use contains all the header of the entry, as rust-bio splits the header into id and desc
        let header = match record.desc() {
            Some(desc) => [record.id(), desc].join(" "), // header was split with a space between id and description
            None => record.id().to_string() // no whitespace in header, only id
        };
        let sequence = record.seq().to_vec();

        if sequence.is_empty() {
            empty_seq_count += 1;
            if args.ignore_empty {
                total_count += 1;
                continue;
            }
        }

        //let part_duration2 = loop_start.elapsed();
        //sum_duration2 += part_duration2;


        // capture taxonomy to use as key
        let (re_match, re_replace) = &selected_regex_line;
        let Some(capture_result) = re_match.captures(&header) else {
            unmatched_record_count += 1;
            if args.ignore_no_match {
                total_count += 1;
                continue;
            }
            else {
                panic!("Could not match given regexp: {}\nProblematic header:\n\n{}\n\n", re_match.as_str(), &header);
            }
        };
        let taxo_key = capture_result.name("taxo")
            .expect(format!("Could not match taxonomy from given regexp: {}\nNeed a capture group named 'taxo'. Problematic header:\n\n{}\n\n", 
                re_match.as_str(), &header).as_str()).as_str().to_string();
        if args.first_seq_only {
            if global_taxo_set.contains(&taxo_key) {
                total_count += 1;
                continue;
            }
            else {
                global_taxo_set.insert(taxo_key.clone());
            }   
        }
        


        //let part_duration3 = loop_start.elapsed();
        //sum_duration3 += part_duration3;

        // replace header according to provided regex
        let new_header = if args.avoid_replace {
            Cow::from(&header)
        }
        else {
            re_match.replace(&header, re_replace.as_ref().expect("No replace string found"))
        };

        //let part_duration4 = loop_start.elapsed();
        //sum_duration4 += part_duration4;

        let aggreg_key = [taxo_key.clone(), String::from_utf8_lossy(&sequence).to_string()].join("");

        // exit after 1st record
        if args.explain {
            print_explain(capture_result, &header, &new_header.to_string());
            break;
        }

        // check for identical sequences
        if !aggreg_hashset.contains(&aggreg_key) {
            aggreg_hashset.insert(aggreg_key.clone());
            let new_record = Record::with_attrs(&new_header, None, &sequence);
            //fastawriter.write_record(&new_record)
             //   .expect(format!("Error while writing record with id: {}", new_header).as_str());
            tx.send(new_record)
                .expect(format!("Error while sending fasta record to writer thread with header: {}", new_header).as_str());
            output_entries_count += 1;
        
        }
        else { // identical seq FOUND
            identical_count += 1;

            // update duplicate seq map
            duplicate_seq_map.entry(taxo_key)
                .and_modify(|e| { *e += 1 })
                .or_insert(1);
        }    

        //let part_duration5 = loop_start.elapsed();
        //sum_duration5 += part_duration5;
        
        total_count += 1;
    }

    // signal that channel is done, actually not sure if necessary
    drop(tx);
 
    // Let last stuff be written, needed to avoid corrupted gzip
    match handle.join() {
        Ok(_) => {},
        Err(err) => panic!("Problem with thread termination\n{:?}", err)
    }

    // write the json report if user chose it
    if let Some(path) = &args.stats_path {
        // impossible to get output redirection
        let full_cmd: Vec<String> = std::env::args().collect();
        let full_cmd_string = full_cmd.join(" ");
        let duration = start.elapsed();
        let stats_content_json = json!({
            "Duplicated entries": identical_count,
            "Empty sequence": empty_seq_count,
            "Unmatched sequence": unmatched_record_count,
            "Output entries": output_entries_count,
            "Total entries": total_count,
            "Time": format!("{:?}", duration),
            "Script": {
                "Version": env!("CARGO_PKG_VERSION"),
                "Full cmd": full_cmd_string,
                "Args": {
                    "--input": args.input_path,
                    "--regex-file": args.regex_path,
                    "--select-line": args.regex_line,
                    "--stats-file": args.stats_path,
                    "--gzip-output": args.is_gzip_output,
                    "--no-replace": args.avoid_replace,
                    "--explain": args.explain,
                    "--ignore-no-match": args.ignore_no_match,
                    "--ignore-empty": args.ignore_empty,
                    "--first-seq-only": args.first_seq_only
                }
            },
            "Taxons with at least 1 duplicate": duplicate_seq_map.len(),
            "Duplicate counts per taxon": duplicate_seq_map
        });
        fs::write(path, serde_json::to_string_pretty(&stats_content_json)
            .expect("Error while prettyfying json"))
            .expect(format!("Error while writing stats file").as_str());
    }

}


fn parse_regex_file(path: &PathBuf, selected_line: usize, avoid_replace:bool) -> (Regex, Option<String>) {
    let data = fs::read_to_string(path.clone())
        .expect(format!("Error while opening regex file {:?}", path).as_str());
    let mut result: (Regex, Option<String>) = (Regex::new("").unwrap(), None);
    let mut line_count = 1;
    let expected_entries_per_line = if avoid_replace {1} else {2};
    let mut line_was_found_and_ok = false; // did the parsing go well

    // sanity check
    if selected_line < 1 { panic!("-l/--regex-line argument should start at 1")}

    for line in data.split("\n").into_iter() {
        // comment or empty, or a line with some data, but not the one selected, skip
        if line.starts_with("#") || line.is_empty() || line_count != selected_line {
            line_count += 1;
            continue;
        }
        // at this point we are at the desired line


        let regs: Vec<&str> = line.split("\t").collect();
        if regs.len() < expected_entries_per_line {
            panic!("Could not parse correctly the following line from the regex file:
\n{}\n
Expected {} entries, found {}.
Ensure that the 2 entries are separated by an actual tab character, not spaces.", line, expected_entries_per_line, regs.len());
        }
        else if regs.len() > 2 {
            panic!("Could not parse correctly the following line from the regex file:
\n{}\n
Expected {} entries, found {}.
The line needs to have 2 entries separated by only 1 tab character.
Ensure that you have used only 1 tab character on the line.
If you need to use tabs in a regex, use '\\t'.", line, expected_entries_per_line, regs.len());
        }

        let regex_match = regs.get(0).unwrap();
        let compiled_re = Regex::new(&regex_match)
            .expect(format!("Error while compiling first regex at the following line from the regex file:\n{}\n",
                                line).as_str());
        
        if expected_entries_per_line == 2 {
            let regex_replace = regs.get(1).unwrap();
            // the replace string isn't used as a compiled regex, this juste allows us to check its validity
            let _compiled_replace = Regex::new(&regex_replace)
                .expect(format!("Error while compiling second regex at the following line from the regex file:\n{}\n",
                                    line).as_str()); 
            result = (compiled_re, Some(regex_replace.to_string()));
        }                        
        else {
            result = (compiled_re, None);
        }
        line_was_found_and_ok = true;
        break; // only 1 line needs parsing

    }

    if line_was_found_and_ok {
        return  result;
    }
    else {
        panic!("Nothing parsed, regex file seems to be empty, or -l/--regex-line argument is invalid or points to a comment or empty line");
    }
}

/// print explanations of the regex on stderr, with nice formatting
fn print_explain(capture_result: Captures, original_header: &String, new_header: &String) {
    let taxo_start = capture_result.name("taxo")
        .expect("Could not retrieve taxo start for explanation").start();
    let taxo_end = capture_result.name("taxo")
        .expect("Could not retrieve taxo end for explanation").end();
    let highlight = "^".repeat(taxo_end - taxo_start);
    let start_spacing = " ".repeat(taxo_start);
    let highlight_line = [start_spacing, highlight].join("");

    
    eprintln!("Original header and taxo capture group (chars {} to {}):", 
        taxo_start, taxo_end);
    // format things on 80 chars wide lines
    let width = 80;
    for i in 0..(original_header.len() / width + 1) {
        //eprintln!("{} {} {}", i, (i * width), ((i + 1) * width));
        let line_start_char_pos = i * width;
        let line_end_char_pos = (i + 1) * width;
        eprintln!("{: >4}-{: >4}|{}",line_start_char_pos, std::cmp::min(original_header.len(), line_end_char_pos),
            &original_header[line_start_char_pos..std::cmp::min(original_header.len(), line_end_char_pos)]);
        if line_start_char_pos < highlight_line.len() {
            eprintln!("          {}", &highlight_line[line_start_char_pos..std::cmp::min(highlight_line.len(), line_end_char_pos)]);
        }
        // highlight line is not filled to match the length of the header,
        // so it is always shorter. Here nothing remains in highlight line.
        else {
            eprintln!("");
        }
    }
    
    eprintln!("\nList of captured groups, enclosed with \":");
    let mut cap_count = 0;
    for cap_group in capture_result.iter() {
        if cap_count == 0 { // ignore 1st group as it is entire match
            cap_count += 1;
            continue;
        }
        if let Some(m) = cap_group {
            eprintln!("|{}: \"{}\"", cap_count, m.as_str());
            cap_count += 1;
        }
    }

    eprintln!("\nNew header:\n{}", &new_header);
}
