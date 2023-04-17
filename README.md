# FastaTaxoFilter

Duplicated sequences within the same taxon are discarded, headers replaced. Outputs to STDOUT.
Several filter options are available.

### Usage

`fastataxofilter -i input.fasta.gz -r regex.tsv -l 1 -s stats.json -g > output.fasta.gz`

For more details: `fastataxofilter -h`

### Input / output

Input can be fasta or gziped fasta, will be autodetected (filename needs to end with 'gz').
Default output is fasta. Use flag `-g, --gzip-output` for gziped fasta.

### What is this script going to do ?

When using `-e, --explain`, the script will only process the first record and output a summary of what happened, what the regex captured and what the new replaced header looks like. Use it to debug your regexes.

### Regex file

The regex configuration file has one entry per line. Each entry has 1 or 2 tab-separated regular expressions (see No replace section below).

!! Need the tab character, not 4 spaces !!

Do not use the tab character in any regex, it is reserved as the delimiter. Use `\t`.
Those regex syntax is what is accepted by the regex rust crate, which is basically standard Perl syntax but without look-ahead and backreferences.
See the [regex crate doc](https://docs.rs/regex/1.7.3/regex/).

Regex 1 needs to match the full header of each fasta sequence. It needs to define the capture groups that are used in regex 2. Regex 1 also needs to define a capture group named 'taxo'. The 'taxo' capture group will be hashed with the sequence to detect duplicates. It needs to be a string that uniquely defines a taxon.
A capture group is defined by enclosing a part of the regex in parentheses `(some_regex)`. A named capture group is defined by `(?P<group_name>some_regex)`.

Regex 2 is the format of the new header. It uses capture groups specified in regex 1. Each capture group is referred to by `${n}` where `n` is a number starting at 1.

Sample line:

`(.+)###(?P<taxo>.+)	${1} ${2}`

You do not need to worry about the > starting character of the header line, as it is automatically discarded on parsing and re-added on writing. It is not considered to be part of the header string for our purpose.

#### Select a line

To chose a regex from the file, use the `-l, --select-line` argument with the line number, starting with 1.
You can use comment lines in this file, the line must start with #. Empty lines are allowed.

#### No replace

If no replacement is necessary, use the `-n, --no-replace` flag.
In this case, the 2nd regex is not used and can be omitted from the regex file. If a 2nd regex is present, it is simply ignored.

#### Ignoring sequences

By default, the program will panic and throw an error if a header cannot match the given regex. `-y, --ignore-no-match` lets you discard sequences whose headers do not match.

`-u, --ignore-empty` will discard empty sequences.

### One sequence per taxon

Option `-f, --first-seq-only` will only process the first sequence it finds for each taxon.

### Misc

The complexity of regex 1 heavily influences the performance of the script.
`-n` and `-f` improve performance.

`-s, --stats-file` will output a json dict at the given path. It includes different stats about the amount of duplicates found, total sequences processed, the arguments used for the script...
Sample stats file:

```
{
  "Duplicate counts per taxon": {
    "root_1;Eukaryota_2759;Arthropoda_6656;Branchiopoda_6658;Diplostraca_84337;Chydoridae_77713;Leydigia_527154;Leydigia_lousi_535388": 1,
    "root_1;Eukaryota_2759;Arthropoda_6656;Hexanauplia_72037;Calanoida_6833;Diaptomidae_236535;Eodiaptomus_903836;Eodiaptomus_wolterecki_1087294": 2,
    "root_1;Eukaryota_2759;Arthropoda_6656;Insecta_50557;Coleoptera_7041;Apionidae_122732;Protapion_202144;Protapion_fulvipes_202145": 2,
    "root_1;Eukaryota_2759;Arthropoda_6656;Insecta_50557;Coleoptera_7041;Buprestidae_50527;Agrilus_195164;Agrilus_cuprescens_324824": 1,
    "root_1;Eukaryota_2759;Arthropoda_6656;Insecta_50557;Coleoptera_7041;Carabidae_41073;Carabus_41074;Carabus_glabratus_1149273": 1,
    
    ...

    "root_1;Eukaryota_2759;Euglenozoa_33682;Euglenida_3035;Euglenales_86650;Euglenaceae_1131320;Euglena_3038;Euglena_gracilis_3039": 1
  },
  "Duplicated entries": 544,
  "Empty sequence": 0,
  "Output entries": 499456,
  "Script": {
    "Args": {
      "--explain": false,
      "--first-seq-only": false,
      "--gzip-output": true,
      "--ignore-empty": false,
      "--ignore-no-match": false,
      "--input": "my_input.fasta.gz",
      "--no-replace": false,
      "--regex-file": "regex.tsv",
      "--select-line": 3,
      "--stats-file": "stats.json"
    },
    "Full cmd": "target\\release\\fastataxofilter.exe -i my_input.fasta.gz -r regex.tsv -l 3 -s stats.json -g",
    "Version": "0.5.0"
  },
  "Taxons with at least 1 duplicate": 285,
  "Time": "12.5462599s",
  "Total entries": 500000,
  "Unmatched sequence": 0
}
```