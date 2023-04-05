### Usage

Duplicated sequences within the same taxon are discarded, headers replaced. Outputs to STDOUT. 
`fastataxofilter -i input.fasta.gz -r regex.tsv -l 1 -s stats.json -g > output.fasta.gz`

For more details:
`fastataxofilter -h`

### Input / output

Input can be fasta or gziped fasta, will be autodetected (filename needs to end with 'gz').
Output can be fasta or gziped fasta, use flag -g/--gzip-output for compressed output.

### Regex file

The regex configuration file has one entry per line. Each entry has 1 or 2 tab-separated regular expressions (see No replace section below).
!! Need the tab character, not 4 spaces !!
Do not use the tab character in any regex, it is reserved as the delimiter. Use `\t`.
Those regex syntax is what is accepted by the regex rust crate, which is basically standard Perl syntax but without look-ahead and backreferences.
See [https://docs.rs/regex/1.7.3/regex/]

Regex 1 needs to match the full header of each fasta sequence. It needs to define the capture groups that are used in regex 2. Regex 1 also needs to define a capture group named 'taxo'. The 'taxo' capture group will be hashed with the sequence to detect duplicates. It needs to be a string that uniquely defines a taxon.
A capture group is defined by enclosing a part of the regex in parentheses `(some_regex)`. A named capture group is defined by `(?P<group_name>some_regex)`.

Regex 2 is the format of the new header. It uses capture groups specified in regex 1. Each capture group is referred to by ${n} where n is a number starting at 1.

Sample line:
(.+)###(?P<taxo>.+)	${1} ${2}

You do not need to worry about the > starting character of the header line, as it is automatically discarded on parsing and re-added on writing. It is not considered to be part of the header string for our purpose.

#### Select a line

To chose a regex from the file, use the `-l, --select-regex` argument with the line number, starting with 1.
You can use comment lines in this file, the line must start with #. Empty lines are allowed.

#### No replace

If no replace is necessary, use the `-n/--no-replace` flag.
In this case, the 2nd regex is not used and can be omitted from the regex file. If a 2nd regex is present, it is simply ignored.

### Misc

Stats include total entries parsed and duplicated entries count, in json format.
The complexity of regex 1 heavily influences the performance of the script.
`-n/--no-replace` improves performance.
