print ([[

The option -keys allows to extract substrings or sequences from the given 
sequence file or fasta index.
The substrings to be extracted are specified in a key file given
as argument to this option. The key file must contain lines of the form

k

or

k i j

where k is a string (the key) and the optional i and j are positive integers 
such that i<=j. k is the key and the optional numbers i and j
specify the first position of the substring and the
last position of the substring to be extracted. The positions are counted
from 1. If k is identical to the string between the first
first and second occurrence of the symbol | in a fasta header, then
the fasta header and the corresponding sequence is output. If i and j  are
both specified, then the corresponding substring is shown in
fasta format. In the latter case the header of the fasta formatted sequence in 
the output begins with

>k i j

followed by the original original fasta header.
Duplicated lines in the input file lead to only one sequence in the output.
The sequences are output according to the order in the original sequence
files.  The formatting of the output can be controlled by the options
-width, -o, -gzip, and -bzip2. If end of the argument list only contains
one filename, say idx, then it is checked if there is a file 
idx.kys. This makes up the fasta index, which is contructed by
calling the gt suffixerator as follows:

gt suffixerator -protein -ssp -tis -des -sds -kys -indexname fastaindex -db inputfile1 [inputfile2 ..]

This reads the protein sequence files given to the option -db and creates
several files:
 - a file fastaindex.ssp specifying the sequence separator positions.
 - a file fastaindex.esq representing the sequence.
 - a file fastaindex.des showing the fasta headers line by line.
 - a file fastaindex.sds giving the sequence header delimter positions.
 - a file fastaindex.kys containing the keys in the fasta files.

Keys of the form 

>k i j

in the file specifying which substrings are extracted can only be specified
when extracting sequences from fasta files, but not from a fasta index.
]])
