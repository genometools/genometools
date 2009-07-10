print ([[

The option -keys allows to extract substrings from the sequence file.
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
-width, -o, -gzip, and -bzip2.]])
