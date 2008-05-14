print ([[

The option -ginum allows to extract substrings from the sequence file.
The substrings to be extracted are specified in a query file given
as argument to this option. The query file must contain lines of
the form

g i j

where g, i, and j are all positive integers or i is 1 and j is 0. g is 
interpreted as the gi number, i is the first position of the substring and j the
last position of the substring to be extracted. The positions are counted
from 1. If g is a gi number which is identical to some gi number
in the input file, then the corresponding substring is shown in
fasta format. In case i is 1 and j is 0, the entire sequence is shown.
The header of the fasta formatted sequence in the output begins with

>g i j

followed by the original original fasta header.
Duplicated lines in the input file lead to only one sequence in the output.
The sequences are output according to the order in the original sequence.
The formatting of the output can be controlled by the options
-width, -o, -gzip, and -bzip2.]])
