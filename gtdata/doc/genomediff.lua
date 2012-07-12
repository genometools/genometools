print ([[

This program only works for DNA!

When used with sequence files or encseq, an enhanced suffix array will be 
build in memory. The ESA will not be created completely, but construction will
use -memlimit as a threshold and build it partwise, calculating kr for each
part.

File format for option -unitfile (lua):

units = {
 genome1 = { "file1", "file2" },
 genome2 = { "file3", "file4" }
}

only give basenames of the files!
Comment lines start with '--' and will be ignored.
See GTDIR/testdata/genomediff/unitfile1.lua for an example.

Options -pl -dc -spmopt -memlimit -dir are options to influence esa
construction. Option -dir is not usable in this context (Sequences should
contain forward and reverse complement, or -mirrored should be used).
]])
-- vim: tw=78
