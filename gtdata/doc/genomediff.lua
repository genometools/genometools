print ([[

The genomediff tool only accepts DNA input.

When used with sequence files or encseq, an enhanced suffix array will be
built in memory. The ESA will not be created completely, but construction will
use -memlimit as a threshold and build it partwise, calculating Kr for each
part.

File format for option -unitfile (lua):

  units = {
   genome1 = { "file1", "file2" },
   genome2 = { "file3", "file4" }
  }

only give basenames of the files!
Comment lines start with '--' and will be ignored.
See GTDIR/testdata/genomediff/unitfile1.lua for an example.

Options -pl -dc -memlimit are options to influence esa construction. 
]])
-- vim: tw=78
