-- This is the GenomeTools annotation viewer config file.
-- All options must be set inside the 'config' table.
-- All items in this table must be tables, called 'sections'.
-- The items in the 'colors' section must be given in
-- the following form:
--      <name> = {red=<val>,green=<val>,blue=<val>}
-- where <val> is a decimal value between 0 and 1.
--
-- All other sections must contain string-string mappings.
--      <name> = "<string>"

config = 
{
  colors = 
  {
    exon ={red=0.1,green=0.2,blue=0.3},
     cds ={red=0.4,green=0.5,blue=0.6}
  },
  collapse = 
  {
    exon = "mRNA"
  }
} 
