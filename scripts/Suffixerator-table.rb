def makesuffixeratorarglisttable(dbdna)
  arglisttable =
  [
   " -dna -db #{dbdna} -bck -pl  -tis -suf -lcp -bwt",
   " -dna -db #{dbdna} -bck -pl  -tis -suf -lcp ",
   " -dna -db #{dbdna} -bck -pl  -tis  -bwt",
   " -dna -db #{dbdna} -bck -pl  -tis  ",
   " -dna -db #{dbdna} -bck -pl   -suf -lcp -bwt",
   " -dna -db #{dbdna} -bck -pl   -suf -lcp ",
   " -dna -db #{dbdna} -bck -pl    -bwt",
   " -dna -db #{dbdna} -bck -pl    ",
   " -dna -db #{dbdna}  -pl  -tis -suf -lcp -bwt",
   " -dna -db #{dbdna}  -pl  -tis -suf -lcp ",
   " -dna -db #{dbdna}  -pl  -tis  -bwt",
   " -dna -db #{dbdna}  -pl  -tis  ",
   " -dna -db #{dbdna}  -pl   -suf -lcp -bwt",
   " -dna -db #{dbdna}  -pl   -suf -lcp ",
   " -dna -db #{dbdna}  -pl    -bwt",
   " -dna -db #{dbdna}  -pl    ",
   " -smap TransDNA -db #{dbdna} -bck -pl  -tis -suf -lcp -bwt",
   " -smap TransDNA -db #{dbdna} -bck -pl  -tis -suf -lcp ",
   " -smap TransDNA -db #{dbdna} -bck -pl  -tis  -bwt",
   " -smap TransDNA -db #{dbdna} -bck -pl  -tis  ",
   " -smap TransDNA -db #{dbdna} -bck -pl   -suf -lcp -bwt",
   " -smap TransDNA -db #{dbdna} -bck -pl   -suf -lcp ",
   " -smap TransDNA -db #{dbdna} -bck -pl    -bwt",
   " -smap TransDNA -db #{dbdna} -bck -pl    ",
   " -smap TransDNA -db #{dbdna}  -pl  -tis -suf -lcp -bwt",
   " -smap TransDNA -db #{dbdna}  -pl  -tis -suf -lcp ",
   " -smap TransDNA -db #{dbdna}  -pl  -tis  -bwt",
   " -smap TransDNA -db #{dbdna}  -pl  -tis  ",
   " -smap TransDNA -db #{dbdna}  -pl   -suf -lcp -bwt",
   " -smap TransDNA -db #{dbdna}  -pl   -suf -lcp ",
   " -smap TransDNA -db #{dbdna}  -pl    -bwt",
   " -smap TransDNA -db #{dbdna}  -pl    ",
  ]
  return arglisttable
end

