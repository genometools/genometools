>module Main where

>import List

>db = [s ++ " -db #{dbdna}" | s<-[" -dna"," -smap TransDNA"]]
>plen = [s ++ " -pl " | s<-[" -bck"," "]]

>calls::[String]
>calls = [
>         "def makesuffixeratorarglisttable(dbdna)",
>         "  arglisttable =",
>         "  ["] ++
>        [gencall (d ++ p ++ tis ++ suf ++ bwt) | 
>         d<-db,
>         p<-plen,
>         tis<-[" -tis"," "],
>         suf<-[" -suf -lcp"] ++ [" "],
>         bwt<-[" -bwt"," "]
>         ] ++
>         ["  ]",
>          "  return arglisttable",
>          "end\n"]

>gencall s = "   \"" ++ s ++ "\","

>main::IO ()
>main = writeFile "Suffixerator-table.rb" (unlines calls)
