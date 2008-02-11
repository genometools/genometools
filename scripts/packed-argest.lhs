>module Main where

>import List

>db = [s ++ " -db #{dbdna}" | s<-[" -dna"," -smap TransDNA"]]
>plen = [s ++ " -pl" | s<-[" -bck"," "]]
>dirs = [" -dir " ++ s | s<-["fwd","cpl","rev","rcl"]]
>parts = [" -parts " ++ show i | i<-[1,2,3]]
>cartesian x y = [s++t | s<-x, t<-y]

>sufbwt::[String]
>sufbwt = cartesian sufargs bwtargs
>         where sufargs = [" -suf -lcp " ++ dir | dir<-dirs] ++ [" "]
>               bwtargs = [" -bwt"," "]

>packedindex::[String]
>packedindex = cartesian (cartesian bsize blbuck) locfreq
>              where bsize = [" -bsize " ++ show i | i<-[8,9,10]]
>                    blbuck = [" -blbuck " ++ show i | i<-[8,9,10]]
>                    locfreq = [" -locfreq " ++ show i | i<-[0,8,16,32]]

>calls::[String]->[String]
>calls xx = [
>            "def makesuffixeratorarglisttable(dbdna)",
>            "  arglisttable =",
>            "  ["] ++
>           [gencall (d ++ p ++ tis ++ x) | 
>            d<-db,
>            p<-plen,
>            part<-parts,
>            tis<-[" -tis"," "],
>            x<-xx
>            ] ++
>            ["  ]",
>             "  return arglisttable",
>             "end\n"]
>            where sufargs = [" -suf -lcp " ++ dir | dir<-dirs] ++ [" "]
>                  bwtargs = [" -bwt"," "]

>gencall::String->String
>gencall s = "   \"" ++ s ++ "\","

>main::IO ()
>main = writeFile "Suffixerator-table.rb" (unlines (calls packedindex))

main = writeFile "Suffixerator-table.rb" (unlines.calls sufbwt)

