def outoptionsnobck
  return "-tis -suf -des -sds -ssp -lcp -bwt"
end

def outoptions
  return outoptionsnobck + " -bck"
end

def checksfx(parts,withsmap,cmp,filelist,alldirs=true)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  if alldirs
    dirlist=["fwd","rev","cpl ","rcl "]
  else
    dirlist=["fwd ","rev"]
  end
  dirlist.each do |dirarg|
    extra=""
    if cmp
      extra=extra + " -cmpcharbychar"
    end
    run_test "#{$bin}gt suffixerator -v -parts #{parts} -pl " +
             "-algbds 10 31 80 #{extra} #{outoptions} " +
             "-indexname esa -dir " + dirarg + " -db " + filearg
    if dirarg == "cpl" or dirarg = "rcl"
      run_test "#{$bin}gt dev sfxmap #{outoptions} -v " +
               "-esa esa",
               :maxtime => 600
    else
      if dirarg == "fwd"
        dirarg_rev = "rev"
      else
        dirarg_rev = "fwd"
      end
      run_test "#{$bin}gt suffixerator -v -parts #{parts} -pl " +
               "-algbds 10 31 80 #{extra} #{outoptions} " +
               "-indexname esa-#{dirarg_rev} -dir " + dirarg_rev +
               " -db " + filearg
      run_test "#{$bin}gt packedindex mkindex -indexname pck -dir " + dirarg +
               " -db " + filearg
      run_test "#{$bin}gt dev sfxmap #{outoptions} -v " +
               "-esa esa -pck pck -cmpsuf",
               :maxtime => 600
      run_test "#{$bin}gt dev sfxmap #{outoptions} -v " +
               "-esa esa-#{dirarg_rev} -pck pck -cmplcp",
               :maxtime => 600
    end
  end
end

def checkdc(filelist)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run_test "#{$bin}gt suffixerator -v -pl -dc 64 -dccheck " +
           "-lcp -suf -ssp -tis -indexname sfx -db " + filearg
  run_test "#{$bin}gt dev sfxmap -suf -tis -ssp -v -esa sfx",
           :maxtime => 600
  run_test "#{$bin}gt suffixerator -v -pl -parts 3 -dc 64 -dccheck " +
           "-lcp -suf -tis -indexname sfx3 -db " + filearg
  run "diff sfx3.suf sfx.suf"
end

def flattenfilelist(filelist)
  s=""
  filelist.each do |f|
    s += "#{$testdata}#{f} "
  end
  return s
end

def checkbwt(filelist)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run_test "#{$bin}gt suffixerator -pl #{outoptions} -indexname sfx -db " +
           flattenfilelist(filelist)
end

allfiles = []
all_fastafiles = ["Arabidopsis-C99826.fna",
                  "Atinsert.fna",
                  "Atinsert_seqrange_13-17_rev.fna",
                  "Atinsert_seqrange_3-7.fna",
                  "Atinsert_single_3.fna",
                  "Atinsert_single_3_rev.fna",
                  "Copysorttest.fna",
                  "Duplicate.fna",
                  "Ecoli-section1.fna",
                  "Ecoli-section2.fna",
                  "Random-Small.fna",
                  "Random.fna",
                  "Random159.fna",
                  "Random160.fna",
                  "RandomN.fna",
                  "Reads1.fna",
                  "Reads2.fna",
                  "Reads3.fna",
                  "Repfind-example.fna",
                  "TTTN.fna",
                  "Small.fna",
                  "Smalldup.fna",
                  "TTT-small.fna",
                  "trna_glutamine.fna",
                  "Verysmall.fna"]

allfiles += all_fastafiles

allmultifiles = []
all_multifastafiles = ["Atinsert.fna",
                       "Duplicate.fna",
                       "Small.fna",
                       "Verysmall.fna",
                       "Random159.fna",
                       "Random160.fna"]

allmultifiles += all_multifastafiles

all_fastqfiles = ["fastq_long.fastq",
                  "test10_multiline.fastq",
                  "test1.fastq",
                  "test5_tricky.fastq"]

allmultifiles += all_fastqfiles
allfiles += all_fastqfiles

alldir = ["fwd","cpl","rev","rcl"]

all_fastafiles.each do |filelist|
  Name "gt suffixerator + sfxmap"
  Keywords "gt_suffixerator tis"
  Test do
    run_test "#{$bin}gt suffixerator -tis -ssp -indexname sfx -db " +
             flattenfilelist(filelist)
    run_test "#{$bin}gt dev sfxmap -ssp -tis -esa sfx"
    filelist.each do |filename|
      run_test "#{$bin}gt suffixerator -tis -ssp -indexname sfx -db " +
               "#{$testdata}#{filename}"
      run_test "#{$bin}gt dev sfxmap -ssp -tis -esa sfx"
    end
  end
end

all_fastafiles.each do |filename|
  Name "gt suffixerator #{filename} all accesstypes"
  Keywords "gt_suffixerator tis"
  Test do
    ["direct", "bit", "uchar", "ushort", "uint32"].each do |sat|
      run_test "#{$bin}gt suffixerator -tis -indexname sfx -sat #{sat} " +
               "-db #{$testdata}#{filename}"
    end
  end
end

Name "gt suffixerator file of reads of equal length"
Keywords "gt_suffixerator reads"
Test do
  run_test "#{$bin}/gt suffixerator -des -tis -ssp -dna " +
           "-db #{$testdata}U89959_genomic.fas -indexname u8idx"
  run_test "#{$bin}/gt simreads -coverage 4 -len 100 -force -o u8.reads u8idx"
  run_test "#{$bin}/gt suffixerator -v #{outoptionsnobck} -dna -db u8.reads"
  run "grep -q '# init character encoding (eqlen ' #{$last_stdout}"
  run_test "#{$bin}/gt dev sfxmap -suf -lcp -des -sds -ssp -esa u8.reads", \
           :maxtime => 200
  run_test "#{$bin}/gt suffixerator -v #{outoptionsnobck} -dir rev -dna -db u8.reads"
  run_test "#{$bin}/gt dev sfxmap -suf -lcp -des -sds -ssp -esa u8.reads", \
           :maxtime => 200
end

all_fastafiles.each do |filename|
  Name "gt suffixerator -dc 64 -dccheck -lcp -parts 1+3 #{filename}"
  Keywords "gt_suffixerator dc"
  Test do
    checkdc([filename])
  end
end

Name "gt suffixerator -dc 64 -dccheck -lcp -parts 1+3 all-fastafiles"
Keywords "gt_suffixerator dc"
Test do
  checkdc(all_fastafiles)
end

alldir.each do |dir|
  all_fastafiles.each do |filelist|
    Name "gt suffixerator #{dir}"
    Keywords "gt_suffixerator"
    Test do
      run_test "#{$bin}gt suffixerator -dir #{dir} -tis -suf -bwt -lcp " +
               "-indexname sfx -pl -db " +
               flattenfilelist(filelist)
      run_test "#{$bin}gt suffixerator -storespecialcodes -dir #{dir} -tis " +
               "-suf -lcp -indexname sfx -pl -db " +
               flattenfilelist(filelist)
      run_test "#{$bin}gt suffixerator -tis -bwt -lcp -pl -ii sfx"
    end
  end
end

faillist = ["-indexname sfx -db /nothing",
            "-indexname /nothing/sfx -db #{$testdata}TTT-small.fna",
            "-smap /nothing -db #{$testdata}TTT-small.fna",
            "-dna -db #{$testdata}sw100K1.fsa",
            "-protein -dir cpl -db #{$testdata}sw100K1.fsa",
            "-dna -db #{$testdata}Random.fna RandomN.fna",
            "-dna -suf -pl 10 -db #{$testdata}Random.fna",
            "-dna -tis -sat plain -db #{$testdata}TTT-small.fna"]

faillist.each do |failcommand|
  Name "gt suffixerator failure"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis " + failcommand,:retval => 1
  end
end

allmultifiles.each do |filename|
  Name "gt suffixerator sfxmap-failure #{filename}"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis -dna -des no -sds no -ssp no " +
             "-indexname localidx  -sat direct " +
             "-db #{$testdata}#{filename}"
    run_test "#{$bin}gt suffixerator -suf -lcp -pl -dir rev -ii localidx"
    run_test "#{$bin}gt dev sfxmap -tis -des -sds -esa localidx",
             :retval => 1
    # In short read files with equal read lengths, ssptabs need not
    # be built explicitly.
    # Thus these tests only fail for non-equal length files.
    if !all_fastqfiles.include?(filename) then
      run_test "#{$bin}gt dev sfxmap -tis -ssp -esa localidx",
               :retval => 1
      run_test "#{$bin}gt dev sfxmap -ssp -esa localidx",
               :retval => 1
    end
    run_test "#{$bin}gt dev sfxmap -des -esa localidx",
             :retval => 1
    run_test "#{$bin}gt dev sfxmap -tis -bck -esa localidx",
             :retval => 1
  end
end

Name "gt suffixerator bwt"
Keywords "gt_suffixerator"
Test do
  checkbwt(all_fastafiles)
end

1.upto(3) do |parts|
  [0,2].each do |withsmap|
    extra=""
    if withsmap == 1
      extra="-protein"
      extraname="protein"
    elsif withsmap == 2
      extra="-smap TransProt11"
      extraname="TransProt11"
    end
    Name "gt suffixerator+map protein filelist #{extraname} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checksfx(parts,extra,true,["sw100K1.fsa","sw100K2.fsa"],false)
    end
  end
end

0.upto(2) do |cmpval|
  1.upto(2) do |parts|
    [0,2].each do |withsmap|
      extra=""
      if withsmap == 1
	extra="-dna"
	extraname=" dna"
      elsif withsmap == 2
	extra="-smap TransDNA"
	extraname=" trans"
      end
      if cmpval == 0
	cmp=false
      elsif cmpval == 1
	cmp=true
      else
	cmp=true
      end
      all_fastafiles.each do |filename|
	Name "gt suffixerator+map #{filename}#{extraname} #{parts} parts"
	Keywords "gt_suffixerator"
	Test do
	  checksfx(parts,extra,cmp,[filename])
	end
      end
      filelist=["RandomN.fna","Random.fna","Atinsert.fna"]
      Name "gt suffixerator+map dna filelist#{extraname} #{parts} parts"
      Keywords "gt_suffixerator"
      Test do
	checksfx(parts,extra,cmp,filelist)
      end
    end
  end
end

def checkmapped(keyword,args)
  Name "gt suffixerator checkmapped #{keyword}"
  Keywords "gt_suffixerator gttestdata"
  Test do
    run_test "#{$bin}gt suffixerator #{outoptions} -algbds 3 34 90 " +
             "-indexname sfxidx #{args}",
             :maxtime => 1200
    run_test "#{$bin}gt dev sfxmap #{outoptions} -v -esa sfxidx",
             :maxtime => 2400
    run_test "#{$bin}gt dev sfxmap #{outoptionsnobck} -stream -v -esa sfxidx",
             :maxtime => 2400
  end
end

def grumbach()
  return "#{$gttestdata}DNA-mix/Grumbach.fna/"
end

if $gttestdata then
  checkmapped("many .fna files","-db " +
              "#{$gttestdata}Iowa/at100K1 " +
              "#{grumbach()}Wildcards.fna " +
              "#{grumbach()}chntxx.fna " +
              "#{grumbach()}hs5hcmvcg.fna " +
              "#{grumbach()}humdystrop.fna " +
              "#{grumbach()}humghcsa.fna " +
              "#{grumbach()}humhdabcd.fna " +
              "#{grumbach()}humhprtb.fna " +
              "#{grumbach()}mipacga.fna " +
              "#{grumbach()}mpocpcg.fna " +
              "#{grumbach()}ychrIII.fna " +
              "-parts 3 -pl")

  checkmapped("swiss with parts=1",
              "-parts 1 -pl -db #{$gttestdata}swissprot/swiss10K " +
              "#{$gttestdata}swissprot/swiss1MB")

  checkmapped("swiss with parts=3",
              "-db #{$gttestdata}swissprot/swiss10K " +
              "#{$gttestdata}swissprot/swiss1MB -parts 3 -pl")

  checkmapped("at100K1",
              "-parts 2 -pl -smap TransDNA -db  #{$gttestdata}Iowa/at100K1")

  checkmapped("swiss with TransProt11",
              "-db #{$gttestdata}swissprot/swiss10K -parts 1 -pl -smap " +
              "TransProt11")
end

SATS = ["direct", "bytecompress", "eqlen", "bit", "uchar", "ushort", "uint32"]

EQLENDNAFILE = {:filename => "#{$testdata}test1.fasta",
                :desc => "equal length DNA",
                :msgs => {
                  "bytecompress" => "cannot use bytecompress on DNA sequences"}}
DNAFILE   = {:filename => "#{$testdata}Atinsert.fna",
             :desc => "non-equal length DNA",
             :msgs => {
                "bytecompress" => "cannot use bytecompress on DNA sequences",
                "eqlen" => "all sequences are of equal length and no " + \
                "sequence contains"}}
EQLENAAFILE = {:filename => "#{$testdata}trembl-eqlen.faa",
                :desc => "equal length AA",
                :msgs => {
                  "eqlen" => "as the sequence is not DNA",
                  "bit" => "as the sequence is not DNA",
                  "uchar" => "as the sequence is not DNA",
                  "ushort" => "as the sequence is not DNA",
                  "uint32" => "as the sequence is not DNA"}}
AAFILE    = {:filename => "#{$testdata}trembl.faa",
                :desc => "non-equal length AA",
                :msgs => {
                  "eqlen" => "as the sequence is not DNA",
                  "bit" => "as the sequence is not DNA",
                  "uchar" => "as the sequence is not DNA",
                  "ushort" => "as the sequence is not DNA",
                  "uint32" => "as the sequence is not DNA"}}

SATTESTFILES = [EQLENDNAFILE, DNAFILE, EQLENAAFILE, AAFILE]

SATTESTFILES.each do |file|
  SATS.each do |sat|
    Name "gt suffixerator sat #{sat} -> #{file[:desc]}"
    Keywords "gt_suffixerator sats"
    Test do
      if !file[:msgs][sat].nil? then
        retval = 1
      else
        retval = 0
      end
      run_test "#{$bin}/gt suffixerator -sat #{sat} -v -suf -lcp -des -sds " + \
               "-ssp -tis -db #{file[:filename]} -indexname myidx", \
               :retval => retval
      if !file[:msgs][sat].nil? then
        grep($last_stderr, /#{file[:msgs][sat]}/)
      end
      run_test "#{$bin}/gt dev sfxmap -suf -lcp -des -sds -ssp -esa myidx", \
               :retval => retval
    end
  end
end

[EQLENDNAFILE, DNAFILE].each do |file|
  (SATS-(file[:msgs].keys)).each do |sat|
    Name "gt suffixerator mirror #{sat} #{file[:desc]}"
    Keywords "gt_suffixerator mirror #{sat}"
    Test do
      alldir.each do |dir|
        run_test "#{$bin}/gt suffixerator -mirrored -sat #{sat} -v -suf " + \
                 "-lcp -bwt -dir #{dir} -db #{file[:filename]} -indexname myidx"
        run      "#{$scriptsdir}/add_revcmp.rb #{file[:filename]} > reverse.fna"
        run_test "#{$bin}/gt suffixerator -sat #{sat} -v -suf -lcp " + \
                 "-bwt -dir #{dir} -db reverse.fna -indexname myidx_r"
        run      "diff myidx.suf myidx_r.suf"
        run      "diff myidx.lcp myidx_r.lcp"
        run      "diff myidx.bwt myidx_r.bwt"
        run_test "#{$bin}/gt dev sfxmap -suf -lcp -des -sds -ssp -esa myidx", \
                 :maxtime => 600
        run_test "#{$bin}/gt dev sfxmap -suf -lcp -des -sds -ssp -esa " + \
                 "myidx_r", :maxtime => 600

      end
    end
  end
end

small_fastafiles = ["Duplicate.fna",
                    "Random-Small.fna",
                    "Copysorttest.fna",
                    "Random159.fna",
                    "Random160.fna",
                    "TTT-small.fna",
                    "trna_glutamine.fna",
                    "Small.fna",
                    "Verysmall.fna",
                    "Arabidopsis-C99826.fna"]
small_fastafiles.each do |file|
  Name "differencecover sortmaxdepth #{file}"
  Keywords "differencecover sortmaxdepth"
  Test do
    file = "#{$testdata}/#{file}"
    run "#{$bin}/gt suffixerator -tis -db #{file} -indexname esa"
    [2, 10].each do |maxinsertionsort|
      [10, 30, 100].each do |maxbltriesort|
        [2, 7, 13].each do |maxdepth|
          run_test "#{$bin}/gt dev sfxmap -sortmaxdepth #{maxdepth} -esa " + \
                   "esa -v -algbds #{maxinsertionsort} #{maxbltriesort} 1000"
        end
        run_test "#{$bin}/gt suffixerator -db #{file} -indexname esa2 " + \
                 "-suf -v -algbds #{maxinsertionsort} #{maxbltriesort} 1000"
        run_test "#{$bin}/gt dev sfxmap -suf -esa esa2"
      end
    end
  end
end

# check whether -memlimit is available as an option to enable
# conditional tests
if !`#{$bin}/gt suffixerator -help`.match(/memlimit/).nil? then
  if $gttestdata then
    ["2L", "2R", "3L", "3R"].each do |chr|
      Name "gt suffixerator -memlimit D. mel #{chr}"
      Keywords "gt_suffixerator memlimit"
      Test do
        size = nil
        run "#{$bin}/gt encseq encode -indexname dmel " + \
            "#{$gttestdata}ltrharvest/d_mel/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz", \
            :maxtime => 600
        out = `env GT_MEM_BOOKKEEPING=on GT_ENV_OPTIONS=-spacepeak #{$bin}/gt suffixerator -v -suf -ii dmel`
        if m = out.match(/space peak in megabytes: ([.0-9]+) /) then
          size = m[1].to_f
        else
            raise "could not get normal space peak"
        end
        [0.25, 0.5, 0.75].each do |frac|
          limit = size * frac
          run_test "env GT_MEM_BOOKKEEPING=on GT_ENV_OPTIONS=-spacepeak " + \
                   "#{$bin}/gt suffixerator " + \
                   "-suf -v " + \
                   "-memlimit #{limit.floor}MB " + \
                   "-ii dmel", :maxtime => 600
          lastout = nil
          File.open($last_stdout) do |f|
            out = f.read
            if m = out.match(/space peak in megabytes: ([.0-9]+) /) then
              mysize = m[1].to_f
            else
              raise "could not get actual space peak!"
            end
            if mysize > limit and mysize - limit > (mysize/10) then
              raise "required size (#{mysize}) was higher than limit (#{limit})"
            end
          end
        end
      end
    end
  end

  Name "gt suffixerator -memlimit excludes -parts"
  Keywords "gt_suffixerator memlimit"
  Test do
    run "#{$bin}/gt suffixerator -db #{$testdata}/at1MB -indexname foo " + \
        "-parts 3 -memlimit 2MB", :retval => 1
    grep($last_stderr, /exclude each other/)
  end

  Name "gt suffixerator -memlimit invalid size"
  Keywords "gt_suffixerator memlimit"
  Test do
    run "#{$bin}/gt suffixerator -db #{$testdata}/at1MB -indexname foo " + \
        "-memlimit 2TB", :retval => 1
    grep($last_stderr, /one of the keywords MB and GB/)
  end

  Name "gt suffixerator -memlimit invalid size"
  Keywords "gt_suffixerator memlimit"
  Test do
    run "#{$bin}/gt suffixerator -db #{$testdata}/at1MB -indexname foo " + \
        "-memlimit 2.2GB", :retval => 1
    grep($last_stderr, /one of the keywords MB and GB/)
  end
end
