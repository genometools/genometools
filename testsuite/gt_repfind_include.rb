repfindtestfiles=["Duplicate.fna",
                  "Wildcards.fna",
                  "hs5hcmvcg.fna",
                  "humhbb.fna",
                  "mpomtcg.fna",
                  "at1MB",
                  "ychrIII.fna"]

def testdatadir(filename)
  if filename == 'Duplicate.fna' or filename == 'at1MB'
    return "#{$testdata}/"
  else
    return "#{$gttestdata}DNA-mix/Grumbach.fna/"
  end
end

def determineminlength(reffile)
  if reffile == 'at1MB'
    return 22
  elsif reffile == 'Duplicate.fna'
    return 8
  else
    return 14
  end
end

def checkrepfind(reffile,withextend = false)
  reffilepath=testdatadir(reffile) + reffile
  if reffile == 'Duplicate.fna'
    testdatadir = $testdata
  else
    testdatadir = $gttestdata
  end
  run_test("#{$bin}gt suffixerator -algbds 3 31 80 -db " +
           "#{reffilepath} -indexname sfxidx -dna -suf -tis -lcp -ssp -pl",
           :maxtime => 320)
  minlength = determineminlength(reffile)
  run_test("#{$bin}gt repfind -l #{minlength} -ii sfxidx", :maxtime => 600)
  resultfile="#{testdatadir}repfind-result/#{reffile}.result"
  run "cmp -s #{last_stdout} #{resultfile}"
  run_test("#{$bin}gt repfind -l #{minlength} -r -ii sfxidx", :maxtime => 600)
  resultfile="#{testdatadir}repfind-result/#{reffile}-r.result"
  run "cmp -s #{last_stdout} #{resultfile}"
  if withextend
    run_test("#{$bin}gt repfind -l #{minlength} -ii sfxidx -extendgreedy " +
             "-minidentity 90 -maxalilendiff 30 -percmathistory 55",
             :maxtime => 600)
    resultfile="#{testdatadir}repfind-result/#{reffile}-gr-ext.result"
    run "cmp -s #{last_stdout} #{resultfile}"
  end
end

def checkrepfindwithquery(reffile,queryfile)
  reffilepath=testdatadir(reffile) + reffile
  queryfilepath=testdatadir(queryfile) + queryfile
  idxname=reffile + "-idx"
  run_test "#{$bin}gt suffixerator -algbds 3 31 80 -db " +
           "#{reffilepath} -indexname #{idxname} -dna -suf -tis -lcp -ssp -pl"
  run_test("#{$bin}gt repfind -l 15 -ii #{idxname} -q #{queryfilepath}",
           :maxtime => 600)
  # run "sort #{last_stdout}"
  #run "/Users/kurtz/bin-ops/i686-apple-darwin/mkvtree.x -indexname mkv-idx " +
  #    "-allout -v -pl -dna -db #{reffilepath}"
  #run "/Users/kurtz/bin-ops/i686-apple-darwin/vmatch-mini.x 15 mkv-idx " +
  #    "#{queryfilepath}"
  #run "sed -e '/^#/d' #{last_stdout} | sort"
  # run "#{$scriptsdir}repfind-cmp.rb #{last_stdout} #{$gttestdata}repfind-result/#{reffile}-#{queryfile}.result"
  run "cmp -s #{last_stdout} #{$gttestdata}repfind-result/#{reffile}-#{queryfile}.result"
end

Name "gt repfind mirror symmetric"
Keywords "gt_repfind"
Test do
  10.times.each do
    run "#{$scriptsdir}gen-randseq.rb --seedlength 200 --length 2200 --mirrored"
    run_test "#{$bin}gt suffixerator -suftabuint -db #{last_stdout} " +
             "-dna -suf -tis -lcp -md5 no -des no -sds no -indexname sfx"
    run_test "#{$bin}gt repfind -minidentity 90 -percmathistory 55 " +
             "-scan -check_extend_symmetry -seedlength 200 " +
             "-extendgreedy -ii sfx -maxalilendiff 30"
  end
end

Name "gt repfind small"
Keywords "gt_repfind"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname sfx -dna -tis -suf -lcp -ssp -pl"
  run_test "#{$bin}gt repfind -l 8 -ii sfx"
  run "grep -v '^#' #{last_stdout}"
  run "diff -w #{last_stdout} #{$testdata}repfind-8-Atinsert.txt"
  run_test "#{$bin}gt repfind -scan -l 8 -ii sfx"
  run "grep -v '^#' #{last_stdout}"
  run "diff -w #{last_stdout} #{$testdata}repfind-8-Atinsert.txt"
  run_test "#{$bin}gt repfind -samples 10 -l 6 -ii sfx",:maxtime => 600
  run "#{$bin}gt repfind -samples 1000 -l 6 -ii sfx",:maxtime => 600
end

Name "gt repfind extend at1MB"
Keywords "gt_repfind extend"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}at1MB " +
           "-indexname sfx -dna -tis -suf -lcp"
  run_test "#{$bin}gt repfind -minidentity 90 -l 20 -xdropbelow 5 -extendxdrop -ii sfx"
  run "diff #{last_stdout} #{$testdata}repfind-20-extend.txt"
  run_test "#{$bin}gt repfind -minidentity 85 -l 20 -extendxdrop -ii sfx -q " +
           "#{$testdata}/U89959_genomic.fas"
  run "diff #{last_stdout} #{$testdata}repfind-20-query-extend.txt"
end

if $gttestdata then
  extendexception = ["hs5hcmvcg.fna","Wildcards.fna","at1MB"]
  repfindtestfiles.each do |reffile|
    Name "gt repfind #{reffile}"
    Keywords "gt_repfind"
    Test do
      withextend = if extendexception.member?(reffile) then false else true end
      checkrepfind(reffile,withextend)
    end
    repfindtestfiles.each do |queryfile|
      if reffile != queryfile
        Name "gt repfind #{reffile} versus #{queryfile}"
        Keywords "gt_repfind"
        Test do
          checkrepfindwithquery(reffile,queryfile)
        end
      end
    end
  end
end
