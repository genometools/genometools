require 'set'
require 'digest/md5'

testdatafiles = Set.new

testdatafiles.add('at1MB')
testdatafiles.add('at100K1')
testdatafiles.add('Random-Small.fna')
testdatafiles.add('Random.fna')
testdatafiles.add('Atinsert.fna')
testdatafiles.add('TTT-small.fna')

def expandfilename(testdatafiles,filename)
  if testdatafiles.member?(filename)
    return "#{$testdata}#{filename}"
  else
    return "#{$gttestdata}DNA-mix/Grumbach.fna/#{filename}"
  end
end

def fromoptions2indexname(options)
  return options.gsub(/[ ]+/,"").gsub(/^-/,"")
end

def fileargs(testdatafiles,options)
  dbarg=['-db']
  doappend = false
  options.split(/ /).each do |opt|
    if opt == '-db'
      doappend = true
    elsif doappend
      dbarg.push(expandfilename(testdatafiles,opt))
    end
  end
  return dbarg.join(" ")
end

def grepprjfile(filename)
begin
  f = File.new(filename)
rescue => error
  STDERR.puts "#{$0}: cannot open \"#{filename}\": #{error}"
  exit 1
end
  return f.readlines.grep(/prefixlength|readmode|realspecialranges/)
end

outoptions="-tis -lcp -suf -bwt"

callargs =   ["#{outoptions} -db Random-Small.fna",
              "#{outoptions} -db Random.fna",
              "#{outoptions} -db Atinsert.fna Random.fna",
              "#{outoptions} -db TTT-small.fna",
              "#{outoptions} -lcp -db at100K1",
              "#{outoptions} -bwt -db at100K1",
              "#{outoptions} -bwt -lcp -db at100K1",
              "#{outoptions} -suf -db at100K1",
              "#{outoptions} -suf -lcp -db at100K1",
              "#{outoptions} -suf -bwt -db at100K1",
              "#{outoptions} -suf -bwt -lcp -db at100K1",
              "#{outoptions} -tis -db at100K1",
              "#{outoptions} -tis -lcp -db at100K1",
              "#{outoptions} -tis -bwt -db at100K1",
              "#{outoptions} -tis -bwt -lcp -db at100K1",
              "#{outoptions} -tis -suf -db at100K1",
              "#{outoptions} -tis -suf -lcp -db at100K1",
              "#{outoptions} -tis -suf -bwt -db at100K1",
              "#{outoptions} -tis -suf -bwt -lcp -db at100K1"]

if $gttestdata then
  callargs.push("#{outoptions} -db at100K1 at1MB Wildcards.fna chntxx.fna " +
               "hs5hcmvcg.fna humdystrop.fna humghcsa.fna " +
               "humhbb.fna humhdabcd.fna humhprtb.fna " +
               "mipacga.fna mpocpcg.fna mpomtcg.fna " +
               "vaccg.fna ychrIII.fna")
end

Name "gt checking project files"
Keywords "gt_checkprjfiles"
Test do
  callargs.each do |args|
    indexfileprefix=fromoptions2indexname(args)
    # long filenames can cause problems on encrypted file systems
    # this will produce a prefix of exactly 100 chars (resulting in 104 byte
    # filenames)
    if indexfileprefix.length >= 100
            indexfileprefix = indexfileprefix[0..67] + 
              Digest::MD5.hexdigest(indexfileprefix)[0..31]
    end
    indexname="#{indexfileprefix}"
    run_test "#{$bin}gt suffixerator -indexname #{indexname} " +
             "-pl -des #{fileargs(testdatafiles,args)}"
    thislist=grepprjfile("#{indexname}.prj")
    reflist=grepprjfile("#{$testdata}prj-files/#{indexname}.prj")
    if thislist != reflist
      STDERR.puts "files #{indexname}.prj and " +
                  "#{$testdata}prj-files/#{indexname}.prj differ"
      exit 1
    end
  end
end
