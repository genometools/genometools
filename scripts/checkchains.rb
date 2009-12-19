#!/usr/bin/env ruby

def runchain2dimvschain2dim(args,matchfile)
  system "cat #{matchfile} | ../scripts/rf2of.rb > tmp.of"
  system "../bin/gt chain2dim " + args + " -m tmp.of > gtchain.out"
  system "chain2dim.x " + args + " tmp.of > vschain.out"
  system "cmp -s gtchain.out vschain.out"
end

params = ["-global",
	  "-silent -global",
	  "-local -wf 1.8",
	  "-local -wf 0.5",
	  "-local -maxgap 20",
	  "-local 2b",
	  "-local 55p -silent",
	  "-global gc",
	  "-global ov",
	  "-global gc -wf 1.5",
	  "-global ov -wf 1.8",
	  "-global -maxgap 10",
	  "-global gc -wf 1.5 -maxgap 10",
	  "-global ov -wf 1.8 -maxgap 10",
	  "-local",
	  "-local 2p",
	  "-local 20",
	  "-local 2p -wf 1.8",
	  "-local 2b -wf 1.8",
	  "-local 20 -wf 1.8",
	  "-local -wf 1.8 -maxgap 20",
	  "-local 2p -wf 1.8 -maxgap 10",
	  "-local 2b -wf 1.8 -maxgap 10",
	  "-local 20 -wf 1.8 -maxgap 10"]


resultdir="#{ENV['GTTESTDATA']}/repfind-result"

Dir.foreach("#{resultdir}") do |matchfile|
  if matchfile.match(/^\S+\.fna-\S+\.fna.result$/)
    params.each do |args|
      runchain2dimvschain2dim(args,"#{resultdir}/#{matchfile}")
    end
  end
end

system "rm -f tmp.of gtchain.out vschain.out"
