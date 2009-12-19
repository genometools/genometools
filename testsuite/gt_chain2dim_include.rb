def runchain2dimfailure(args,matchfile='matchfile')
  Name "gt chain2dim failure"
  Keywords "gt_chain2dim"
  Test do
    run_test "#{$bin}gt chain2dim -m #{matchfile} " + args,:retval => 1
  end
end

def runchain2dim(args,matchfile)
  Name "gt chain2dim #{args}"
  Keywords "gt_chain2dim"
  Test do
    run_test "#{$bin}gt chain2dim -m #{matchfile} " + args
    argstring = args.gsub(/[ ]/,"")
    run "cmp -s #{$last_stdout} #{$testdata}/chaindata/chain#{argstring}"
  end
end

def runchain2dimvschain2dim(args,matchfile)
  Name "gt chain2dim #{args}"
  Keywords "gt_chain2dim"
  Test do
    run "rf2
    run_test "#{$bin}gt chain2dim -m #{matchfile} " + args
    run "cp #{$last_stdout} gtchain.out"
    run "/Users/stefan/bin-ops/i686-apple-darwin/chain2dim.x " + args + 
        " #{matchfile}"
    run "cp #{$last_stdout} vschain.out"
    run "cmp -s gtchain.out vschain.out"
  end
end

def runchain2dimall(params,matchfile)
  params.each do |args|
    runchain2dim(args,matchfile)
  end
end

runchain2dimfailure("-maxgap 0")
runchain2dimfailure("-maxgap -1")
runchain2dimfailure("-wf 0.0")
runchain2dimfailure("-wf -1.0")
runchain2dimfailure("-wf 1.0")
runchain2dimfailure("-global xv")
runchain2dimfailure("-global ov h")
runchain2dimfailure("-local 2p h")
runchain2dimfailure("-local -global")
runchain2dimfailure("-global","#{$testdata}/ecolicmp-neg.of")
runchain2dimfailure("-local","#{$testdata}/ecolicmp-4.of")
runchain2dimfailure("-global","#{$testdata}/ecolicmp-seÂ.of")

# runchain2dimall("#{$testdata}/ecolicmp.of")

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

runchain2dimall(params,"#{$testdata}/ecolicmp250.of")

resultdir="#{$gttestdata}/repfind-result"
Dir.foreach("#{resultdir}") do |matchfile|
  if matchfile.match(/^\S+\.fna-\S+\.fna.result$/)
    params.each do |args|
      runchain2dimvschain2dim(args,"#{resultdir}/#{matchfile}")
    end
  end
end
