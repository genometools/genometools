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
    run "cmp -s #{last_stdout} #{$testdata}chaindata/chain#{argstring}"
  end
end

def runchain2dimall(params,matchfile)
  params.each do |args|
    runchain2dim(args,matchfile)
  end
end

Name "gt chain2dim small all"
Keywords "gt_chain2dim all"
Test do
  run_test "#{$bin}gt chain2dim -global all -m #{$testdata}chaindata/matches-nd.txt"
  run "cmp -s #{last_stdout} #{$testdata}chaindata/matches-nd.chains"
end

Name "gt chain2dim ecoli all"
Keywords "gt_chain2dim all"
Test do
  [283,298,304,753].each do |lines|
    run "head -n #{lines} #{$testdata}ecolicmp250.of | sed -e 's/[0-9]*$/1/'"
    run_test "#{$bin}gt chain2dim -global all -m #{last_stdout}"
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
runchain2dimfailure("-global","#{$testdata}ecolicmp-neg.of")
runchain2dimfailure("-global","#{$testdata}ecolicmp-se.of")

# runchain2dimall("#{$testdata}ecolicmp.of")

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

runchain2dimall(params,"#{$testdata}ecolicmp250.of")
