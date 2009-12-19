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
    run_test "#{$bin}gt chain2dim -m #{matchfile} " + args
    run "cp #{$last_stdout} gtchain.out"
    run "/Users/stefan/bin-ops/i686-apple-darwin/chain2dim.x " + args + 
        " #{matchfile}"
    run "cp #{$last_stdout} vschain.out"
    run "cmp -s gtchain.out vschain.out"
    argstring = args.gsub(/[ ]/,"")
    run "cp gtchain.out chain#{argstring}"
  end
end

def runchain2dimall(matchfile)
  runchain2dim("-global",matchfile)
  runchain2dim("-silent -global",matchfile)
  runchain2dim("-local -wf 1.8",matchfile)
  runchain2dim("-local -wf 0.5",matchfile)
  runchain2dim("-local -maxgap 20",matchfile)
  runchain2dim("-local 2b ",matchfile)
  runchain2dim("-local 55p -silent",matchfile)
  runchain2dim("-global gc",matchfile)
  runchain2dim("-global ov",matchfile)
  runchain2dim("-global gc -wf 1.5",matchfile)
  runchain2dim("-global ov -wf 1.8",matchfile)
  runchain2dim("-global -maxgap 10",matchfile)
  runchain2dim("-global gc -wf 1.5 -maxgap 10",matchfile)
  runchain2dim("-global ov -wf 1.8 -maxgap 10",matchfile)
  runchain2dim("-local",matchfile)
  runchain2dim("-local 2p",matchfile)
  runchain2dim("-local 20",matchfile)
  runchain2dim("-local 2p -wf 1.8",matchfile)
  runchain2dim("-local 2b -wf 1.8",matchfile)
  runchain2dim("-local 20 -wf 1.8",matchfile)
  runchain2dim("-local -wf 1.8 -maxgap 20",matchfile)
  runchain2dim("-local 2p -wf 1.8 -maxgap 10",matchfile)
  runchain2dim("-local 2b -wf 1.8 -maxgap 10",matchfile)
  runchain2dim("-local 20 -wf 1.8 -maxgap 10",matchfile)
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
runchain2dimall("#{$testdata}/ecolicmp250.of")
