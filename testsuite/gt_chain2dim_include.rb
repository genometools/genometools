def runchain2dimfailure(args,matchfile='matchfile')
  Name "gt chain2dim failure"
  Keywords "gt_chain2dim"
  Test do
    run_test "#{$bin}gt chain2dim -m #{matchfile} " + args,:retval => 1
  end
end

def runchain2dim(args,matchfile="#{$testdata}/ecolicmp.of")
  Name "gt chain2dim failure"
  Keywords "gt_chain2dim"
  Test do
    run_test "#{$bin}gt chain2dim -m #{matchfile} " + args
    run "cp #{$last_stdout} gtchain.out"
    run "/Users/kurtz/bin-ops/i686-apple-darwin/chain2dim.x " + args + 
        " #{matchfile}"
    run "cp #{$last_stdout} vschain.out"
    run "cmp -s gtchain.out vschain.out"
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
runchain2dimfailure("-local","#{$testdata}/ecolicmp-4Â.of")
runchain2dimfailure("-global","#{$testdata}/ecolicmp-seÂ.of")

runchain2dim("-global")
runchain2dim("-silent -global")
runchain2dim("-local -wf 1.8")
runchain2dim("-local -wf 0.5")
runchain2dim("-local -maxgap 20")
runchain2dim("-local 2b ")
runchain2dim("-local 55p -silent")

runchain2dim("-global gc")
runchain2dim("-global ov")
runchain2dim("-global gc -wf 1.5")
runchain2dim("-global ov -wf 1.8")
runchain2dim("-global -maxgap 10")
runchain2dim("-global gc -wf 1.5 -maxgap 10")
runchain2dim("-global ov -wf 1.8 -maxgap 10")
runchain2dim("-local")
runchain2dim("-local 2p")
runchain2dim("-local 2b")
runchain2dim("-local 20")
runchain2dim("-local -wf 1.8")
runchain2dim("-local 2p -wf 1.8")
runchain2dim("-local 2b -wf 1.8")
runchain2dim("-local 20 -wf 1.8")
runchain2dim("-local -wf 1.8 -maxgap 20")
runchain2dim("-local 2p -wf 1.8 -maxgap 10")
runchain2dim("-local 2b -wf 1.8 -maxgap 10")
runchain2dim("-local 20 -wf 1.8 -maxgap 10")
