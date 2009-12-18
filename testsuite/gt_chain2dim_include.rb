def runchain2dimfailure(args)
  Name "gt chain2dim failure"
  Keywords "gt_chain2dim"
  Test do
    run_test "#{$bin}gt chain2dim -m matchfile " + args,:retval => 1
  end
end

def runchain2dim(args)
  Name "gt chain2dim failure"
  Keywords "gt_chain2dim"
  Test do
    run_test "#{$bin}gt chain2dim -m #{$testdata}/ecolicmp.of " + args
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

runchain2dim("-global gc")
runchain2dim("-global ov")
runchain2dim("-global")
runchain2dim("-global gc -wf 1.5")
runchain2dim("-global ov -wf 1.8")
runchain2dim("-global -maxgap 10")
runchain2dim("-global gc -wf 1.5 -maxgap 10")
runchain2dim("-global ov -wf 1.8 -maxgap 10")
runchain2dim("-local")
runchain2dim("-local 2p")
runchain2dim("-local 2b")
runchain2dim("-local 20")
runchain2dim("-local -wf 1.6")
runchain2dim("-local 2p -wf 1.6")
runchain2dim("-local 2b -wf 1.6")
runchain2dim("-local 20 -wf 1.6")
runchain2dim("-local -wf 1.6 -maxgap 10")
runchain2dim("-local 2p -wf 1.6 -maxgap 10")
runchain2dim("-local 2b -wf 1.6 -maxgap 10")
runchain2dim("-local 20 -wf 1.6 -maxgap 10")
