if $gttestdata then

Name "gt mgth testdata mixed options 1"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Pyrococcus_NTo1_horikoshii -g yes -r 3 -e 3 -x yes -t yes #{$gttestdata}mgth/Pyrococcus_NTo1_horikoshii.xml.gz #{$gttestdata}mgth/Pyrococcus_horikoshii.txt.gz #{$gttestdata}mgth/Hits_NTo1_Pyrococcus_horikoshii.txt.gz", :maxtime => 100
  run "diff Pyrococcus_NTo1_horikoshii.xml #{$gttestdata}mgth/Pyrococcus_NTo1_horikoshii_ori_1.xml"
end

Name "gt mgth testdata mixed options 2"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Pyrococcus_NTo1_horikoshii -s 1.17 -n 3 -b -9.10 -q -3.0 -h -3.98 -l -1.3 -p 368 -f 217.0 -t yes -g yes -m yes -x yes -r 2 -e 3 -d 0.73 #{$gttestdata}mgth/Pyrococcus_NTo1_horikoshii.xml.gz #{$gttestdata}mgth/Pyrococcus_horikoshii.txt.gz #{$gttestdata}mgth/Hits_NTo1_Pyrococcus_horikoshii.txt.gz", :maxtime => 100
  run "diff Pyrococcus_NTo1_horikoshii.html #{$gttestdata}mgth/Pyrococcus_NTo1_horikoshii_ori_2.html"
end

Name "gt mgth testdata mixed options 3"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Pyrococcus_NTo1_horikoshii -s 5 -n 3.5 -b -14.23 -q -2.88 -h -5 -l -0.75 -p 295.75 -f 300.20 -t yes -g yes -m no -x no -r 3 -e 2 -a 25 -d 0.2 #{$gttestdata}mgth/Pyrococcus_NTo1_horikoshii.xml.gz #{$gttestdata}mgth/Pyrococcus_horikoshii.txt.gz #{$gttestdata}mgth/Hits_NTo1_Pyrococcus_horikoshii.txt.gz", :maxtime => 100
  run "diff Pyrococcus_NTo1_horikoshii.xml #{$gttestdata}mgth/Pyrococcus_NTo1_horikoshii_ori_3.xml"
end

Name "gt mgth testdata mixed options 4"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Pyrococcus_horikoshii -s 3.17 -n 2.32 -b -13.64 -q -3.45 -h -4.77 -l -1.5 -p 325.50 -f 225.25 -t yes -g yes -m no -x yes -r 1 -e 1 -d 0.21  #{$gttestdata}mgth/Pyrococcus_horikoshii.xml.gz #{$gttestdata}mgth/Pyrococcus_horikoshii.txt.gz #{$gttestdata}mgth/Hits_Pyrococcus_horikoshii.txt.gz", :maxtime => 100
  run "diff Pyrococcus_horikoshii.txt #{$gttestdata}mgth/Pyrococcus_horikoshii_ori_4.txt"
end

Name "gt mgth testdata mixed options 5"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Pyrococcus_horikoshii -s 2.55 -n 4.77 -b -15 -q -3 -h -3.27 -l -1.2 -p 177.13 -f 123 -t yes -g yes -m no -x yes -r 2 -e 2 -d 0.05 #{$gttestdata}mgth/Pyrococcus_horikoshii.xml.gz #{$gttestdata}mgth/Pyrococcus_horikoshii.txt.gz #{$gttestdata}mgth/Hits_Pyrococcus_horikoshii.txt.gz", :maxtime => 100
  run "diff Pyrococcus_horikoshii.html #{$gttestdata}mgth/Pyrococcus_horikoshii_ori_5.html"
end


Name "gt mgth testdata mixed options 6"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Pyrococcus_horikoshii -s 4.12 -n 2.9 -b -10 -q -7.5 -h -2 -l -0.5 -p 300 -f 211 -t yes -g yes -m no -x no -r 1 -e 3 -a 17 -d 0.13 #{$gttestdata}mgth/Pyrococcus_horikoshii.xml.gz #{$gttestdata}mgth/Pyrococcus_horikoshii.txt.gz #{$gttestdata}mgth/Hits_Pyrococcus_horikoshii.txt.gz", :maxtime => 100
  run "diff Pyrococcus_horikoshii.txt #{$gttestdata}mgth/Pyrococcus_horikoshii_ori_6.txt"
end

Name "gt mgth testdata mixed options 7"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Metagenome -s 3.17 -n 2.32 -b -13.64 -q -3.45 -h -4.77 -l -1.5 -p 325.50 -f 225.25 -t yes -g yes -m no -x yes -r 1 -e 1 -d 0.01 #{$gttestdata}mgth/Metagenome_NT.xml.gz #{$gttestdata}mgth/Metagenome.txt.gz #{$gttestdata}mgth/Hits_NT_Metagenome.txt.gz", :maxtime => 800
  run "diff Metagenome.txt #{$gttestdata}mgth/Metagenome_ori_1.txt"
end

Name "gt mgth testdata mixed options 8"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Metagenome -s 1.17 -n 3 -b -9.10 -q -3.0 -h -3.98 -l -1.3 -p 368 -f 217.0 -t yes -g yes -m no -x yes -r 2 -e 3 -d 0.05  #{$gttestdata}mgth/Metagenome_NTo1.xml.gz #{$gttestdata}mgth/Metagenome.txt.gz #{$gttestdata}mgth/Hits_NTo1_Metagenome.txt.gz", :maxtime => 800
  run "diff Metagenome.html #{$gttestdata}mgth/Metagenome_ori_2.html"
end

Name "gt mgth testdata mixed options 9"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Metagenome_454 -s 3.17 -n 2.32 -b -13.64 -q -3.45 -h -4.77 -l -1.5 -p 325.50 -f 225.25 -t yes -g yes -m no -x yes -r 1 -e 1 -d 0.01  #{$gttestdata}mgth/Metagenome_NT_454.xml.gz #{$gttestdata}mgth/Metagenome_454.txt.gz #{$gttestdata}mgth/Hits_NT_Metagenome_454.txt.gz", :maxtime => 800
  run "diff Metagenome_454.txt #{$gttestdata}mgth/Metagenome_454_ori_3.txt"
end

Name "gt mgth testdata mixed options 10"
Keywords "gt_mgth"
Test do
  run_test "#{$bin}gt mgth -o Metagenome_454 -s 1.17 -n 3 -b -9.10 -q -3.0 -h -3.98 -l -1.3 -p 368 -f 217.0 -t yes -g yes -m no -x yes -r 2 -e 3 -d 0.05  #{$gttestdata}mgth/Metagenome_NTo1_454.xml.gz #{$gttestdata}mgth/Metagenome_454.txt.gz #{$gttestdata}mgth/Hits_NTo1_Metagenome_454.txt.gz", :maxtime => 800
  run "diff Metagenome_454.html #{$gttestdata}mgth/Metagenome_454_ori_4.html"
  run "rm #{$gttestdata}mgth/*.gt_bsi #{$gttestdata}mgth/*.gt_bsr"
end

end
