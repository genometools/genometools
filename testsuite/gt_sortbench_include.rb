methods = ["thomas","system","inlinedptr","inlinedarr",
           "direct","radixlinsmall","radixlinlarge","radixrec","radixdiv"]

lenlist=[10,20,30,1000,2000,4000,1000000,2000000]

methods.each do |met|
  lenlist.each do |len|
    Name "gt sortbench #{met}"
    Keywords "gt_sortbench"
    Test do
      run "#{$bin}gt dev sortbench -impl #{met} -size #{len}"
      run "#{$bin}gt dev sortbench -impl #{met} -size #{len} -maxval 10000"
    end
  end
end
