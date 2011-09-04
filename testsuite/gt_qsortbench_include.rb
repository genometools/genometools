methods = ["thomas","system","inlinedptr","inlinedarr",
           "direct","radix","radixrec","radixit"]

lenlist=[1000,2000,4000,1000000,2000000]

methods.each do |met|
  lenlist.each do |len|
    Name "gt qsortbench #{met}"
    Keywords "gt_qsortbench"
    Test do
      run "#{$bin}gt dev qsortbench -impl #{met} -size #{len}"
    end
  end
end
