methods = ["radixinplace","radixlsb",
           "thomas","system","inlinedptr","inlinedarr","direct","dual-pivot",
           "radixkeypair"]

lenlist=[10,20,30,1000,2000,4000,1000000,2000000]

methods.each do |met|
  Name "gt sortbench #{met}"
  Keywords "gt_sortbench"
  Test do
    lenlist.each do |len|
      if met.match(/^radixinplace/)
        ["","-j 4"].each do |opt|
          run "#{$bin}gt #{opt} dev sortbench -verify -impl #{met} -size #{len} -maxval 1000"
          run "#{$bin}gt #{opt} dev sortbench -verify -impl #{met} -size #{len}"
        end
      else
        run "#{$bin}gt dev sortbench -verify -impl #{met} -size #{len}"
        run "#{$bin}gt dev sortbench -verify -impl #{met} -size #{len} -maxval 10000"
      end
    end
  end
end
