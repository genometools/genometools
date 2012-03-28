methods = ["radixinplace","thomas",
           "system","inlinedptr","inlinedarr","direct",
           "radixlinsmall","radixlinlarge",
           "radixlinsmall -parts 2","radixlinlarge -parts 2",
           "radixlinsmall -parts 3","radixlinlarge -parts 3",
           "radixlinsmall -parts 5","radixlinlarge -parts 5",
           "radixrec"]

lenlist=[10,20,30,1000,2000,4000,1000000,2000000]

methods.each do |met|
  Name "gt sortbench #{met}"
  Keywords "gt_sortbench"
  Test do
    lenlist.each do |len|
      if met.match(/^radixlin/) and not met.match(/-parts/)
        run "#{$bin}gt -j 2 dev sortbench -impl #{met} -size #{len}"
        run "#{$bin}gt -j 3 dev sortbench -impl #{met} -size #{len}"
      end
      if met.match(/^radixinplace/)
        ["","-j 4"].each do |opt|
          run "#{$bin}gt #{opt} dev sortbench -impl #{met} -size #{len} -maxval 1000"
          run "#{$bin}gt #{opt} dev sortbench -impl #{met} -size #{len}"
        end
      else
        run "#{$bin}gt dev sortbench -impl #{met} -size #{len}"
        run "#{$bin}gt dev sortbench -impl #{met} -size #{len} -maxval 10000"
      end
    end
  end
end
