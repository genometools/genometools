methods = ["thomas","system","inlinedptr","inlinedarr",
           "direct","radixlinsmall","radixlinlarge","radixrec","radixdiv"]

lenlist=[10,20,30,1000,2000,4000,1000000,2000000]

maxvallist = nil
if Kernel.system("#{$bin}gt -64bit") then
  maxvallist = [922332036854775807,10000]
else
  maxvallist = [3432223432,10000]
end

methods.each do |met|
  lenlist.each do |len|
    maxvallist.each do |maxv|
      Name "gt sortbench #{met}"
      Keywords "gt_sortbench"
      Test do
        run "#{$bin}gt dev sortbench -impl #{met} -size #{len} -maxval #{maxv}"
      end
    end
  end
end
