Name "gmap2gff3"
Keywords "scripts"
Test do
  run "#{$scriptsdir}gmap2gff3 #{$testdata}gmap2gff3_prob.gmap"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gmap2gff3_prob.out"
end
