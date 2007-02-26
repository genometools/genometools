Name "$GT_ENV_OPTIONS parsing (-help)"
Keywords "gt_env_options"
Test do
  run("env GT_ENV_OPTIONS=-help #{$bin}gt", :retval => 1)
  grep $last_stdout, /spacepeak/
end

Name "$GT_ENV_OPTIONS parsing (-spacepeak)"
Keywords "gt_env_options"
Test do
  run "env GT_ENV_OPTIONS=-spacepeak #{$bin}gt gff3 #{$testdata}standard_gene_as_tree.gff3"
  grep $last_stdout, /space peak in megabytes/
end
