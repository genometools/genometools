SCRIPT_FILTER_FIELDS = ["name", "author", "version", "email", "short_description", "description"]

Name "gt scriptfilter metadata"
Keywords "gt_scriptfilter"
Test do
  run "#{$bin}gt scriptfilter -scriptname false " + \
      " #{$testdata}/gtscripts/filter_metadata_test_all_strings.lua"
  run_test "diff #{$testdata}/script_filter_output.txt #{last_stdout}"
end

["", "-oneline"].each do |par|
  Name "gt scriptfilter metadata as functions #{par}"
  Keywords "gt_scriptfilter"
  Test do
    run "#{$bin}gt scriptfilter -scriptname false #{par} " + \
        " #{$testdata}/gtscripts/filter_metadata_test_all_strings.lua " + \
        " > strings.txt"
    SCRIPT_FILTER_FIELDS.each do |field|
      run_test "#{$bin}gt scriptfilter -scriptname false #{par} " + \
               "#{$testdata}/gtscripts/filter_metadata_test_#{field}_function.lua"
      run_test "diff strings.txt #{last_stdout}"
    end
  end
end
  
