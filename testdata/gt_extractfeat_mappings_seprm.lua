function mapping(sequence_region)
  return os.getenv("GT_TESTDATA").."gt_extractfeat_mappings_seprm_"..sequence_region..".fas"
end

