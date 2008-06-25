def aPrefix(prefix, list)
  return (list.collect {|x| prefix + x})
end

def prependTestdata(filelist)
  return aPrefix("#{$testdata}", filelist)
end

def paramList(paramsHash)
  if paramsHash != nil
    return paramsHash.to_a.flatten.compact
  else
    return []
  end
end

def runAndCheckPackedIndex(indexName,dbFiles, extraParams=Hash.new)
  params = {
    :create => { '-tis' => nil, '-des' => nil },
    :timeOuts => { :bdxcreat => 100, :suffixerator => 100,
      :chkintegrity => 400, :chksearch => 400, :trsuftab => 100,
      :mkctxmap => 100 },
    :bdx => {},
    :chksearch => { '-chksfxarray' => nil, '-nsamples' => '100' },
    :mkctxmap => { '-ctxilog' => -1 }
  }
  extraParams.keys.each do |key|
    params[key].merge!(extraParams[key]) if params.has_key?(key)
  end
  if indexName != nil
    params[:create]['-indexname'] = indexName
  else
    indexName = dbFiles.compact[0].sub(/.*\//,'')
  end
  #  puts('timeout: ', params[:timeOuts][:bdxcreat])
  if !extraParams.has_key?(:useSuftabTranslation) ||
      !extraParams[:useSuftabTranslation]
    run_test((["#{$bin}gt", '-debug', 'packedindex', 'mkindex'] +
              paramList(params[:create]) + paramList(params[:bdx]) +
              ['-db'] + dbFiles).join(' '),
             :maxtime => params[:timeOuts][:bdxcreat])
  end
  run_test((["#{$bin}gt", '-debug', 'suffixerator'] +
            paramList(params[:create]) +
            ['-bwt', '-suf', '-db'] + dbFiles).join(' '),
           :maxtime => params[:timeOuts][:suffixerator])
  if extraParams.has_key?(:useSuftabTranslation) &&
      extraParams[:useSuftabTranslation]
    run_test((["#{$bin}gt", '-debug', 'packedindex', 'trsuftab'] +
              paramList(params[:bdx]) + [indexName]).join(' '),
             :maxtime => params[:timeOuts][:trsuftab])
  end
  if extraParams.has_key?(:mkctxmap) &&
      extraParams[:mkctxmap]
    run_test((["#{$bin}gt", '-debug', 'packedindex', 'mkctxmap'] +
              paramList(params[:mkctxmap]) + [indexName]).join(' '),
             :maxtime => params[:timeOuts][:mkctxmap])    
  end
  run_test(["#{$bin}gt", '-debug', 'packedindex', 'chkintegrity',
            '-ticks', '1000', indexName].join(' '),
           :maxtime => params[:timeOuts][:chkintegrity])
  run_test((["#{$bin}gt", '-debug', 'packedindex', 'chksearch'] +
            paramList(params[:chksearch]) + [indexName]).join(' '),
           :maxtime => params[:timeOuts][:chksearch])
end

Name "gt packedindex check tools for simple sequences"
Keywords "gt_packedindex"
Test do
  allfiles = prependTestdata(["RandomN.fna","Random.fna","Atinsert.fna",
                              "TTT-small.fna","trna_glutamine.fna",
                              "Random-Small.fna","Duplicate.fna"])
  runAndCheckPackedIndex('miniindex', allfiles)
end

Name "gt packedindex check tools for simple sequences w/o locate"
Keywords "gt_packedindex"
Test do
  allfiles = prependTestdata(["RandomN.fna","Random.fna","Atinsert.fna",
                              "TTT-small.fna","trna_glutamine.fna",
                              "Random-Small.fna","Duplicate.fna"])
  runAndCheckPackedIndex('miniindex', allfiles,
                         :bdx => { '-locfreq' => 0 },
                         :chksearch => { '-chksfxarray' => 'no' })
end

Name "gt packedindex check tools for simple sequences with sprank"
Keywords "gt_packedindex"
Test do
  allfiles = prependTestdata(["RandomN.fna","Random.fna","Atinsert.fna",
                              "TTT-small.fna","trna_glutamine.fna",
                              "Random-Small.fna","Duplicate.fna"])
  runAndCheckPackedIndex('miniindex', allfiles,
                         :bdx => { '-sprank' => nil },
                         :chksearch => { '-full-lfmap' => nil },
                         :timeOuts => { :chksearch => 800 })
end

Name "gt packedindex check tools for revcom sequence with sprank"
Keywords "gt_packedindex"
Test do
  allfiles = prependTestdata(["Atinsert.fna"])
  runAndCheckPackedIndex('miniindex', allfiles,
                         :create => { '-dir' => 'rev' },
                         :bdx => { '-sprank' => nil },
                         :chksearch => { '-full-lfmap' => nil },
                         :timeOuts => { :chksearch => 800 })
end

Name "gt packedindex check tools for simple sequences with context"
Keywords "gt_packedindex"
Test do
  allfiles = prependTestdata(["RandomN.fna","Random.fna","Atinsert.fna",
                              "TTT-small.fna","trna_glutamine.fna",
                              "Random-Small.fna","Duplicate.fna"])
  runAndCheckPackedIndex('miniindex', allfiles,
                         :bdx => { '-sprank' => nil, '-ctxilog' => -1 },
                         :chksearch => { '-chkcontext' => nil },
                         :timeOuts => { :chksearch => 800 })
end

Name "gt packedindex check tools for protein sample"
Keywords "gt_packedindex"
Test do
  runAndCheckPackedIndex(nil, prependTestdata(['sw100K2.fsa']),
                         :bdx => { '-bsize' => 1 })
end

Name "gt packedindex check tools for simple sequences, tr-mode"
Keywords "gt_packedindex"
Test do
  allfiles = prependTestdata(["RandomN.fna","Random.fna","Atinsert.fna",
                              "TTT-small.fna","trna_glutamine.fna",
                              "Random-Small.fna","Duplicate.fna"])
  runAndCheckPackedIndex('miniindex', allfiles,
                         :useSuftabTranslation => true,
                         :bdx => { '-sprank' => nil },
                         :mkctxmap => {})
end

Name "gt packedindex check tools for boundary-case sequences"
Keywords "gt_packedindex"
Test do
  allfiles = prependTestdata(['Random160.fna', 'Random159.fna'])
  allfiles.each do |file|
    runAndCheckPackedIndex(nil, [file], :bdx => { '-blbuck' => 20 })
  end
  runAndCheckPackedIndex(nil, prependTestdata(['Random80.fna']),
                         :bdx => { '-bsize' => 10 })
end

if $gttestdata then
  Name "gt packedindex check tools for chr01 yeast"
  Keywords "gt_packedindex"
  Test do
    runAndCheckPackedIndex('chr01.19960731',
      ["#{$gttestdata}ltrharvest/s_cer/chr01.19960731.fsa.gz"])
  end

  Name "gt packedindex check tools for at1MB"
  Keywords "gt_packedindex at1MB"
  Test do
    runAndCheckPackedIndex('at1MB', ["#{$gttestdata}Iowa/at1MB"],
                           :timeOuts => { :bdxcreat => 400,
                             :suffixerator => 400,
                             :chkintegrity => 800, :chksearch => 400 })
  end

  dmelFiles = aPrefix("#{$gttestdata}ltrharvest/d_mel/",
                      [
                       'X_genomic_dmel_RELEASE3-1.FASTA.gz',
                       '2L_genomic_dmel_RELEASE3-1.FASTA.gz',
                       '2R_genomic_dmel_RELEASE3-1.FASTA.gz',
                       '3L_genomic_dmel_RELEASE3-1.FASTA.gz',
                       '3R_genomic_dmel_RELEASE3-1.FASTA.gz',
                       '4_genomic_dmel_RELEASE3-1.FASTA.gz',
                       ])
  dmelFiles.each do |file|
    shortFileName = file.match('[^/]*$')[0]
    maxDisplay = 12
    if shortFileName.length > maxDisplay
      shortFileName = shortFileName[0,maxDisplay - 4] + '...'
    end
    Name 'gt packedindex check tools for d.melanogaster (' +
      shortFileName + ')'
    Keywords 'gt_packedindex'
    Test do
      runAndCheckPackedIndex('dmel',
                             [file],
                             :timeOuts =>
                             { :trsuftab => 7200, :suffixerator => 28800,
                               :chkintegrity => 3200, :chksearch => 800 },
                             :useSuftabTranslation => true)
    end
  end
end
