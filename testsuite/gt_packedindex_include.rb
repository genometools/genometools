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

def runAndCheckPackedIndex(indexName,dbFiles,createExtraParams,bdxParams,
                           timeOutExtras=nil)
  createParams = { '-dna' => nil, '-tis' => nil, '-des' => nil }
  createParams.merge!(createExtraParams) if (createExtraParams != nil)
  if indexName != nil
    createParams['-indexname'] = indexName
  else
    indexName = dbFiles.compact[0].sub(/.*\//,'')
  end
  timeOuts = { 'bdxcreat' => 100, 'suffixerator' => 100,
    'chkintegrity' => 400, 'chksearch' => 400 }
  timeOuts.merge!(timeOutExtras) if (timeOutExtras != nil)
#  puts('timeout: ', timeOuts['bdxcreat'])
  run_test((["#{$bin}gt", 'packedindex', 'mkindex'] +
            paramList(createParams) + paramList(bdxParams) +
            ['-db'] + dbFiles).join(' '), :maxtime => timeOuts['bdxcreat'])
  run_test((["#{$bin}gt", 'suffixerator'] +
            paramList(createParams) +
            ['-bwt', '-suf', '-db'] + dbFiles).join(' '),
           :maxtime => timeOuts['suffixerator'])
  run_test(["#{$bin}gt", 'packedindex', 'chkintegrity', '-ticks', '1000',
            indexName].join(' '), :maxtime => timeOuts['chkintegrity'])
  run_test(["#{$bin}gt", 'packedindex', 'chksearch', '-chksfxarray',
            '-nsamples', '100', indexName].join(' '),
           :maxtime => timeOuts['chksearch'])
end

Name "gt packedindex check tools for simple sequences"
Keywords "gt_packedindex"
Test do
  allfiles = prependTestdata(["RandomN.fna","Random.fna","Atinsert.fna",
                              "TTT-small.fna","trna_glutamine.fna",
                              "Random-Small.fna","Duplicate.fna"])
  runAndCheckPackedIndex('miniindex', allfiles, nil, nil)
end

Name "gt packedindex check tools for boundary-case sequences"
Keywords "gt_packedindex"
Test do
  allfiles = prependTestdata(['Random160.fna', 'Random159.fna'])
  allfiles.each do |file|
    runAndCheckPackedIndex(nil, [file], nil, { '-blbuck' => 20 })
  end
  runAndCheckPackedIndex(nil, prependTestdata(['Random80.fna']),
                         nil, { '-bsize' => 10 })
end

if $gttestdata then
  Name "gt packedindex check tools for chr01 yeast"
  Keywords "gt_packedindex"
  Test do
    runAndCheckPackedIndex('chr01.19960731',
      ["#{$gttestdata}ltrharvest/s_cer/chr01.19960731.fsa.gz"], nil, nil)
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
    Keywords 'gt_packedindex gt_notworking'
    Test do
      runAndCheckPackedIndex('dmel',
                             [file], nil, nil,
                             { 'bdxcreat' => 7200, 'suffixerator' => 7200,
                               'chkintegrity' => 1600, 'chksearch' => 800 })
    end
  end
end
