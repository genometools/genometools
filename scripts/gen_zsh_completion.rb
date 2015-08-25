#!/usr/bin/env ruby

scriptdir = File.dirname(__FILE__)
help2comp = File.join(scriptdir, 'help2comp.py')

tools = %w{
  bed_to_gff3
  cds
  chain2dim
  chseqids
  clean
  compreads
  condenseq
  congruence
  convertseq
  csa
  dot
  dupfeat
  encseq
  encseq2spm
  eval
  extractfeat
  extractseq
  fastq_sample
  featureindex
  fingerprint
  genomediff
  gff3
  gff3_to_gtf
  gff3validator
  gtf_to_gff3
  hop
  id_to_md5
  interfeat
  loccheck
  ltrclustering
  ltrdigest
  ltrharvest
  matchtool
  matstat
  md5_to_id
  merge
  mergefeat
  mgth
  mkfeatureindex
  mkfmindex
  mmapandread
  orffinder
  packedindex
  prebwt
  readjoiner
  repfind
  scriptfilter
  select
  seq
  seqfilter
  seqids
  seqmutate
  seqorder
  seqstat
  seqtransform
  seqtranslate
  sequniq
  shredder
  shulengthdist
  simreads
  snpper
  speck
  splicesiteinfo
  splitfasta
  stat
  suffixerator
  tagerator
  tallymer
  tirvish
  uniq
  uniquesub
  wtree
}

completion =<<OUTER
#compdef gt

typeset -A opt_args

_arguments -C \\
  '1:cmd:->cmds' \\
  '*:: :->args' \\

case "$state" in
  (cmds)
    local arguments
    arguments=(
#{
`gt -help | #{help2comp} gt`.rstrip
}
    )

    local commands
    commands=(
    #{tools.join("\n    ")}
    )

    _describe -t commands 'command' commands && ret=0
    _arguments -s $arguments
  ;;
  (args)
    case $line[1] in

#{tools.map { |tool|
  opts = `gt #{tool} -help | #{help2comp} gt_#{tool}`.split("\n")

<<INNER
      (#{tool})
        local arguments
        arguments=(
    #{opts.join("\n    ")}
          '*:filename:_files'
        )
        _arguments -s $arguments
        ret=0
      ;;
INNER
}.join}

    esac
esac

return 1
OUTER

puts completion
