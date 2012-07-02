#
# Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
# common functions for genomediff scripts
#

require 'tempfile'
require 'fileutils'

module Genomediff
  def Genomediff.reduceN(file,removeN)
    alphabet = [?a, ?A, ?c, ?C, ?g, ?G, ?t, ?T]
    newline = [?\n]
    charcount = 0
    header = false
    wildcard = false
    first = true
    newfile = ""
    dirname = File.dirname(file)
    File.open(file, 'r') {|fp|
      basename = File.basename(file, ".*")
      newfile = File.join(dirname, basename + "_redN.fas")
      #return newfile if File.exist?(newfile)
      File.open(newfile, 'w') {|nfp|
        while cc = fp.getc
          if charcount == 80
            nfp.putc ?\n
            charcount = 0
          end
          if not header and cc == ?>
            header = true
            nfp.putc ?\n unless first
            first = false
            nfp.putc cc
            next
          elsif header
            if cc == ?\n
              header = false
              charcount = 0
            end
            nfp.putc cc
            next
          end
          if cc == ?\n
            next
          end
          if not wildcard and not alphabet.include?(cc)
            wildcard = true
            unless removeN
              nfp.putc ?N 
              charcount += 1
            end
            next
          elsif wildcard
            if alphabet.include?(cc)
              wildcard = false
              nfp.putc cc
              charcount += 1
              next
            end
          else
            nfp.putc cc
            charcount += 1
          end
        end
        nfp.putc ?\n
      }
    }
    exit 1 unless File.exist?(newfile)
    return newfile
  end

  def Genomediff.reverse_and_concat(file)
    tf = Tempfile.new("gt_rev")
    tmpfile = tf.path
    tf.close
    dirname = File.dirname(file)
    `$GTDIR/bin/gt            \
     convertseq -o #{tmpfile} \
     -force -r #{file}`
    basename = File.basename(file, ".*")
    newfile = File.join(dirname, basename + "_plus_rev.fas")
    `cat #{file} #{tmpfile} > #{newfile}`
    exit 1 unless File.exist?(newfile)
    File.unlink(tmpfile)
    return newfile
  end

  def Genomediff.pck_index(files,idxname,parameter)
    return `$GTDIR/bin/gt                  \
            packedindex mkindex            \
            -mirrored                      \
            -dna -dir rev -ssp -sprank -pl \
            -db #{files}                   \
            -indexname #{idxname}          \
            #{parameter}`
  end

  def Genomediff.esa_index(files,idxname,parameter)
    return `$GTDIR/bin/gt suffixerator \
            -mirrored                  \
            -dna -suf -tis -lcp -ssp   \
            -db #{files}               \
            -indexname #{idxname}      \
            #{parameter}`
  end

  def Genomediff.pck_genomediff(idxname,parameter)
    return `$GTDIR/bin/gt genomediff \
                     -indextype pck #{parameter} #{idxname}`
  end

  def Genomediff.esa_genomediff(idxname,parameter)
    return `$GTDIR/bin/gt genomediff \
                     -indextype esa #{parameter} #{idxname}`
  end
end

module Eval

  def Eval.pck_index(code,
                revfiles,
                log_fp,
                verbose,
                bsize,
                parts)
    puts "# INDEX PCK" if verbose

    idxname = File.join(File.dirname(revfiles[0]), code+"_pck")

    log_fp.puts "# RUN #{code} -pck -parts #{parts} -bsize #{bsize}"
    puts "# RUN #{code} -pck -parts"+
         " #{parts} -bsize #{bsize}" if verbose

    output = Genomediff.pck_index(revfiles.join(" "),bsize,parts,idxname)

    /^# TIME overall \d+.\d+.*$/.match output
    log_fp.puts $&
    /^# space peak .*$/.match output
    log_fp.puts $&
    puts output if verbose

    return idxname
  end

  def Eval.esa_index(code,
                revfiles,
                log_fp,
                verbose,
                parts)
    puts "# INDEX ESA" if verbose
    m = code.match(/(\d+_\d+)_(\d+\.\d+)_(\d+)/)
    filecode, div, n = m[1], m[2], m[3]

    idxname = File.join(File.dirname(revfiles[0]), code+"_esa")

    log_fp.puts "# RUN #{filecode} d:#{div} n:#{n} -parts #{parts}"
    puts "# RUN #{filecode} d:#{div} n:#{n}"+
         " -esa" if verbose

    output = Genomediff.esa_index(revfiles.join(" "),parts,idxname)

    /^# TIME overall \d+.\d+.*$/.match output
    log_fp.puts $&
    /^# space peak .*$/.match output
    log_fp.puts $&
    puts output if verbose

    return idxname
  end

  def Eval.pck_genomediff(code,
                     idxname,
                     parameter,
                     log_fp,
                     verbose)
    puts "# GENOMEDIFF PCK" if verbose
    m = code.match(/(\d+_\d+)_(\d+\.\d+)_(\d+)/)
    filecode, div, n = m[1], m[2], m[3]
    log_fp.puts "# RUN #{filecode} d:#{div} n:#{n} #{parameter}"
    puts "# RUN #{filecode} d:#{div} n:#{n} #{parameter}" if verbose

    output = Genomediff.pck_genomediff(idxname,parameter)

    /^# TIME overall \d+.\d+.*$/.match output
    log_fp.puts $&
    /^# space peak .*$/.match output
    log_fp.puts $&
    /^# mmap .*$/.match output
    log_fp.puts $&
    puts output if verbose
  end

  def Eval.esa_genomediff(code, idxname, parameter, log_fp, verbose)
    puts "# GENOMEDIFF ESA" if verbose
    m = code.match(/(\d+_\d+)_(\d+\.\d+)_(\d+)/)
    filecode, div, n = m[1], m[2], m[3]
    log_fp.puts "# RUN #{filecode} d:#{div} n:#{n} #{parameter}"
    puts "# RUN #{filecode} d:#{div} n:#{n} #{parameter}" if verbose

    output = Genomediff.esa_genomediff(idxname,parameter)

    /^# TIME overall \d+.\d+.*$/.match output
    log_fp.puts $&
    /^# space peak .*$/.match output
    log_fp.puts $&
    /^# mmap .*$/.match output
    log_fp.puts $&
    puts output if verbose
  end

  def Eval.run_pck_index_creation(code, revfiles, log_fp, verbose)
    self.pck_index(code,
              revfiles,
              log_fp,
              verbose,
              8,
              1)
  end

  def Eval.run_pck_index_bsize(code, revfiles, log_fp, verbose, bsize)
    self.pck_index(code,
              revfiles,
              log_fp,
              verbose,
              bsize,
              1)
  end

  def Eval.run_pck_index_parts(code, revfiles, log_fp, verbose, parts)
    self.pck_index(code,
              revfiles,
              log_fp,
              verbose,
              8,
              parts)
  end


  def Eval.run_esa_index_creation(code, revfiles, log_fp, verbose)
    self.esa_index(code,
              revfiles,
              log_fp,
              verbose,
              1)
  end

  def Eval.run_esa_index_parts(code, revfiles, log_fp, verbose,parts)
    self.esa_index(code,
              revfiles,
              log_fp,
              verbose,
              parts)
  end

  def Eval.run_genomediff_pck(code, idxname, parameter, log_fp, verbose)
    parameter += " -v " if verbose
    self.pck_genomediff(code,
                   idxname,
                   parameter,
                   log_fp,
                   verbose)
  end

  def Eval.run_genomediff_esa(code, idxname, parameter, log_fp, verbose)
    parameter += " -v " if verbose
    self.esa_genomediff(code,
                   idxname,
                   parameter,
                   log_fp,
                   verbose)
  end
end
