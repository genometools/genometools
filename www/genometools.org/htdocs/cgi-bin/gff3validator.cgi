#!/usr/bin/env ruby
#
# Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2012 Center for Bioinformatics, University of Hamburg
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

GENOMETOOLS_PATH = "/home/satta/genometools_for_web"
GTRUBY_PATH      = "#{GENOMETOOLS_PATH}/gtruby"
# the LD_LIBRARY_PATH has to be set externally to "#{GENOMETOOLS_PATH}/lib"!
SCRIPT_PATH      = "/var/www/servers/genometools.org/htdocs/cgi-bin"
UPLOAD_PATH      = "/tmp"
MAXSIZE          = 52428800

$: << (GTRUBY_PATH)
require "gtruby"
require "cgi"
require 'tempfile'

puts "Content-type: text/html"
puts ""

cgi = CGI.new("html4")

HTML_HEADER = <<END
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
<title>GFF3 Online Validator</title>
<link rel="stylesheet" type="text/css" href="../style.css">
<style type="text/css">
table.padded-table td {
  padding: 0px 7px;
	}
</style>
<script type="text/javascript">
 function disable(field) {
   field.disabled = true;
   field.readonly = true;
 }

 function enable(field) {
   field.disabled = false;
   field.readonly = false;
 }
</script>
</head>
<body>
<div id="menu">
<ul>
<li><a href="../index.html">Overview</a></li>
<li><a href="../pub/">Download</a></li>
<li><a href="gitweb.cgi?p=genometools.git;a=summary">Browse source</a></li>
<li><a href="../mailman/listinfo/gt-users">Mailing list</a></li>
<li><a href="http://genometools.lighthouseapp.com/">Issue tracker</a></li>
<li><a href="../design.html">Design</a></li>
<li><a href="../libgenometools.html">C API</a></li>
<li><a href="../docs.html"><tt>gtscript</tt> docs</a></li>
<li><a href="../annotationsketch.html"><tt>AnnotationSketch</tt></a></li>
<li><a id="current" href="/cgi-bin/gff3validator.cgi">GFF3 validator</a></li>
<li><a href="../manuals.html">Manuals</a></li>
<li><a href="../license.html">License</a></li>
</ul>
</div>
<div id="main">
  <h1><em>GFF3 online validator</h1>
  <p>Use this form to upload a GFF3 annotation file (up to 50 MB, can be .gz or
.bz2 compressed) which is then validated against the <a href="http://www.sequenceontology.org/gff3.shtml">GFF3 specification</a> using the current <a href="http://song.cvs.sourceforge.net/song/ontology/so.obo?view=log">Sequence Ontology OBO file</a>.</p>
END

HTML_FOOTER = <<END
<div id="footer">
Copyright &copy; 2012 Sascha Steinbiss. Last update: 2012-04-11
</div>
</div>
<!-- Piwik -->
<script type="text/javascript">
var pkBaseURL = (("https:" == document.location.protocol) ?  "https://gremme.org/piwik/" : "http://gremme.org/piwik/");
document.write(unescape("%3Cscript src='" + pkBaseURL + "piwik.js' type='text/javascript'%3E%3C/script%3E"));
</script><script type="text/javascript">
try {
var piwikTracker = Piwik.getTracker(pkBaseURL + "piwik.php", 5);
piwikTracker.trackPageView();
piwikTracker.enableLinkTracking();
} catch( err ) {}
</script><noscript><p><img src="http://gremme.org/piwik/piwik.php?idsite=5" style="border:0" alt="" /></p></noscript>
<!-- End Piwik Tracking Tag -->
</body>
</html>
END

UPLOAD_FORM = <<END
    <form action="gff3validator.cgi" method="POST" enctype="multipart/form-data">
    <table>
      <input type="hidden" name="submitted" value="true">
      <tr><td>Annotation file:</td>
        <td>
          <input name="file" type="file">
        </td>
      </tr>
      <tr><td></td>
        <td>
          <input type="checkbox" name="tidy" %s>enable &quot;tidy&quot; mode (tries to fix errors) <br />
          <input type="checkbox" name="show" %s>also report content of affected GFF lines <br />
        </td>
      </tr>
      <tr>
        <td colspan=2><input value="Validate this file!" type="submit"></td>
      </tr>
      </table>
    </form>
    <p>This GFF3 validator is part of the <em>GenomeTools</em> distribution which you
can <a href="http://genometools.org/pub">download</a> to your computer.
Use the <tt>gff3validator</tt> tool to validate your own &ndash; possibly larger &ndash; GFF3 files and the <tt>gff3</tt> tool with option <tt>-tidy</tt> to tidy them up (<tt>-help</tt> shows further options).</p>
END

puts HTML_HEADER

def print_gff3line(line, lineno, red = false)
  puts "<tr>"
  puts "<td style='font-family: Verdana, Geneva, Arial, sans-serif;'>#{lineno}</td>"
  puts "<td>#{if red then "&rarr;" end}</td>"
  line.split("\t").each do |v|
    print "<td #{if red then "style='color:red;'" end}>#{v}</td>"
  end
  puts
  puts "</tr>"
end

if cgi.params.has_key?('submitted') then
  tidy = (cgi.params["tidy"].length > 0)
  show = (cgi.params["show"].length > 0)
  begin
    show_checked = ''
    tidy_checked = ''
    if tidy then
      tidy_checked = 'checked="checked"'
    end
    if show then
      show_checked = 'checked="checked"'
    end

    puts UPLOAD_FORM % [tidy_checked, show_checked]

    ufile = cgi.params['file']
    if ufile.first.length == 0 then
      GT::gterror("No file was uploaded!")
    elsif ufile.first.length > MAXSIZE then
      GT::gterror("Your uploaded file was too large! This service
                   supports annotation files up to #{(MAXSIZE/1024)}KB. Please
                   upload a smaller file.")
    end
    # mangle file path to avoid directory traversal attacks
    originalfilename = ufile.first.original_filename
    truncated_filename = File.basename(File.expand_path(originalfilename))
    targetfilename = "#{UPLOAD_PATH}/#{truncated_filename}"
    File.open(targetfilename, "w+") do |file|
      file.write(ufile.first.read)
    end

    errfile = Tempfile.new('validator')
    $stderr.reopen(errfile)

    checker = GT::TypeCheckerOBO.new("#{GENOMETOOLS_PATH}/gtdata/obo_files/so.obo")

    stream = GT::GFF3InStream.new(targetfilename, false)
    stream.set_type_checker(checker)
    if tidy then
      stream.enable_tidy_mode
    end

    gn = stream.next_tree()
    while (gn) do
      gn = stream.next_tree()
    end
    errfile.flush

    puts "<h2 style='color:green;'>Validation successful!</h2>"
    errfile.rewind
    entries = errfile.readlines
    if entries.length > 0 then
      puts "<p>#{entries.length} issue#{"s" unless entries.length == 1} remaining:</p>"
      file = File.open(targetfilename).readlines
      entries.each do |entry|
        puts "<div style='background-color:#EFEFEF; padding:0px 5px;'>"
        puts "<p>#{entry}</p>"
        if show and (m = /line ([0-9]+)/.match(entry)) then
          linenumber = m[1].to_i
          puts "<p style='background-color:#EEEEEE;'><pre><table class='padded-table'>"
          if linenumber-1 > 0 then
            print_gff3line(file[linenumber-2].chomp,linenumber-1)
          end
          print_gff3line(file[linenumber-1].chomp, linenumber, true)
          if linenumber <= (file.length)-1 then
            print_gff3line(file[linenumber].chomp, linenumber+1)
          end
          puts "</table></pre></p>"
        end
        puts "</div>"
      end
    end
    errfile.close
    errfile.unlink

  rescue Exception => err
    puts "<h2 style='color:red;'>Validation unsuccessful!</h2><p>#{err}</p>"
  end
else
  puts UPLOAD_FORM % ['', 'checked="checked"']
end

print HTML_FOOTER
