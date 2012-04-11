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
DEFAULT_ANNOTATION_FILE = "#{GENOMETOOLS_PATH}/testdata/standard_gene_as_tree.gff3"
SCRIPT_PATH      = "/var/www/servers/genometools.org/htdocs/cgi-bin"
UPLOAD_PATH      = "/tmp"
MAXSIZE          = 2097152

$: << (GTRUBY_PATH)
require "gtruby"
require "cgi"
require 'tempfile'

puts "Content-type: text/html"
puts ""

cgi = CGI.new("html4")

class String
# String#strip_html - Removes HTML tags from a string.
# Author:: Rob Pitt
# Removes HTML tags from a string. Allows you to specify some tags to be kept.
  def strip_html( allowed = [] )
    re = if allowed.any?
      Regexp.new(
        %(<(?!(\\s|\\/)*(#{
          allowed.map {|tag| Regexp.escape( tag )}.join( "|" )
        })( |>|\\/|'|"|<|\\s*\\z))[^>]*(>+|\\s*\\z)),
        Regexp::IGNORECASE | Regexp::MULTILINE, 'u'
      )
    else
      /<[^>]*(>+|\s*\z)/m
    end
    gsub(re,'')
  end
end

HTML_HEADER = <<END
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
<title>The AnnotationSketch module</title>
<link rel="stylesheet" type="text/css" href="../style.css">
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
<li><a href="../license.html">License</a></li>
</ul>
</div>
<div id="main">
  <h1><em>GFF3 online validator</h1>
  <p>Use this form to upload a GFF3 annotation file (up to 2MB) which is then
     validated against the current <a href="#">Sequence Ontology OBO file</a>.</p>
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
    <form action="validator.cgi" method="POST" enctype="multipart/form-data">
    <table>
      <input type="hidden" name="submitted" value="true">
      <tr><td>Annotation file:</td>
        <td>
          <input type="radio" name="example" value="example" onclick="disable(this.form.file);" %s>Example file<br>
          <input type="radio" name="example" value="file" onclick="enable(this.form.file);" %s>Custom file: <input name="file" type="file" %s>
        </td>
      </tr>
      <tr><td>Validation options:</td>
        <td>
          <input type="radio" name="mode" value="default" %s>Default (SO specification compliant)<br>
          <input type="radio" name="mode" value="strict" %s>Strict mode (stricter than SO specification)<br>
          <input type="radio" name="mode" value="tidy" %s>Tidy mode (tries to fix errors)
        </td>
      </tr>
      <tr>
        <td colspan=2><input value="Validate this file!" type="submit"></td>
      </tr>
      </table>
    </form>
END

HTML_IMAGE = <<END
 <h2>Validation results for %s</h2>
 <div>
  <p>%s</p>
 </div>
END

puts HTML_HEADER

if cgi.params.has_key?('submitted') then
  # CGI parameters behave differently with uploaded file size.
  # StringIO.read does not seem to work right either.
  # account for that by checking for the type of the parameters
  if cgi["example"].kind_of?(StringIO) then
    read_method = :string
  else
    read_method = :read
  end
  begin
    if cgi["example"].nil? or cgi["example"].send(read_method) == "example" then
      e_str = 'checked="checked"'
      d_str1 = ""
      d_str2 = 'disabled="disabled"'
    else
      e_str = ""
      d_str1 = 'checked="checked"'
      d_str2 = ""
    end
    if cgi["mode"].nil? or cgi["mode"].send(read_method) == "default" then
      m_str1 = 'checked="checked"'
      m_str2 = ""
      m_str3 = ""
    end
    if cgi["mode"].send(read_method) == "tidy" then
      m_str1 = ""
      m_str2 = ""
      m_str3 = 'checked="checked"'
    end
    if cgi["mode"].send(read_method) == "strict" then
      m_str1 = ""
      m_str2 = 'checked="checked"'
      m_str3 = ""
    end

    puts UPLOAD_FORM % [e_str,  \
                        d_str1, \
                        d_str2, \
                        m_str1, \
                        m_str2, \
                        m_str3]

    if cgi["example"].nil? or cgi["example"].send(read_method) == "example" then
      targetfilename = DEFAULT_ANNOTATION_FILE
      originalfilename = File.basename(DEFAULT_ANNOTATION_FILE)
    else
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
    end

    errfile = Tempfile.new('validator')
    $stderr.reopen(errfile)

    stream = GT::GFF3InStream.new(targetfilename)
    if cgi["mode"].send(read_method) == "tidy" then
      stream.enable_tidy_mode
    elsif cgi["mode"].send(read_method) == "strict" then
      stream.enable_strict_mode
    end
    gn = stream.next_tree()

    while (gn) do
      gn = stream.next_tree()
    end
    puts "<h2>Validation successful!</h2>"
    errfile.rewind
    puts "<p>#{errfile.read}</p>"
    errfile.close
    errfile.unlink

  rescue Exception => err:
    puts "<h2>Validation unsuccessful!</h2><p>#{err}</p>"
  end
else
  puts UPLOAD_FORM % ['checked="checked"', '', 'disabled="disabled"', 'checked="checked"', '' ,'']
end

print HTML_FOOTER
