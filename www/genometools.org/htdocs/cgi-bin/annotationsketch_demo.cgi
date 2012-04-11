#!/usr/bin/env ruby
#
# Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg
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
STYLE_FILE       = "#{GENOMETOOLS_PATH}/gtdata/sketch/default.style"
DEFAULT_ANNOTATION_FILE = "#{GENOMETOOLS_PATH}/testdata/standard_gene_as_tree.gff3"
SCRIPT_PATH      = "/var/www/servers/genometools.org/htdocs/cgi-bin"
UPLOAD_PATH      = "/tmp"
IMAGE_DIR        = "imgs"   # relative paths from SCRIPT_PATH please
MAXSIZE          = 2097152

$: << (GTRUBY_PATH)
require "gtruby"
require "cgi"

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
  <ul class="submenu">
    <li><a href="../annotationsketch.html#collapsing">Collapsing</a></li>
    <li><a href="../annotationsketch.html#styles">Styles</a></li>
    <li><a href="../trackselectors.html">Track assignment</a></li>
    <li><a href="../customtracks.html">Custom tracks</a></li>
    <li><a href="../annotationsketch.html#gtsketch">The <tt>gt sketch</tt> tool</a></li>
    <li><a href="../examples.html">Code examples</a></li>
    <li><a id="current" href="cgi-bin/annotationsketch_demo.cgi">Try it online</a></li>
    <li><a href="../libgenometools.html">API reference</a></li>
  </ul>
<li><a href="/cgi-bin/gff3validator.cgi">GFF3 validator</a></li>
<li><a href="../license.html">License</a></li>
</ul>
</div>
<div id="main">
  <h1><em>AnnotationSketch</em> online demo</h1>
  <p>Use this form to upload a GFF3 annotation file (up to 2MB) which is then
     drawn using an <em>AnnotationSketch</em>-based Ruby script and output to
     your browser. Some basic options, such as displayed sequence region and
     range, or the generated image width, can be set. For more options, the
     stand-alone <tt>gt sketch</tt> tool can be used.</p>
END

HTML_FOOTER = <<END
<div id="footer">
Copyright &copy; 2007-2011 Sascha Steinbiss. Last update: 2011-02-11
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
    <form action="annotationsketch_demo.cgi" method="POST" enctype="multipart/form-data">
    <table>
      <tr><td>Annotation file:</td>
        <td>
          <input type="radio" name="example" value="example" onclick="disable(this.form.file);" %s>Example file<br>
          <input type="radio" name="example" value="file" onclick="enable(this.form.file);" %s>Custom file: <input name="file" type="file" %s>
          <input type="hidden" name="submitted" value="true">
        </td>
      </tr>
      <tr>
        <td>Sequence region:</td>
        <td><input name="seqid" type="text" value="%s"></td>
      </tr>
      <tr><td></td><td style="font-size:small;">
                   (leave blank to use first sequence region in file)</td>
      </tr>
      <tr>
        <td>Range to display:</td>
        <td><input name="rangestart" size="10" type="text" value="%s">bp
         &ndash;
        <input name="rangeend" size="10" type="text" value="%s">bp</td>
        <tr><td></td><td style="font-size:small;">
                     (leave blank to show complete sequence region)</td>
      </tr>
      <tr>
        <td>Image width:</td>
        <td><input name="width" type="text" value="%s"> pixels</td>
      </tr>
      <tr>
        <td colspan=2><input value="Sketch this file!" type="submit"></td>
      </tr>
      </table>
    </form>
END

HTML_IMAGE = <<END
 <h2><em>AnnotationSketch</em> diagram of %s</h2>
 <div>
  <p><img src="%s" alt="AnnotationSketch diagram"></p>
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
    puts UPLOAD_FORM % [e_str, \
                        d_str1, \
                        d_str2, \
                        cgi['seqid'].string.strip_html, \
                        cgi['rangestart'].string.strip_html, \
                        cgi['rangeend'].string.strip_html, \
                        cgi['width'].string.strip_html]

    if cgi["example"].nil? or cgi["example"].send(read_method) == "example" then
      targetfilename = DEFAULT_ANNOTATION_FILE
      originalfilename = File.basename(DEFAULT_ANNOTATION_FILE)
    else
      ufile = cgi.params['file']
      if ufile.first.length == 0 then
        GT::gterror("No file was uploaded!")
      elsif ufile.first.length > MAXSIZE then
        GT::gterror("Your uploaded file was too large! This demo service
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
    feature_index = GT::FeatureIndexMemory.new()
    feature_index.add_gff3file(targetfilename)

    if cgi['width'].string != "" then
      width = cgi['width'].send(read_method).to_i
      if width < 70 then
        GT::gterror("Please set a width of more than 70 pixels!")
      end
    else
      width = 800
    end

    if cgi['seqid'].string.to_s != "" then
      seqid = cgi['seqid'].string
      seqids = feature_index.get_seqids
      if !seqids.include?(seqid) then
        if not cgi["example"].send(read_method) == "example" then
          File.unlink(targetfilename)
        end
        GT::gterror("Invalid sequence region '#{seqid}':
                     must be #{"one of" unless seqids.length == 1}
                     #{seqids.collect{|v|"&quot;#{v}&quot;"}.join(" or ")}!")
      end
    else
      seqid = feature_index.get_first_seqid
    end

    if cgi['rangestart'] and cgi['rangestart'].string != "" \
     and cgi['rangeend'] and  cgi['rangeend'].string != "" then
      range = GT::Range.malloc
      range.start = cgi['rangestart'].string.to_i
      range.end   = cgi['rangeend'].string.to_i
      if range.start >= range.end then
        GT::gterror("Invalid range, must be numeric and
                     <i>start</i> must be &lt; <i>end</i>!")
      end
    else
      range = feature_index.get_range_for_seqid(seqid)
    end

    style = GT::Style.new()
    style.load_file(STYLE_FILE)
    d = GT::Diagram.from_index(feature_index, seqid, range, style)
    l = GT::Layout.new(d, width, style)
    c = GT::CanvasCairoFile.new(style, width, l.get_height())
    l.sketch(c)
    c.to_file("#{SCRIPT_PATH}/#{IMAGE_DIR}/#{originalfilename}.png")
    puts HTML_IMAGE % [originalfilename.strip_html, \
                       "#{IMAGE_DIR}/#{originalfilename}.png"]
    if not cgi["example"].send(read_method) == "example" then
      File.unlink(targetfilename)
    end
  rescue Exception => err:
    puts "<h2>An error has occurred</h2><p>#{err}</p>"
  end
else
  puts UPLOAD_FORM % ['checked="checked"', '', 'disabled="disabled"', '', '', \
                      '', 800]
end

print HTML_FOOTER
