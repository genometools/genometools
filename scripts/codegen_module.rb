$:.unshift File.join(File.dirname(__FILE__), ".")
require "code_templates.rb"

$GITCONFIG = ENV["GT_CODEGEN_CONFIG"] || "~/.gitconfig"

module CodeGen
  def CodeGen.perror(msg)
    STDERR.puts "#$0: error: #{msg}"
    STDERR.puts "Use '#$0 help' to print usage instructions"
    exit(1)
  end

  def CodeGen.find_user_data
    configfile = File.expand_path($GITCONFIG)
    if File.exist?(configfile)
      name = nil
      email = nil
      File.open(configfile) do |f|
        while (line = f.gets) && (!name||!email) do
          name = $1  if !name  && line =~ /\s*name\s*=\s*(.*)\s*\n/
          email = $1 if !email && line =~ /\s*email\s*=\s*(.*)\s*\n/
        end
      end
      return name, email
    else
      self.perror("#$GITCONFIG file not found")
      exit
    end
  end

  def CodeGen.create_license
    name, email = CodeGen.find_user_data
    year = Time.now.year
    ERB.new($license).result(binding)
  end

  def CodeGen.filename2classname(filename)
    "Gt" + File.basename(filename).gsub(/(^|_)(.)/){$2.upcase}
  end

  def CodeGen.extract_functions_from_rep(classN, repfile)
    get_parameters=false
    functions = Hash.new
    repfile.each do |line|
      if line.match /^typedef (\S+( \S+)* \*?)\(\*#{classN}(\w+)\)/
        functions[$3] = [$1]
        get_parameters = $3
      end
      if get_parameters
        line.scan /(\w+(\s\w+)*[*]?)(,|\);)/ do |parameter|
          functions[get_parameters].push parameter[0]
        end
        get_parameters = false if line.match /\);/
      end
      break if line.match /^struct #{classN}/
    end
    return functions
  end

end
