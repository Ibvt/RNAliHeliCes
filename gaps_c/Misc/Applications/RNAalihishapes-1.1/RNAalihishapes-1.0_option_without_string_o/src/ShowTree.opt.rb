require 'optparse'
require 'optparse/time'
require 'ostruct'
require 'pp'

class RNAHeliCesBarrierOptparser

  CODES = %w[iso-2022-jp shift_jis euc-jp utf8 binary]
  CODE_ALIASES = { "jis" => "iso-2022-jp", "sjis" => "shift_jis" }

  #
  # Return a structure describing the options.
  #
  def self.parse(args)
    # The options specified on the command line will be collected in *options*.
    # We set default values here.
    options = OpenStruct.new
    options.inputfiles = []
    options.nvalue = 40
    options.kvalue = 50
    options.hishape_type = 2
    #options.version = false

    opts = OptionParser.new do |opts|
      opts.banner = "Usage: main.rb [options]"

      opts.separator ""
      opts.separator "Specific options:"

      
      # TODO1: how can a positional argument be added?
      # TODO2: why inputfiles only store the first element?
      # Mandatory argument.
      opts.on("-f", "--file FILE", Array,
              "Read input from specified file") do |input|
        options.inputfiles << input
      end

      opts.on("-n", "--nvalue [ANCHOR NUMBER]",
              "Specify n value for HiPath calculation") do |nvalue|
          options.nvalue = nvalue
      end

      opts.on("-k", "--kbest [NUMBER]",
              "Choose the k-best classes as leaves") do |kbest|
          options.kvalue = kbest
      end
      
      opts.on("-t", "--type [HISHAPE TYPE]",
              "Specify hishape abstraction type") do |hishape_type|
          options.hishape_type = hishape_type
      end      


      opts.separator ""
      opts.separator "Common options:"

      # No argument, shows at tail.  This will print an options summary.
      # Try it and see!
      opts.on_tail("-h", "--help", "Produce help message") do
        puts opts
        exit
      end

    end

    opts.parse!(args)
    options
  end  # parse()

end  # class OptparseExample

#options = OptparseExample.parse(ARGV)
#pp options
