#!/usr/bin/env ruby

require "pp"
require "open3"

# usage: split the first two lines into .fa and the following block as .hi
if ARGV.length != 3 then
  STDERR.puts "Usage: #{$0} dis-file fa-file hi-file"
  exit(-1)
end
faafilename = ARGV[0]
fafilename = ARGV[1]
hifilename = ARGV[2]


##############################################################
####  read filename from faafile  #### 
if (File.exist?(faafilename))
  begin 
    faafile = File.new(faafilename, "r")    
  rescue
    STDERR.print "Could not open file #{faafilename}!\n"
    exit 1
  end
else
  STDERR.print "File #{faafilename} does not exist!\n"
  exit 1
end

###############################
####  output filename .fa  #### 
begin 
  fafile = File.new(fafilename, "w")    
rescue
  STDERR.print "Could not open file #{fafilename}!\n"
  exit 1
end


###############################
####  output filename .hi  #### 
begin 
  hifile = File.new(hifilename, "w")    
rescue
  STDERR.print "Could not open file #{hifilename}!\n"
  exit 1
end


line_no = 0
faafile.each do |line|
    if line_no == 0 then
        fafile.puts(">#{line}")
    elsif line_no == 1 then
        fafile.puts(line)
    elsif !line.start_with?("S") then
        hifile.puts(line)
    end
    line_no = line_no+1
end