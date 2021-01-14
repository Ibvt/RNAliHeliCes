#!/usr/bin/env ruby

require "pp"
require "open3"

# usage: split the first two lines into .fa and the following block as .hi
if ARGV.length != 4 then
  STDERR.puts "Usage: #{$0} hi.dis-file hi.dis.1-file rates.out kin.plt"  
  exit(-1)
end
hidisfilename = ARGV[0]
hidis1filename = ARGV[1]
ratesoutfilename = ARGV[2]
kinpltfilename = ARGV[3]


##############################################################
####  read filename from hi.dis  #### 
if (File.exist?(hidisfilename))
  begin 
    hidisfile = File.new(hidisfilename, "r")    
  rescue
    STDERR.print "Could not open file #{hidisfilename}!\n"
    exit 1
  end
else
  STDERR.print "File #{hidisfilename} does not exist!\n"
  exit 1
end

#####################################
####  output filename .hi.dis.1  #### 
begin 
  hidis1file = File.new(hidis1filename, "w")    
rescue
  STDERR.print "Could not open file #{hidis1filename}!\n"
  exit 1
end

#####################################
####  output filename rates.out  #### 
begin 
  ratesoutfile = File.new(ratesoutfilename, "w")    
rescue
  STDERR.print "Could not open file #{ratesoutfilename}!\n"
  exit 1
end

####################################
####  output filename .kin.plt  #### 
begin 
  kinpltfile = File.new(kinpltfilename, "w")    
rescue
  STDERR.print "Could not open file #{kinpltfilename}!\n"
  exit 1
end


line_no = 0
node_no = 0
hidisfile.each do |line|
    if line_no == 0 then

    elsif !line.start_with?("S") then
        hidis1file.puts(line) 
    else
        #hidis2file.puts(line)
	line_data = line.split(%r{\s+})
	2.upto(line_data.length-1) do |i|
	    ratesoutfile.printf("%10.4g", line_data[i])
        end
	ratesoutfile.printf("\n");
	node_no = node_no+1
    end
    line_no = line_no+1
end

rootname = kinpltfilename.split(".")[0]
kinpltfile.puts("set title 'Hishape based kinetic analysis'")
kinpltfile.puts("set xlabel 'arbitrary units'")
#kinpltfile.puts("set xrange [0.1:100000]")
kinpltfile.puts("set logscale x")
kinpltfile.puts("set ylabel 'Population density'")
kinpltfile.puts("set yrange [0:1]")
kinpltfile.puts("set term postscript enhanced font 'Time-roman,12'")
kinpltfile.puts("set term pdf")
kinpltfile.puts("set output '#{rootname}.kin.pdf'")
kinpltfile.puts("set key left top spacing 2.4") #title 'Legend' box 1
kinpltfile.puts("set key width 5")
kinpltfile.puts("set key height 5")

kinpltfile.puts("plot '#{rootname}.hi.dis.kin' using 1:2 title 'mfe' with lines linewidth 3, \\")
2.upto(node_no-1) do |i|
     kinpltfile.puts("'#{rootname}.hi.dis.kin' using 1:#{i+1} title '#{i}' with lines linewidth 3, \\")
end
kinpltfile.puts("'#{rootname}.hi.dis.kin' using 1:#{node_no+1} title '#{node_no}' with lines linewidth 3")

kinpltfile.puts("pause -1 'Hit any key to continue'")
#exit