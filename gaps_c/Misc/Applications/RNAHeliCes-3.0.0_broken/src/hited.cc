#include "hited_util.hh"



int main(int argc, char **argv)
{  
  
    // ########################
    // ## options            ##
    // ########################
        /* options declaration */
    //bool window_mode = false;
    //#ifdef WINDOW_MODE
        //    bool window_mode = true;
    //#endif


    unsigned int hishape_type;
    float ratio = 5.0f;  // TODO: describe the ratio is correspondent to \lambda in the paper in option description
    unsigned int convert = 0;

    try {
      
        // idea is the option is not compulsory and it will be given a relexed code  
        //po::options_description desc("Allowed options", 120, 60);
        po::options_description desc("Allowed options");
	
	desc.add_options()
	("file,f", po::value< std::vector<std::string> >(), "Read sequence(s) from specified file")
	("type,t", po::value<unsigned int>(&hishape_type)->default_value(2), "Specify hishape abstract type")
        ("ratio,r", po::value<float>(&ratio)->implicit_value(5.0f), "Specify weight between the score from helix index distance and score from helix index type conversion")
	("convert,c", po::value<unsigned int>(&convert)->implicit_value(1), "Convert secondary structure to hishape")
	("help,h", "Produce help message")
	("version,v", "Show version");


        /*
         * it can add at most additional 2 records, that means,
         * it is still probable that the sum of the records >= 2
         *
         * ./main --help -f aaa dddd ccc -s bbb
         * will get the result
         * input files are: aaa
         * input sequences are: dddd ccc bbb
         * positional options are: dddd ccc bbb
         */
        po::positional_options_description p;
        p.add("file", -1);


        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).positional(p).run(), vm);
        //po::store(parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        // option exclusions and implications
	conflicting_options(vm, "help", "version");
	conflicting_options(vm, "file", "sequence");
	conflicting_options(vm, "align", "convert");


        //std::cout << "vm.size()=" << vm.size() << std::endl;
        if (vm.count("help") || vm.size()==1) {
            std::cout << "Usage: HiTed [<first hishape>] [<second hishape>] [options]\n";
	    std::cout << "Examples:" << std::endl;
	    std::cout << "  ./HiTed '[27]' '[36.5m,(,27,38,)]' -r 1" << std::endl;
	    std::cout << "  ./HiTed ../examples/riboswitches.fas -t 1 -r 1" << std::endl;  // compare the 1st and 2nd, 3rd and 4nd
	    std::cout << "  ./HiTed ../examples/riboswitches.fas -c" << std::endl;            // convert every record 
            std::cout << desc; //<< '\n';
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "HiTed 2.0.12 (Feb. 09, 2012)\n";  // TODO: use the C++ buildin date function 
            return 0;
        }
        

       if (vm.count("file"))
       {
	    std::vector<std::string> files = vm["file"].as< std::vector<std::string> >();
            if (!files.empty()) {
	        // ##################
	        // ## convert mode ##
	        if (convert==1 && files.size()==1) 
		{
		  	if (hishape_type < 1 || hishape_type > 4) 
			{
				std::cout << "The hishape abstract type can only be specified between 1 and 4." << std::endl;  // For multiple hishape alignment, 
				return 1;
			}

			std::string line;
			std::fstream fiInputFile;
			if (false == fiInputFile.is_open())
			{       
					// get first input file
					fiInputFile.open(files[0].c_str(), std::ios::in);
					// checks if the file was opened
					if ((false == fiInputFile) || (false == fiInputFile.is_open()))
					{
						// not opened, throws an exception
						throw "ERROR: loading input file (\"" + files[0] + "\")!";
					}
					
					std::string header = "";
					std::string seq = "";
					std::string ss = "";
					int line_no = 0;
					while ( getline(fiInputFile, line) )
					{

						if ( (line_no%3 == 0) ) {
						    if (line[0] == '>')
						    {
							// get header only
							//header = fasta.substr(1);
							header = line;
						    }
						    else
						    {
							throw "ERROR: wrong faa-format in input file (\"" + files[0] + "\")!";
						    }
						}
						else 
						{
						    if (line_no%3 == 1)
							seq = line;
						    if (line_no%3 == 2)
							ss = line;
						}


						if (header != "" && seq != "" && ss != "") {
							std::cout << header << std::endl;
							std::cout << seq << std::endl;
							std::cout << ss << std::endl;
                                                        std::cout << sSToHishape(header+"\n"+seq+"\n"+ss+"\n", hishape_type) /*<< " (hishape type = " << hishape_type << ")" */<< std::endl;
							header = "";
							seq = "";
							ss = "";
						}
						line_no ++;
					}
			}
			else
			{
				std::cout << "ERROR: unable to open input file (\"" + files[0] + "\")!\n";
			}
			
			if (true == fiInputFile.is_open())
			{
				fiInputFile.close();
			}
		}
                else
	        // ##################
	        // #### align mode ##
		if (files.size()==2) 
		{   
		        std::string strInput1,strInput2;
		        if (!files[0].empty() && files[0].at(0) == '[' &&  !files[1].empty() && files[1].at(0) == '[')
			{
				strInput1 = files[0];
				strInput2 = files[1];
				pairAlign(strInput1, strInput2, hishape_type, ratio);
			}
		}
		else
		if (files.size()==1) 
		{
			std::string line;
			std::fstream fiInputFile;
			if (false == fiInputFile.is_open())
			{       
					// get first input file
					fiInputFile.open(files[0].c_str(), std::ios::in);
					// checks if the file was opened
					if ((false == fiInputFile) || (false == fiInputFile.is_open()))
					{
						// not opened, throws an exception
						throw "ERROR: loading input file (\"" + files[0] + "\")!";
					}

					std::string header1 = "";
					std::string seq1 = "";
					std::string ss1 = "";
					std::string header2 = "";
					std::string seq2 = "";
					std::string ss2 = "";
					int line_no = 0;
					while ( getline(fiInputFile, line) )
					{

						if ( (line_no%6 == 0) ) {
						    if (line[0] == '>')
						    {
							// get header only
							//header1 = line.substr(1);
							header1 = line;
						    }
						    else
						    {
							throw "ERROR: wrong faa-format in input file (\"" + files[0] + "\")!";
						    }
						}
						else if ( (line_no%6 == 3) ) {
						    if (line[0] == '>')
						    {
							// get header only
							//header2 = line.substr(1);
							header2 = line;
						    }
						    else
						    {
							throw "ERROR: wrong faa-format in input file (\"" + files[0] + "\")!";
						    }
						}
						else 
						{
						    if (line_no%6 == 1)
							seq1 = line;
						    if (line_no%6 == 2)
							ss1 = line;
						    if (line_no%6 == 4)
							seq2 = line;
						    if (line_no%6 == 5)
							ss2 = line;
						}


						if (header1 != "" && header2 != "" && seq1 != "" && seq2 != "" && ss1 != "" && ss2 != "") {
						        std::cout << "#### " << header1 + " vs. " + header2 + " ####" << std::endl;  
                                                        pairAlign(sSToHishape(header1+"\n"+seq1+"\n"+ss1+"\n", hishape_type), sSToHishape(header2+"\n"+seq2+"\n"+ss2+"\n", hishape_type), hishape_type, ratio);
							header1 = "";
							header2 = "";
							seq1 = "";
							seq2 = "";
							ss1 = "";
							ss2 = "";
						}
						line_no ++;
					}
			}
			else
			{
				std::cout << "ERROR: unable to open input file (\"" + files[0] + "\")!\n";
			}
			
			if (true == fiInputFile.is_open())
			{
				fiInputFile.close();
			}
		}
		else
		{
		    std::cout << "Status: files.size()=" << files.size() << std::endl;
		}
	    }
        }
    } catch (std::exception &e) {
        std::cerr << "Exception: " << e.what() << '\n';
        std::exit(1);
    }


  return 0;
}