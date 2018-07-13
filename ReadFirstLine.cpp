#include <iostream>
#include "ReadFirstLine.h"

bool get_first_line(std::string *str, const char *filename) {
	std::ifstream ifpfile;
	boost::iostreams::filtering_stream<boost::iostreams::input> infile;

	if( !openfile(filename, &infile, &ifpfile)) {
		std::cout<< "Error opening file\n";
		return(false);
	};
	
	getline(infile,*str);

	return(true);
}
