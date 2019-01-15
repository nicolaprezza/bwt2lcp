// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * induce_lcp.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: nico
 */


#include <iostream>
#include "internal/lcp.hpp"
#include <unistd.h>
#include "internal/dna_bwt_n.hpp"

using namespace std;

string input_bwt;
string output_file;

bool out_da = false;
uint8_t lcp_size = 1;

bool containsN = false;

char TERM = '#';

void help(){

	cout << "bwt2lcp [options]" << endl <<
	"Input: BWT of a collection of reads. Output: LCP array of the collection." << endl <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-i <arg>    Input BWT (REQUIRED)" << endl <<
	"-o <arg>    Output file name (REQUIRED)" << endl <<
	"-l <arg>    Number of Bytes used to represent LCP values. <arg>=1,2,4,8 Bytes. Default: 1." << endl <<
	//"-n          Alphabet is {A,C,G,N,T," << TERM << "}. Default: alphabet is {A,C,G,T," << TERM << "}."<< endl <<
	"-t          ASCII code of the terminator. Default:" << int('#') << " (#). Cannot be the code for A,C,G,T,N." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 3) help();

	int opt;
	while ((opt = getopt(argc, argv, "hi:o:l:t:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'i':
				input_bwt = string(optarg);
			break;
			case 'o':
				output_file = string(optarg);
			break;
			case 'l':
				lcp_size = atoi(optarg);
			break;
			case 't':
				TERM = atoi(optarg);
			break;
			/*case 'n':
				containsN=true;
			break;*/
			default:
				help();
			return -1;
		}
	}

	if(lcp_size != 1 and lcp_size != 2 and lcp_size != 4 and lcp_size != 8) help();

	if(TERM == 'A' or TERM == 'C' or TERM == 'G' or TERM == 'T' or TERM == 'N'){

		cout << "Error: invalid terminator '" << TERM << "'" << endl;
		help();

	}

	if(input_bwt.size()==0) help();
	if(output_file.size()==0) help();

	cout << "Input bwt file: " << input_bwt << endl;
	cout << "Output LCP file: " << output_file << endl;

	containsN = hasN(input_bwt);

	if(not containsN){

		cout << "Alphabet: A,C,G,T,'" << TERM << "'" << endl;

		cout << "Loading and indexing BWT ... " << endl;

		dna_bwt_t BWT = dna_bwt_t(input_bwt, TERM);

		cout << "Done. Size of BWT: " << BWT.size()<< endl;

		switch(lcp_size){

			case 1: { lcp<dna_bwt_t, uint8_t> M1(&BWT);
					cout << "Storing output to file ... " << endl;
					M1.save_to_file(output_file);
					break; }
			case 2: { lcp<dna_bwt_t, uint16_t> M2(&BWT);
					cout << "Storing output to file ... " << endl;
					M2.save_to_file(output_file);
					break; }
			case 4: { lcp<dna_bwt_t, uint32_t> M4(&BWT);
					cout << "Storing output to file ... " << endl;
					M4.save_to_file(output_file);
					break; }
			case 8: { lcp<dna_bwt_t, uint64_t> M8(&BWT);
					cout << "Storing output to file ... " << endl;
					M8.save_to_file(output_file);
					break; }
			default:break;

		}

	}else{

		cout << "Alphabet: A,C,G,N,T,'" << TERM << "'" << endl;

		cout << "Loading and indexing BWT ... " << endl;

		dna_bwt_n_t BWT = dna_bwt_n_t(input_bwt, TERM);

		cout << "Done. Size of BWT: " << BWT.size()<< endl;

		switch(lcp_size){

			case 1: { lcp<dna_bwt_n_t, uint8_t> M1(&BWT);
					cout << "Storing output to file ... " << endl;
					M1.save_to_file(output_file);
					break; }
			case 2: { lcp<dna_bwt_n_t, uint16_t> M2(&BWT);
					cout << "Storing output to file ... " << endl;
					M2.save_to_file(output_file);
					break; }
			case 4: { lcp<dna_bwt_n_t, uint32_t> M4(&BWT);
					cout << "Storing output to file ... " << endl;
					M4.save_to_file(output_file);
					break; }
			case 8: { lcp<dna_bwt_n_t, uint64_t> M8(&BWT);
					cout << "Storing output to file ... " << endl;
					M8.save_to_file(output_file);
					break; }
			default:break;

		}

	}

	cout << "Done. " << endl;

}

