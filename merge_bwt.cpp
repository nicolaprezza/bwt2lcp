// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * merge_bwt.cpp
 *
 *  Created on: Nov 26, 2018
 *      Author: nico
 */

#include <iostream>
#include "internal/bwt_merger.hpp"
#include <unistd.h>
#include "internal/dna_bwt.hpp"
#include "internal/dna_bwt_n.hpp"

using namespace std;

string input_bwt1;
string input_bwt2;
string output_file;

bool out_da = false;
uint8_t lcp_size = 0;

bool containsN = false;

char TERM = '#';

void help(){

	cout << "merge_bwt [options]" << endl <<
	"Merges the eBWTs of two collections of reads by navigating the (compressed) generalized suffix tree of their union." << endl <<
	"Options:" << endl <<
	"-h          Print this help" << endl <<
	"-1 <arg>    Input BWT 1 (REQUIRED)" << endl <<
	"-2 <arg>    Input BWT 2 (REQUIRED)" << endl <<
	"-o <arg>    Output prefix (REQUIRED)" << endl <<
	"-d          Output document array as an ASCII file of 0/1. Default: do not output." << endl <<
	"-l <arg>    Output LCP of the merged BWT using <arg>=0,1,2,4,8 Bytes" << endl <<
	"            per integer. If arg=0, LCP is not computed (faster). Default: 0." << endl <<
	//"-n          Alphabet is {A,C,G,N,T," << TERM << "}. Default: alphabet is {A,C,G,T," << TERM << "}."<< endl <<
	"-t          Ascii code of the terminator. Default:" << int('#') << " (#). Cannot be the code for A,C,G,T,N." << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 4) help();

	int opt;
	while ((opt = getopt(argc, argv, "h1:2:o:l:dt:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case '1':
				input_bwt1 = string(optarg);
			break;
			case '2':
				input_bwt2 = string(optarg);
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
			case 'd':
				out_da = true;
			break;
			/*case 'n':
				containsN = true;
			break;*/
			default:
				help();
			return -1;
		}
	}

	if(TERM == 'A' or TERM == 'C' or TERM == 'G' or TERM == 'T' or TERM == 'N'){

		cout << "Error: invalid terminator '" << TERM << "'" << endl;
		help();

	}

	if(input_bwt1.size()==0) help();
	if(input_bwt2.size()==0) help();
	if(output_file.size()==0) help();

	cout << "Input bwt 1: " << input_bwt1 << endl;
	cout << "Input bwt 2: " << input_bwt2 << endl;
	cout << "Output prefix: " << output_file << endl;

	containsN = hasN(input_bwt1) or hasN(input_bwt2);

	if(not containsN){

		cout << "Alphabet: A,C,G,T,'" << TERM << "'" << endl;

		cout << "Loading and indexing BWTs ... " << endl;

		dna_bwt_t BWT1 = dna_bwt_t(input_bwt1, TERM);
		dna_bwt_t BWT2 = dna_bwt_t(input_bwt2, TERM);

		cout << "Done. Size of BWTs: " << BWT1.size() << " and " << BWT2.size() << endl;

		switch(lcp_size){

			case 0: { bwt_merger<dna_bwt_t, uint8_t> M0(&BWT1, &BWT2, false, out_da);
					cout << "Storing output to file ... " << endl;
					M0.save_to_file(output_file);
					break; }
			case 1: { bwt_merger<dna_bwt_t, uint8_t> M1(&BWT1, &BWT2, true, out_da);
					cout << "Storing output to file ... " << endl;
					M1.save_to_file(output_file);
					break;}
			case 2: {bwt_merger<dna_bwt_t, uint16_t> M2(&BWT1, &BWT2, true, out_da);
					cout << "Storing output to file ... " << endl;
					M2.save_to_file(output_file);
					break;}
			case 4: {bwt_merger<dna_bwt_t, uint32_t> M4(&BWT1, &BWT2, true, out_da);
					cout << "Storing output to file ... " << endl;
					M4.save_to_file(output_file);
					break;}
			case 8: {bwt_merger<dna_bwt_t, uint64_t> M8(&BWT1, &BWT2, true, out_da);
					cout << "Storing output to file ... " << endl;
					M8.save_to_file(output_file);
					break;}
			default:break;

		}

	}else{

		cout << "Alphabet: A,C,G,N,T,'" << TERM << "'" << endl;

		cout << "Loading and indexing BWTs ... " << endl;

		dna_bwt_n_t BWT1 = dna_bwt_n_t(input_bwt1, TERM);
		dna_bwt_n_t BWT2 = dna_bwt_n_t(input_bwt2, TERM);

		cout << "Done. Size of BWTs: " << BWT1.size() << " and " << BWT2.size() << endl;

		switch(lcp_size){

			case 0: { bwt_merger<dna_bwt_n_t, uint8_t> M0(&BWT1, &BWT2, false, out_da);
					cout << "Storing output to file ... " << endl;
					M0.save_to_file(output_file);
					break; }
			case 1: { bwt_merger<dna_bwt_n_t, uint8_t> M1(&BWT1, &BWT2, true, out_da);
					cout << "Storing output to file ... " << endl;
					M1.save_to_file(output_file);
					break;}
			case 2: {bwt_merger<dna_bwt_n_t, uint16_t> M2(&BWT1, &BWT2, true, out_da);
					cout << "Storing output to file ... " << endl;
					M2.save_to_file(output_file);
					break;}
			case 4: {bwt_merger<dna_bwt_n_t, uint32_t> M4(&BWT1, &BWT2, true, out_da);
					cout << "Storing output to file ... " << endl;
					M4.save_to_file(output_file);
					break;}
			case 8: {bwt_merger<dna_bwt_n_t, uint64_t> M8(&BWT1, &BWT2, true, out_da);
					cout << "Storing output to file ... " << endl;
					M8.save_to_file(output_file);
					break;}
			default:break;

		}

	}

	cout << "Done. " << endl;

}

