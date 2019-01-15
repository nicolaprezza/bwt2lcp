// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * lcp.hpp
 *
 *  Created on: Dec 5, 2018
 *      Author: nico
 *
 * Induces the LCP array of a collection of reads from its compact BWT.
 *
 * Assumption: terminators are all represented with special character TERM=#
 *
 * RAM n*(B+0.5) Bytes, where B = sizeof(lcp_int_t) is the number of Bytes needed to represent an LCP value.
 *
 * Based on an extension (to BWTs of collections) of the suffix-tree navigation algorithm described in
 *
 * "Linear time construction of compressed text indices in compact space" by Djamal Belazzougui.
 *
 */

#ifndef INTERNAL_LCP_HPP_
#define INTERNAL_LCP_HPP_

#include "dna_bwt.hpp"
#include "include.hpp"
#include <stack>
#include <algorithm>

using namespace std;

template<class bwt_t, typename lcp_int_t>
class lcp{

public:

	/*
	 * Build LCP from BWT. Parameters:
	 *
	 */
	lcp(bwt_t * bwt){

		this->bwt;

		n = bwt->size();

		LCP = vector<lcp_int_t>(n, nil);
		LCP[0] = 0;

		/*
		 * FIRST PASS: LEAVES NAVIGATION. COMPUTE LCP VALUES INSIDE SUFFIX TREE LEAVES.
		 */

		cout << "\nNow navigating suffix tree leaves of size >= 2 to compute internal LCP values." << endl;

		uint64_t m = 0;//portion of text covered by visited leaves
		uint64_t leaves = 0;//number of visited leaves
		uint64_t max_stack = 0;
		uint64_t lcp_values = 1;//number of filled LCP values

		{

			auto TMP_LEAVES = vector<sa_leaf>(5);

			stack<sa_leaf> S;

			sa_leaf root = bwt->first_leaf();
			S.push(root);

			int last_perc = -1;
			int perc = 0;

			while(not S.empty()){

				sa_leaf L = S.top();
				S.pop();
				leaves++;

				max_stack = S.size() > max_stack ? S.size() : max_stack;

				assert(L.rn.second > L.rn.first);

				for(uint64_t i = L.rn.first+1; i<L.rn.second; ++i){

					assert(LCP[i]==nil);

					LCP[i] = L.depth;

					lcp_values++;
					m++;

				}

				m++;

				assert(m<=n);

				int t = 0;

				bwt->next_leaves(L, TMP_LEAVES, t, 2);

				for(int i=t-1;i>=0;--i) S.push(TMP_LEAVES[i]);

				perc = (100*lcp_values)/n;

				if(perc > last_perc){

					cout << "LCP: " << perc << "%." << endl;
					last_perc = perc;

				}

			}
		}

		cout << "Visited leaves cover " << m << "/" << n << " input characters." << endl;
		cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;

		cout << "Max stack size = " << max_stack << endl;
		cout << "Processed " << leaves << " suffix-tree leaves of size >= 2." << endl;

		cout << "\nNow navigating suffix tree nodes to compute remaining LCP values." << endl;

		{

			auto TMP_NODES = vector<typename bwt_t::sa_node_t>(5);

			uint64_t nodes = 0;//visited ST nodes
			max_stack = 0;

			stack<typename bwt_t::sa_node_t> S;

			typename bwt_t::sa_node_t root = bwt->root();

			S.push(root);

			int last_perc = -1;
			int perc = 0;

			while(not S.empty()){

				max_stack = S.size() > max_stack ? S.size() : max_stack;

				typename bwt_t::sa_node_t N = S.top();
				S.pop();
				nodes++;

				update_lcp<lcp_int_t>(N,LCP,lcp_values);

				int t = 0;

				bwt->next_nodes(N, TMP_NODES, t);

				for(int i=t-1;i>=0;--i){

					S.push(TMP_NODES[i]);

				}

				perc = (100*lcp_values)/n;

				if(perc > last_perc){

					cout << "LCP: " << perc << "%." << endl;
					last_perc = perc;

				}

			}

			cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;
			cout << "Max stack size = " << max_stack << endl;
			cout << "Processed " << nodes << " suffix-tree nodes." << endl;

		}


	}

	/*
	 * store LCP to file
	 */
	void save_to_file(string lcp_path){

		std::ofstream out(lcp_path);

		out.write((char*)LCP.data(),LCP.size() * sizeof(lcp_int_t));

		out.close();

	}

	vector<lcp_int_t> LCP;

private:

	uint64_t n = 0;//total size

	bwt_t * bwt = NULL;

	lcp_int_t nil = ~lcp_int_t(0);

};



#endif /* INTERNAL_LCP_HPP_ */
