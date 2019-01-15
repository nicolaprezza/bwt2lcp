// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * bwt_merger.hpp
 *
 *  Created on: Dec 5, 2018
 *      Author: nico
 *
 * Merges the compact BWTs of length n1 and n2 of two string collections. Always computes the document array DA (n1+n2 bits).
 *
 * Assumption: terminators are all represented with special character TERM=#
 *
 * RAM:
 *
 *  - If LCP is not computed: n*0.625 Bytes, where n = combined length of the two input BWTs.
 *  - If LCP is computed: n*(B+0.625) Bytes, where B = sizeof(lcp_int_t) is the number of Bytes needed
 *    to represent an LCP value.
 *
 *  If alphabet is {A,C,G,T,N,#} (dna_bwt_n used instead of dna_bwt) then space usage is increased by 0.05*n Bytes.
 *
 * Based on an extension (to BWTs of collections) of the suffix-tree navigation algorithm described in
 *
 * "Linear time construction of compressed text indices in compact space" by Djamal Belazzougui.
 *
 */

#ifndef INTERNAL_BWT_MERGER_HPP_
#define INTERNAL_BWT_MERGER_HPP_

#include "dna_bwt.hpp"
#include "dna_bwt_n.hpp"
#include "include.hpp"
#include <stack>
#include <algorithm>

using namespace std;

template<class bwt_t, typename lcp_int_t>
class bwt_merger{

public:

	/*
	 * merge bwt1 and bwt2. Parameters:
	 *
	 * bwt1, bwt2: the two BWTs.
	 * compute_lcp = if true, compute the LCP array of the merged BWT.
	 * out_da = store to file document array.
	 *
	 */
	bwt_merger(bwt_t * bwt1, bwt_t * bwt2, bool compute_lcp = false, bool out_da = false){

		this->out_da = out_da;
		this->bwt1 = bwt1;
		this->bwt2 = bwt2;

		n = bwt1->size() + bwt2->size();

		DA = vector<bool>(n);

		if(compute_lcp){

			LCP = vector<lcp_int_t>(n, nil);
			LCP[0] = 0;
		}

		/*
		 * FIRST PASS: NAVIGATE LEAVES AND MERGE BWTs (i.e. build DA). If enabled, compute LCP inside leaves.
		 */

		cout << "\nNow navigating suffix tree leaves to compute Document Array";
		if(compute_lcp) cout << " and compute internal LCP values";
		cout << "." << endl;

		uint64_t da_values = 0;//number computed DA values
		uint64_t leaves = 0;//number of visited leaves
		uint64_t max_stack = 0;
		uint64_t lcp_values = 1;//number of computed LCP values

		{

			auto TMP_LEAVES = vector<pair<sa_leaf, sa_leaf> >(5);

			stack<pair<sa_leaf, sa_leaf> > S;
			S.push({bwt1->first_leaf(), bwt2->first_leaf()});

			int last_perc_lcp = -1;
			int last_perc_da = -1;
			int perc_lcp = 0;
			int perc_da = 0;

			while(not S.empty()){

				pair<sa_leaf, sa_leaf> L = S.top();
				S.pop();
				leaves++;

				assert(leaf_size(L)>0);
				max_stack = S.size() > max_stack ? S.size() : max_stack;

				sa_leaf L1 = L.first;
				sa_leaf L2 = L.second;

				update_DA(L1,L2,compute_lcp,lcp_values,da_values);

				//insert leaf in stack iff size(L1) + size(L2) >= min_size
				//optimization: if we are computing LCP and if size(L1) + size(L2) = 1,
				//then we will find that leaf during the internal nodes traversal (no need to visit leaf here)
				int t = 0;//number of children leaves
				next_leaves(bwt1, bwt2, L1, L2, TMP_LEAVES, t, compute_lcp ? 2 : 1);
				for(int i=t-1;i>=0;--i) S.push(TMP_LEAVES[i]);

				perc_lcp = (100*lcp_values)/n;
				perc_da = (100*da_values)/n;

				if(perc_da > last_perc_da){

					cout << "DA: " << perc_da << "%. ";
					if(compute_lcp) cout << "LCP: " << perc_lcp << "%.";
					cout << endl;

					last_perc_lcp = perc_lcp;
					last_perc_da = perc_da;

				}

			}
		}

		cout << "Computed " << da_values << "/" << n << " DA values." << endl;

		if(compute_lcp)
		cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;

		cout << "Max stack depth = " << max_stack << endl;
		cout << "Processed " << leaves << " suffix-tree leaves." << endl;

		if(compute_lcp){

			cout << "\nNow navigating suffix tree nodes to compute remaining LCP and DA values." << endl;

			auto TMP_NODES = vector<pair<typename bwt_t::sa_node_t, typename bwt_t::sa_node_t> >(5);

			uint64_t nodes = 0;//visited ST nodes
			max_stack = 0;

			stack<pair<typename bwt_t::sa_node_t, typename bwt_t::sa_node_t> > S;
			S.push({bwt1->root(), bwt2->root()});

			int last_perc_lcp = -1;
			int perc_lcp = 0;
			int last_perc_da = -1;
			int perc_da = 0;

			while(not S.empty()){

				max_stack = S.size() > max_stack ? S.size() : max_stack;

				pair<typename bwt_t::sa_node_t, typename bwt_t::sa_node_t> N = S.top();
				S.pop();
				nodes++;

				typename bwt_t::sa_node_t N1 = N.first;
				typename bwt_t::sa_node_t N2 = N.second;
				typename bwt_t::sa_node_t merged = merge_nodes(N1, N2);

				//find leaves in the children of N1 and N2 that were
				//skipped in the first pass, and update DA accordingly
				find_leaves(N1, N2, da_values);

				//compute LCP values at the borders of merged's children
				update_lcp<lcp_int_t>(merged, LCP, lcp_values);

				//follow Weiner links
				int t = 0;
				next_nodes(bwt1, bwt2, N1, N2, TMP_NODES, t);
				for(int i=t-1;i>=0;--i) S.push(TMP_NODES[i]);

				perc_lcp = (100*lcp_values)/n;
				perc_da = (100*da_values)/n;

				if(perc_da > last_perc_da){

					cout << "DA: " << perc_da << "%. ";
					cout << "LCP: " << perc_lcp << "%.";
					cout << endl;

					last_perc_lcp = perc_lcp;
					last_perc_da = perc_da;

				}

			}

			cout << "Computed " << da_values << "/" << n << " DA values." << endl;
			cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;
			cout << "Max stack depth = " << max_stack << endl;
			cout << "Processed " << nodes << " suffix-tree nodes." << endl;

		}

	}

	/*
	 * store to file the components that have been built (BWT, DA for sure; optionally, LCP). Adds extensions .bwt, .da, .lcp
	 */
	void save_to_file(string base_path){

		string bwt_path = base_path;
		bwt_path.append(".bwt");

		string da_path = base_path;
		da_path.append(".da");

		string lcp_path = base_path;
		lcp_path.append(".lcp");

		{
			std::ofstream out(bwt_path);

			uint64_t rank1 = 0;
			for(uint64_t i=0;i<n;++i){

				char c = DA[i]==0 ? bwt1->operator[](i-rank1) : bwt2->operator[](rank1);

				out.write((char*)&c,sizeof(c));

				rank1 += DA[i];

			}
			out.close();
		}

		if(out_da){
			std::ofstream out(da_path);
			for(uint64_t i=0;i<n;++i){

				uint8_t c = DA[i] ? '1' : '0';
				out.write((char*)&c,sizeof(c));

			}
			out.close();
		}

		if(LCP.size()>0){

			std::ofstream out(lcp_path);

			out.write((char*)LCP.data(),LCP.size() * sizeof(lcp_int_t));

			out.close();

		}

	}

private:

	vector<bool> DA; //document array
	vector<lcp_int_t> LCP;

	void update_DA(sa_leaf L1,sa_leaf L2, bool compute_lcp, uint64_t & lcp_values, uint64_t & m){

		uint64_t start1 = L1.rn.first + L2.rn.first;//start position of first interval in merged intervals
		uint64_t start2 = L2.rn.first + L1.rn.second;//start position of second interval in merged intervals
		uint64_t end = L1.rn.second + L2.rn.second;//end position (excluded) of merged intervals

		assert(end>start1);

		for(uint64_t i = start1; i<start2; ++i){
			DA[i] = 0;
			m++;
		}

		for(uint64_t i = start2; i<end; ++i){
			DA[i] = 1;
			m++;
		}

		assert(L1.depth==L2.depth);

		if(compute_lcp){

			for(uint64_t i = start1+1; i<end; ++i){

				assert(LCP[i]==nil);

				LCP[i] = L1.depth;

				lcp_values++;

			}

		}

		assert(m<=n);

	}

	void find_leaves(sa_node x1, sa_node x2, uint64_t & m){

		//find leaves that were ignored in the first pass
		if(range_length(child_TERM(x1))+range_length(child_TERM(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_TERM(x1),0};
			sa_leaf L2 = {child_TERM(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

		if(range_length(child_A(x1))+range_length(child_A(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_A(x1),0};
			sa_leaf L2 = {child_A(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

		if(range_length(child_C(x1))+range_length(child_C(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_C(x1),0};
			sa_leaf L2 = {child_C(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

		if(range_length(child_G(x1))+range_length(child_G(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_G(x1),0};
			sa_leaf L2 = {child_G(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

		if(range_length(child_T(x1))+range_length(child_T(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_T(x1),0};
			sa_leaf L2 = {child_T(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

	}


	void find_leaves(sa_node_n x1, sa_node_n x2, uint64_t & m){

		//find leaves that were ignored in the first pass
		if(range_length(child_TERM(x1))+range_length(child_TERM(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_TERM(x1),0};
			sa_leaf L2 = {child_TERM(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

		if(range_length(child_A(x1))+range_length(child_A(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_A(x1),0};
			sa_leaf L2 = {child_A(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

		if(range_length(child_C(x1))+range_length(child_C(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_C(x1),0};
			sa_leaf L2 = {child_C(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

		if(range_length(child_G(x1))+range_length(child_G(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_G(x1),0};
			sa_leaf L2 = {child_G(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

		if(range_length(child_N(x1))+range_length(child_N(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_N(x1),0};
			sa_leaf L2 = {child_N(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

		if(range_length(child_T(x1))+range_length(child_T(x2)) == 1){

			//symbolic depth = 0. It will not be used in update_DA
			sa_leaf L1 = {child_T(x1),0};
			sa_leaf L2 = {child_T(x2),0};

			uint64_t dummy = 0;
			update_DA(L1,L2,false,dummy,m);

		}

	}

	uint64_t n = 0;//total size

	bool out_da = false;

	bwt_t * bwt1 = NULL;
	bwt_t * bwt2 = NULL;

	lcp_int_t nil = ~lcp_int_t(0);

};


void next_leaves(dna_bwt_t * bwt1, dna_bwt_t * bwt2, sa_leaf & L1, sa_leaf & L2, vector<pair<sa_leaf, sa_leaf> > & TMP_LEAVES, int & t, int min_size){

	p_range ext_1 = bwt1->LF(L1.rn);
	p_range ext_2 = bwt2->LF(L2.rn);

	//push non-empty leaves on stack in decreasing size order

	t = 0;

	if(range_length(ext_1.A) + range_length(ext_2.A) >= min_size) TMP_LEAVES[t++] = {{ext_1.A, L1.depth+1},{ext_2.A, L2.depth+1}};
	if(range_length(ext_1.C) + range_length(ext_2.C) >= min_size) TMP_LEAVES[t++] = {{ext_1.C, L1.depth+1},{ext_2.C, L2.depth+1}};
	if(range_length(ext_1.G) + range_length(ext_2.G) >= min_size) TMP_LEAVES[t++] = {{ext_1.G, L1.depth+1},{ext_2.G, L2.depth+1}};
	if(range_length(ext_1.T) + range_length(ext_2.T) >= min_size) TMP_LEAVES[t++] = {{ext_1.T, L1.depth+1},{ext_2.T, L2.depth+1}};

	std::sort( TMP_LEAVES.begin(), TMP_LEAVES.begin()+t, [ ]( const pair<sa_leaf, sa_leaf>& lhs, const pair<sa_leaf, sa_leaf>& rhs )
	{
		return leaf_size(lhs) < leaf_size(rhs);
	});

}

void next_leaves(dna_bwt_n_t * bwt1, dna_bwt_n_t * bwt2, sa_leaf & L1, sa_leaf & L2, vector<pair<sa_leaf, sa_leaf> > & TMP_LEAVES, int & t, int min_size){

	p_range_n ext_1 = bwt1->LF(L1.rn);
	p_range_n ext_2 = bwt2->LF(L2.rn);

	//push non-empty leaves on stack in decreasing size order

	t = 0;

	if(range_length(ext_1.A) + range_length(ext_2.A) >= min_size) TMP_LEAVES[t++] = {{ext_1.A, L1.depth+1},{ext_2.A, L2.depth+1}};
	if(range_length(ext_1.C) + range_length(ext_2.C) >= min_size) TMP_LEAVES[t++] = {{ext_1.C, L1.depth+1},{ext_2.C, L2.depth+1}};
	if(range_length(ext_1.G) + range_length(ext_2.G) >= min_size) TMP_LEAVES[t++] = {{ext_1.G, L1.depth+1},{ext_2.G, L2.depth+1}};
	if(range_length(ext_1.N) + range_length(ext_2.N) >= min_size) TMP_LEAVES[t++] = {{ext_1.N, L1.depth+1},{ext_2.N, L2.depth+1}};
	if(range_length(ext_1.T) + range_length(ext_2.T) >= min_size) TMP_LEAVES[t++] = {{ext_1.T, L1.depth+1},{ext_2.T, L2.depth+1}};

	std::sort( TMP_LEAVES.begin(), TMP_LEAVES.begin()+t, [ ]( const pair<sa_leaf, sa_leaf>& lhs, const pair<sa_leaf, sa_leaf>& rhs )
	{
		return leaf_size(lhs) < leaf_size(rhs);
	});

}


void next_nodes(dna_bwt_t * bwt1, dna_bwt_t * bwt2, sa_node & N1, sa_node & N2, vector<pair<sa_node, sa_node> > & TMP_NODES, int & t){

	p_node left_exts1 = bwt1->LF(N1);
	p_node left_exts2 = bwt2->LF(N2);

	pair<sa_node, sa_node> A = {left_exts1.A, left_exts2.A};
	pair<sa_node, sa_node> C = {left_exts1.C, left_exts2.C};
	pair<sa_node, sa_node> G = {left_exts1.G, left_exts2.G};
	pair<sa_node, sa_node> T = {left_exts1.T, left_exts2.T};

	t = 0;

	if(number_of_children(A) >= 2) TMP_NODES[t++] = A;
	if(number_of_children(C) >= 2) TMP_NODES[t++] = C;
	if(number_of_children(G) >= 2) TMP_NODES[t++] = G;
	if(number_of_children(T) >= 2) TMP_NODES[t++] = T;

	//push right-maximal nodes on stack in decreasing size (i.e. interval length) order

	std::sort( TMP_NODES.begin(), TMP_NODES.begin()+t, [ ]( const pair<sa_node, sa_node>& lhs, const pair<sa_node, sa_node>& rhs )
	{
		return node_size(lhs) < node_size(rhs);
	});

}

void next_nodes(dna_bwt_n_t * bwt1, dna_bwt_n_t * bwt2, sa_node_n & N1, sa_node_n & N2, vector<pair<sa_node_n, sa_node_n> > & TMP_NODES, int & t){

	p_node_n left_exts1 = bwt1->LF(N1);
	p_node_n left_exts2 = bwt2->LF(N2);

	pair<sa_node_n, sa_node_n> A = {left_exts1.A, left_exts2.A};
	pair<sa_node_n, sa_node_n> C = {left_exts1.C, left_exts2.C};
	pair<sa_node_n, sa_node_n> G = {left_exts1.G, left_exts2.G};
	pair<sa_node_n, sa_node_n> N = {left_exts1.N, left_exts2.N};
	pair<sa_node_n, sa_node_n> T = {left_exts1.T, left_exts2.T};

	t = 0;

	if(number_of_children(A) >= 2) TMP_NODES[t++] = A;
	if(number_of_children(C) >= 2) TMP_NODES[t++] = C;
	if(number_of_children(G) >= 2) TMP_NODES[t++] = G;
	if(number_of_children(N) >= 2) TMP_NODES[t++] = N;
	if(number_of_children(T) >= 2) TMP_NODES[t++] = T;

	//push right-maximal nodes on stack in decreasing size (i.e. interval length) order

	std::sort( TMP_NODES.begin(), TMP_NODES.begin()+t, [ ]( const pair<sa_node_n, sa_node_n>& lhs, const pair<sa_node_n, sa_node_n>& rhs )
	{
		return node_size(lhs) < node_size(rhs);
	});

}

#endif /* INTERNAL_BWT_MERGER_HPP_ */
