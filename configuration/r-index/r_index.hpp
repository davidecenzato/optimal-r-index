// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.
//
// Modified by Davide Cenzato, 2026
/*
 * r_index.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: nico
 *
 * Small version of the r-index: O(r) words of space, O(log(n/r)) locate time per occurrence
 *
 */

#ifndef R_INDEX_S_H_
#define R_INDEX_S_H_

#include <definitions.hpp>
#include <rle_string.hpp>
#include "sparse_sd_vector.hpp"
#include "sparse_hyb_vector.hpp"
#include "utils.hpp"

template<typename T>
void read_file(const char *filename, std::vector<T>& ptr){
    struct stat filestat;
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        std::cerr << "open() file " + std::string(filename) + " failed\n";

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        std::cerr << "stat() file " + std::string(filename) + " failed\n";

    if(filestat.st_size % sizeof(T) != 0)
        std::cerr << "invilid file " + std::string(filename) + "\n";

    size_t length = filestat.st_size / sizeof(T);
    ptr.resize(length);

    if ((fread(&ptr[0], sizeof(T), length, fd)) != length) 
        std::cerr << "fread() file " + std::string(filename) + " failed\n";

    fclose(fd);
}

template <typename T>
vector<ulint> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<ulint> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

using namespace sdsl;

namespace ri{

template	<	class sparse_bv_type = sparse_sd_vector,
				class rle_string_t = rle_string_sd
			>
class r_index{

public:

	using triple = std::tuple<range_t, ulint, ulint>;

	r_index(){}

	r_index(string &input, bool concat = true)
	{
		if(concat)
			//r_index_concat(input);
			r_index_concat_bigbwt(input);
		else
			//r_index_multidollar_wrong(input);
			r_index_multidollar(input);
	}

	void r_index_concat(string &input_file){

		std::ifstream fs(input_file);
		std::stringstream buffer;
		buffer << fs.rdbuf();

		std::string input = buffer.str();

		this->sais = true;

		if(contains_reserved_chars(input)){

			cout << "Error: input string contains one of the reserved characters 0x0, 0x1" << endl;
			exit(1);

		}

		cout << "Text length = " << input.size() << endl << endl;

		cout << "(1/3) Building BWT and computing SA samples";
		//if(sais) cout << " (SE-SAIS) ... " << flush;
		//else cout << "(DIVSUFSORT) ... " << flush;

		//build run-length encoded BWT
		std::cout << "Processing = " << input << std::endl;
		auto bwt_and_samples = sufsort(input);

		string& bwt_s = get<0>(bwt_and_samples);
		vector<pair<ulint,ulint> >& samples_first_vec = get<1>(bwt_and_samples);
		vector<ulint>& samples_last_vec = get<2>(bwt_and_samples);

		for(int i=0;i<samples_first_vec.size();++i)
			std::cout << samples_first_vec[i].first << " - " << samples_last_vec[i] << std::endl;
		std::cout << bwt_s.size() << std::endl;
		std::cout << bwt_s << std::endl;

		cout << "done.\n(2/3) RLE encoding BWT ... " << flush;

		bwt = rle_string_t(bwt_s);

		//build F column
		F = vector<ulint>(256,0);
		for(uchar c : bwt_s)
			F[c]++;

		for(ulint i=255;i>0;--i)
			F[i] = F[i-1];

		F[0] = 0;

		for(ulint i=1;i<256;++i)
			F[i] += F[i-1];

		for(ulint i=0;i<bwt_s.size();++i)
			if(bwt_s[i]==TERMINATOR)
				terminator_position = i;

		assert(input.size()+1 == bwt.size());

		cout << "done. " << endl<<endl;

		r = bwt.number_of_runs();

		assert(samples_first_vec.size() == r);
		assert(samples_last_vec.size() == r);

		int log_r = bitsize(uint64_t(r));
		int log_n = bitsize(uint64_t(bwt.size()));

		cout << "Number of BWT equal-letter runs: r = " << r << endl;
		cout << "Rate n/r = " << double(bwt.size())/r << endl;
		cout << "log2(r) = " << log2(double(r)) << endl;
		cout << "log2(n/r) = " << log2(double(bwt.size())/r) << endl << endl;

		cout << "(3/3) Building phi function ..." << flush;

		//sort samples of first positions in runs according to text position
		std::sort(samples_first_vec.begin(), samples_first_vec.end());

		//build Elias-Fano predecessor
		{

			auto pred_bv = vector<bool>(bwt_s.size(),false);

			for(auto p : samples_first_vec){

				assert(p.first < pred_bv.size());
				std::cout << p.first << std::endl;
				pred_bv[p.first] = true;

			}

			pred = sparse_bv_type(pred_bv);

		}

		assert(pred.rank(pred.size()) == r);

		//last text position must be sampled
		assert(pred[pred.size()-1]);

		samples_last = int_vector<>(r,0,log_n); //text positions corresponding to last characters in BWT runs, in BWT order
		pred_to_run = int_vector<>(r,0,log_r); //stores the BWT run (0...R-1) corresponding to each position in pred, in text order

		for(ulint i=0;i<samples_last_vec.size();++i){

			assert(bitsize(uint64_t(samples_last_vec[i])) <= log_n);
			samples_last[i] = samples_last_vec[i];

		}

		for(ulint i=0;i<samples_first_vec.size();++i){

			assert(bitsize(uint64_t(samples_first_vec[i].second)) <= log_r);
			pred_to_run[i] = samples_first_vec[i].second;

		}

		cout << " done. " << endl<<endl;

	}

	/*
	 * Build index from Big-BWT output
	 */
	void r_index_concat_bigbwt(string &input){

		cout << "(1/2) RLE encoding BWT (and compute C vector) ... " << flush;

		string bwt_file = input + ".bwt";
		std::ifstream bwt_s(bwt_file);

		bwt_s.seekg(0, std::ios::end);
		ulint bwt_size = bwt_s.tellg(), text_size = 0;
		bwt_s.seekg(0, std::ios::beg);

		vector<char> heads; vector<ulint> lens;
		char prev, curr;
		ulint l = 1; this->r = 1;
		// F column
		F = vector<ulint>(256,0);

		//cout << "bwt size: " << bwt_size << endl;
		bwt_s.get(prev); F[prev]++;
		//cout << int(prev) << endl;
		for(ulint i=1;i<bwt_size;++i)
		{
			bwt_s.get(curr);
			curr = curr==0 ? 1 : curr; F[curr]++;
			//cout << int(curr) << endl;
			if(curr != prev)
			{
				heads.push_back(prev);
				lens.push_back(l); text_size += l;
				prev = curr;
				r++; l = 1;
				if(curr == 1)
					terminator_position = i;
			}
			else
				l++;
		}
		// insert the last run
		heads.push_back(prev);
		lens.push_back(l); text_size += l;

		assert(heads.size() == lens.size() == r);

		for(ulint i=255;i>0;--i)
			F[i] = F[i-1];

		F[0] = 0;

		for(ulint i=1;i<256;++i)
			F[i] += F[i-1];

		bwt = rle_string_t(heads,lens,bwt_size);
		assert(r == bwt.number_of_runs(););

		bwt_s.close();

		int log_r = bitsize(uint64_t(r));
		int log_n = bitsize(uint64_t(bwt_size));

		cout << "BWT size: " << bwt_size << endl;
		cout << "Number of BWT equal-letter runs: r = " << r << endl;
		cout << "Rate n/r = " << double(bwt_size/(double)r) << endl;
		cout << "log2(r) = " << log2(double(r)) << endl;
		cout << "log2(n/r) = " << log2(double(bwt_size/(double)r)) << endl << endl;

		cout << "(2/2) Building phi function ...\n" << flush;

		string ssa_file = input + ".ssa";
		std::ifstream ssa_s(ssa_file, std::ios::binary);
		string esa_file = input + ".esa";
		std::ifstream esa_s(esa_file, std::ios::binary);

		samples_last = int_vector<>(r,0,log_n);
		vector<ulint> Samples_first = vector<ulint>(r,0);
		uint64_t sa; ulint clen = 0;

		for(ulint i=0;i<heads.size();++i)
		{
			esa_s.read(reinterpret_cast<char*>(&sa), 5);
			assert(sa == i);
			esa_s.read(reinterpret_cast<char*>(&sa), 5);
			//std::cout << sa << " ";
			if(sa == 0)
				sa = text_size-1;
			else
				sa--;
			samples_last[i] = sa;
			ssa_s.read(reinterpret_cast<char*>(&sa), 5);
			assert(sa == i);
			ssa_s.read(reinterpret_cast<char*>(&sa), 5);
			//std::cout << sa << std::endl;
			if(sa == 0)
				sa = text_size-1;
			else
				sa--;
			Samples_first[i] = sa;
		}
		ssa_s.close(); esa_s.close();


		//sort samples at the beginning of a BWT run
		vector<ulint> idx_vec = sort_indexes(Samples_first);
		//build Elias-Fano predecessor
		{

			auto pred_bv = vector<bool>(bwt_size,false);

			for(auto p : Samples_first){

				assert(p < pred_bv.size());
				pred_bv[p] = true;

			}

			pred = sparse_bv_type(pred_bv);
		}

		assert(pred.rank(pred.size()) == r);

		//last text position must be sampled
		assert(pred[pred.size()-1]);

		pred_to_run = int_vector<>(r,0,log_r); //stores the BWT run (0...R-1) corresponding to each position in pred, in text order

		for(ulint i=0;i<idx_vec.size();++i){

			assert(bitsize(uint64_t(idx_vec[i])) <= log_r);
			pred_to_run[i] = idx_vec[i];

		}

		cout << " done. " << endl<<endl;
	}

	/*
	 * Build index for BWT with different dollars
	 */
	void r_index_multidollar(string &input){

		cout << endl << "(1/8) computing the RLE multidollar BWT ..." << endl << flush;

		string bwt_file = input + ".ebwt";
		std::ifstream bwt_s(bwt_file);

		bwt_s.seekg(0, std::ios::end);
		ulint bwt_size = bwt_s.tellg();
		bwt_s.seekg(0, std::ios::beg);

		vector<char> heads; vector<ulint> lens;
		char prev, curr;
		ulint l = 1; this->r = 1;
		// F column
		F = vector<ulint>(256,0);

		bwt_s >> prev; F[prev]++;
		//cout << "prev: " << prev << endl;
		for(ulint i=1;i<bwt_size;++i)
		{
			bwt_s >> curr; F[curr]++;
			if(curr != prev or curr == '$')
			{
				heads.push_back(prev);
				lens.push_back(l);
				prev = curr;
				r++; l = 1;
			}
			else
				l++;
		}
		if(curr != prev or curr == '$')
		{
			heads.push_back(prev);
			lens.push_back(l);
		}

		cout << "(2/8) computing C vector ..." << endl << flush;

		assert(heads.size() == lens.size());
		assert(heads.size() == r);

		for(ulint i=255;i>0;--i)
			F[i] = F[i-1];

		F[0] = 0;

		for(ulint i=1;i<256;++i)
			F[i] += F[i-1];

		cout << "(3/8) Computing rank select RLE BWT data structure ..." << endl << flush;

		bwt = rle_string_t(heads,lens,bwt_size);

		bwt_s.close();

		cout << "(4/8) Find terminator_position ..." << endl << flush;

		ulint curr_index = 0;
		while(bwt[curr_index] != '$')
		{
			curr = bwt[curr_index];
			curr_index = bwt.rank(curr_index+1,curr) + F[curr] - 1;
		}
		terminator_position = curr_index;

		int log_r = bitsize(uint64_t(r));
		int log_n = bitsize(uint64_t(bwt_size));

		cout << endl << "BWT size = " << bwt_size << endl;
		cout << "Number of BWT equal-letter runs: r = " << r << endl;
		cout << "Rate n/r = " << double(bwt_size/(double)r) << endl;
		cout << "log2(r) = " << log2(double(r)) << endl;
		cout << "log2(n/r) = " << log2(double(bwt_size/(double)r)) << endl << endl;

		cout << "(5/8) Building phi function ..." << endl << flush;

		string len_file = input + ".len";
		std::ifstream len_s(len_file);
		len_s.seekg(0, std::ios::end);
		ulint no_strings = len_s.tellg();
		len_s.seekg(0, std::ios::beg);

		std::vector<ulint> offsets(no_strings,0);
		uint8_t cl;
		for(ulint i=1;i<no_strings;++i)
		{
			len_s.read(reinterpret_cast<char*>(&cl), sizeof(uint8_t));
			offsets[i] = offsets[i-1] + cl + 1;
		}

		string da_file = input + ".da";
		std::ifstream da_s(da_file, std::ios::binary);
		string sa_file = input + ".posSA";
		std::ifstream sa_s(sa_file);

		cout << "(6/8) Read sample files ..." << endl << flush;

		samples_last = int_vector<>(r,0,log_n);
		vector<ulint> Samples_first = vector<ulint>(r,0);
		uint32_t da; uint8_t sa; ulint clen = 0;

		sa_s.read(reinterpret_cast<char*>(&sa), sizeof(uint8_t));
		da_s.read(reinterpret_cast<char*>(&da), sizeof(uint32_t));

		Samples_first[0] =	offsets[da] + sa;
		for(ulint i=0;i<heads.size();++i)
		{
			clen += lens[i];
			// set file pointer in the correct position
			da_s.seekg((clen-1)*sizeof(uint32_t), std::ios::beg);
			sa_s.seekg((clen-1)*sizeof(uint8_t), std::ios::beg);

			sa_s.read(reinterpret_cast<char*>(&sa), sizeof(uint8_t));
			da_s.read(reinterpret_cast<char*>(&da), sizeof(uint32_t));
			samples_last[i] = offsets[da] + sa;

			sa_s.read(reinterpret_cast<char*>(&sa), sizeof(uint8_t));
			da_s.read(reinterpret_cast<char*>(&da), sizeof(uint32_t));
			Samples_first[i+1] = offsets[da] + sa;
		}

		clen += lens[lens.size()-1];
		da_s.seekg((clen-1)*sizeof(uint32_t), std::ios::beg);
		sa_s.seekg((clen-1)*sizeof(uint8_t), std::ios::beg);

		sa_s.read(reinterpret_cast<char*>(&sa), sizeof(uint8_t));
		da_s.read(reinterpret_cast<char*>(&da), sizeof(uint32_t));
		samples_last[samples_last.size()-1] = offsets[da] + sa;

		len_s.close(); da_s.close(); sa_s.close();

		cout << "(7/8) Sorting first samples ..." << endl << flush;

		//sort samples at the beginning of a BWT run
		vector<ulint> idx_vec = sort_indexes(Samples_first);

		cout << "(8/8) Elias-Fano predecessor data structure ..." << endl << flush;

		//build Elias-Fano predecessor
		{

			auto pred_bv = vector<bool>(bwt_size,false);

			for(auto p : Samples_first)
			{
				assert(p < pred_bv.size());
				pred_bv[p] = true;
			}

			pred = sparse_bv_type(pred_bv);
		}

		assert(pred.rank(pred.size()) == r);

		//last text position must be sampled
		assert(pred[pred.size()-1]);

		//stores the BWT run (0...R-1) corresponding to each position in pred, in text order
		pred_to_run = int_vector<>(r,0,log_r); 

		for(ulint i=0;i<idx_vec.size();++i)
		{
			assert(bitsize(uint64_t(idx_vec[i])) <= log_r);
			pred_to_run[i] = idx_vec[i];
		}

		cout << " done. " << endl << endl;
	}

	/*
	 * Build index for BWT with different dollars
	 */
	void r_index_multidollar_old(string &input){

		cout << "(1/2) RLE encoding BWT (and compute C vector) ... " << flush;

		string bwt_file = input + ".ebwt";
		std::ifstream bwt_s(bwt_file);

		bwt_s.seekg(0, std::ios::end);
		ulint bwt_size = bwt_s.tellg();
		bwt_s.seekg(0, std::ios::beg);

		vector<char> heads; vector<ulint> lens;
		char prev, curr;
		ulint l = 1; this->r = 1;
		// F column
		F = vector<ulint>(256,0);

		bwt_s >> prev; F[prev]++;
		//cout << "prev: " << prev << endl;
		for(ulint i=1;i<bwt_size;++i)
		{
			bwt_s >> curr; F[curr]++;
			if(curr != prev or curr == '$')
			{
				heads.push_back(prev);
				lens.push_back(l);
				prev = curr;
				r++; l = 1;
				//if(curr == '$')
					//terminator_position = i;
			}
			else
				l++;
		}
		if(curr != prev or curr == '$')
		{
			heads.push_back(prev);
			lens.push_back(l);
		}

		assert(heads.size() == lens.size());

		for(ulint i=255;i>0;--i)
			F[i] = F[i-1];

		F[0] = 0;

		for(ulint i=1;i<256;++i)
			F[i] += F[i-1];

		bwt = rle_string_t(heads,lens,bwt_size);

		bwt_s.close();

		ulint curr_index = 0;
		while(bwt[curr_index] != '$')
		{
			curr = bwt[curr_index];
			curr_index = bwt.rank(curr_index+1,curr) + F[curr] - 1;
		}
		terminator_position = curr_index;
		//cout << "term pos: " << terminator_position << endl;
		//exit(1);

		int log_r = bitsize(uint64_t(r));
		int log_n = bitsize(uint64_t(bwt_size));

		cout << "BWT size: " << bwt_size << endl;
		cout << "Number of BWT equal-letter runs: r = " << r << endl;
		cout << "Rate n/r = " << double(bwt_size/(double)r) << endl;
		cout << "log2(r) = " << log2(double(r)) << endl;
		cout << "log2(n/r) = " << log2(double(bwt_size/(double)r)) << endl << endl;

		cout << "(2/2) Building phi function ..." << flush;

		string len_file = input + ".len";
		std::ifstream len_s(len_file);
		len_s.seekg(0, std::ios::end);
		ulint no_strings = len_s.tellg();
		len_s.seekg(0, std::ios::beg);

		int_vector<> offsets = int_vector<>(no_strings,0,log_n); 
		uint8_t cl;
		for(ulint i=1;i<no_strings;++i)
		{
			len_s >> cl;
			offsets[i] = offsets[i-1] + cl + 1;
		}

		string da_file = input + ".da";
		std::ifstream da_s(da_file, std::ios::binary);
		string sa_file = input + ".posSA";
		std::ifstream sa_s(sa_file);

		samples_last = int_vector<>(r,0,log_n);
		vector<ulint> Samples_first = vector<ulint>(r,0);
		//pred_to_run = int_vector<>(r,0,log_r);
		uint32_t da; uint8_t sa; ulint clen = 0;

		sa_s >> sa;
		da_s.read(reinterpret_cast<char*>(&da), sizeof(uint32_t));
		Samples_first[0] =	offsets[da] + sa;
		std::cout << Samples_first[0] << " ";
		for(ulint i=0;i<heads.size();++i)
		{
			clen += lens[i];
			// get end
			da_s.seekg((clen-1)*sizeof(uint32_t), std::ios::beg);
			sa_s.seekg(clen-1, std::ios::beg);

			sa_s >> sa;
			da_s.read(reinterpret_cast<char*>(&da), sizeof(uint32_t));
			samples_last[i] = offsets[da] + sa;

			sa_s >> sa;
			da_s.read(reinterpret_cast<char*>(&da), sizeof(uint32_t));
			Samples_first[i+1] = offsets[da] + sa;
			std::cout << Samples_first[i+1] << " ";
		}

		clen += lens[lens.size()-1];
		da_s.seekg((clen-1)*sizeof(uint32_t), std::ios::beg);
		sa_s.seekg(clen-1, std::ios::beg);
		sa_s >> sa;
		da_s.read(reinterpret_cast<char*>(&da), sizeof(uint32_t));
		samples_last[samples_last.size()-1] = offsets[da] + sa;

		len_s.close(); da_s.close(); sa_s.close();

		//sort samples at the beginning of a BWT run
		vector<ulint> idx_vec = sort_indexes(Samples_first);


		cout << "bwt_size: " << bwt_size << endl;
		ulint no_bv = 0;
		//build Elias-Fano predecessor
		{

			auto pred_bv = vector<bool>(bwt_size,false);

			for(auto p : Samples_first){

				assert(p < pred_bv.size());
				no_bv++;
				pred_bv[p] = true;

			}
			std::cout << "no_bv: " << no_bv << std::endl;
			pred = sparse_bv_type(pred_bv);
			for(ulint i=0;i<pred_bv.size();++i)
				std::cout << pred_bv[i];
			std::cout << std::endl;
		}

		assert(pred.rank(pred.size()) == r);

		//last text position must be sampled
		assert(pred[pred.size()-1]);

		pred_to_run = int_vector<>(r,0,log_r); //stores the BWT run (0...R-1) corresponding to each position in pred, in text order

		for(ulint i=0;i<idx_vec.size();++i){

			assert(bitsize(uint64_t(idx_vec[i])) <= log_r);
			pred_to_run[i] = idx_vec[i];

		}

		cout << " done. " << endl<<endl;
	}

	/*
	 * get full BWT range
	 */
	range_t full_range(){

		//inclusive range
		return {0,bwt_size()-1};

	}

	uchar operator[](ulint i ){
		return bwt[i];
	}

	/*
	 * \param r inclusive range of a string w
	 * \param c character
	 * \return inclusive range of cw
	 */
	range_t LF(range_t rn, uchar c){

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==bwt_size()) || F[c]>=F[c+1])
			return {1,0};

		//number of c before the interval
		ulint c_before = bwt.rank(rn.first,c);
		//std::cout << "c before " << bwt.rank(rn.first,c) << std::endl;

		//number of c inside the interval rn
		ulint c_inside = bwt.rank(rn.second+1,c) - c_before;
		//std::cout << "c inside " << bwt.rank(rn.second,c) << std::endl;

		//if there are no c in the interval, return empty range
		if(c_inside==0) return {1,0};

		ulint l = F[c] + c_before;

		return {l,l+c_inside-1};

	}

	/*
	 * Phi function. Phi(SA[0]) is undefined
	 */
	ulint Phi(ulint i){

		assert(i != bwt.size()-1);

		//jr is the rank of the predecessor of i (circular)
		ulint jr = pred.predecessor_rank_circular(i);

		assert(jr<=r-1);

		//the actual predecessor
		ulint j = pred.select(jr);

		assert(jr<r-1 or j == bwt.size()-1);

		//distance from predecessor
		ulint delta = j<i ? i-j : i+1;

		//cannot fall on first run: this can happen only if I call Phi(SA[0])
		assert(pred_to_run[jr]>0);

		//sample at the end of previous run
		assert(pred_to_run[jr]-1 < samples_last.size());
		ulint prev_sample = samples_last[ pred_to_run[jr]-1 ];

		return (prev_sample + delta) % bwt.size();

	}

	//backward navigation of the BWT
	ulint LF(ulint  i){

		auto c = bwt[i];
		return F[c] + bwt.rank(i,c);

	}

	//forward navigation of the BWT
	ulint FL(ulint  i){

		//i-th character in first BWT column
		auto c = F_at(i);

		//this c is the j-th (counting from 0)
		ulint j = i - F[c];

		return bwt.select(j,uchar(c));

	}

	//forward navigation of the BWT, where for efficiency we give c=F[i] as input
	ulint FL(ulint  i, uchar c){

		//i-th character in first BWT column
		assert(c == F_at(i));

		//this c is the j-th (counting from 0)
		ulint j = i - F[c];

		return bwt.select(j,uchar(c));

	}

	/*
	 * access column F at position i
	 */
	uchar F_at(ulint i){

		ulint c = (upper_bound(F.begin(),F.end(),i) - F.begin()) - 1;
		assert(c<256);
		assert(i>=F[c]);

		return uchar(c);

	}

	/*
	 * Return BWT range of character c
	 */
	range_t get_char_range(uchar c){

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==bwt_size()) || F[c]>=F[c+1])
			return {1,0};

		ulint l = F[c];
		ulint r = bwt_size()-1;

		if(c<255)
			r = F[c+1]-1;

		return {l,r};

	}

	/*
	 * Return BWT range of pattern P
	 */
	range_t count(string &P){

		auto range = full_range();
		ulint m = P.size();
		//std::cout << range.first << " - " << range.second << std::endl;

		for(ulint i=0;i<m and range.second>=range.first;++i)
		{
			range = LF(range,P[m-i-1]);
			//std::cout << range.first << " - " << range.second << std::endl;
		}

		return range;

	}

	/*
	 * Return number of occurrences of P in the text
	 */
	ulint occ(string &P){

		auto rn = count(P);

		return rn.second>=rn.first ? (rn.second-rn.first)+1 : 0;

	}

	/*
	 * iterator locate(string &P){
	 *
	 * 	return iterator to iterate over all occurrences without storing them
	 * 	in memory
	 *
	 * }
	 */

	/*
	 * locate all occurrences of P and return them in an array
	 * (space consuming if result is big).
	 */
	vector<ulint> locate_all(string& P){

		vector<ulint> OCC;

		pair<range_t, ulint> res = count_and_get_occ(P);

		ulint L = std::get<0>(res).first;
		ulint R = std::get<0>(res).second;
		ulint k = std::get<1>(res);	//SA[R]

		ulint n_occ = R>=L ? (R-L)+1 : 0;

		if(n_occ>0){

			OCC.push_back(k);

			for(ulint i=1;i<n_occ;++i){

				k = Phi(k);
				OCC.push_back(k);

			}

		}

		return OCC;

	}


	/*
	 * get number of runs in the BWT (terminator character included)
	 */
	ulint number_of_runs(){
		return bwt.number_of_runs();
	}

	/*
	 * get terminator (0x1) position in the BWT
	 */
	ulint get_terminator_position(){
		return terminator_position;
	}

	/*
	 * get BWT in string format
	 */
	string get_bwt(){
		return bwt.toString();
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out)
	{
		ulint w_bytes = 0;

		assert(F.size()>0);
		assert(bwt.size()>0);

		out.write((char*)&terminator_position,sizeof(terminator_position));
		out.write((char*)F.data(),256*sizeof(ulint));

		w_bytes += sizeof(terminator_position) + 256*sizeof(ulint);
		//std::cout << w_bytes << " ";

		w_bytes += bwt.serialize(out);
		//std::cout << bwt.serialize(out) << " ";

		w_bytes += pred.serialize(out);
		//std::cout << pred.serialize(out) << " ";
		w_bytes += samples_last.serialize(out);
		//std::cout << samples_last.serialize(out) << " ";
		w_bytes += pred_to_run.serialize(out);
		//std::cout << pred_to_run.serialize(out) << std::endl;

		//std::cout << "end serialize: " << w_bytes << std::endl;
		return w_bytes;
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&terminator_position,sizeof(terminator_position));

		F = vector<ulint>(256);
		in.read((char*)F.data(),256*sizeof(ulint));

		bwt.load(in);
		r = bwt.number_of_runs();
		pred.load(in);
		samples_last.load(in);
		pred_to_run.load(in);

	}

	/*
	 * save the structure to the path specified.
	 * \param path_prefix prefix of the index files. suffix ".ri" will be automatically added
	 */
	void save_to_file(string path_prefix){

		string path = string(path_prefix).append(".ri");

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	/*
	 * load the structure from the path specified.
	 * \param path: full file name
	 */
	void load_from_file(string path){
		std::cout << "path= " << path << std::endl;
		std::ifstream in(path);
		cout << "start loading" << endl;
		load(in);
		cout << "end loading" << endl;
		in.close();

	}

	ulint text_size(){
		return bwt.size()-1;
	}

	ulint bwt_size(){
		return bwt.size();
	}

	uchar get_terminator(){
		return TERMINATOR;
	}

	ulint print_space(){

		cout << "Number of runs = " << bwt.number_of_runs() << endl<<endl;

		ulint tot_bytes = bwt.print_space();

		cout << "\nTOT BWT space: " << tot_bytes << " Bytes" <<endl<<endl;

		return tot_bytes;

	}

private:

	/*
	 * returns <<l,r>, SA[r] >, where l,r are the inclusive ranges of the pattern P. If P does not occur, then l>r
	 *
	 * returns <range, j,k>
	 *
	 */
	pair<range_t, ulint> count_and_get_occ(string &P){

		//k = SA[r]
		ulint k = 0;

		range_t range = full_range();
		assert(r-1 < samples_last.size());
		k = (samples_last[r-1]+1) % bwt.size();

		range_t range1;

		ulint m = P.size();

		for(ulint i=0;i<m and range.second>=range.first;++i){

			uchar c = P[m-i-1];

			range1 = LF(range,c);

			//if suffix can be left-extended with char
			if(range1.first <= range1.second){


				if(bwt[range.second] == c){

					// last c is at the end of range. Then, we have this sample by induction!
					assert(k>0);
					k--;

				}else{

					//find last c in range (there must be one because range1 is not empty)
					//and get its sample (must be sampled because it is at the end of a run)
					//note: by previous check, bwt[range.second] != c, so we can use argument range.second
					ulint rnk = bwt.rank(range.second,c);

					//there must be at least one c before range.second
					assert(rnk>0);

					//this is the rank of the last c
					rnk--;

					//jump to the corresponding BWT position
					ulint j = bwt.select(rnk,c);

					//the c must be in the range
					assert(j>=range.first and j < range.second);

					//run of position j
					ulint run_of_j = bwt.run_of_position(j);

					k = samples_last[run_of_j];

				}

			}

			range = range1;

		}

		return {range, k};

	}

	/*
	 * returns a triple containing BWT of input string
	 * (uses 0x1 character as terminator), text positions corresponding
	 * to first letters in BWT runs (plus their ranks from 0 to R-1), and text positions corresponding
	 * to last letters in BWT runs (in BWT order)
	 */
	tuple<string, vector<pair<ulint, ulint> >, vector<ulint> > sufsort(string &s){

		string bwt_s;

	    cache_config cc;

	    int_vector<8> text(s.size());
	    assert(text.size()==s.size());

	    for(ulint i=0;i<s.size();++i)
	    	text[i] = (uchar)s[i];

	    assert(text.size()==s.size());

	    append_zero_symbol(text);

	    store_to_cache(text, conf::KEY_TEXT, cc);

	    construct_config::byte_algo_sa = sais ? SE_SAIS : LIBDIVSUFSORT;
	    construct_sa<8>(cc);

	    //now build BWT from SA
	    int_vector_buffer<> sa(cache_file_name(conf::KEY_SA, cc));

	    vector<pair<ulint, ulint> > samples_first; 	//text positions corresponding to first characters in BWT runs, and their ranks 0...R-1
	    vector<ulint> samples_last;	 				//text positions corresponding to last characters in BWT runs

	    {

	        for (ulint i=0; i<sa.size(); i++){
	            auto x = sa[i];

	            assert(x<=text.size());

	            if ( x > 0 )
	            	bwt_s.push_back((uchar)text[x-1]);
	            else
	            	bwt_s.push_back(TERMINATOR);

	            //Insert samples at begin of runs
	            if(i>0){

	            	if(	i==1 ||									//case 1: i-1 == 0 is at run begin
						(i>1 && bwt_s[i-1] != bwt_s[i-2])		//case 2: i-1 is at the begin of a run
	            	){

	            		samples_first.push_back( {sa[i-1]>0?sa[i-1]-1:sa.size()-1, samples_first.size()} );

	            	}

	            	//check last BWT letter
	            	if(i==sa.size()-1 && bwt_s[i]!=bwt_s[i-1]) samples_first.push_back( {sa[i]>0?sa[i]-1:sa.size()-1, samples_first.size()} );

				}

	            //Insert samples at end of runs
	            if(i>0){

	            	if(	bwt_s[i-1] != bwt_s[i]		//i-1 is at the end of a run
	            	){

	            		samples_last.push_back( sa[i-1]>0?sa[i-1]-1:sa.size()-1 );

	            	}

	            	//last BWT letter is always at end of a run and is never checked in the previous if
	            	if(i==sa.size()-1) samples_last.push_back( sa[i]>0?sa[i]-1:sa.size()-1 );

	            }

	        }

	    }

	    assert(samples_first.size() == samples_last.size());

	    sdsl::remove(cache_file_name(conf::KEY_TEXT, cc));
	    sdsl::remove(cache_file_name(conf::KEY_SA, cc));

	    return std::make_tuple(std::move(bwt_s), std::move(samples_first), std::move(samples_last));

	}

	static bool contains_reserved_chars(string &s){

		for(auto c : s)
			if(c == 0 or c == 1)
				return true;

		return false;

	}

	static const uchar TERMINATOR = 1;

	bool sais = true;

	/*
	 * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
	 */

	//F column of the BWT (vector of 256 elements)
	vector<ulint> F;
	//L column of the BWT, run-length compressed
	rle_string_t bwt;
	ulint terminator_position = 0;
	ulint r = 0;//number of BWT runs


	//the predecessor structure on positions corresponding to first chars in BWT runs
	sparse_bv_type pred;
	int_vector<> samples_last; //text positions corresponding to last characters in BWT runs, in BWT order
	int_vector<> pred_to_run; //stores the BWT run (0...R-1) corresponding to each position in pred, in text order

};

}

#endif /* R_INDEX_S_H_ */