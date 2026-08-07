/* ******************************************************************
 * A simple function that takes in input the optimal BWT and writes the 
 * original string collections in Bentely et al. order.
 * ****************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <assert.h>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <algorithm>

#include <sdsl/wavelet_trees.hpp>

#define TERMINATE_CHAR '$'     //it is the symbol used as "end of strings"

using namespace std;

const int alph_size = 128;

void remapToOptimalOrder(string wt_input, string output_file)
{
    cout << "Building the Wavelet tree of the optimal BWT." << endl;
    cout << wt_input << endl;
    sdsl::wt_blcd<> ebwt; sdsl::construct(ebwt,wt_input,1);

    vector<uint32_t> C(alph_size,0); vector<uint32_t> bkt(alph_size,0);
    cout << "Building the C vector of the optimal BWT." << endl;
    for(int i=0;i<ebwt.size();i++){ bkt[ebwt[i]]++; }
    for(int i=1;i<C.size();i++){ C[i]=C[i-1]+bkt[i-1]; }

    uint32_t noStrings = bkt[TERMINATE_CHAR];
    if(noStrings==0){
        cerr << "WARNING: TERMINATE_CHAR must be "<< TERMINATE_CHAR << endl;
        exit(1);
    }
    FILE* fp = fopen(output_file.c_str(), "w+");

    for(size_t i=0; i<noStrings; ++i)
    {
        std::vector<uint8_t> RP;  
        size_t index = i; 
        uint8_t p = ebwt[index]; 
        RP.push_back(p);
        size_t starting = index;
        index = C[p] + ebwt.rank(index,p);

        while(true)
        {
            p = ebwt[index];
            if(p == TERMINATE_CHAR)
                break;
            RP.push_back(p);
            index = C[p] + ebwt.rank(index,p);
        }

        std::string header = std::string("> SEQUENCE ") + std::to_string(i) + "\n";
        if((fwrite(&header[0], sizeof(uint8_t), header.size(), fp))!=header.size()) {cerr << "fwrite failed" << endl;}
        reverse(RP.begin(),RP.end()); RP.push_back('\n');
        if((fwrite(&RP[0], sizeof(uint8_t), RP.size(), fp))!=RP.size()) {cerr << "fwrite failed" << endl;}
    }
}

int main(int argc, char *argv[])
{  
  time_t start_wc = time(NULL);
  
    // check input data
  if(argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " input output" << std::endl;
    std::cout << "Invert the optBWT of the string collection and output the original string collection to name.rfasta" << std::endl;
    exit(1);
  }
  std::cout << "==== Command line:" << std::endl;
  for(int i=0;i<argc;i++)
    std::cout << " " << argv[i];
  std::cout << std::endl;
  
  string input_file = argv[1];
  string output_file = argv[2];

  remapToOptimalOrder(input_file, output_file);
  
  return 0;
}
