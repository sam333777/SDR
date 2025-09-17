

#ifndef DY4_IOFUNC_H
#define DY4_IOFUNC_H

// Add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

// Declaration of function prototypes
void printRealVector(const std::vector<real> &);

void printComplexVector(const std::vector<std::complex<real>> &);

void readBinData(const std::string, std::vector<real> &);

void writeBinData(const std::string, const std::vector<real> &);

void readSdrRaw(const std::string in_fname, std::vector<char> &bin_data);

void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<real> &block_data);

void writeStoutBlockData(int block_size, std::vector<real> &processed_data);

void writeBinData(const std::string &out_fname, const std::vector<float> &bin_data);


void writeStoutBlockData(int block_size, std::vector<real> &processed_data);


//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void generate_random_values(std::vector<real>& x, const real& lower_bound, const real& upper_bound);

#endif // DY4_IOFUNC_H
