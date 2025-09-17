

#include "dy4.h"
#include "iofunc.h"
#include <stdint.h>

// Some basic functions for printing information from vectors
// or to read from/write to binary files in real format

void printRealVector(const std::vector<real> &x) {
	std::cout << "Printing real vector of size " << x.size() << "\n";
	for (int i = 0; i < (int)x.size(); i++) {
		std::cout << x[i] << " ";
	}
	std::cout << "\n";
}

void printComplexVector(const std::vector<std::complex<real>> &X) {
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (int i = 0; i < (int)X.size(); i++) {
		std::cout << X[i] << " ";
	}
	std::cout << "\n";
}

// Assumes data in the raw binary file is in the format corresponding to 'real'
void readBinData(const std::string in_fname, std::vector<real> &bin_data)
{
	std::ifstream fdin(in_fname, std::ios::binary);
	if (!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cout << "Reading raw binary from \"" << in_fname << "\"\n";
	}

	fdin.seekg(0, std::ios::end);
	const unsigned int num_samples = fdin.tellg() / sizeof(real);

	bin_data.resize(num_samples);
	fdin.seekg(0, std::ios::beg);
	fdin.read(reinterpret_cast<char *>(&bin_data[0]), num_samples * sizeof(real));
	fdin.close();
}

// Assumes data in the raw binary file is in the format corresponding to 'real'
void writeBinData(const std::string out_fname, const std::vector<real> &bin_data)
{
	std::cout << "Writing raw binary to \"" << out_fname << "\"\n";
	std::ofstream fdout(out_fname, std::ios::binary);
	for (int i = 0; i < (int)bin_data.size(); i++) {
		fdout.write(reinterpret_cast<const char *>(&bin_data[i]), sizeof(bin_data[i]));
	}
	fdout.close();
}

//Modified from project iofunc
//SDR will output 8-bit unsigned int
void readSdrRaw(const std::string in_fname, std::vector<char> &bin_data)
{
	std::ifstream fdin(in_fname, std::ios::binary);
	if(!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cout << "Reading raw binary from \"" << in_fname << "\"\n";
	}
	fdin.seekg(0, std::ios::end);
	const unsigned int num_samples = fdin.tellg() / sizeof(char);

	bin_data.resize(num_samples);
	fdin.seekg(0, std::ios::beg);
	fdin.read(reinterpret_cast<char*>(&bin_data[0]), num_samples*sizeof(float));
	fdin.close();

}

void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<real> &block_data)
{
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));

	for (int k=0; k<(int)num_samples; k++)
	{
		block_data[k] = real(((unsigned char)raw_data[k]-128)/128.0);
	}
}

void writeStoutBlockData(int block_size, std::vector<real> &processed_data)
{
	std::vector<short int> audio_data(block_size);

	for (int k=0; k<processed_data.size(); k++)
	{
		if (std::isnan(processed_data[k]))
		{
			audio_data[k] = 0;
		}
		else
		{
			audio_data[k] = static_cast<short int>(processed_data[k] * 16384);
		}
	}

	fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
}

//Writing the data (16-bit little-endian) into binary file for aplay
void writeBinData(const std::string &out_fname, const std::vector<float> &bin_data)
{
    std::cout << "Writing 16-bit little-endian binary to \"" << out_fname << "\"\n";
    std::ofstream fdout(out_fname, std::ios::binary);
    
    if (!fdout) {
        std::cerr << "Error: Unable to open file for writing!" << std::endl;
        return;
    }

    for (float value : bin_data) {
        // Scale float (-1.0 to 1.0) to int16_t (-32768 to 32767)
        int16_t int16_value = static_cast<int16_t>(value * 32767.0f);


        fdout.write(reinterpret_cast<const char*>(&int16_value), sizeof(int16_value));
    }

    fdout.close();
}

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void generate_random_values(std::vector<real>& x, const real& lower_bound = -1., const real& upper_bound = 1.) {

	for (int i = 0; i < (int)x.size(); i++) {
		x[i] = lower_bound + (upper_bound - lower_bound) * ((real)(rand()) / RAND_MAX);
	}
}

