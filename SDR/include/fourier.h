

#ifndef DY4_FOURIER_H
#define DY4_FOURIER_H

// Add headers as needed
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

// Declaration of function prototypes
void DFT(const std::vector<real> &,
	std::vector<std::complex<real>> &);

// You should add your own IDFT
// time-permitting you can build your own function for FFT

void computeVectorMagnitude(const std::vector<std::complex<real>> &,
	std::vector<real> &);

// Provide the prototype to estimate PSD
// ...

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void DFT_reference(const std::vector<real>& x, std::vector<std::complex<real>>& Xf);
void DFT_init_bins(const std::vector<real>& x, std::vector<std::complex<real>>& Xf);

void generate_DFT_twiddles(const int& N, std::vector<std::complex<real>> &Twiddle1D);
void generate_DFT_matrix(const int& N, std::vector<std::vector<std::complex<real>>>& Twiddle2D);
void generate_DFT_twiddles(const int& N, std::vector<std::complex<real>> &Twiddle1D);
void generate_DFT_matrix(const int& N, std::vector<std::vector<std::complex<real>>> &Twiddle2D);
void matrixPSD_precomputes(int nfft, std::vector<real> &hann_window, std::vector<std::vector<std::complex<real>>> &dft_matrix);
void matrixPSD_optimized_precompute(const std::vector<real> &samples, std::vector<real> &psd_est, std::vector<real> &freq, int nfft, real Fs, std::vector<real> &hann_window, std::vector<std::vector<std::complex<real>>> &dft_matrix);
void estimatePSD(std::vector<real> &samples, std::vector<real> &psd_result, std::vector<real> &freq_result, unsigned int nfft, unsigned int Fs);
void lightweight_fmDemodArctan(std::vector<real> &fm_demod, real &phase_diff, const std::vector<real> &I, const std::vector<real> &Q);

#endif // DY4_FOURIER_H
