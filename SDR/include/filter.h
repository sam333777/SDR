

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// Add headers as needed
#include <iostream>
#include <vector>

// Declaration of function prototypes
void impulseResponseLPF(real, real, unsigned short int, std::vector<real> &);
void impulseResponseLPFwithResamplerFactor(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h, const int up_sample_coeff);
void convolveFIR(std::vector<real> &, const std::vector<real> &, const std::vector<real> &);
void impulseResponseBPF(real Fs, real Fb,real Fe, unsigned short int num_taps, std::vector<real> &h);
void impulseResponseBPFwithResamplerFactor(real Fs, real Fb,real Fe, unsigned short int num_taps, std::vector<real> &h);
void allPassFilter(std::vector<real> &y, std::vector<real> &x, std::vector<real> &state);
void block_conv(std::vector<real> &y, const std::vector<real> &h, const std::vector<real> &x, std::vector<real> &state);
void fmPllwithState(std::vector<real> &pllIn, real freq, real Fs, real ncoScale, real phaseAdjust, real normBandwidth, real &saveI, real &saveQ, real &integrator, real &eastmatePhase, real &trigOffset, std::vector<real> &ncoOut, real &ncoFirst);
//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void convolveFIR_inefficient(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void convolveFIR_reference(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h);
void single_pass_convoul_with_precal_range(const std::vector<real> &x, const std::vector<real> &h, std::vector<real> &Xf);
void single_pass_convoul_with_unrolling(const std::vector<real> &x, const std::vector<real> &h, std::vector<real> &Xf);
void single_pass_convoul_with_unrolling(const std::vector<real> &x, const std::vector<real> &h, std::vector<real> &Xf);
void blockConvolution(std::vector<real>& y, const std::vector<real>& x, const std::vector<real>& h, int blockSize);
void resampler(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h, const int down_sample_coeff, const int up_sample_coeff, std::vector<real> &state);
void resampler_optimized(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h, const int down_sample_coeff, const int up_sample_coeff, std::vector<real> &state);
void downSampler(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h, const int down_sample_coeff, std::vector<real> &state);
void fmPll(const std::vector<real>& pllIn, real freq, real Fs, real& nextncoOut, std::vector<real>& ncoOut, real ncoScale, real phaseAdjust, real normBandwidth);
void impulseResponseRootRaisedCosine(std::vector<real> &y, real Fs, unsigned short int nTaps);


#endif // DY4_FILTER_H



