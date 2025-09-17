

#include "dy4.h"
#include "filter.h"
#include <cmath>

// Function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h)
{
	// Bring your own functionality
	h.clear();
	h.resize(num_taps, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab

	real Normcutoff = Fc / (Fs / 2);
	int middleindex = (num_taps - 1) / 2;

	for (int i = 0; i < num_taps; i++) 
	{
		if (i == middleindex)
		{
			h[i] = Normcutoff;
		}
		else
		{
			h[i] = Normcutoff * std::sin(PI * Normcutoff * (i - middleindex)) / (PI * Normcutoff * (i - middleindex));
		}
		h[i] *= (0.5 - 0.5 * std::cos(2* PI * i / (num_taps - 1)));
	}
}

void impulseResponseLPFwithResamplerFactor(real Fs, real Fc, unsigned short int num_taps, std::vector<real> &h, const int up_sample_coeff)
{
    real Fs_Resam = Fs * up_sample_coeff;
    unsigned short int num_taps_Resam = num_taps * up_sample_coeff;


    h.clear();
    h.resize(num_taps_Resam, 0.0);

    real Normcutoff = Fc/ (Fs_Resam/2);
    int middleindex = (num_taps_Resam - 1)/2;

	for (int i = 0; i < num_taps_Resam; i++) {
		if (i == (num_taps_Resam-1)/2) {
			h[i] = Normcutoff;
		}
		else {
			h[i] = Normcutoff * std::sin(PI * Normcutoff * (i-middleindex)) / (PI * Normcutoff * (i-middleindex));
		}
    h[i] *= up_sample_coeff * (0.5 - 0.5 * std::cos(2 * PI * i / (num_taps_Resam - 1)));
	}
}
// Function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h)
{
	y.clear();
	y.resize(x.size() + h.size() - 1, 0.0);

	// The rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < h.size(); ++j) {
			if (i+j < y.size())
            	y[i + j] += x[i] * h[j];
        }
    }
}

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void convolveFIR_inefficient(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);
    for (auto n = 0; n < (int)y.size(); n++) {
        for (auto k = 0; k < (int)x.size(); k++) {
            if ((n - k >= 0) && (n - k) < (int)h.size())
                y[n] += x[k] * h[n - k];
        }
    }
}

void convolveFIR_reference(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h) {
    y.clear();
    y.resize(int(x.size() + h.size()-1), 0.0);

    for (auto n = 0; n < (int)y.size(); n++) {
        for (auto k = 0; k < (int)h.size(); k++) {
            if ((n - k >= 0) && (n - k) < (int)x.size())
                y[n] += h[k] * x[n - k];
        }
    }
}

void blockProcessing(std::vector<real>& y, const std::vector<real>& x, const std::vector<real>& h, std::vector<real>& state, size_t position, size_t thisBlockSize, size_t hLength) {
    for (size_t n = 0; n < thisBlockSize; ++n) {
        real acc = 0.0;
        
        for (size_t k = 0; k < h.size(); ++k) {
            int idx = n - k;
            if (idx >= 0) {
                acc += x[position + idx] * h[k];
            } else {
                // Use the saved state
                int sIndex = (state.size() + idx);
                if (sIndex >= 0 && sIndex < (int)state.size()) {
                    acc += state[sIndex] * h[k];
                }
            }
        }
        y[position + n] = acc;
    }

    int overlapLength = (int)hLength - 1;
    for (int i = 0; i < overlapLength; i++) {
        int readIndex = (int)thisBlockSize - overlapLength + i;
        if (readIndex >= 0) {
            state[i] = x[position + readIndex];
        } else {
            state[i] = state[thisBlockSize + readIndex];
        }
    }
}

void blockConvolution(std::vector<real>& y, const std::vector<real>& x, const std::vector<real>& h, int blockSize) {
    y.resize(x.size() + h.size() - 1, 0.0);
    std::vector<real> state(h.size() - 1, 0.0);

    size_t position = 0;
    while (position < x.size()) {
        size_t thisBlockSize = std::min((size_t)blockSize, x.size() - position);

        blockProcessing(y, x, h, state, position, thisBlockSize, h.size());

        position += thisBlockSize;
    }
}

//resampler for block processing
void resampler(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h, const int down_sample_coeff, const int up_sample_coeff, std::vector<real> &state) 
{

    y.clear();
    y.resize(x.size() * up_sample_coeff / down_sample_coeff, 0.0);

    int index;
    int y_block_size = static_cast<int>(y.size());
    int h_size = static_cast<int>(h.size());


    for (int n = 0; n < y_block_size; n++) {
        for (int k = (n * down_sample_coeff) % up_sample_coeff; k < h_size; k += up_sample_coeff) {
            index = (n * down_sample_coeff - k) / up_sample_coeff;
            if (index >= 0) {
                y[n] += h[k] * x[index];
            } else {
                y[n] += h[k] * state[state.size() + index];
            }
        }
    }

    std::vector<real> temp(&x[x.size() - state.size()], &x[x.size()]);
    state = temp;
}

//resampler removed modulo from the calculation

void resampler_optimized(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h, const int down_sample_coeff, const int up_sample_coeff, std::vector<real> &state) {
    int ntaps = h.size();
    int sig_len = x.size();
    int output_len = sig_len * up_sample_coeff / down_sample_coeff;
    y.clear();
    y.resize(output_len, 0.0);

    int index;

    for (int i = 0; i < output_len; i++) //Y_index
    {
        for (int k = ((i * down_sample_coeff) % up_sample_coeff); k < ntaps; k += up_sample_coeff) //h_index
        {
            index = (i * down_sample_coeff - k)/up_sample_coeff;
            if (index >= 0)
            {
                y[i] += h[k] * x[index];
            }
            else 
            {
                y[i] += h[k] * state[state.size() + index];
            }
        }
    }

    std::vector<real> temp(&x[x.size()-state.size()], &x[x.size()]);
    state = temp;
}

//resampler for block processing
void downSampler(std::vector<real> &y, const std::vector<real> &x, const std::vector<real> &h, const int down_sample_coeff, std::vector<real> &state) 
{

    y.clear();
    y.resize(x.size() / down_sample_coeff, 0.0);

    int index = 0;

    for (int n = 0; n < x.size(); n+= down_sample_coeff) {
        for (int k = 0; k < h.size(); k++) {
            if (n - k >= 0) {
                y[index] += h[k] * x[n - k];
            } else {
                y[index] += h[k] * state[state.size() + n - k];
            }
        }
        index++;
    }

    for (int n = (x.size() - h.size()); n < x.size(); n++) 
    {
       state[n - (x.size() - h.size())] = x[n];
    }

}

void impulseResponseBPF(real Fs, real Fb,real Fe, unsigned short int num_taps, std::vector<real> &h)
{
    h.clear(); 
    h.resize(num_taps, 0.0);

    const real center = (num_taps - 1) / 2;
    const real fMid = (Fb + Fe) / (2.0);
    real scaleFactor = 0;
    for (unsigned short int k = 0; k < num_taps; k++){
	int n = k - center;
        if(n == 0){
            h[k] = (2*Fe)/Fs - (2*Fb)/Fs;
        } else{
            h[k] = (std::sin(2* PI * Fe * (n / Fs) ) / (PI*n)) - (std::sin(2* PI * Fb * (n / Fs) ) / (PI*n));
        }
	h[k] *= (0.5 - 0.5 * std::cos((2*PI*k)/(num_taps-1)));
	scaleFactor += h[k] * std::cos((2*PI*n*fMid)/Fs);
    }
    for (unsigned short int k = 0; k < num_taps - 1; k++){
	h[k] = h[k]/scaleFactor;
    }
}


void impulseResponseBPFwithResamplerFactor(real Fs, real Fc, real BW, unsigned short int num_taps, std::vector<real> &h, const int up_sample_coeff)
{
    real Fs_Resam = Fs * up_sample_coeff;
    unsigned short int num_taps_Resam = num_taps * up_sample_coeff;

    h.clear();
    h.resize(num_taps_Resam, 0.0);

    const real Fb = Fc - BW / 2.0;
    const real Fe = Fc + BW / 2.0;

    const real normCenter = ((Fb + Fe) / 2.0) / (Fs_Resam / 2.0);
    const real normPass = (Fe - Fb) / (Fs_Resam / 2.0);
    int middleindex = (num_taps_Resam - 1) / 2;

    for (int i = 0; i < num_taps_Resam; i++) {
        if (i == middleindex) {
            h[i] = normPass;
        } else {
            h[i] = normPass * (std::sin(PI * normPass / 2.0 * (i - middleindex))) / (PI * normPass / 2.0 * (i - middleindex));
        }
        h[i] *= std::cos(i * PI * normCenter) * std::pow(std::sin(i * PI / static_cast<real>(num_taps_Resam)), 2.0);
        h[i] *= up_sample_coeff * (0.5 - 0.5 * std::cos(2 * PI * i / (num_taps_Resam - 1)));
    }
}
void allPassFilter(std::vector<real> &y, std::vector<real> &x, std::vector<real> &state)
{
	// allocate memory for the impulse response
	y.clear();
	y = state;
	y.reserve(x.size());
	for(unsigned int n = 0; n<(x.size()-state.size());n++){
		y.push_back(x[n]);
	}

	std::vector<real> temp(&x[x.size()-state.size()], &x[x.size()]);
	state = temp;

}


void fmPll(const std::vector<real>& pllIn, real freq, real Fs, real& nextncoOut, std::vector<real>& ncoOut, real ncoScale, real phaseAdjust, real normBandwidth) {

    const real Cp = 2.666;
    const real Ci = 3.555;
    const real Kp = normBandwidth * Cp;
    const real Ki = normBandwidth * normBandwidth * Ci;
    ncoOut.resize(pllIn.size() + 1);

    real integrator = 0.0;
    real eastmatePhase = 0.0;
    real saveI = 1.0;
    real saveQ = 0.0;
    ncoOut[0] = nextncoOut;
    int trigOffset = 0;

    for (size_t k = 0; k < pllIn.size(); k++) {
    
        real errorI = pllIn[k] * saveI;  
        real errorQ = pllIn[k] * (-saveQ); 
        real errorD = std::atan2(errorQ, errorI);
        integrator += Ki * errorD;
        eastmatePhase += Kp * errorD + integrator;

        trigOffset++;
        real trigArg = 2 * M_PI * (freq / Fs) * trigOffset + eastmatePhase;
        saveI = std::cos(trigArg);
        saveQ = std::sin(trigArg);

        if (k < pllIn.size() - 1) {
            ncoOut[k + 1] = std::cos(trigArg * ncoScale + phaseAdjust);
        } else {
            nextncoOut = std::cos(trigArg * ncoScale + phaseAdjust);
        }
    }
}

void fmPllwithState(std::vector<real> &pllIn, real freq, real Fs, real ncoScale, real phaseAdjust, real normBandwidth, real &saveI, real &saveQ, real &integrator, real &eastmatePhase, real &trigOffset, std::vector<real> &ncoOut, real &ncoFirst)
{
    real Cp = 2.666;
    real Ci = 3.555;

    real Kp = normBandwidth * Cp;
    real Ki = normBandwidth * normBandwidth * Ci;

    ncoOut.clear();
    ncoOut.resize(pllIn.size(), 0.0);
    ncoOut[0] = ncoFirst;
    
    real trigArg, errorI, errorQ, errorD;

    for (int i = 0; i < pllIn.size(); i++)
    {
        errorI = pllIn[i] * (+saveI);
        errorQ = pllIn[i] * (-saveQ);
        errorD = std::atan2(errorQ, errorI);

        integrator = integrator + Ki * errorD;

        eastmatePhase = eastmatePhase + Kp * errorD + integrator;

        trigOffset += 1;
        trigArg = 2 * PI * (freq/Fs) * trigOffset + eastmatePhase;
        saveI = std::cos(trigArg);
        saveQ = std::sin(trigArg);

        if ((i + 1) == pllIn.size())
        {
            ncoFirst = std::cos(trigArg * ncoScale + phaseAdjust);
        }
        else
        {
            ncoOut[i + 1] = std::cos(trigArg * ncoScale + phaseAdjust);
        }
    }

}

void block_conv(std::vector<real> &y, const std::vector<real> &h, const std::vector<real> &x, std::vector<real> &state)
{
    int ntaps = h.size();
    int sig_len = x.size();
    y.clear();
    y.resize(sig_len, 0.0);

    for (int i = 0; i < y.size(); i++)
    {
        for (int k = 0; k < ntaps; k++)
        {
            int delt = i - k;
            if (delt >= 0)
            {
                y[i] += h[k] * x[delt];
            }
            else
            {
                y[i] += h[k] * state[state.size()+delt];
            }
        }
    }

    std::vector<real> temp(&x[x.size()-state.size()], &x[x.size()]);
    state = temp;
}

void impulseResponseRootRaisedCosine(std::vector<real> &y, real Fs, unsigned short int nTaps){
    real tSymbol = 1/2375.0;
    real beta = 0.90;

    for(int k = 0; k < nTaps-1; k++){
        real t = ((k-nTaps)/2)/Fs;
        if (t == 0.0){
            y.push_back(1 + beta*((4/PI)-1));
        } else if (t == -tSymbol/(4*beta) || t == tSymbol/(4*beta)){
            y.push_back((beta/std::pow(2,1/2))*(((1+2/PI)*(std::sin(PI/(4*beta)))) + ((1-2/PI)*(std::cos(PI/(4*beta))))));
        }else {
            y.push_back((std::sin(PI*t*(1-beta)/tSymbol) + 4*beta*(t/tSymbol)*std::cos(PI*t*(1+beta)/tSymbol))/(PI*t*(1-(4*beta*t/tSymbol)*(4*beta*t/tSymbol))/tSymbol));
        }
    }
}

