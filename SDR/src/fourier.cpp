
// Source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"
#include <cmath>

// Just DFT function (no FFT)
void DFT(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {
	Xf.clear();
	Xf.resize(x.size(), std::complex<real>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
			std::complex<real> expval(0, -2 * PI * (k * m) / x.size());
			Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// Function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<real>> &Xf, std::vector<real> &Xmag)
{
	Xmag.clear();
	Xmag.resize(Xf.size(), real(0));
	for (int i = 0; i < (int)Xf.size(); i++) {
		Xmag[i] = std::abs(Xf[i]) / Xf.size();
	}
}

// Add your own code to estimate the PSD
// ...

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void DFT_reference(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	Xf.clear();
	Xf.resize(x.size(), std::complex<real>(0));
	for (int m = 0; m < (int)Xf.size(); m++) {
		for (int k = 0; k < (int)x.size(); k++) {
			std::complex<real> expval(0, -2 * M_PI * (k * m) / x.size());
			Xf[m] +=  + x[k] * std::exp(expval);
		}
	}
}

void DFT_init_bins(const std::vector<real> &x, std::vector<std::complex<real>> &Xf) {

	int N = (int)x.size();
	std::fill(Xf.begin(), Xf.end(), std::complex<real>(0., 0.));
	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
			std::complex<real> expval(0, -2 * M_PI * (k * m) / N);
			Xf[m] += x[k] * std::exp(expval);
		}
	}
}

void generate_DFT_twiddles(const int& N, std::vector<std::complex<real>> &Twiddle1D) {

	Twiddle1D.resize(N);
	for (int k = 0; k < N; k++) {
		std::complex<real> expval(0, -2 * M_PI * k / N);
		Twiddle1D[k] = std::exp(expval);
	}
}

void generate_DFT_matrix(const int& N, std::vector<std::vector<std::complex<real>>> &Twiddle2D) {

	Twiddle2D.resize(N, std::vector<std::complex<real>>(N));
    std::vector<std::complex<real>> Twiddle1D;
	generate_DFT_twiddles(N, Twiddle1D);

	for (int m = 0; m < N; m++) {
		for (int k = 0; k < N; k++) {
			Twiddle2D[m][k] = Twiddle1D[(k * m) % N];
		}
	}
}

void matrixPSD_precomputes(int nfft, std::vector<real> &hann_window, std::vector<std::vector<std::complex<real>>> &dft_matrix)
{
    int freq_bins = nfft;
    dft_matrix.resize(freq_bins, std::vector<std::complex<real>>(freq_bins, {0.0, 0.0}));
    hann_window.resize(freq_bins, 0.0);
    //dft_matrix.resize(freq_bins, std::vector<std::complex<real>>(freq_bins, {0.0, 0.0}));
    const std::complex<real> j(0.0, 1.0);

    for(int m = 0; m < freq_bins; m++){
        for(int k = 0; k < freq_bins; k++){
            real exponent = 2.0 * M_PI * (-static_cast<real>(k) * m) / freq_bins;
            dft_matrix[m][k] = std::exp(j * exponent);
        }
        hann_window[m] = 0.5 * (1.0 - std::cos(2.0 * M_PI * m / (freq_bins - 1)));
    }
}

void matrixPSD_optimized_precompute(const std::vector<real> &samples, std::vector<real> &psd_est, std::vector<real> &freq, int nfft, real Fs, std::vector<real> &hann_window, std::vector<std::vector<std::complex<real>>> &dft_matrix) {

	
    int freq_bins = nfft;
    real df = Fs / static_cast<real>(freq_bins);
    freq.clear();
    for(int i = 0; i < freq_bins / 2; i++) {
        freq.push_back(i * df);
    }
    int no_segments = static_cast<int>(samples.size()) / freq_bins;

	std::vector<std::vector<std::complex<real>>> Xf(no_segments, std::vector<std::complex<real>>(freq_bins, {0.0, 0.0}));

    // Applying Hann window
    for(int seg = 0; seg < no_segments; seg++){
        for(int m = 0; m < freq_bins; m++){
            // Adding samples * hann_window * dft_matrix
            for(int k = 0; k < freq_bins; k++){
                int index = seg * freq_bins + k;
                real val = (index < static_cast<int>(samples.size())) ? samples[index] : 0.0;
                Xf[seg][m] += val * hann_window[k] * dft_matrix[m][k];
            }
        }
    }


    psd_est.resize(freq_bins / 2, 0.0);
    for(int m = 0; m < freq_bins / 2; m++){
        real sum_power = 0.0;
        for(int seg = 0; seg < no_segments; seg++){
            real pwr = std::norm(Xf[seg][m]);
            // Normalization as in the Python implementation
            sum_power += (1.0 / ((Fs / 2.0) * (freq_bins / 2.0))) * pwr;
        }
        // Average PSD for each frequency bin across all segments, convert to dB
        real avg_pwr = sum_power / static_cast<real>(no_segments);
        psd_est[m] = 10.0 * std::log10(avg_pwr);
    }
}

void estimatePSD(std::vector<real> &samples, std::vector<real> &psd_result, std::vector<real> &freq_result, unsigned int nfft, unsigned int Fs){

    int freq_bins = nfft;
    real df = (Fs * 1000) / freq_bins;
    std::vector<std::complex<real>> Xf(freq_bins);
    
    // Create the frequency vector with only positive frequencies
    //numpy.arange([start, ]stop, [step, ], dtype=None) -> numpy.ndarray

    for (int i = 0; i < freq_bins; i++) {
        freq_result.push_back((i * df) / 1000);
    }

    // Design the Hann window to smooth the data
    std::vector<real> hann(freq_bins);
    for (int i = 0; i < freq_bins; i++) {
        hann[i] = 0.5 * (1 - std::cos(2 * PI * i / (freq_bins - 1)));
    }

    int no_segments = static_cast<int>(samples.size() / freq_bins);
    std::vector<real> psd_list;
    //psd_list.reserve(no_segments * (freq_bins / 2));

    // Compute PSD for each segment
    for (int k = 0; k < no_segments; k++) {
        std::vector<real> windowed_samples(freq_bins);
        for (int i = 0; i < freq_bins; i++) {
            windowed_samples[i] = samples[k * freq_bins + i] * hann[i];
        }
        
        Xf.clear();
        DFT(windowed_samples, Xf);

        Xf.resize(freq_bins / 2);

        for (int i = 0; i < freq_bins / 2; i++) {
            real psd_val = 2 * (1.0f / (Fs * (freq_bins / 2.0f))) * abs(pow(Xf[i],2));
            psd_list.push_back(psd_val);
        }
    }

    // Average PSD for each frequency bin across all segments
    std::vector<real> psd_avg(freq_bins / 2, 0.0f);
    for (int k = 0; k < freq_bins / 2; k++) {
        for (int l = 0; l < no_segments; l++) {
            psd_avg[k] += psd_list[k + l * (freq_bins / 2)];
        }
        psd_avg[k] /= no_segments;
    }

    // Convert to dB
    psd_result.resize(freq_bins / 2, 0);
    for (int k = 0; k < freq_bins / 2; k++) {
        psd_result[k] = 10.0f * std::log10(psd_avg[k]);
    }
}

void lightweight_fmDemodArctan(std::vector<real> &fm_demod, real &prev_phase, const std::vector<real> &I, const std::vector<real> &Q){
	//std::cerr << I.size() << " ";
	//std::cerr << fm_demod.size() << " ";
	real temp;
	real currentPhase;
	for(int i = 1; i < I.size(); i++){
	    real num = (I[i] * (Q[i] - Q[i-1]) - Q[i] * (I[i] - I[i-1]));
	    real den = (pow(I[i], 2) + pow(Q[i],2));
	    if (den != 0){
		    temp = num/den;
	    }else{
		    temp = 0;
	    }
	    currentPhase = prev_phase + temp;
	    fm_demod[i] = temp;
	    prev_phase = currentPhase;
	    //std::cerr << fm_demod[i] << " " << i << " ";
	}

}