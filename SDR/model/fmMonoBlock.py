

"""
For the project release, we recommend working in 64-bit double format, as the
memory overhead is acceptable. Therefore, the normalization of 8-bit unsigned
I/Q samples has been performed in 64-bit double format.

By default, the block size must match the amount of I/Q data acquired in 40 ms
"""

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

from fmSupportLib import fmDemodArctan, fmPlotPSD

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_decim = 10
rf_taps = 101

audio_Fs = 48e3
audio_decim = 5
audio_up_decim = 1
audio_Fc = 16e3
audio_taps = 101

def filterBlock(h, x, audio_state, audio_decim = 1, audio_up_decim = 1):
    y = np.zeros(int(len(x)*audio_up_decim/audio_decim))
    print(len(h)//audio_up_decim)
    for n in range(len(y)):
        for k in range(len(h)//audio_up_decim):
            if (n*audio_decim-(k*audio_up_decim+n*audio_decim%audio_up_decim))/audio_up_decim >= 0:
                y[n] += h[k*audio_up_decim+(n*audio_decim)%audio_up_decim] * x[(n*audio_decim-(k*audio_up_decim+(n*audio_decim)%audio_up_decim))//audio_up_decim]
            else:
                y[n] += h[k*audio_up_decim+(n*audio_decim)%audio_up_decim] * audio_state[(n*audio_decim-(k*audio_up_decim+(n*audio_decim)%audio_up_decim))//audio_up_decim]
    audio_state = x[-(len(h)-1):]
    return y, audio_state
	
def Demod(I, Q, prev_I, prev_Q):
	fm_demod = np.empty(len(I))
	for k in range(len(I)):
		if ((I[k]**2+Q[k]**2)==0):
			fm_demod[k] = 0.0
		else:
			fm_demod[k] = (I[k]*(Q[k]-prev_Q) - Q[k]*(I[k]-prev_I))/(I[k]**2+Q[k]**2)
		prev_I = I[k]
		prev_Q = Q[k]

	return fm_demod, prev_I, prev_Q

def LowPassFilter(Fc, Fs, Ntaps):
	h = np.zeros(Ntaps)
	norm_f = Fc/(Fs/2)
	for i in range(Ntaps):
		if (i == (Ntaps -1)/2):
			h[i] = norm_f
		else:
			h[i] = norm_f*np.sin(np.pi*norm_f*(i-(Ntaps-1)/2))/(np.pi*norm_f*(i-(Ntaps-1)/2))
		h[i] = h[i] *(np.sin(i*np.pi/Ntaps))**2
	return h

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/samples1.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	'''
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")
	'''

	# IQ data is normalized between -1 and +1 in 64-bit double format
	iq_data = (np.float64(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 64-bit double format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	rf_coeff = LowPassFilter(rf_Fc, rf_Fs, rf_taps)
	audio_coeff = audio_up_decim*LowPassFilter(audio_Fc, rf_Fs/rf_decim*audio_up_decim, audio_taps*audio_up_decim)

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# define block_size 
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(audio_taps - 1)
	state_q_lpf_100k = np.zeros(audio_taps - 1)
	state_audio = np.zeros(audio_taps-1)

	state_phase = 0
	
	# add state as needed for the mono channel filter
	state_I = 0
	state_Q = 0
	# audio buffer that stores all the audio blocks

	audio_data = np.array([])

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count + 1) * block_size < len(iq_data):

		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[block_count * block_size:(block_count + 1) * block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[block_count * block_size + 1:(block_count + 1) * block_size:2],
				zi=state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		#fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
		fm_demod, state_I, state_Q = Demod(i_ds, q_ds, state_I, state_Q)
		audio_filtered, state_audio = filterBlock(audio_coeff, fm_demod, state_audio, audio_decim, audio_up_decim)
		audio_block = audio_filtered

		audio_data = np.concatenate((audio_data, audio_block))

		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			# (for easier visualization purposes we divide Fs by 1e3 to imply the kHz units on the x-axis)
			# (this scales the y axis of the PSD, but not the relative strength of different frequencies)
			ax0.clear()
			fmPlotPSD(ax0, fm_demod, (rf_Fs / rf_decim) / 1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			'''
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)
			'''

			# create binary file where each sample is a 64-bit double
			fm_demod.astype('float64').tofile(fm_demod_fname)

			# save figure to file
			fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmMonoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data / 2) * 32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	# plt.show()
 
