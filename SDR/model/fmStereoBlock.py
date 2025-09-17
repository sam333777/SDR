

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD
from fmPll import fmPll
# for take-home add your functions

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 101
rf_decim = 10

intermediate_Fs = 240e3

audio_Fs = 48e3
audio_decim = 5
audio_taps = 101
audio_Fc = 16e3
# add other settings for audio, like filter taps, ...

def my_firwin(n_taps, fc, fs):
    norm_cutoff = fc / (fs / 2)
    middle_index = (n_taps - 1) / 2
    taps = np.zeros(n_taps)

    for i in range(n_taps):
        if i == middle_index:
            taps[i] = norm_cutoff
        else:
            taps[i] = (
                norm_cutoff
                * np.sin(np.pi * norm_cutoff * (i - middle_index))
                / (np.pi * norm_cutoff * (i - middle_index))
            )
        # Apply Hann window
        taps[i] *= 0.5 - 0.5 * np.cos(2 * np.pi * i / (n_taps - 1))

    return taps

def my_block_process(block, impulse_response, state):
    block_filtered = np.zeros(len(block))
    for i in range(len(block)):
        for j in range(len(impulse_response)):
            if i - j >= 0:
                block_filtered[i] += block[i-j] * impulse_response[j]
            else:
                #use the state from the previous block
                block_filtered[i] += state[i - j + len(state)] * impulse_response[j]
    #update the state
    updated_state = np.concatenate((state[len(block):], block))[-len(state):]
    return block_filtered, updated_state

def delayBlock(input_block, state_block):
	output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
	state_block = input_block[-len(state_block):]

	return output_block, state_block

def fastArctan(I, Q, prev_phase = 0.0):
	#
	# the default prev_phase phase is assumed to be zero, however
	# take note in block processing it must be explicitly controlled

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))
	
	for k in range(len(I)):
		if k > 0:
			dQ = Q[k] - Q[k - 1] 
			dI = I[k] - I[k - 1]

		else: 
			dQ = 0
			dI = 0

		num = (I[k] * dQ) - (Q[k] * dI)
		den = (I[k]**2) + (Q[k]**2)
		if den != 0:
			temp = num/den
		else:
			temp = 0

		current_phase = prev_phase + temp

		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		[prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		fm_demod[k] = temp

		# save the state of the current phase
		# to compute the next derivative
		prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/samples8.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	'''
	# IQ data is normalized between -1 and +1 in 64-bit double format
	iq_data = (np.float64(raw_data) - 128.0) / 128.0
	print("Reformatted raw RF data to 64-bit double format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")
	'''

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc / (rf_Fs / 2), window=('hann'))

	# coefficients for the filter to extract mono audio
	# to be updated by you during the in-lab session based on firwin
	# same principle as for rf_coeff (but different arguments, of course)
	audio_coeff = signal.firwin(audio_taps, audio_Fc / (intermediate_Fs / 2), window=('hann'))
	stereo_pilot_coeff = signal.firwin(audio_taps, [(18.5e3 / (intermediate_Fs / 2)), (19.5e3 / (intermediate_Fs / 2)) ], window=('hann'), pass_zero=False)
	stereo_audio_coeff = signal.firwin(audio_taps, [(22e3 / (intermediate_Fs / 2)), (54e3 / (intermediate_Fs / 2)) ], window=('hann'), pass_zero=False)

	# to be updated by you for the take-home exercise
	# with your own code for impulse response generation
	# audio_coeff = my_firwin(audio_taps, audio_Fc, intermediate_Fs)
        

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = (int)(rf_Fs/25)
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps - 1)
	state_q_lpf_100k = np.zeros(rf_taps - 1)
	state_phase = 0
	# add state as needed for the mono channel filter
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)

	# audio buffer that stores all the audio blocks
	audio_filt_state = np.zeros(audio_taps - 1)
	stereo_filt_state = np.zeros(audio_taps - 1)
	ncoOutState = 0

	stereo_pilot_state = np.zeros(audio_taps - 1)
	stereo_audio_state = np.zeros(audio_taps - 1)
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

		fm_demod, state_phase = fastArctan(i_ds, q_ds, state_phase)
		#Stereo Pilot
		stereo_pilot_filt, stereo_pilot_state = signal.lfilter(stereo_pilot_coeff, 1.0, fm_demod, zi = stereo_pilot_state)
		#Stereo Audio
		stereo_audio_filt, stereo_audio_state = signal.lfilter(stereo_audio_coeff, 1.0, fm_demod, zi = stereo_audio_state)

		#All-Pass Filter
		audio_filt_state_delayed = np.zeros(len(fm_demod))
		audio_filt_delayed, audio_filt_state_delayed = delayBlock(fm_demod, audio_filt_state_delayed)

		#Run through PLL
		ncoOut, ncoOutState = fmPll(stereo_pilot_filt, 19e3, intermediate_Fs, ncoOutState)
			
		#Mixer Section
		for i in range(len(ncoOut)-1):
			stereo_audio_filt[i] = 2 * ncoOut[i] * stereo_audio_filt[i]

		#final stereo filter
		stereo_audio_filt, stereo_filt_state = signal.lfilter(audio_coeff, 1.0, fm_demod, zi = stereo_filt_state)

		#Mono Filtering
		audio_filt, audio_filt_state = signal.lfilter(audio_coeff, 1.0, fm_demod, zi = audio_filt_state)

		#Signal Processing
		mono_audio_block = audio_filt[::audio_decim]
		stereo_audio_block = stereo_audio_filt[::audio_decim]
			
		stereo_left = np.zeros(len(mono_audio_block) - 1)
		stereo_right = np.zeros(len(mono_audio_block) - 1)
			
		#Seperate the channels
		for i in range(len(mono_audio_block)-1):
			stereo_left[i] = stereo_audio_block[i] + mono_audio_block[i]
			stereo_right[i] = stereo_audio_block[i] - mono_audio_block[i]
		
		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data

		for i in range(len(stereo_left)-1):
			audio_data = np.concatenate((audio_data, [stereo_left[i]]))
			audio_data = np.concatenate((audio_data, [stereo_right[i]]))
		
		print(len(audio_data))

		# to save runtime, select the range of blocks to log data
		# this includes both saving binary files and plotting PSD
		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			# (for easier visualization purposes we divide Fs by 1e3 to imply the kHz units on the x-axis)
			# (this scales the y axis of the PSD, but not the relative strength of different frequencies)
			ax0.clear()
			fmPlotPSD(ax0, fm_demod, (rf_Fs / rf_decim) / 1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			'''
			# create binary file where each sample is a 64-bit double
			fm_demod.astype('float64').tofile(fm_demod_fname)
			'''

			fmPlotPSD(ax1, audio_filt, (rf_Fs / rf_decim) / 1e3, subfig_height[1], 'Extracted Mono')

			fmPlotPSD(ax2, audio_data, audio_Fs / 1e3, subfig_height[2], 'Downsampled Mono Audio')

			# save figure to file
			fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmSBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data / 2) * 32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	# plt.show()
