

import matplotlib.pyplot as plt
import numpy as np
import math

def fmPll(pllIn, freq, Fs, nextncoOut, ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01):

	"""

	pllIn 	 		array of floats
					input signal to the PLL (assume known frequency)

	freq 			float
					reference frequency to which the PLL locks

	Fs  			float
					sampling rate for the input/output signals

	ncoScale		float
					frequency scale factor for the NCO output

	phaseAdjust		float
					phase adjust to be added to the NCO output only

	normBandwidth	float
					normalized bandwidth for the loop filter
					(relative to the sampling rate)

	state 			to be added

	"""

	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp

	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output array for the NCO
	ncoOut = np.empty(len(pllIn)+1)

	# initialize internal state
	integrator = 0.0
	eastmatePhase = 0.0
	saveI = 1.0
	saveQ = 0.0
	ncoOut[0] = nextncoOut
	trigOffset = 0

	# note: state saving will be needed for block processing
	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+saveI)  # complex conjugate of the
		errorQ = pllIn[k] * (-saveQ)  # feedback complex exponential

		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)

		# loop filter
		integrator = integrator + Ki*errorD

		# update phase estimate
		eastmatePhase = eastmatePhase + Kp*errorD + integrator

		# internal oscillator
		trigOffset += 1
		trigArg = 2*math.pi*(freq/Fs)*(trigOffset) + eastmatePhase
		saveI = math.cos(trigArg)
		saveQ = math.sin(trigArg)

		# for stereo only the in-phase NCO component should be returned
		# for block processing you should also return the state
		if (k < len(pllIn)-1):
			ncoOut[k+1] = math.cos(trigArg * ncoScale + phaseAdjust)
		else:
			nextncoOut = math.cos(trigArg * ncoScale + phaseAdjust)

	return ncoOut, nextncoOut

	# for RDS add also the quadrature NCO component to the output


def fmPllRds(pllIn, freq, Fs, state, ncoScale, phaseAdjust, normBandwidth):

	"""

	pllIn 	 		array of floats
					input signal to the PLL (assume known frequency)

	freq 			float
					reference frequency to which the PLL locks

	Fs  			float
					sampling rate for the input/output signals

	ncoScale		float
					frequency scale factor for the NCO output

	phaseAdjust		float
					phase adjust to be added to the NCO output only

	normBandwidth	float
					normalized bandwidth for the loop filter
					(relative to the sampling rate)

	state 			to be added

	"""

	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp

	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output array for the NCO for RDS boht In-phase and Quadrature are required
	ncoOutInPhase = np.empty(len(pllIn)+1)
	ncoOutQuadrature = np.empty(len(pllIn)+1)
	# initialize internal state
	integrator = state[0]
	eastimatePhase = state[1]
	saveI = state[2]
	saveQ = state[3]
	ncoOutInPhase[0] = state[4]
	ncoOutQuadrature[0] = state[5]
	trigOffset = state[6]

	# note: state saving will be needed for block processing
	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+saveI)  # complex conjugate of the
		errorQ = pllIn[k] * (-saveQ)  # feedback complex exponential

		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)

		# loop filter
		integrator = integrator + Ki*errorD

		# update phase estimate
		eastimatePhase = eastimatePhase + Kp*errorD + integrator

		# internal oscillator
		trigOffset += 1
		trigArg = 2*math.pi*(freq/Fs)*(trigOffset) + eastimatePhase
		saveI = math.cos(trigArg)
		saveQ = math.sin(trigArg)

		# for stereo only the in-phase NCO component should be returned
		# for block processing you should also return the state
		if (k < len(pllIn)-1):
			ncoOutInPhase[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
			ncoOutQuadrature[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)
		else:
			nextncoOutIn = math.cos(trigArg * ncoScale + phaseAdjust)
			nextncoOutQuad = math.sin(trigArg * ncoScale + phaseAdjust)
		
	state[0] = integrator
	state[1] = eastimatePhase
	state[2] = saveI
	state[3] = saveQ
	state[4] = ncoOutInPhase[-1]
	state[5] = ncoOutQuadrature[-1]
	state[6] = trigOffset


#	fig, (p_adjust1) = plt.subplots(nrows=1)
#	p_adjust1.set_ylim(-1.25, 1.25)
#	plt.plot(100*pllIn[0:512], c='r')
#	plt.plot(ncoOutInPhase[0:512], c='b')
#	fig.subplots_adjust(hspace = 1.0)
#	plt.show()
#
	return ncoOutInPhase, ncoOutQuadrature, state


if __name__ == "__main__":

	pass
