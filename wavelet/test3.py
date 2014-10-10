
"""
read real HRM samples from a csv and use a wavelet transform with overlapping windows
to extract heart rate data
"""

# import numpy as np
from numpy import *
from math import pi
from random import sample, random
import matplotlib.pyplot as plt
from pywt import dwt, Wavelet
from pylab import plot, show, title, xlabel, ylabel, subplot, xlim, ylim
from scipy import fft, arange
import csv

def mainLoop(sampleList):
	lastWindow = -1
	windowDiff = 5
	# windowLength = 20
	windowLength = 256
	windowCount = 0
	for sample in sampleList:
		time = sample[0]
		# if we're at the start of a window make an array to be processed
		if lastWindow == -1 or time-lastWindow > windowDiff:
			windowCount += 1
			lastWindow = time
			windowList = list()
			for sample in sampleList:
				# if sample[0] >= time and sample[0] < lastWindow + windowLength:
				if sample[0] >= time and len(windowList) < windowLength:
					windowList.append(sample)
			# if windowList[len(windowList)-1][0] - windowList[0][0] < windowLength - 1:
			if len(windowList) < windowLength:
				return
			windowTransform(windowList)
			# print 'in window ' + str(lastWindow) + ' to ' + str(windowList[len(windowList)-1][0]) + ', heart rate: ' + str(60*windowTransform(windowList))

def windowTransform(windowList):
	n = len(windowList)			# samples
	samp_rate = getRate(windowList)	# Hz
	# print "rate: " + str(samp_rate)
	levels = 8
	wave = getWave(windowList)

	# Can't detect frequencies that could not oscillate >= 2 times in the window
	# low_cutoff = 1 / ((len(windowList)/samp_rate)/2)
	low_cutoff = .4
	# print "low cutoff: " + str(low_cutoff)

	# cA, cD = dwt(wave, 'db2')
	coeffs = wavedec(wave, 'haar', level=levels)
	# coeffs = my_wavedec(wave, level=levels)

	i = 0
	# print len(coeffs)
	nyquist = samp_rate / 2.

	freqs = zeros(levels)
	freq_widths = zeros(levels)
	powers = zeros(levels)
	freq_range = nyquist - low_cutoff

	for a in coeffs:
		if (i <> 0):
			lower_limit = nyquist / 2.**(levels-i+1)
			upper_limit = nyquist / 2.**(levels-i)
			# print i, len(a), mean(abs(a)), std(a), lower_limit, upper_limit
			# I don't know why the 2 to 1 ratio works better...
			# freqs[i-1] = (2 * lower_limit + upper_limit) / 3.
			freqs[i-1] = (lower_limit + upper_limit) / 2.
			# freqs[i-1] = (lower_limit + 2 * upper_limit) / 3.

			# weight the frequencies by the widths of their bins
			# this could be a problem since the largest bin is 512 times the smallest
			# and later on the sum will always be skewed that way
			freq_widths[i-1] = (upper_limit - lower_limit) / nyquist
			# if ( i <> 1):
			# 	freq_widths[i-1] = freqs[i-1] - freqs[i-2]
			if (freqs[i-1] <= 8.0 and freqs[i-1] >= low_cutoff):
				powers[i-1] = mean(abs(a))
			else:
				powers[i-1] = 0
			 # powers[i-1] = mean(abs(a)) if (freqs[i-1] <= 8.0 and freqs[i=1] >= .5) else 0	# ignore frequencies higher than 4 Hz
		i += 1

	# find three biggest consecutive sums
	big_sum = 0.
	core_freq = 0.
	weight_factor = 0.

	# why is 6 so much better than 5?
	# probably still something to do with grabbing one of the lowest bins
	# with low frequency interference
	num_bins = 4
	start_bin = 0
	bins = zeros(levels)

	# use a weighted average of the num_bins largest consecutive bins.
	# create a weighted frequency by adding up the frequencies multiplied by the 'area'
	# of their respective bins, then divide to get the core_freq
	# also try weighting the frequencies in the middle of the set
	# for i in arange(levels-num_bins):
	# 	sum = 0.
	# 	freq_weight = 0.
	# 	sum_widths = 0.
	# 	for j in arange(num_bins):
	# 		weight_factor = getFactor(num_bins, j)
	# 		skew = freqSkew(num_bins, j)
	# 		freq_weight += freqs[i+j] * powers[i+j]  * weight_factor / freq_widths[i+j]# * skew
	# 		sum_widths += powers[i+j] * weight_factor / freq_widths[i+j] 
	# 		sum +=  powers[i+j] * weight_factor # freq_widths[i+j]
	# 	if (sum > big_sum):
	# 		start_bin = i
	# 		big_sum = sum
	# 		core_freq = freq_weight / sum_widths

	# alternatively find the centroid of the curve as if it was a shape
	numerator = 0.
	denominator = 0.
	for i in arange(levels):
		if (i <> 0):
			numerator += (((freqs[i]-freqs[i-1])**2)/6)*(2*powers[i] + powers[i-1]) * freq_widths[i]
			denominator += ((freqs[i]-freqs[i-1])/2)*(powers[i] + powers[i-1]) * freq_widths[i]
	core_freq = numerator/denominator

	# for i in arange(len(freqs)):
	# 	bins[i] = powers[i] / freq_widths[i]#freqs[i] * 60 * powers[i]
	# 	if powers[i] == 0:
	# 		bins[i] = 0
	# print "here da bins " + str(bins) + ", start bin " + str(start_bin)

	print 'in window ' + str(windowList[0][0]) + ' to ' + str(windowList[len(windowList)-1][0]) + ', heart rate: ' + str(60*core_freq)
	# plotWindow(windowList, freqs, powers)
	
	# print 'averaging in the range ' + str(nyquist / 2.**(levels-first_freq + num_bins)) + ' to ' + str(nyquist / 2.**(levels-first_freq))

	# if 60*core_freq < 50:
	# print 'in window ' + str(windowList[0][0]) + ', weird heart rate: ' + str(60*core_freq)
	return core_freq	

def getRate(windowList):
	numSamples = len(windowList)
	timeElapsed = windowList[numSamples-1][0]-windowList[0][0]
	rate = numSamples/timeElapsed
	return rate

def getWave(windowList):
	wave = list()
	for sample in windowList:
		wave.append(sample[1])
	return wave

def getTime(windowList):
	time = list()
	for sample in windowList:
		time.append(sample[0])
	return time

def getFactor(n, i):
	# factor = (n+1)/2 + abs((n-1)/2 - i)
	factor = 1 + abs((n-1)/2 - i)
	scale = 1
	return 1/float(scale*factor)

def freqSkew(n, i):
	scale = .25
	if(i < (n-1)/2):
		factor = 1 + scale*abs((n-1)/2 - i)
	else:
		factor = 1 - scale*abs((n-1)/2 - i)
	return factor


def plotWindow(windowList, freqs, powers):

	# freqs[:] = [f * 60 for f in freqs]
	# plt.plot(freqs, powers)
	# plt.xlim(0, 400)
	# plt.grid(True)
	# plt.show()
	subplot(3,1,1)
	plot(getTime(windowList),getWave(windowList))
	xlabel('Time')
	ylabel('Amplitude')
	subplot(3,1,2)
	plotFFTSpectrum(getWave(windowList), getRate(windowList))
	subplot(3,1,3)
	xlim(0, 4)
	ylim(0, 3)
	plot(freqs, powers)
	ylabel('Power')

	show()
	# 	print 'in window ' + str(windowList[0][0]) + ', weird heart rate: ' + str(60*core_freq)
	# 	i = 0
	# 	for a in coeffs:
	# 		if (i <> 0):
	# 			print i, len(a), mean(abs(a)), std(a), nyquist / 2.**(levels-i+1), nyquist / 2.**(levels-i)
	# 		i += 1
	# 	plt.plot(getTime(windowList), getWave(windowList))
	# 	plt.show()
	# freqs[:] = [f * 60 for f in freqs]
	# plt.plot(freqs, powers)
	# plt.xlim(0, 400)
	# plt.grid(True)
	# plt.show()

def plotFFTSpectrum(y,Fs):
	"""
	Plots a Single-Sided Amplitude Spectrum of y(t)
	"""
	n = len(y) # length of the signal
	k = arange(n)
	T = n/Fs
	frq = k/T # two sides frequency range
	frq = frq[range(n/2)] # one side frequency range

	Y = fft(y)/n # fft computing and normalization
	Y = Y[range(n/2)]

	plot(frq,abs(Y),'r') # plotting the spectrum
	xlim(0, 4)
	# set_xticks(arange(0, 5, .25))
	xlabel('Freq (Hz)')
	ylabel('|Y(freq)|')


def readFile():
	sampleList = list()
	with open('samples_finger_02.csv','rU') as samplesFile:
		reader = csv.reader(samplesFile, dialect=csv.excel_tab)
		for row in reader:
			valuePair = row[0].split(",", 1)
			sampleList.append(map(float, valuePair))
	return sampleList

def waveGen():
	n = 4096			# samples
	freq0 = 0 	# Hz
	samp_rate = 64	# Hz
	levels = 8

	start_freq = 1	# Hz
	end_freq = 2	# Hz
	if (start_freq != end_freq):
		freq0 = arange(start_freq, end_freq, (end_freq - start_freq) / (n * 1.0))
	else:
		freq0 = start_freq


	factor0 = samp_rate / freq0
	time = arange(n)/float(samp_rate)
	wave0 = sin(2 * pi * freq0 * time)

	# errors = [random() - 0.5 for _ in range(n)]
	# wave0 += errors

	sampleList = list()
	for t in arange(len(time)):
		sample = [time[t], wave0[t]]
		sampleList.append(sample)

	return sampleList

def my_wavedec(data, level=None):
	coeffs_list = []
	length = len(data) >> 1
	a = array(zeros(length), float) #list(zeros(length))
	d = array(zeros(length), float) #list(zeros(length))
	for j in xrange(level):
		a = array(zeros(length), float) #list(zeros(length))
		d = array(zeros(length), float) #list(zeros(length))
		for i in xrange(length):
			total = data[i*2] + data[i*2+1]
			diff = data[i*2] - data[i*2+1]
			a[i] = total
			d[i] = diff

		data = a
		coeffs_list.append(d)

		# print "length = " + str(length) + " a: " + str(a) + ", d: " + str(d)
		length = length >> 1
	coeffs_list.append(a)
	coeffs_list.reverse()
	return coeffs_list

def wavedec(data, wavelet, mode='sym', level=None):
    """
    Multilevel 1D Discrete Wavelet Transform of data.
    Returns coefficients list - [cAn, cDn, cDn-1, ..., cD2, cD1]

    data    - input data
    wavelet - wavelet to use (Wavelet object or name string)
    mode    - signal extension mode, see MODES
    level   - decomposition level. If level is None then it will be
              calculated using `dwt_max_level` function.
    """

    if not isinstance(wavelet, Wavelet):
        wavelet = Wavelet(wavelet)

    # if level is None:
    #     level = dwt_max_level(len(data), wavelet.dec_len)
    # elif level < 0:
    #     raise ValueError(
    #         "Level value of %d is too low . Minimum level is 0." % level)

    coeffs_list = []

    a = data
    for i in xrange(level):
        a, d = dwt(a, wavelet, mode)
        # a, d = my_dwt(a)
        # print "length of a: " + str(len(a))
        coeffs_list.append(d)

    # print "type: " + str(type(a))
    coeffs_list.append(a)
    coeffs_list.reverse()

    return coeffs_list

def my_dwt(data):
	length = len(data) >> 1;
	a = ndarray() #list(zeros(length))
	d = ndarray() #list(zeros(length))
	# while(length >= 1):
	# 	for i in xrange(length):
	# 		total = data[i * 2] + data[i * 2 + 1]
	# 		difference = data[i * 2] - data[i * 2 + 1]
	# 		coeffs[i] = total
	# 		coeffs[length + i] = difference

	# 	if(length == 1):
	# 		return coeffs

	# 	for j in xrange(length << 1):
	# 		data[i] = coeffs[i]
	# 	length >>= 1

	for i in xrange(length):
		total = data[i * 2] + data[i * 2 + 1]
		difference = data[i * 2] - data[i * 2 + 1]
		append(a, total) #a.append(total)#a[i] = total
		append(d, difference) #d.append(difference)#d[i] = difference
	return a, d

if __name__ == '__main__':
	sampleList = readFile()
	# sampleList = waveGen()
	mainLoop(sampleList)
	plt.plot(getTime(sampleList), getWave(sampleList))
	plt.show()
