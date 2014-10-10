
"""
read real HRM samples from a csv and use a wavelet transform with overlapping windows
to extract heart rate data
"""

# import numpy as np
from numpy import *
from math import pi
from random import sample, random
import matplotlib.pyplot as plt
from pylab import plot, show, title, xlabel, ylabel, subplot, xlim, ylim
from scipy import fft, arange
import csv
import sys

def mainLoop(sampleList):
	lastWindow = -1
	windowDiff = 64
	# windowLength = 20
	windowLength = 256
	windowCount = 0
	index = 0
	for sample in sampleList:
		time = sample[0]
		# if we're at the start of a window make an array to be processed
		# if lastWindow == -1 or time-lastWindow > windowDiff:
		if lastWindow == -1 or index-lastWindow > windowDiff:
			windowCount += 1
			# lastWindow = time
			lastWindow = index
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
		index += 1

def windowTransform(windowList):
	n = len(windowList)			# samples
	samp_rate = 64 # getRate(windowList)	# Hz
	print "rate: " + str(samp_rate)
	levels = 7
	wave = getWave(windowList)

	# Can't detect frequencies that could not oscillate >= 2 times in the window
	# low_cutoff = 1 / ((len(windowList)/samp_rate)/2)
	low_cutoff = .4
	# print "low cutoff: " + str(low_cutoff)

	# cA, cD = dwt(wave, 'db2')
	coeffs = my_wavedec(wave, level=levels)

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
			# freq_widths[i-1] = (upper_limit - lower_limit) / nyquist
			if ( i <> 1):
				freq_widths[i-1] = freqs[i-1] - freqs[i-2]
			if (freqs[i-1] <= 8.0 and freqs[i-1] >= low_cutoff):
				powers[i-1] = mean(abs(a))
			else:
				powers[i-1] = 0
			print "power " + str(powers[i-1])
			 # powers[i-1] = mean(abs(a)) if (freqs[i-1] <= 8.0 and freqs[i=1] >= .5) else 0	# ignore frequencies higher than 4 Hz
		i += 1

	# find three biggest consecutive sums
	core_freq = 0.
	numerator = 0.
	denominator = 0.

	for i in arange(levels):
		if (i <> 0):
			numerator += (((freqs[i]-freqs[i-1])**2)/6)*(2*powers[i] + powers[i-1]) * freq_widths[i]
			denominator += ((freqs[i]-freqs[i-1])/2)*(powers[i] + powers[i-1]) * freq_widths[i]
			print "N: " + str(numerator) + "\nD: " + str(denominator)
	core_freq = numerator/denominator

	print 'in window ' + str(windowList[0][0]) + ' to ' + str(windowList[len(windowList)-1][0]) + ', heart rate: ' + str(60*core_freq)
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

def printTransformArray(sampleList):
	i = 0
	transformArray = []
	sys.stdout.write("{")
	for sample in sampleList:
		if (i <> 256 and i <> 0):
			sys.stdout.write(", ")
		sys.stdout.write(str(sample[1]))
		transformArray.append(sample[1])
		i += 1 
		if (i == 256):
			sys.stdout.write("}; \n {")
		if(i == 326):
			sys.stdout.write("};")
			return transformArray


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

def rawData():
	# sampleList = list();
	sampleList = [[22.95, 64], [22.95, 74], [22.97, 101], [22.99, 140], [23.00, 187], [23.02, 224], [23.03, 262], [23.05, 315], [23.06, 376], [23.08, 428], [23.09, 482], [23.11, 538], [23.13, 596], [23.14, 646], [23.16, 693], [23.17, 757], [23.19, 760], [23.20, 760], [23.22, 760], [23.24, 760], [23.25, 760], [23.27, 760], [23.28, 760], [23.30, 760], [23.31, 760], [23.33, 767], [23.34, 759], [23.36, 758], [23.38, 671], [23.39, 525], [23.41, 381], [23.42, 283], [23.44, 214], [23.45, 167], [23.47, 134], [23.49, 114], [23.50, 101], [23.52, 95], [23.53, 95], [23.55, 94], [23.56, 89], [23.58, 87], [23.59, 89], [23.61, 90], [23.63, 87], [23.64, 84], [23.66, 75], [23.67, 63], [23.69, 54], [23.70, 48], [23.72, 42], [23.74, 39], [23.75, 36], [23.77, 33], [23.78, 32], [23.80, 31], [23.81, 32], [23.83, 32], [23.84, 255], [23.86, 40], [23.88, 45], [23.89, 57], [23.91, 73], [23.92, 93], [23.94, 127], [23.95, 172], [23.97, 224], [23.99, 281], [24.00, 336], [24.02, 393], [24.03, 448], [24.05, 501], [24.06, 551], [24.08, 592], [24.09, 636], [24.11, 680], [24.13, 710], [24.14, 724], [24.16, 701], [24.17, 628], [24.19, 505], [24.20, 367], [24.22, 265], [24.24, 194], [24.25, 144], [24.27, 109], [24.28, 83], [24.30, 65], [24.31, 52], [24.33, 43], [24.34, 37], [24.36, 32], [24.38, 29], [24.39, 26], [24.41, 25], [24.42, 23], [24.44, 22], [24.45, 22], [24.47, 21], [24.48, 21], [24.50, 21], [24.52, 21], [24.53, 21], [24.55, 21], [24.56, 20], [24.58, 20], [24.59, 20], [24.61, 20], [24.63, 20], [24.64, 20], [24.66, 21], [24.67, 21], [24.69, 23], [24.70, 24], [24.72, 29], [24.73, 255], [24.75, 80], [24.77, 122], [24.78, 172], [24.80, 228], [24.81, 287], [24.83, 349], [24.84, 413], [24.86, 468], [24.88, 520], [24.89, 573], [24.91, 622], [24.92, 677], [24.94, 732], [24.95, 760], [24.97, 760], [24.98, 760], [25.00, 760], [25.02, 677], [25.03, 531], [25.05, 382], [25.06, 277], [25.08, 255], [25.09, 150], [25.11, 113], [25.12, 86], [25.14, 67], [25.16, 54], [25.17, 45], [25.19, 38], [25.20, 33], [25.22, 29], [25.23, 27], [25.25, 25], [25.27, 25], [25.28, 24], [25.30, 24], [25.31, 24], [25.33, 24], [25.34, 27], [25.36, 26], [25.37, 27], [25.39, 28], [25.41, 31], [25.42, 36], [25.44, 45], [25.45, 64], [25.47, 91], [25.48, 125], [25.50, 162], [25.52, 205], [25.53, 249], [25.55, 286], [25.56, 324], [25.58, 370], [25.59, 511], [25.61, 475], [25.62, 531], [25.64, 580], [25.66, 623], [25.67, 667], [25.69, 708], [25.70, 749], [25.72, 760], [25.73, 760], [25.75, 760], [25.76, 759], [25.78, 760], [25.80, 760], [25.81, 760], [25.83, 722], [25.84, 603], [25.86, 439], [25.87, 317], [25.89, 231], [25.91, 170], [25.92, 127], [25.94, 96], [25.95, 74], [25.97, 59], [25.98, 48], [26.00, 40], [26.01, 34], [26.03, 30], [26.05, 28], [26.06, 26], [26.08, 24], [26.09, 23], [26.11, 22], [26.12, 21], [26.14, 21], [26.16, 21], [26.17, 21], [26.19, 21], [26.20, 20], [26.22, 20], [26.23, 20], [26.25, 21], [26.26, 23], [26.28, 28], [26.30, 39], [26.31, 65], [26.33, 106], [26.34, 163], [26.36, 234], [26.37, 315], [26.39, 399], [26.40, 486], [26.42, 572], [26.44, 655], [26.45, 736], [26.47, 760], [26.48, 767], [26.50, 760], [26.51, 760], [26.53, 759], [26.55, 759], [26.56, 760], [26.58, 760], [26.59, 760], [26.61, 760], [26.62, 760], [26.64, 760], [26.65, 760], [26.67, 759], [26.69, 741], [26.70, 602], [26.72, 435], [26.73, 314], [26.75, 229], [26.76, 169], [26.78, 126], [26.80, 96], [26.81, 74], [26.83, 255], [26.84, 48], [26.86, 40], [26.87, 34], [26.89, 30], [26.90, 27], [26.92, 25], [26.94, 24], [26.95, 23], [26.97, 23], [26.98, 22], [27.00, 22], [27.01, 21], [27.03, 22], [27.05, 22], [27.06, 23], [27.08, 27], [27.09, 35], [27.11, 51], [27.12, 81], [27.14, 122], [27.15, 175], [27.17, 238], [27.19, 309], [27.20, 386], [27.22, 467], [27.23, 552], [27.25, 638], [27.26, 726], [27.28, 760], [27.30, 760], [27.31, 760], [27.33, 760], [27.34, 767], [27.36, 760], [27.37, 760], [27.39, 760], [27.40, 760], [27.42, 760], [27.44, 759], [27.45, 760], [27.47, 760], [27.48, 760], [27.50, 760], [27.51, 760], [27.53, 760], [27.55, 694], [27.56, 527], [27.58, 379], [27.59, 275], [27.61, 201], [27.62, 149], [27.64, 112], [27.65, 86], [27.67, 67], [27.69, 54], [27.70, 44], [27.72, 37], [27.73, 33], [27.75, 29], [27.76, 26], [27.78, 24], [27.80, 23], [27.81, 23], [27.83, 22], [27.84, 23], [27.86, 28], [27.87, 33], [27.89, 37], [27.90, 38], [27.92, 38], [27.94, 38], [27.95, 38], [27.97, 42], [27.98, 51], [28.00, 71], [28.01, 107], [28.03, 160], [28.05, 226], [28.06, 299], [28.08, 382], [28.09, 464], [28.11, 546], [28.12, 627], [28.14, 705], [28.15, 760], [28.17, 760], [28.19, 760], [28.20, 760], [28.22, 760], [28.23, 767], [28.25, 760], [28.26, 760], [28.28, 760], [28.30, 760], [28.31, 760], [28.33, 760], [28.34, 760], [28.36, 760], [28.37, 759], [28.39, 760], [28.40, 760], [28.42, 760], [28.44, 760], [28.45, 760], [28.47, 696], [28.48, 522], [28.50, 376], [28.51, 272], [28.53, 199], [28.55, 148], [28.56, 112], [28.58, 255], [28.59, 67], [28.61, 53], [28.62, 44], [28.64, 38], [28.65, 33], [28.67, 30], [28.69, 30], [28.70, 36], [28.72, 50], [28.73, 71], [28.75, 91], [28.76, 107], [28.78, 119], [28.79, 127], [28.81, 139], [28.83, 156], [28.84, 182], [28.86, 220], [28.87, 270], [28.89, 331], [28.90, 402], [28.92, 482], [28.94, 569], [28.95, 662], [28.97, 759], [28.98, 760], [29.00, 760], [29.01, 760]]
	while len(sampleList) > 256:
		sampleList.pop()
	return sampleList

if __name__ == '__main__':
	# sampleList = readFile()
	# sampleList = waveGen()
	sampleList = rawData();
	mainLoop(sampleList)
	# printTransformArray(sampleList)
	# plt.plot(getTime(sampleList), getWave(sampleList))
	# plt.show()
	# plt.plot(printTransformArray(sampleList))
	# plt.show()
