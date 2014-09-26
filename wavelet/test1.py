#! /opt/local/bin/python2.7
from numpy import *
from math import pi
from random import sample, random
import matplotlib.pyplot as plt
from pywt import dwt, wavedec

n = 512			# samples
freq0 = 2.5 	# Hz
samp_rate = 64	# Hz
levels = 8

start_freq = 1	# Hz
end_freq = 1	# Hz
if (start_freq != end_freq):
	freq0 = arange(start_freq, end_freq, (end_freq - start_freq) / (n * 1.0))
else:
	freq0 = start_freq


factor0 = samp_rate / freq0
time = arange(n)/float(samp_rate)
wave0 = sin(2 * pi * freq0 * time)

# errors = [random() - 0.5 for _ in range(n)]
# wave0 += errors

cA, cD = dwt(wave0, 'db2')
coeffs = wavedec(wave0, 'db2', level=levels)

i = 0
print len(coeffs)
nyquist = samp_rate / 2.

freqs = zeros(levels)
powers = zeros(levels)
freq_widths = zeros(levels)

for a in coeffs:
	if (i <> 0):
		print i, len(a), mean(abs(a)), std(a), nyquist / 2.**(levels-i+1), nyquist / 2.**(levels-i)
		# I don't know why the 2 to 1 ratio works better...
		freqs[i-1] = (2 * nyquist / 2.**(levels-i+1) + nyquist / 2.**(levels-i))/3
		# freqs[i-1] = (nyquist / 2.**(levels-i+1) + nyquist / 2.**(levels-i)) / 2.
		freq_widths[i-1] = nyquist / 2.**(levels-i) - nyquist / 2.**(levels-i+1) 
		powers[i-1] = mean(abs(a)) # if freqs[i-1] <= 8.0 else 0	# ignore frequencies higher than 4 Hz
	i += 1

# find three biggest consecutive sums
big_sum = 0.
start_pos = 0
core_freq = 0.

for i in arange(levels-5):
	sum = 0.
	freq_weight = 0.
	for j in arange(5):
		freq_weight += freqs[i+j] * powers[i+j] * freq_widths[i+j]
		sum +=  powers[i+j]
	if (sum > big_sum):
		big_sum = sum
		core_freq = freq_weight / sum

# for i in arange(levels-3):
# 	sum = 0.
# 	freq_weight = 0.
# 	for j in arange(3):
# 		freq_weight += freqs[i+j] * powers[i+j]
# 		sum +=  powers[i+j]
# 	if (sum > big_sum):
# 		big_sum = sum
# 		start_pos = i
# 		core_freq = freq_weight / sum

print core_freq

plt.plot(time, wave0)
plt.show()

"""
plt.plot(cA)
plt.show()

plt.plot(cD)
plt.show()
"""

# print coeffs