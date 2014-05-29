import time
import numpy, sys, os
from matplotlib import pyplot as plt
from collections import deque
import serial

class HRMplot:

	def __init__ (self):
		self.maxLen = 1500
		self.time = deque([0.0]*self.maxLen, self.maxLen)
		self.data = deque([0.0]*self.maxLen, self.maxLen)
		self.beats = deque([0.0]*self.maxLen, self.maxLen)
		self.bpm = deque([0.0]*self.maxLen, self.maxLen)
		self.time_stamps = []

		self.scroll_width = 20
		self.threshold = 400
		self.pos = 0
		self.stamps_tail = 0
		self.stamps_head = 0
		self.buffer_pulses = 0

		plt.ion()
		self.plot_raw = plt.plot(self.time,self.data,color=	'r',linewidth=2)
		self.plot_beats = plt.plot(self.time,self.beats,linewidth=1) 
		self.title = plt.title('HRM Output',fontsize=15)
		plt.ylabel('Arduino',fontsize=10)
		plt.xlabel("Time (s)",fontsize=10)
		self.local_min = 0
		self.local_max = 1000
		plt.axis([0, self.scroll_width, self.local_min, self.local_max])
		plt.grid(True)

	def update(self, data_pair):
		bpm_str = ''

		self.time.append(data_pair[0])
		self.data.append(data_pair[1])
		self.pos = len(self.time)-1 			# array index of last stored value
		self.setAxes()
		plt.setp(self.plot_raw,'xdata',self.time,'ydata',self.data);
		if(self.data[self.pos] > self.threshold):
			self.beats.append(900)
			self.beatFinder(self.beats, self.pos, self.bpm, self.time_stamps, self.stamps_head, self.time[self.pos], self.buffer_pulses)
		else:
			self.beats.append(0)
		self.calcBPM(self.pos, self.time_stamps, self.stamps_head, self.stamps_tail, self.buffer_pulses, self.bpm, bpm_str)
		plt.setp(self.plot_beats, 'xdata', self.time, 'ydata', self.beats)
		plt.draw()

	def setAxes(self):
		# scale the local_min and local_max based on the range of recent data values
		buffer_time = 2
		finish_pos = self.pos
		start_pos = self.pos
		current_time = self.time[self.pos]
		start_time = current_time-buffer_time
		while(start_pos > 0 and self.time[start_pos] > start_time):
			start_pos -= 1

		self.local_min = self.data[start_pos]
		self.local_max = self.data[finish_pos]

		for i in range (start_pos, finish_pos):
			data_point = self.data[i]
			self.local_min = min(data_point, self.local_min)
			self.local_max = max(data_point, self.local_max)

		self.local_min = max(self.local_min-5, 0)
		self.local_max = min(self.local_max+5, 1000)
		self.threshold = self.local_max - (self.local_max - self.local_min)/2
		# print 'threshold = ' + str(self.threshold)

		if(self.time[self.pos] > self.scroll_width):
			plt.axis([self.time[self.pos] - self.scroll_width, self.time[self.pos], self.local_min, self.local_max])
		else:
			plt.axis([0, self.scroll_width, self.local_min, self.local_max])

	def beatFinder(self, beats, pos, bpm, time_stamps, stamps_head, current_time, buffer_pulses):
		# filter the beats deque to determine which values are actually heartbeats.
		# delete false positives and add in missed beats as necessary

		# allow the heartbeat to vary by (at most) this amount each time a beat
  		# is found
		variation = 1.4
		diff = 0
		lower_bound = 0
		upper_bound = 250
		missed_1 = 0
		missed_2 = 0
		inst_bpm = 0

  		# only accept a new beat if it's close enough to the established
    	# bpm. EX. when you divide by 1.75 this beat could have come in at
    	# <= 1.75 times the rate of the previously established heart rate
		if(pos > 0 and bpm[pos-1] > 0):
			upper_bound = bpm[pos-1]*variation
			missed_1 = bpm[pos-1]/variation
			missed_2 = bpm[pos-1]/(2 * variation)
			lower_bound = bpm[pos-1]/(3*variation)

		if(stamps_head > 0):
			diff = current_time-time_stamps[stamps_head]
			inst_bpm = 60/diff

		if(pos > 0 and beats[pos-1] == 0 and (bpm[pos-1] == 0 or inst_bpm < upper_bound)):
			print 'beat found, stamps_head = ' + str(self.stamps_head)
			# detect if one beat was missed
			if(bpm[pos-1] > 0 and inst_bpm < missed_1 and inst_bpm > missed_2):
				print 'add 1, instantaneous: ' + str(inst_bpm) + ', last: ' + str(bpm[pos-1])
				
				time_stamps.append(diff/2 + time_stamps[stamps_head])
				buffer_pulses += 1

			# detect if two were missed
			elif(bpm[pos-1] > 0 and inst_bpm < missed_2 and inst_bpm > lower_bound):
				print 'add 2, instantaneous: ' + str(inst_bpm) + ', last: ' + str(bpm[pos-1])
				time_stamps.append(diff/3 + time_stamps[stamps_head])
				time_stamps.append(2*diff/3 + time_stamps[stamps_head])
				buffer_pulses += 2
			time_stamps.append(current_time)
			stamps_head = len(time_stamps) - 1
			buffer_pulses += 1
		# check if a beat shouild be ignored
		else:
			if(pos > 0 and beats[pos-1] == 0):
				print 'ignoring beat, instantaneous: ' + str(inst_bpm) + ', last: ' + str(bpm[pos-1])
		
		self.time_stamps = time_stamps
		self.stamps_head = stamps_head
		self.buffer_pulses = buffer_pulses

	def calcBPM(self, pos, time_stamps, stamps_head, stamps_tail, buffer_pulses, bpm, bpm_str=''):
		# Calculate beats per minute from the last 'buffer' seconds of data. 
		# If there is not enough new data keep the last calculated heart rate (or 0,
		# the sentinel). If the calculated heart rate is too extreme set it to the
		# last legitimate rate (or 0)

		buffer_time = 10
		lower_bound = 45
		upper_bound = 250

		if (pos > 0 and bpm[pos-1] > 0):
			bpm_str = str(bpm[pos-1])
		else:
			bpm_str = '...'

		if(stamps_head > 0 and time_stamps[stamps_head] - time_stamps[stamps_tail] > buffer_time):
			while(time_stamps[stamps_head] - time_stamps[stamps_tail] > buffer_time):
				buffer_pulses -= 1
				stamps_tail += 1
			diff = time_stamps[stamps_head] - time_stamps[stamps_tail]

			calc_bpm = 60*((buffer_pulses-1)/diff)

			if(calc_bpm > upper_bound or calc_bpm < lower_bound):
				bpm_str = '...'
				if(pos > 0):
					bpm.append(bpm[pos-1])
				else:
					bpm.append(0)
			elif(pos > 0 and bpm[pos-1] > 0):
				print 'calc_bpm: ' + str(calc_bpm) + ' and ' + str(bpm[pos-1])
				if(calc_bpm < bpm[pos-1]):
					bpm.append(bpm[pos-1]-1)
				elif(calc_bpm > bpm[pos-1]):
					bpm.append(bpm[pos-1]+1)
				else:
					bpm.append(bpm[pos-1])
				bpm_str = str(bpm[pos])
			else:
				bpm.append(int(calc_bpm))
				bpm_str = str(bpm[pos])
			if(pos > 0 and (bpm[pos] != bpm[pos-1])):
				print 'buffer pulses: ' + str(buffer_pulses) + ', diff: ' + str(diff)
		elif(pos > 0):
			bpm.append(bpm[pos-1])
		else:
			bpm.append(0)

		bpm_str = 'HRM Output: BPM = ' + bpm_str
		if(pos > 0 and bpm[pos] != bpm[pos-1]):
			print bpm_str

		self.bpm = bpm
		self.bpm_str = bpm_str
		self.buffer_pulses = buffer_pulses
		self.stamps_tail = stamps_tail






