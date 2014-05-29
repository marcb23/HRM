import time
import numpy as np
from collections import deque
import serial
# from matplotlib import pyplot as plt
import HRMplot

# credit for following two classes goes to https://gist.github.com/electronut/5641933
# class that holds analog data for N samples
class AnalogData:
  # constr
	def __init__(self, maxLen):
		self.time = deque([0.0]*maxLen, maxLen)
		self.data = deque([0.0]*maxLen, maxLen)
		self.maxLen = maxLen

  # ring buffer
 	def addToBuf(self, buf, val):
 		buf.append(val)
		# if len(buf) < self.maxLen:
		# 	buf.append(val)
		# else:
		# 	buf.pop()
		# 	buf.append(val)
 
  # add data
	def add(self, data_pair):
		# print 'adding ' + str(data[1]) + ' at time ' + str(data[0])
		assert(len(data) == 2)
		self.addToBuf(self.time, data_pair[0])
		self.addToBuf(self.data, data_pair[1])
    
# plot class
class AnalogPlot:
  # constr
	def __init__(self, analogData):
    # set plot to animated
		plt.ion() 
		self.axline, = plt.plot(analogData.ax, analogData.ay)
		self.ayline, = plt.plot(analogData.ay)
		plt.ylim([0, 1023])
 
  # update plot
	def update(self, analogData):
		self.axline.set_ydata(analogData.ay)
		self.axline.set_xdata(analogData.ax)
		plt.draw()

def dataLoop(ser):
	# plot parameters
  	# analogData = AnalogData(1000)
  	# analogPlot = AnalogPlot(analogData)

  	plot = HRMplot.HRMplot()

	start = time.time()
	print 'start time: ' + str(start)
	while True:
		try:
			ser.write("1")
			# print "waiting for response.."
			raw_in = ser.readline()
			# print raw_in
			# try:
			data_pair = [(time.time()-start), float(raw_in)]
			# analogData.add(data)
			# analogPlot.update(analogData)
			# print "float works"

			plot.update(data_pair)

			# except:
			# 	print "rogue string: " + raw_in
			time.sleep(.005)
		except KeyboardInterrupt:
			print 'exiting'
			break
		
	ser.flush()
	ser.close()


if __name__ == '__main__':
	ser = serial.Serial('/dev/tty.usbmodem1411', 9600)
	connected = False
	while not connected:
		print ser.readline() # Read the newest output from the Arduino	
		connected = True
	ser.write("go")	
	dataLoop(ser)
	
