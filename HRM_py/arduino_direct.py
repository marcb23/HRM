import time
import numpy as np
from collections import deque
import serial
from matplotlib import pyplot as plt

# credit for following two classes goes to https://gist.github.com/electronut/5641933
# class that holds analog data for N samples
class AnalogData:
  # constr
	def __init__(self, maxLen):
		self.ax = deque([0.0]*maxLen, maxLen)
		self.ay = deque([0.0]*maxLen, maxLen)
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
	def add(self, data):
		print 'adding to buffer...'
		assert(len(data) == 2)
		self.addToBuf(self.ax, data[0])
		self.addToBuf(self.ay, data[1])
    
# plot class
class AnalogPlot:
  # constr
	def __init__(self, analogData):
    # set plot to animated
		plt.ion() 
		self.axline, = plt.plot(analogData.ax)
		self.ayline, = plt.plot(analogData.ay)
		plt.ylim([0, 1023])
 
  # update plot
	def update(self, analogData):
		self.axline.set_ydata(analogData.ax)
		self.ayline.set_ydata(analogData.ay)
		plt.draw()

def dataLoop(ser):
	# plot parameters
  	analogData = AnalogData(1000)
  	analogPlot = AnalogPlot(analogData)

	print 'plotting...'
	start = time.clock()
	while True:
		try:
			raw_in = ser.readline()
			print raw_in
			data = [time.clock(), float(ser.readline())]
			analogData.add(data)
			analogPlot.update(analogData)
			time.sleep(.01)
		except KeyboardInterrupt:
			print 'exiting'
			break
		
	ser.flush()
	ser.close()


if __name__ == '__main__':
	ser = serial.Serial('/dev/tty.usbmodem1421', 9600)
	print ser.readline() # Read the newest output from the Arduino	
	ser.write("go")	
	dataLoop(ser)
	
