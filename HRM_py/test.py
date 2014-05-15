import time

# for Mac OSX
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pylab as plt
import random

# from mpltools import style
# style.use('ggplot')

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

def test_fps(use_blit=True):

    ax1.cla()
    ax1.set_title('Sensor Input vs. Time -' + 'Blit [{0:3s}]'.format("On" if use_blit else "Off"))
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Sensor Input (mV)')

    plt.ion()  # Set interactive mode ON, so matplotlib will not be blocking the window
    plt.show(False)  # Set to false so that the code doesn't stop here

    cur_time = time.time()
    ax1.hold(True)

    x, y = [], []
    times = [time.time() - cur_time]  # Create blank array to hold time values
    y.append(0)

    line1, = ax1.plot(times, y, '.-', alpha=0.8, color="gray", markerfacecolor="red")

    fig.show()
    fig.canvas.draw()

    if use_blit:
        background = fig.canvas.copy_from_bbox(ax1.bbox) # cache the background

    tic = time.time()

    niter = 200
    i = 0
    while i < niter:

        fields = random.random() * 100

        times.append(time.time() - cur_time)
        y.append(fields)

        # this removes the tail of the data so you can run for long hours. You can cache this
        # and store it in a pickle variable in parallel.

        if len(times) > 50:
            del y[0]
            del times[0]

        xmin, xmax, ymin, ymax = [min(times) / 1.05, max(times) * 1.1, -5,110]

        # feed the new data to the plot and set the axis limits again
        line1.set_xdata(times)
        line1.set_ydata(y)

        plt.axis([xmin, xmax, ymin, ymax])

        if use_blit:
            fig.canvas.restore_region(background)    # restore background
            ax1.draw_artist(line1)                   # redraw just the points
            fig.canvas.blit(ax1.bbox)                # fill in the axes rectangle
        else:
            fig.canvas.draw()

        i += 1

    fps = niter / (time.time() - tic)
    print "Blit [{0:3s}] -- FPS: {1:.1f}, time resolution: {2:.4f}s".format("On" if use_blit else "Off", fps, 1/fps)
    return fps

fps1 = test_fps(use_blit=True)
fps2 = test_fps(use_blit=False)

print "-"*50
print "With Blit ON plotting is {0:.2f} times faster.".format(fps1/fps2)