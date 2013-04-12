import sys
import string
import math
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

total_x_ticks = 5
total_y_ticks = 5

def SaveImage(num, filename, x_axis, y_axis, values, hist2d, axis_range, title, axis_names):
  plt.figure(num)
  plt.subplot(121)
  plt.title(title)
  plt.hexbin(x_axis, y_axis, values, bins='log', cmap=cm.jet)
  plt.axis( [ axis_range[0], axis_range[1], axis_range[2], axis_range[3] ] )
  plt.xlabel(axis_names[0])
  plt.ylabel(axis_names[1])

  plt.subplot(122)
  plt.contour(hist2d,20)
  plt.axis(axis_range)
  plt.colorbar()
  plt.savefig(filename)

for num, filename in enumerate(sys.argv[1:]): # loop over input files
  print(filename)
  hist2d_file=open(filename)                                                   # open file
  first_line = hist2d_file.readline()                                          # Read first line
  print string.split(first_line, "  ")
  [chains, starts, ends] = string.split(first_line, "  ")   # first line
  chain = [n for n in chains.split(',')]
  start = [int(n) for n in starts.split(',')]
  end   = [int(n) for n in ends.split(',')]

  size = 0.;
  for i in range(len(chain)):
    size += end[i]-start[i] + 1
  hist2d = np.zeros([ int(size), int(size)], float)                # histogram 2d

  x_axis=[]        # it contains x-coordinates of non-null values on the hist
  y_axis=[]        # it contains y-coordinates of non-null values on the hist
  values=[]        # values[i] contains the value of the histogram at coordinates (x_axis[i],y_axis[i])
  hist2d = np.zeros([ int(size), int(size)], float)                # histogram 2d

  axis_range = ( 0, size, 0,size )
  axis_labels = ("residue","residue")

  for line in hist2d_file:
    [x_ch, x_am, y_ch, y_am, value]=string.split(line, "  ")                               # Read data
    x_am = int(x_am)
    y_am = int(y_am)
    for i in range(len(chain)):
      if (chain[i]==x_ch): 
        x_am -= start[i]
        break
      else: x_am += end[i]-start[i] + 1
    for i in range(len(chain)):
      if (chain[i]==y_ch): 
        y_am -= start[i]
        break
      else: y_am += end[i]-start[i] + 1
    x_axis.append(x_am)                                                  # store x
    y_axis.append(y_am)                                                  # store y
    values.append( float(value) )                                     # store value at (x,y)
    hist2d[x_am][y_am]=float(value)                      # store histogram values at (x_bin,y_bin)
  SaveImage(num, filename[:-3]+'png', x_axis, y_axis, values, hist2d, axis_range, "", axis_labels) # save image
  hist2d_file.close()                                                                                  # close file
