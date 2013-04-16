import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import sys

def imsave(filename, X, xtics, xtics_labels, ytics, ytics_labels, title, **kwargs):
  fig = plt.figure()    # Create a figure
  ax = fig.add_subplot(111)
  ax.grid()
  plt.xticks(xtics, xtics_labels, rotation=45)
  plt.yticks(ytics, ytics_labels)
  plt.title( title )
  xylabel = r"$chainId-resSeq$"
  plt.xlabel(xylabel)
  plt.ylabel(xylabel)
  imgplot = plt.imshow(X , origin='lower', interpolation='nearest', **kwargs) # Plot figure
  plt.colorbar(drawedges=False, **kwargs)  # Add the colorbar
  plt.savefig(filename, dpi=1000, bbox_inches='tight')  # Save to file
  plt.close(fig)

for num, filename in enumerate(sys.argv[1:]): # Save a matrix of contacts and print to a png file
  print "Now producing a graph for ", filename, " file."
  contact_file=open(filename)                                                   # open file
  title = contact_file.readline()
  first_line = contact_file.readline()                                          # Read first line
  [x_name, x_start, x_end] = first_line.split( "  ")                           # Read first line
  second_line = contact_file.readline()                                         # Read second line
  [y_name, y_start, y_end] = second_line.split("  ")                          # Read second line

  x_names = [ n for n in x_name.split(',')]
  x_starts = [int(n) for n in x_start.split(',')]
  x_ends = [int(n) for n in x_end.split(',')]
  y_names = [ n for n in y_name.split(',')]
  y_starts = [int(n) for n in y_start.split(',')]
  y_ends = [int(n) for n in y_end.split(',')]
  xoffsets = {x_names[0]:0}
  yoffsets = {y_names[0]:0}
  xstarts = {}
  ystarts = {}
  
  if( len(x_names)!=len(x_starts) or len(x_names)!=len(x_ends) ):
    print "file formatting error"
  if( len(y_names)!=len(y_starts) or len(y_names)!=len(y_ends) ):
    print "file formatting error"

  x_sizes = [ 0 for n in range(len(x_names)) ]
  y_sizes = [ 0 for n in range(len(y_names)) ]
  cols = 0;
  rows = 0;
  for i in range(len(x_names)):
    x_sizes[i] = x_ends[i] - x_starts[i] + 1
    cols = cols + x_sizes[i]
    xstarts[x_names[i]] = x_starts[i]
    if(i>0):
      xoffsets[x_names[i]] = xoffsets[x_names[i-1]]+x_sizes[i]
  for i in range(len(y_names)):
    y_sizes[i] = y_ends[i] - y_starts[i] + 1
    rows = rows + y_sizes[i]
    ystarts[y_names[i]] = y_starts[i]
    if(i>0):
      yoffsets[y_names[i]] = yoffsets[y_names[i-1]]+y_sizes[i]
      
  contactmap = np.zeros([int(cols),int(rows)], float)                              # Initialize contact map
  for line in contact_file:
    [xch,xam,ych,yam,value]=line.split("  ")                               # Read contact map from file
    x = xoffsets[xch]+int(xam)-xstarts[xch]
    y = yoffsets[ych]+int(yam)-ystarts[ych]
    contactmap[int(y)][int(x)]=float(value) if float(value)>0. else np.nan                                   # Update contact map
  
  xtics = [ xoffsets[n] for n in x_names ]
  ytics = [ yoffsets[n] for n in y_names ]
  xtics.append(cols-1)
  ytics.append(rows-1)
  
  xtics_labels = []
  ytics_labels = []
  for i in range(len(x_names)):
    xtics_labels.append(r"$"+ x_names[i]+str(x_starts[i])+"$" )
  xtics_labels.append(r"$"+ x_names[len(x_names)-1]+str(x_ends[len(x_names)-1])+"$"  )
  for i in range(len(y_names)):
    ytics_labels.append(r"$"+ y_names[i]+str(y_starts[i])+"$"  )
  ytics_labels.append(r"$"+ y_names[len(y_names)-1]+str(y_ends[len(y_names)-1])+"$"  )
  
  imsave(filename[:-3]+'png', contactmap, xtics, xtics_labels, ytics, ytics_labels, title, cmap=plt.cm.jet )            # Print contact map to file
  contact_file.close()                                                          # close file
