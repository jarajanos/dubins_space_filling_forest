import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Polygon
import sys

#COLORS = ['black', 'grey', 'lightgray', 'rosybrown', 'brown', 'darkred', 'red', 'sienna', 'saddlebrown', 'peru', 'darkorange', 'tan', 'gold', 'olivedrab', 'chartreuse', 'seagreen', 'darkgreen', 'darkslategrey', 'deepskyblue', 'slategray', 'royalblue', 'navy', 'indigo', 'plum', 'purple', 'deeppink', 'crimson']
COLORS = ['tan', 'plum', 'brown', 'gold', 'purple', 'deepskyblue', 'crimson', 'deeppink', 'rosybrown', 'darkred', 'peru', 'sienna', 'seagreen', 'indigo', 'chartreuse', 'darkorange', 'olivedrab', 'slategray', 'darkslategrey', 'navy', 'grey', 'saddlebrown', 'red', 'darkgreen']

pts = []
p_col = []
lines = []
tris = []
lim_x = [float('inf'), float('-inf')]
lim_y = [float('inf'), float('-inf')]

color_index = 0

ax = plt.gca()
#ax.set_title(sys.argv[1])
for file in sys.argv[1:]:
  nodes = open(file, 'r')

  #first line
  line = nodes.readline()
  arr = line[1:].split()
  
  file_type = arr[0].strip()
  dimension = arr[1].strip()

  num_point = 0
  if dimension == "2D":
    num_point = 2
  elif dimension == "2DDubins":
    num_point = 3
  elif dimension == "2DPolynom":
    num_point = 3
  else:
    num_point = 6

  num_point_offset = 0
  if file_type == "Trees":
    num_point_offset = 2

  for line in nodes:
    if line.strip() == "" or line[0] == '#':
      continue

    ar = list(map(float, line.strip().split(' ')))
    iter = 0
    for pt in ar:
      if iter % num_point == 0:
        lim_x[0] = min(lim_x[0], pt)
        lim_x[1] = max(lim_x[1], pt)
      elif iter % num_point == 1:
        lim_y[0] = min(lim_y[0], pt)
        lim_y[1] = max(lim_y[1], pt)
      
      iter += 1
      if iter == len(ar) - num_point_offset:
        break

    if len(ar) == 3 * num_point + num_point_offset:
      # triangle
      loc_pts = np.array([[ar[0], ar[1]], [ar[num_point], ar[num_point + 1]], [ar[2 * num_point], ar[2 * num_point + 1]]])
      tris.append(Polygon(loc_pts, fill=True, color='royalblue'))
    elif len(ar) == 2 * num_point + num_point_offset:
      # line
      lines.append(mlines.Line2D([ar[0], ar[num_point]], [ar[1], ar[num_point + 1]], lw=1.5, color=(COLORS[int(ar[-num_point_offset]) % len(COLORS)] if num_point_offset != 0 else 'black'))) # COLORS[color_index])))
      #lines.append(mlines.Line2D([ar[0], ar[num_point]], [ar[1], ar[num_point + 1]], lw=1.5, color=(COLORS[0] if num_point_offset != 0 else 'black'))) # COLORS[color_index])))
    elif len(ar) == num_point + num_point_offset:
      # point
      pts.append([ar[0], ar[1]])
      #if num_point_offset == 0:
        #p_col.append([COLORS[ar[2] % len(COLORS)]])
      #  p_col.append('lightgray')
      #else:
      p_col.append('black')
  
  color_index += 1
  color_index = color_index % len(COLORS)

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.set_xlim(lim_x[0], lim_x[1])
#print(lim_x)
#ax.set_xlim(-10, 2010)
ax.set_ylim(lim_y[0], lim_y[1])
#print(lim_y)
for pol in tris:
  ax.add_patch(pol)

for line in lines:
  ax.add_line(line)

np_pts = np.asarray(pts)
np_col = np.asarray(p_col)
if len(np_pts) != 0:
  colors = np.unique(np_col)
  for c in colors:
    #print(np_col)
    points = np_pts[np_col == c, :]
    size = 6 ** 2 if c == 'black' else 3 ** 2
    ax.scatter(points[:, 0], points[:, 1], c=c, zorder=5)

plt.show()

