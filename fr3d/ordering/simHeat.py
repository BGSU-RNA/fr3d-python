import matplotlib.pyplot as plt
import numpy as np
from metric import Metric


def simHeat(m, path):
  data = np.zeros((len(m.points), len(m.points)))
  for i in range(0,len(m.points)):
    for j in range(0, len(m.points)):
      data[i][j] = m.d(m.points[i], m.points[j])
  column_labels = m.points
  row_labels = m.points
  fig, ax = plt.subplots()
  heatmap = ax.pcolor(data, cmap=plt.cm.bwr)

  # put the major ticks at the middle of each cell
  ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
  ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

  # want a more natural, table-like display
  ax.invert_yaxis()
  ax.xaxis.tick_top()

  ax.set_xticklabels(row_labels, minor=False)
  ax.set_yticklabels(column_labels, minor=False)
  plt.show()
