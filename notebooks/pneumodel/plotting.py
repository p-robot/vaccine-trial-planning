import matplotlib.pyplot as plt
import os

def remove_spines(axes=None, top=False, right=False, left=True, bottom=True):
  """ Minimize chartjunk by stripping out unnecessary plot borders and axis ticks.

  :param axes: If None, gets the current axis through matplotlib.pyplot.gca().
  :param top/right/left/bottom: These toggle whether the corresponding plot border is drawn.

  .. _link: https://github.com/cs109/content/blob/caffc21c8f7c758c1884852ed023d29dccea063f/HW2.ipynb

  """
  ax = axes or plt.gca()
  ax.spines['top'].set_visible(top)
  ax.spines['right'].set_visible(right)
  ax.spines['left'].set_visible(left)
  ax.spines['bottom'].set_visible(bottom)
    
  # turn off all ticks
  ax.yaxis.set_ticks_position('none')
  ax.xaxis.set_ticks_position('none')
  
  # now re-enable visibles
  if top:
    ax.xaxis.tick_top(True)
  if bottom:
    ax.xaxis.tick_bottom()
  if left:
    ax.yaxis.tick_left()
  if right:
    ax.yaxis.tick_right()

def save_figure(fig, filename, folder='../figures', exts=['pdf', 'png'], **kwargs):
  """ Save a matplotlib figure.

  :param fig: The matplotlib figure to save.
  :param filename: The name of the saved file. "-[n]" will be appended to this name, where n is the smallest number that makes the name unique.
  :param folder: Save the file to this folder.
  :param exts: Save the file as these file types. Default is 'pdf'.
  :param kwargs: Keyword arguments for plt.savefig(), e.g. additional_artists.

  """
  if not os.path.exists(folder):
    os.makedirs(folder)
  paths = []
  for ext in exts:
    i = 0
    while True:
      path = '{}/{}-{:d}.{}'.format(folder, filename, i, ext)
      if not os.path.exists(path):
        break
      i += 1
    plt.savefig(path, **kwargs)
    paths.append(path)
  return paths

def get_boxplot_style(color, lw=2, alpha=1.0, width=0.5):
  """ Get keyword arguments for a boxplot style.
  
  :param color: Color of the boxes, whiskers, median, caps, and fliers
  :param lw: Line width
  :param alpha: Opacity
  :param width: Width of the box


  :return: A dictionary of keyword arguments for use with matplotlib.pyplot.boxplot().

  """
  return dict(
    sym=None, 
    widths=width, 
    whis='range',
    boxprops=dict(color=color, alpha=alpha, lw=lw),
    whiskerprops=dict(color=color, alpha=alpha, lw=lw),
    medianprops=dict(color=color, alpha=alpha, lw=lw),
    capprops=dict(color=color,alpha=alpha, lw=lw),
    flierprops=dict(color=color, alpha=alpha, marker='o', markersize=7)
  )

def get_errorbar_style(color, lw=1.5):
  """ Get keyword arguments for a boxplot style.
  
  :param color: Color of the error bars
  :param lw: Line width

  :return: A dictionary of keyword arguments for use with matplotlib.pyplot.errorbar().

  """
  return dict(
    marker='s', 
    ls='none', 
    lw=lw, 
    capsize=3, 
    capthick=1, 
    ecolor=color, 
    color=color
  )

def style_violin(violin_parts, violin_color, stem_color):
  for key, part in violin_parts.iteritems():
    if key == 'bodies':
      for subpart in violin_parts['bodies']:
        subpart.set_facecolor(c=violin_color)
        subpart.set_edgecolor(c=violin_color)
    else:
      part.set_facecolor(stem_color)
      part.set_edgecolor(stem_color)