from BayesianClustering import *
import matplotlib.pyplot as plt
import numpy as np


# -------------- Key Press Handler --------------
def on_press(event):
  global Scale , RoIs , RoiPlots , txt

  if event is None: pass
  elif event.key == '+': Scale += 0.1
  elif event.key == '-': Scale -= 0.1
  else: return

  txt.set_text( f"Scale={Scale}nm | Press + or - to modify" )
  for Plot in RoiPlots: Plot.remove()
  RoiPlots = [ plt.plot( np.array( x ) * Scale * nanometer , np.array( y ) * Scale * nanometer )[0] for x,y in RoIs ]
  plt.gcf().canvas.draw()
# -----------------------------------------------


# ---------------- Main process -----------------
Scale , RoiPlots = 25 , []
txt = plt.title( "" )

plt.scatter( *GetLocalizations( "stem.csv" ) , s=0.1 , marker='.' , c="black" )

RoIs = GetRoIs( "ROI/stem.zip" )
for RoI in RoIs:
  RoI[0].append( RoI[0][0] )
  RoI[1].append( RoI[1][0] )

on_press( None )

plt.gcf().canvas.mpl_connect( 'key_press_event' , on_press )
plt.axis('off')
plt.margins(x=0,y=0)
plt.tight_layout() 
plt.show()
# -----------------------------------------------
