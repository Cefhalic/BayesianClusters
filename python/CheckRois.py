from BayesianClusteringTools import *
import matplotlib.pyplot as plt
import numpy as np

# ---------------- Main process -----------------
Locs , RoIs = CheckRoIs( "stem.csv"  , "ROI/stem.zip" , 25 * nanometer )

plt.scatter( *Locs , s=0.1 , marker='.' , c="black" )

for RoI in RoIs: plt.scatter( *RoI , s=0.1 , marker='.' )

plt.axis('off')
plt.margins(x=0,y=0)
plt.tight_layout() 
plt.show()
# -----------------------------------------------
