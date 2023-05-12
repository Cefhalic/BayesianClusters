from BayesianClustering import *
from sys import argv
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
def Callback( clusters , background  ):
	print( "Displaying")
	pts = [ len( cluster ) for cluster in clusters ]
	plt.xlabel('Cluster size')	
	plt.ylabel('Count')
	plt.hist( pts , 100, facecolor='green', alpha=0.75)
	plt.show()
# ----------------------------------------------------------------------------

Configuration.FromVector( argv[1:] );
OneStopGetClusters( Callback )