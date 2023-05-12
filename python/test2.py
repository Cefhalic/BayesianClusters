from BayesianClustering import *
from sys import argv
import matplotlib.pyplot as plt
import numpy as np

def Callback( clusters , background  ):
	print( "Displaying")
	pts = [ len( cluster ) for cluster in clusters ]
	plt.boxplot( pts , labels= ["$Dataset A$"] , widths= 0.8 , showfliers=False , showmeans=True , meanprops=dict(color="grey"), meanline=True, medianprops=dict(color="black") )
	plt.scatter( np.random.normal( 1 , 0.05 , len( pts ) ) , pts , color="r" , alpha=0.5 , s=1 )
	plt.xlabel('x-label')
	plt.ylabel('Cluster size')
	plt.yscale( "log" )
	plt.show()

Configuration.FromVector( argv[1:] );
OneStopGetClusters( Callback )