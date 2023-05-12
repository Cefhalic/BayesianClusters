from BayesianClustering import *
from sys import argv
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
def Callback( clusters , background  ):
	print( "Displaying")

	X , Y = [] , []
	for point in background:
		X.append( point.x )
		Y.append( point.y )

	x , y , colours = [] , [] , []
	for colour , cluster in enumerate( clusters ) :
		for point in cluster:
			x.append( point.x )
			y.append( point.y )
			colours.append( colour )

	plt.scatter( X, Y, s=0.1 , marker='.' , c="black" )
	plt.scatter( x, y, s=0.1 , marker='.' , c=colours , cmap="prism" )
	plt.show()
# ----------------------------------------------------------------------------

Configuration.FromVector( argv[1:] );
OneStopGetClusters( Callback )