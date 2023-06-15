from BayesianClustering import *
# from sys import argv
# import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
def Callback( ScanResults ):
	#print( "Python displaying" , ScanResults )
	for i in ScanResults:
		print( i.r , i.t , i.score );

	exit()

	# X , Y = [] , []
	# for point in background:
	# 	X.append( 1e6 * point.x )
	# 	Y.append( 1e6 * point.y )

	# x , y , colours = [] , [] , []
	# for colour , cluster in enumerate( clusters ) :
	# 	for point in cluster:
	# 		x.append( 1e6 * point.x )
	# 		y.append( 1e6 * point.y )
	# 		colours.append( colour )

	# plt.scatter( X, Y, s=0.1 , marker='.' , c="black" )
	# plt.scatter( x, y, s=0.1 , marker='.' , c=colours , cmap="prism" )
	# plt.xlabel( r'$\mu$m' )
	# plt.ylabel( r'$\mu$m' )
	# plt.show()
# ----------------------------------------------------------------------------

Cfg = ScanConfiguration( "example-configs/config.txt" )
AutoRoi_Scan_SimpleCallback( "1_un_red.csv" , Cfg , Callback )
