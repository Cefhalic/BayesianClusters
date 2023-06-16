from BayesianClustering import *

# ----------------------------------------------------------------------------
def Callback( Clusters ):
	print( "Python callback reporting scan scores:" )
	print( "+----------------+----------------+----------------+----------------+----------------+" )
	print( "| Localizations  |      Area      |    Perimeter   |   Centroid x   |   Centroid y   |" )
	print( "+----------------+----------------+----------------+----------------+----------------+" )
	for i in Clusters: print( f"| {i.localizations:14} | {i.area:14.6} | {i.area:14.6} | {i.centroid_x:+14.6} | {i.centroid_y:+14.6} |" )
	print( "+----------------+----------------+----------------+----------------+----------------+" )
# ----------------------------------------------------------------------------

AutoRoi_Cluster_SimpleCallback( "1_un_red.csv" , 50*1e-9 , 170*1e-9 , Callback )