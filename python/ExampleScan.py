from BayesianClustering import *

# ----------------------------------------------------------------------------
def Callback( ScanResults ):
	for i in ScanResults: print( f"Python callback reporting scan scores: r={i.r:16.6e}, t={i.t:16.6e} => {i.score:+16}" );
# ----------------------------------------------------------------------------

Cfg = ScanConfiguration( "example-configs/config.txt" )
AutoRoi_Scan_SimpleCallback( "1_un_red.csv" , Cfg , Callback )
