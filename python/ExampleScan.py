from BayesianClustering import *

# ----------------------------------------------------------------------------
def Callback( ScanResults ):
	for i in ScanResults: print( f"Python callback reporting scan scores: r={i.r:16.6e}, t={i.t:16.6e} => {i.score:+16}" );
# ----------------------------------------------------------------------------

Cfg = ScanConfiguration( "example-configs/config.txt" )

# def Interpolator( value ):
# 	return 1

#Cfg = ScanConfiguration( aSigmaBins=100 , aSigmaMin=5*nanometer , aSigmaMax=100*nanometer , aInterpolator=Interpolator ,	aRbins=35 , aMinScanR=0*nanometer , aMaxScanR=200*nanometer , aTbins=35, aMinScanT=0*nanometer, aMaxScanT=500*nanometer ,  aPB=0.2 , aAlpha=20)
# Cfg = ScanConfiguration( 100 , 5*nanometer , 100*nanometer , Interpolator , 35 , 0*nanometer , 200*nanometer , 35, 0*nanometer, 500*nanometer , 0.2 , 20.0 )

#Cfg = ScanConfiguration( aSigmaBins=100 , aSigmaMin=5*nanometer , aSigmaMax=100*nanometer , aInterpolator={ 0*nanometer:0.03631079 , 20*nanometer:0.110302441 , 30*nanometer:0.214839819 , 40*nanometer:0.268302465 , 50*nanometer:0.214839819 , 60*nanometer:0.110302441 , 70*nanometer:0.03631079 , 80*nanometer:0.007664194 , 90*nanometer:0.001037236 , 100*nanometer:9.00054E-05 } ,	aRbins=35 , aMinScanR=0*nanometer , aMaxScanR=200*nanometer , aTbins=35, aMinScanT=0*nanometer, aMaxScanT=500*nanometer ,  aPB=0.2 , aAlpha=20)
# Cfg = ScanConfiguration( 100 , 5*nanometer , 100*nanometer , 
# 												{ 0*nanometer:0.03631079 , 20*nanometer:0.110302441 , 30*nanometer:0.214839819 , 40*nanometer:0.268302465 , 50*nanometer:0.214839819 , 60*nanometer:0.110302441 , 70*nanometer:0.03631079 , 80*nanometer:0.007664194 , 90*nanometer:0.001037236 , 100*nanometer:9.00054E-05 } ,	
# 												35 , 0*nanometer , 200*nanometer , 35, 0*nanometer, 500*nanometer , 0.2 , 20.0 )


# AutoRoi_Scan_SimpleCallback( "1_un_red.csv" , Cfg , Callback )
# AutoRoi_Scan_ToJson( "1_un_red.csv" , Cfg ,  "./{input}/RoI-{roi}/scan.json" )

RunScan( "stem.csv" , ImageJRoI( "ROI/stem.zip" , 25*nanometer ) , Cfg , "./{input}/RoI-{roi}/scan.json" )
