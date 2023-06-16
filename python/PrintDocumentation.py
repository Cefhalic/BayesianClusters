from inspect import getmembers
import BayesianClustering 

members = getmembers( BayesianClustering )

print()
print( "Functions:" )
for i,j in members: 
  if str( type( j ) ) == "<class 'Boost.Python.function'>":
    doc = [ x.strip() for x in j.__doc__.split( "\n" ) ]
    print( f"  {doc[1][:-2]}\n    {doc[2]}\n" ) 

print( "Classes:" )
for i,j in members: 
  if isinstance( j , type ):
    print( f"  {i}\n    {j.__doc__}" )
    for k in [ x.strip()[:-2] for x in j.__init__.__doc__.split( "\n" ) if "None" in x ]:
      print( f"      {k}" )

print()
