import os,re,sys

KnownFiles = {}
def FindFile( name ):
	global KnownFiles
	if name in KnownFiles: return KnownFiles[name]

	for (dirpath, dirnames, filenames) in os.walk( '.' ):
		for filename in filenames:
			if filename == name:
				filename = os.path.join( dirpath , filename )
				KnownFiles[name] = filename
				return filename


_ , URL , branch = sys.argv
URL = URL.strip("/")

IsTex = re.compile( r'.*\.tex' )
Pattern = re.compile( r'Definition at line ([0-9]*) of file (.*)\.' )

for (dirpath, dirnames, filenames) in os.walk( '.doxygen/latex/' ):
	for file in filenames:
		if IsTex.match( file ):
			file = os.path.join( dirpath , file )
			with open( file , "r") as file1 , open( "temp.tex" , "w" ) as file2:
				for line in file1:
					m = Pattern.search( line )
					if not m is None : 
						filename =  FindFile( m.group(2).replace( "\\-" , "" ) )
						filename = "{}/{}/{}#L{}".format( URL , branch , filename[2:] ,  m.group(1) )
						line = "\\href{{{}}}{{{}}}".format( filename , m.group(0) )
						print( line )
					file2.write( line )
		
			os.rename( "temp.tex" , file )
		
