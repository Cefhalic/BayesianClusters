

/* ===== Cluster sources ===== */
#include "Cluster_GlobalVars.hpp"
#include "Cluster_ConfigFile.hpp"


void LoadConfigFile( const std::string& aFilename  )
{
  auto f = fopen( aFilename.c_str() , "r");
  if ( f == NULL ) throw std::runtime_error( "File is not available" );

  // fseek(f, 0, SEEK_END); // seek to end of file
  // auto lSize = ftell(f); // get current file pointer
  // fseek(f, 0, SEEK_SET); // seek back to beginning of file

  // {
  //   char ch[256];
  //   char* lPtr( ch );
  //   ProgressBar lProgressBar( "Reading File" , lSize );

  //   auto ReadUntil = [ &ch , &f , &lPtr , &lProgressBar ]( const char& aChar ){
  //     lPtr = ch;
  //     while ( ( *lPtr = fgetc(f)) != EOF )
  //     {
  //       lProgressBar++;
  //       if( *lPtr == aChar ) return;
  //       lPtr++;
  //     }
  //   };

  //   while( true )
  //   {
  //     ReadUntil( '=' ); 

  //     ReadUntil( '\n' );
  //     double s = m * strtod( ch , &lPtr );      

  //     if( fabs(x) < 1 and fabs(y) < 1 ) lData.emplace_back( x , y , s );
  //   }
  // }
  fclose(f);

  // std::cout << "Read " << lData.size() << " points" << std::endl;

  // return lData;
}
