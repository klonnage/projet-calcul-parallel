#include "Input.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

#define READ_INT_VAR( nameInStruct, name, val, out ) \
    if ( name == string( #nameInStruct ) ) {         \
        out.nameInStruct = atoi( val.c_str() );      \
        goto while_end;                              \
    }

#define READ_FLOAT_VAR( nameInStruct, name, val, out ) \
    if ( name == string( #nameInStruct ) ) {           \
        out.nameInStruct = atof( val.c_str() );        \
        goto while_end;                                \
    }

InputData ReadInput( const char *file )
{
    InputData data;
    ifstream  inFile;

    inFile.open( file );
    if ( inFile.fail() ) {
        cerr << "Failed to read input from file " << file << endl;
        return data;
    }

    string       line;
    auto         my_isspace = []( char c ) { return isspace( c ); };
    unsigned int lineNumber = 0;

    while ( getline( inFile, line ) ) {
        auto separatorIndex = line.find_first_of( ':' );
        if ( separatorIndex == string::npos ) {
            cerr << "No separator character (':') at line " << lineNumber << endl;
            continue;
        }
        auto   name              = line.substr( 0, separatorIndex );
        auto   commentStartIndex = line.find_first_of( '#' );
        string val;
        if ( commentStartIndex == string::npos ) { val = line.substr( separatorIndex + 1 ); }
        else {
            val = line.substr( separatorIndex + 1, commentStartIndex - separatorIndex - 1 );
        }

        val.erase( remove_if( val.begin(), val.end(), my_isspace ), val.end() );

        READ_INT_VAR( rowCount, name, val, data );
        READ_INT_VAR( colCount, name, val, data );
        READ_FLOAT_VAR( Lrow, name, val, data );
        READ_FLOAT_VAR( Lcol, name, val, data );
        READ_FLOAT_VAR( D, name, val, data );
        READ_FLOAT_VAR( dt, name, val, data );
        READ_FLOAT_VAR( tMax, name, val, data );
        READ_FLOAT_VAR( beta, name, val, data );
        READ_INT_VAR( kmax, name, val, data );
        READ_FLOAT_VAR( eps, name, val, data );
        READ_INT_VAR( mode, name, val, data );
        READ_INT_VAR( coverage, name, val, data );

        cerr << "Line " << lineNumber << " was not parsed. Read parameter \"" << name << "\" and val \"" << val << "\""
             << endl;

    while_end:
        ++lineNumber;
    }

    inFile.close();

    return data;
}