#include "Input.hpp"

#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>

using namespace std;

#define READ_INT_VAR(nameInStruct, name, val, out) \
if (name == string(#nameInStruct)) \
{\
    out##.##nameInStruct = atoi(val.c_str()); \
    goto while_end; \
}

#define READ_FLOAT_VAR(nameInStruct, name, val, out) \
if (name == string(#nameInStruct)) \
{\
    out##.##nameInStruct = atof(val.c_str()); \
    goto while_end; \
}

InputData ReadInput(const char *file)
{

    InputData data;
    ifstream inFile;

    inFile.open(file);
    if (inFile.fail())
    {
        cerr << "Failed to read input from file " << file << endl;
        return data;
    }

    string line;

    unsigned int lineNumber = 0;
    while( getline(inFile, line) )
    {
        auto separatorIndex = line.find_first_of(':');
        if (separatorIndex == string::npos)
        {
            cerr << "No separator character (':') at line " << lineNumber << endl;
            continue;
        }
        auto name = line.substr(0, separatorIndex + 1);
        auto val = line.substr(separatorIndex + 1);

        READ_INT_VAR(Nx, name, val, data);
        READ_INT_VAR(Ny, name, val, data);
        READ_FLOAT_VAR(Lx, name, val, data);
        READ_FLOAT_VAR(Ly, name, val, data);
        READ_FLOAT_VAR(D, name, val, data);
        READ_FLOAT_VAR(dt, name, val, data);
        READ_FLOAT_VAR(tMax, name, val, data);
        READ_FLOAT_VAR(beta, name, val, data);
        READ_INT_VAR(kmax, name, val, data);

while_end:
        ++lineNumber;
    }

    inFile.close();

    return data;
}