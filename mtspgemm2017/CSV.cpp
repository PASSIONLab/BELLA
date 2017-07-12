#include "CSV.hpp"

#include <fstream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

bool CSVwriter<IT,FT>::Open()
{
    if(stream.is_open())
	    return true;
	else std::cout << "error file not open" << endl;
}

bool CSVwriter<IT,FT>::Close()
{
    if(stream.is_open())
	    stream.close();

	if(!stream.is_open())
		return true;
	else std::cout << "error file still open" << endl;
}

template <class IT, class FT>
void CSVwriter<IT,FT>::WriteBlock()
{
	for(IT i = skiprows; i < maxrows; ++i)
		stream << i << _delimiter << batchC[i] << _delimiter << batchR[i] << _delimiter << batchV[i] << endl;

	Close();
}

bool CSVreader<IT,FT>::Open()
{
    if(stream.is_open())
	    return true;
	else std::cout << "error file not open" << endl;
}

bool CSVreader<IT,FT>::Close()
{
    if(stream.is_open())
	    stream.close();

	if(!stream.is_open())
		return true;
	else std::cout << "error file still open" << endl;
}

template <class IT, class FT>
void CSVreader<IT,FT>::ReadBlock()
{

}

