#ifndef _LINEAR_BUFFERED_READER_H_
#define _LINEAR_BUFFERED_READER_H_

#include <fstream>
#include <sstream>
#include <iostream>


namespace dipl
{

typedef unsigned char byte;

template <unsigned int bufsize>
class LinearBufferedReader
{
	char buffer[bufsize];
	std::string filename;
	std::ios_base::openmode mode;
	std::ifstream infile;
	unsigned int pos;
	unsigned int end; // end within buffer, might be != buf_size at end of file

	std::streampos file_length;

	inline void refill_buffer();

public:
	LinearBufferedReader(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in );

	inline char get();
	void getline(std::string &line);
	inline bool good();
	bool is_open(){ return infile.is_open();}

	std::streampos tellg();

	void seekg(std::streamoff pos);

};


template <unsigned int bufsize>
LinearBufferedReader<bufsize>::LinearBufferedReader(const std::string &filename, std::ios_base::openmode mode)
: filename(filename), mode(mode), infile(filename.c_str(), mode), pos(0), end(0)
{
	infile.seekg(0, std::ios::end);
	file_length = infile.tellg();
	infile.seekg(0, std::ios::beg);
}

template <unsigned int bufsize>
inline void LinearBufferedReader<bufsize>::refill_buffer()
{
	if(pos<end || !good()) return ; // we only refill the buffer if it's been completely read

	infile.read(buffer, bufsize); // try to read buf_size many characters/bytes
	pos = 0;
	end = infile.gcount();
}


template <unsigned int bufsize>
inline char LinearBufferedReader<bufsize>::get()
{
	if(pos>=end) refill_buffer();
	if(!good()) return -1;

	return buffer[pos++];
}

template <unsigned int bufsize>
inline bool LinearBufferedReader<bufsize>::good()
{
	return pos<end || infile.good();
}


template <unsigned int bufsize>
inline void LinearBufferedReader<bufsize>::getline(std::string &line)
{
	std::stringstream ss;
	char c;

	while((c = get()) != '\n' && c!=-1)
	{
		ss << c;
	}

	line = ss.str();
}


template <unsigned int bufsize>
inline std::streampos LinearBufferedReader<bufsize>::tellg()
{
	std::streampos result;

	result = infile.tellg();

	if(result==-1)
	{
		result = file_length;
	}

	result -= (end-pos);

	return result;
}


template <unsigned int bufsize>
inline void LinearBufferedReader<bufsize>::seekg(std::streamoff seek_pos)
{
	if(!infile.good())
	{// reopen the inputfile
		infile.close();
		infile.open(filename.c_str(), mode);
	}

	infile.seekg(seek_pos, std::ios::beg);
	// make sure that the buffer is filled during the next read step
	pos = 0;
	end = 0;
}


} // end of namespace dipl

#endif
