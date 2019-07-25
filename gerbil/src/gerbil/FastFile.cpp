/*
 * FastFile.cpp
 *
 *  Created on: 17.09.2015
 *      Author: marius
 */

#include "../../include/gerbil/FastFile.h"

gerbil::FastFile::FastFile(const bfs::path &path) :
		_path(path), _type(ft_unknown), _compr(fc_none), _size(0) {
	std::string fileExtension = _path.extension().string();
	if (fileExtension == ".bz2") {
		_compr = fc_bz2;
		fileExtension = _path.stem().extension().string();
	} else if (fileExtension == ".gz") {
		_compr = fc_gzip;
		fileExtension = _path.stem().extension().string();
	}
	if (fileExtension == ".fastq" || fileExtension == ".fq")
		_type = ft_fastq;
	else if (fileExtension == ".fasta" || fileExtension == ".fa")
		_type = ft_fasta;
	else if (fileExtension == ".ml")
		_type = ft_multiline;
	else
		_type = ft_unknown;

	FILE *p_file = NULL;
	p_file = fopen(_path.c_str(), "rb");
	fseek(p_file, 0, SEEK_END);
	_size = ftell(p_file);
	fclose(p_file);
}

gerbil::FastFile::~FastFile() {

}

const bfs::path& gerbil::FastFile::getPath() const {
	return _path;
}

const gerbil::TFileCompr& gerbil::FastFile::getCompr() const {
	return _compr;
}

const gerbil::TFileType& gerbil::FastFile::getType() const {
	return _type;
}

const uint64_t& gerbil::FastFile::getSize() const {
	return _size;
}
