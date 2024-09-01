//----------------------------------------------------------
// Copyright 2019 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef DETECT_H
#define DETECT_H

#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include "reads.h"
#include "tensor.h"
#include "error_handling.h"
#include "../htslib/htslib/hts.h"
#include "../htslib/htslib/sam.h"

struct HMMdetection{

	public:
		std::map<unsigned int, std::pair<double,double>> refposToLikelihood;
		std::string stdout;
};


class OutputWriter {
	public:
		virtual ~OutputWriter() {}
		virtual void open(const std::string& filename) = 0;
		virtual void write(const DNAscent::read &r) = 0;
		virtual void close() = 0;
		virtual void writeHeader_HR(const std::string& header_str) = 0;
		virtual void writeHeader_sam(bam_hdr_t *sam_hdr) = 0;		
};


class HumanReadableWriter : public OutputWriter {

	private:
		std::ofstream file;
		
	public:
		void open(const std::string& filename) override {
			file.open(filename);
			if (not file.is_open()) throw IOerror(filename);
		}

		void write(const DNAscent::read &r) override {
			if (file.is_open()) {
				file << r.humanReadable_detectOut;
			}
		}

		void close() override {
			if (file.is_open()) {
				file.close();
			}
		}
		void writeHeader_HR(const std::string& header_str){
			if (file.is_open()) {
				file << header_str;
			}
		}
		void writeHeader_sam(bam_hdr_t *sam_hdr){}
};


class SamWriter : public OutputWriter {

	private:
		htsFile *file;
		bam_hdr_t *header;
		std::string fn;

	public:
		void open(const std::string& filename) override {
			file = sam_open(filename.c_str(), "wb");
			if (file == nullptr) throw IOerror(filename);
			fn = filename;
		}

		void write(const DNAscent::read &r) override {
			if (file){
				
				int status = sam_write1(file, header, r.record);
				if (status < 0) throw BamWriteError(fn);				
			}
		}
		void close() override {
			if (file){
				sam_close(file);
			}
		}
		void writeHeader_HR(const std::string& header_str){ }
		void writeHeader_sam(bam_hdr_t *sam_hdr){
			header = sam_hdr;
			if (file) {
				int status = sam_hdr_write(file, sam_hdr);
				if (status < 0) throw BamWriteError(fn);
			}
		}		
};


enum class OutputFormat { HumanReadable, Sam };

class OutputWriterFactory {
	public:
		static std::unique_ptr<OutputWriter> createWriter(OutputFormat format) {
			switch (format) {
				case OutputFormat::HumanReadable:
					return std::make_unique<HumanReadableWriter>();
				case OutputFormat::Sam:
					return std::make_unique<SamWriter>();
				default:
					return nullptr;
			}
		}
};


int detect_main( int argc, char** argv );
std::vector< unsigned int > getPOIs( std::string &, int );
double sequenceProbability( std::vector <double> &, std::string &, size_t, bool, PoreParameters, size_t, size_t );
void runCNN(DNAscent::read & , std::shared_ptr<ModelSession> , std::vector<TF_Output>, bool );
HMMdetection llAcrossRead( DNAscent::read &, unsigned int );

#endif
