/*
Copyright (c) 2008-2010 Daniel Marbach & Thomas Schaffter

We release this software open source under an MIT license (see below). If this
software was useful for your scientific work, please cite our paper(s) listed
on http://gnw.sourceforge.net.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

$Id: Perturbation.cpp 16 2011-03-06 14:22:36Z jbao $
*/

#include <fstream>
#include <sstream>
#include "Perturbation.h"
#include "Exception.h"
#include "logging/logging.h"
using namespace ::logging;
	
// ============================================================================
// PUBLIC METHODS

/**
 * Constructor
 */
Perturbation::Perturbation() {
	numGenes_ = 0;
    //double dt = GnwSettings::Instance()->getDt();
    //maxt = GnwSettings::Instance()->getMaxtTimeSeries();
    //numTimePoints_ = static_cast<int>(maxt_/dt) + 1;
}

 
Perturbation::Perturbation(GeneNetwork& grn) {
	grn_ = &grn;
	numGenes_ = grn_->getSize();
	numPerturbations_ = -1;
	//perturbations_ = null;		
}

Perturbation::~Perturbation() {

}

// ----------------------------------------------------------------------------

/**
 * Print the perturbations to a file 
 */
void Perturbation::printPerturbations(std::string postfix) {
	
	try { 
		// Filename
		std::string filename = GnwSettings::Instance()->getOutputDirectory() + grn_->getId() + "_" + postfix + "_perturbations.tsv";
		ofstream data_file(filename.c_str()); 
		
		::logging::log::emit<Info>() << "Writing file " << filename.c_str() << ::logging::log::endl;
		
		// Header
		data_file << grn_->getHeader(false);

		// perturbations
		for (int p=0; p<numPerturbations_; p++) {
			for (int i=0; i<numGenes_-1; i++)
				data_file << perturbations_[p][i] << "\t";
            data_file << perturbations_[p][numGenes_-1];
			data_file << std::endl;
		}

		data_file.close();
		

	} catch (exception& e) {
		::logging::log::emit<Info>() << "Perturbation::printPerturbations(): " << e.what() << ::logging::log::endl;
		throw Exception();
	}
}


// ----------------------------------------------------------------------------

/**
 * Load perturbations from the given file.
 */
void Perturbation::loadPerturbations(std::string label) {
	
	try { 
		// Filename
		std::string filename = GnwSettings::Instance()->getOutputDirectory() + grn_->getId() + "_" + label + "_perturbations.tsv";
		ifstream data_file(filename.c_str()); 
		
		::logging::log::emit<Info>() << "Reading file " << filename.c_str() << ::logging::log::endl;
		
		std::string line, entry;
    
    	// header
		getline(data_file, line);
		
		// data matrix
		std::vector< std::vector<double> > data;
		while (getline(data_file, line)){
        	std::istringstream ls(line);
        	std::vector<double> ln;
        	//std::string entry;
        	//getline(ls, entry, '\t');
        	//std::cerr << entry << " ";
        	while (std::getline(ls, entry, '\t')) {
            	ln.push_back(atof(entry.c_str()));
            	//std::cerr << entry << " ";
        	}
        	data.push_back(ln);
        	//std::cerr << std::endl;
    	}
		
		numPerturbations_ = data.size();
		perturbations_ = Mat_DP(numPerturbations_, numGenes_);

		// perturbations
		for (int p=0; p<numPerturbations_; p++) {
			for (int i=0; i<numGenes_; i++)
				perturbations_[p][i] = data[p][i];
		}

		data_file.close();
		

	} catch (exception& e) {
		::logging::log::emit<Info>() << "Perturbation::loadPerturbations(): " << e.what() << ::logging::log::endl;
		throw Exception();
	}
}
	

	
