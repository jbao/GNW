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

$Id: RegulatoryModule.cpp 13 2011-02-28 14:29:38Z jbao $
*/
    
#include <cmath>
#include <cassert>
#include <sstream>
#include "GnwSettings.h"
#include "RegulatoryModule.h"
#include "logging/logging.h"
using namespace ::logging;
    
// ============================================================================
// PUBLIC METHODS

/**
 * Constructor
 */
RegulatoryModule::RegulatoryModule() {
	isEnhancer_ = true;
	bindsAsComplex_ = false;
	numActivators_ = -1;
	numDeactivators_ = -1;
}

// ----------------------------------------------------------------------------

RegulatoryModule::RegulatoryModule(const std::vector<std::string>& labels) {
	isEnhancer_ = true;
	bindsAsComplex_ = false;
	numActivators_ = -1;
	numDeactivators_ = -1;
	labels_ = labels;
}

// ----------------------------------------------------------------------------

RegulatoryModule::RegulatoryModule(const RegulatoryModule& rm) {
	isEnhancer_ = rm.isEnhancer_;
	bindsAsComplex_ = rm.bindsAsComplex_;
	numActivators_ = rm.numActivators_;
	numDeactivators_ = rm.numDeactivators_;
	k_ = rm.k_;
	n_ = rm.n_;
}

// ----------------------------------------------------------------------------

RegulatoryModule& RegulatoryModule::operator=(const RegulatoryModule& rhs) {
	// do the copy
	isEnhancer_ = rhs.isEnhancer_;
	bindsAsComplex_ = rhs.bindsAsComplex_;
	numActivators_ = rhs.numActivators_;
	numDeactivators_ = rhs.numDeactivators_;
	k_ = rhs.k_;
	n_ = rhs.n_;
	
	// return the existing object
	return *this;
}

// ----------------------------------------------------------------------------
/*
void RegulatoryModule::getK(double *res, int len) { 
	res = new double[len]; 
	for (int i = 0; i < len; ++i)
		res[i] = k_[i]; 
}
	
void RegulatoryModule::setK(double k[], int len) { 
	double k_[len];
	for (int i = 0; i < len; ++i)
		k_[i] = k[i]; 
}
	
void RegulatoryModule::getN(double *res, int len) { 
	res = new double[len];
	for (int i = 0; i < len; ++i)
		res[i] = n_[i]; 
}
	
void RegulatoryModule::setN(double n[], int len) { 
	double n_[len];
	for (int i = 0; i < len; ++i)
		n_[i] = n[i]; 
}
*/

// ----------------------------------------------------------------------------

/** 
 * Compute the activation of this module as a function of the regulator concentrations x.
 * x must be ordered, first come the numActivators_ activators, then the deactivators.
 */
//inline
double RegulatoryModule::computeActivation(Vec_I_DP& x) {
	
	int numInputs = numActivators_ + numDeactivators_;
	assert(x.size() == numInputs);
	
	// define xi_i := (x_i/k_i)^n_i
	Vec_DP xi(numInputs);
	for (int i=0; i<numInputs; i++) {
		assert(x[i] >= 0.0);
		xi[i] = pow(x[i] / k_[i], n_[i]);
		//std::stringstream ss;
		//ss << x[i];
		//::logging::log::emit<Debug>() << "xi = " << ss.str().c_str() <<  ::logging::log::endl;
	}
	
	// compute the numerator
	double multiplyActivators = 1;
	for (int i=0; i<numActivators_; i++)
		multiplyActivators *= xi[i];
	double numerator = multiplyActivators;
	
	// compute the partition function
	double denominator = 1;
	
	if (bindsAsComplex_) {
		// if it is a complex, there are three states of the module,
		// either the activated complex is bound, the deactivated complex is bound, or it is not bound
		
		// activated complex bound
		denominator += multiplyActivators;
		
		if (numDeactivators_ > 0) { // BUG FIXED: this if was not here in v1
			double multiplyAllInputs = multiplyActivators;
			for (int j=numActivators_; j<numInputs; j++)
				multiplyAllInputs *=  xi[j];
			denominator += multiplyAllInputs;
		}
		
	} else {
		// I was actually computing (x0+1)(x1+1) ... = 1 + x0 + x1 + x0x1 + ... in the latter form!
		/*int numStates = (int)Math.round(Math.pow(2, numInputs));
		for (int i=1; i<numStates; i++) {
			String s = Integer.toBinaryString(i); // note, leading 0's are not included
			int slength = s.length();
			
			// if input j is active in this state i, add it to the term
			double term = 1;
			for (int j=0; j<numInputs; j++) {
				// if s.length()-j-1 smaller than zero, the corresponding entry is one of the leading zeros
				if (slength-j-1 >= 0 && s.charAt(slength-j-1) == '1')	
					term *=  xi[j];
			}
			denominator += term;
		}*/
		// ok, this *is* arguably faster
		for (int i=0; i<numInputs; i++)
			denominator *=  (xi[i] + 1);
	}
	double activation = numerator / denominator;
	
	std::stringstream num, denom;
	num << numerator;
	denom << denominator;
	//::logging::log::emit<Debug>() << "numerator = " << num.str().c_str() << " denominator = " << denom.str().c_str() << ::logging::log::endl;
	
	assert(activation >= 0 && activation <= 1);
	return activation;
}

// ----------------------------------------------------------------------------

/**
 * Random initialization of parameter values. isEnhancer_, numActivators_, and
 * numDeactivators_ must already be set. If the module is a complex, the k_'s 
 * are interpreted as follows:
 * The activation is half-maximal if x1 = k1, x2 = k2, ..., kN = kN, i.e.,
 * k_complex = k1^n1 * k2^n2 * ... * kN^nN. In the computation of the activation,
 * there will be only two states, the complex is bound or not.
 */
void RegulatoryModule::randomInitializationOfParameters() {
	
	GnwSettings *uni = GnwSettings::Instance();
	int numInputs = numActivators_ + numDeactivators_;
	
	k_ = Vec_DP(numInputs);
	n_ = Vec_DP(numInputs);
	for (int i=0; i<numInputs; i++) {
		k_[i] = uni->getRandomK();
		n_[i] = uni->getRandomN();
	}
}

// ----------------------------------------------------------------------------
	
/**
 * Returns an array of booleans, where element i is true if input i is an enhancer
 * and false if it is an inhibitor of the gene (not of the regulatory module).
 * Note, an enhancing input can be either an activator of an enhancing module, or
 * a deactivator of a repressor module.
 */
std::vector<bool>& RegulatoryModule::getEdgeSigns() {

	std::vector<bool> edgeSigns;
	if (isEnhancer_) {
		for (int i=0; i<numActivators_; i++)
			edgeSigns.push_back(true);
		for (int i=0; i<numDeactivators_; i++)
			edgeSigns.push_back(false);
	} else {
		for (int i=0; i<numActivators_; i++)
			edgeSigns.push_back(false);
		for (int i=0; i<numDeactivators_; i++)
			edgeSigns.push_back(true);	
	}
	
	return edgeSigns;
}
