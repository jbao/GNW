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

$Id: HillGene.cpp 29 2012-01-04 17:06:55Z jbao $
*/

#include <assert.h>
#include <math.h>
#include <sstream>
#include <bitset>
#include <algorithm>
#include "GnwSettings.h"
#include "Edge.h"
#include "HillGene.h"
#include "nr/nr.h"
#include "Exception.h"
#include "logging/logging.h"
//using namespace ::logging;

// ============================================================================
// PUBLIC METHODS

/** Default constructor */
HillGene::HillGene() {
	//alpha_ = NULL;
	//alphaWildType_ = NULL;
	//regulatoryModules_ = null;
}

HillGene::HillGene(std::string label) : Node(label) {
	
}

// ----------------------------------------------------------------------------

HillGene::HillGene(const HillGene& hg) {
	max_ = hg.max_;
	delta_ = hg.delta_;
	maxTranslation_ = hg.maxTranslation_;
	deltaProtein_ = hg.deltaProtein_;
	label_ = hg.label_;
	//regulatoryModules_ = hg.regulatoryModules_;
	alpha_ = hg.alpha_;
	inputGenes_ = hg.inputGenes_;
	outputGenes_ = hg.outputGenes_;
	inputSigns_ = hg.inputSigns_;
	//grn_ = hg.grn_;
    numActivators_ = hg.numActivators_;
    numInhibitors_ = hg.numInhibitors_;
	k_ = hg.k_;
	n_ = hg.n_;
}

// ----------------------------------------------------------------------------

HillGene& HillGene::operator=(const HillGene& rhs) {
	// do the copy
	max_ = rhs.max_;
	delta_ = rhs.delta_;
	maxTranslation_ = rhs.maxTranslation_;
	deltaProtein_ = rhs.deltaProtein_;
	label_ = rhs.label_;
	//regulatoryModules_ = rhs.regulatoryModules_;
	alpha_ = rhs.alpha_;
	inputGenes_ = rhs.inputGenes_;
	outputGenes_ = rhs.outputGenes_;
	inputSigns_ = rhs.inputSigns_;
	//grn_ = hg.grn_;
    numActivators_ = rhs.numActivators_;
    numInhibitors_ = rhs.numInhibitors_;
	k_ = rhs.k_;
	n_ = rhs.n_;
	
	// return the existing object
	return *this;
}

// ----------------------------------------------------------------------------

HillGene::~HillGene() {
	
}

// ----------------------------------------------------------------------------
	
/** Initialization of the gene with the given list of parameters and inputs */
void HillGene::initialization(std::map<std::string, double>& params, const std::vector<std::string>& inputGenes) {
	
	inputGenes_ = inputGenes;
	// Initialize all Gene parameters (max_, delta_, maxTranslation_, and deltaProtein_)
	initializeHalfLife(params);
	// Initialize the gene regulation function and its parameters in the subclass
	initializeRegulation(params);
}

// ----------------------------------------------------------------------------
	
/** Initialization of delta_, max_, deltaProtein_, and maxTranslation_ from the given list of parameters */
void HillGene::initializeHalfLife(std::map<std::string, double>& params) {

	//char buf[256];
	//sprintf(buf, "%f", params["max"]);
	//::logging::log::emit<Info>() << buf << ::logging::log::endl;

	delta_ = params["delta"];
	max_ = params["max"];
	if (GnwSettings::Instance()->getModelTranslation()) {
		deltaProtein_ = params["deltaProtein"];
		maxTranslation_ = params["maxTranslation"];
	}
}

// ----------------------------------------------------------------------------

/** 
 * Initialize the gene regulation function from the given parameters 
 */	
void HillGene::initializeRegulation(std::map<std::string, double>& params) {
	
    // numActivators and numInhibitors
	std::stringstream ss;
    ss.str("");
    ss << "numActivators";
    double numActivators = params[ss.str()];
    
    ss.str("");
    ss << "numInhibitors";
    double numInhibitors = params[ss.str()];
    
    setNumActivators(static_cast<int>(numActivators));
    setNumInhibitors(static_cast<int>(numInhibitors));
		
	// create the modules
	//regulatoryModules_ = new ArrayList<RegulatoryModule>();
	int numInputs = numActivators_ + numInhibitors_;
    Vec_DP k(numInputs);
    Vec_DP n(numInputs);
        
	//std::map<std::string, double>::iterator index = params.find("k_1");
	//int inputCounter = 1;
	/*
	char buf[256];
	for (std::map<std::string, double>::iterator it = params.begin(); it != params.end(); ++it) {
		sprintf(buf, "%s\t%f", it->first.c_str(), it->second);
		::logging::log::emit<Debug>() << buf <<
				::logging::log::endl;
	}
	*/
	for (int i = 0; i < numInputs; ++i) {
		
		//std::vector<std::string> label;
		//label.push_back(inputGenes_.at(numInputs));
		//RegulatoryModule *module = new RegulatoryModule(label);  
		
		// k_ and n_
		//int sumInputs = static_cast<int>(numActivators + numInhibitors);
        ss.str("");
        ss << "k_" << i + 1;
        k[i] = params[ss.str()];
        
        //char buf[256];
        //sprintf(buf, "k_%d\t%f", inputCounter, params[ss.str()]);
        //::logging::log::emit<Info>() << buf <<	::logging::log::endl;
        
        ss.str("");
        ss << "n_" << i + 1;
        n[i] = params[ss.str()];
		
		//regulatoryModules_.push_back(*module);
	}
	assert(numInputs == inputGenes_.size());

    setK(k);
    setN(n);
	
    // states
    int numStates = static_cast<int>(pow(static_cast<double>(2), 
                static_cast<double>(numInputs)));
	alpha_ = Vec_DP(numStates);

	// Set the alpha for all possible states
	// State s is interpreted as a binary number, bit k indicates whether module k
	// is active (1) or inactive (0) in this state.
	//System.out.println(numStates+" "+alpha_.length);
	for (int i=0; i<numStates; i++) {
		ss.str("");
		ss << "a_" << i;
		alpha_[i] = params[ss.str()];
		
		//char buf[256];
		//sprintf(buf, "a_%d\t%f", i, alpha_[i]);
		//::logging::log::emit<Info>() << buf << ::logging::log::endl;
	}
	
}

// ----------------------------------------------------------------------------

/** Return the names and values of the gene's parameters (the first one is always the degradation rate) */
void HillGene::compileParameters(std::map<std::string, double>& params) {
	
	std::stringstream ss;
	
	params.clear();
	params.insert(std::make_pair("delta", delta_));
	params.insert(std::make_pair("max", max_));
	
	if (GnwSettings::Instance()->getModelTranslation()) {
		params.insert(std::make_pair("deltaProtein", deltaProtein_));
		params.insert(std::make_pair("maxTranslation", maxTranslation_));
	}

    // numActivators and numInhibitors
    ss.str("");
    ss << "numActivators";
    params[ss.str()] = static_cast<double>(getNumActivators());
    ss.str("");
    ss << "numInhibitors";
    params[ss.str()] = static_cast<double>(getNumInhibitors());
		
	//::logging::log::emit<Info>() << ::logging::log::dec << regulatoryModules_.size() <<
	//			::logging::log::endl;
	
	// alpha_
	//if (alpha_ != null) {
	//int numStates = static_cast<int>(pow(static_cast<double>(2), static_cast<double>(getRegulatoryModules().size())));
	for (int i=0; i<alpha_.size(); i++) {
		//char buf[256];
		//sprintf(buf, "%f", alpha_[i]);
		//::logging::log::emit<Info>() << buf << ::logging::log::endl;
	
		ss.str("");
		ss << "a_" << i;
		params[ss.str()] = alpha_[i];
	}
	
	// k_
	int counter = 1;
    Vec_DP k = getK();
    Vec_DP n = getN();
    
    for (int j=0; j<k.size(); j++) {
        ss.str("");
        ss << "k_" << counter++;
        params[ss.str()] = k[j];
    }
	
	// n_
	counter = 1;
    for (int j=0; j<n.size(); j++) {
        ss.str("");
        ss << "n_" << counter++;
        params.insert(std::make_pair(ss.str(), n[j]));
    }
	
	//char buf[256];
	//sprintf(buf, "%f", params["max"]);
	//::logging::log::emit<Info>() << buf << ::logging::log::endl;
	//::logging::log::emit<Info>() << ::logging::log::dec << params.size() <<
	//			::logging::log::endl;
}

// ----------------------------------------------------------------------------

/** Perturb the basal activation with this value (subclass defines what exactly that means) */
void HillGene::perturbBasalActivation(double deltaBasalActivation) {
	
	// first, backup the wild-type
	//alphaWildType_ = Vec_DP(alpha_.size());
	//for (int i=0; i<alpha_.size(); i++)
	//	alphaWildType_[i] = alpha_[i];
	
	// Note: what we want to perturb is the basal transcription rate alpha_[0].
	// Remember that alpha_i = alpha_0 + something. Thus, if we perturb alpha_0,
	// it appears in all terms of the vector alpha_
	
	// first, adapt deltaBasalActivation so that alpha_0 is in [0 1]
	if (alpha_[0] + deltaBasalActivation > 1)
		deltaBasalActivation = 1 - alpha_[0];
	else if (alpha_[0] + deltaBasalActivation < 0)
		deltaBasalActivation = 0 - alpha_[0];		
	
	for (int i=0; i<alpha_.size(); i++) {
		alpha_[i] += deltaBasalActivation;
		// truncate to [0 1]
		if (alpha_[i] < 0)
			alpha_[i] = 0;
		else if (alpha_[i] > 1)
			alpha_[i] = 1;
	}
}

// ----------------------------------------------------------------------------

/** Restore the wild-type basal activation */
void HillGene::restoreWildTypeBasalActivation() {
	
	for (int i=0; i<alpha_.size(); i++)
		alpha_[i] = alphaWildType_[i];
}

// ----------------------------------------------------------------------------

/**
 * Compute the production rate of this gene
 * @param geneIndex Index of this gene in the network that it belongs to
 * @param c Current expression levels of all gene of the network
 * @return Production rate of this gene
 */
//inline
double HillGene::computeMRnaProductionRate(Vec_INT& geneIndex, Vec_I_DP& c) {
	
	//std::stringstream mm, ii;
	//mm << regulatoryModules_.size();
	//ii << inputGenes_.size();
	//::logging::log::emit< ::logging::Debug >() << this->getLabel().c_str() << " " << mm.str().c_str() << " " << ii.str().c_str() << ::logging::log::endl;
	
	int numInputs = inputGenes_.size();
	//assert(numModules = inputGenes_.size());
	// The mean activations of the modules
	Vec_DP m(numInputs);
	// Index of the next input
	int nextInput = 0;
	
	// Compute the mean activations
    Vec_DP x(numInputs);
    for (int k=0; k<x.size(); k++) {
        std::stringstream ss, ssx;
        //ss << geneIndex.size();
        //s2 << numModules;
    
        x[k] = c[geneIndex[nextInput++]];
        ss << x[k];
        ssx << numInputs;
        //x[k] = c[ grn_.getIndexOfNode(inputGenes_.at(nextInput++).getLabel()) ];
        //::logging::log::emit< ::logging::Debug >() << getLabel().c_str() << 
        //    " " << inputGenes_[k].c_str() << " " << ss.str().c_str()
        //    <<  ::logging::log::endl;
        
        m[k] = computeActivation(x[k], k);
    }

	// The relative activation of the gene
	double alpha = 0;
	double sum = 0;
	
	for (int i=0; i<alpha_.size(); i++) {
		// binary representation of state i
		// note, leading 0's are not included
		char bits[numInputs];
		dec2bin(i, bits);
		std::stringstream ss;
		ss << bits;
		std::string s = ss.str();
		
		double p = 1; // the probability of being in state i
		
		// if module j is active in this state, multiply with m_j, otherwise with (1-m_j)
		for (int j=0; j<numInputs; j++) {
			std::stringstream ss, s2;
			ss << i;
			s2 << j;
			//::logging::log::emit<Debug>() << "alpha " << ss.str().c_str() << " module " << s2.str().c_str() <<  ::logging::log::endl;
		
			// if s.length()-j-1 smaller than zero, the corresponding module is off (it's one of the leading zeros)
			if (static_cast<int>(s.length())-j-1 >= 0 && s.at(static_cast<int>(s.length())-j-1) == '1')	
				p *= m[j];
			else
				p *= 1 - m[j];
		}
		assert(p >= 0 && p <= 1);
		assert((sum += p) >= 0); // always true, just to compute the sum
		
		std::stringstream s1,s2;
		s1 << i;
		s2 << sum;
		//::logging::log::emit< ::logging::Debug >() << "alpha " << s1.str().c_str() << " sum " << s2.str().c_str() << " bitstring " << s.c_str() << ::logging::log::endl;
		
		alpha += alpha_[i] * p;
	}
	
	assert(sum < 1+1e-6 && sum > 1-1e-6);
	
	return max_ * alpha;
}

// ----------------------------------------------------------------------------
	
/**¬†
 * Compute the mRNA¬†degradation rate of this gene (>= 0)
 * @param geneIndex Index of this gene in the network that it belongs to
 * @param c Current expression levels of all gene of the network
 */
double HillGene::computeMRnaDegradationRate(double x) {

	return delta_*x;
}

// ---------------------------------------------------------------------------

//
// accepts a decimal integer and returns a binary coded string
//
void HillGene::dec2bin(long decimal, char *binary) {
	int k = 0, n = 0;
	int neg_flag = 0;
	int remain;
	int old_decimal; // for test
	char temp[80];

 	// take care of negative input
	if (decimal < 0) {
		decimal = -decimal;
		neg_flag = 1;
	}

	do {
		old_decimal = decimal; // for test
		remain = decimal % 2;
		// whittle down the decimal number
		decimal = decimal / 2;
		// this is a test to show the action
		//printf("%d/2 = %d remainder = %d\n", old_decimal, decimal, remain);
		// converts digit 0 or 1 to character '0' or '1'
		temp[k++] = remain + '0';
	} while (decimal > 0);

 	if (neg_flag)
		temp[k++] = '-'; // add - sign
	else
		//temp[k++] = ' '; // space

 	// reverse the spelling
	while (k >= 0)
		binary[n++] = temp[--k];

 	binary[n-1] = 0; // end with NULL

}

// ----------------------------------------------------------------------------
	
/**
 * Compute the protein production rate of this gene
 */
double HillGene::computeProteinProductionRate(double x) {
	
	return maxTranslation_*x;
}


// ----------------------------------------------------------------------------

/**¬†
 * Compute the protein degradation rate of this gene (>= 0)
 */
double HillGene::computeProteinDegradationRate(double y) {
	
	return deltaProtein_*y;
}

// ============================================================================

/**
 * Set random half-lives for this gene's products. In the non-dimensionalized
 * model, max_ = delta_ and maxTranslation_ = deltaProtein_. We have exponential
 * decay, the half-life and the decay rate are thus related by:
 *    t_1/2 = ln(2) / delta
 *    delta = ln(2) / (t_1/2)
 */
void HillGene::randomKineticInitialization() {
	
	GnwSettings *uni = GnwSettings::Instance();
	delta_ = log(2) / uni->getRandomHalfLife();
	max_ = delta_;
	
	if (uni->getModelTranslation()) {
		deltaProtein_ = log(2) / uni->getRandomHalfLife();
		maxTranslation_ = deltaProtein_;
	}
}

// ============================================================================

/** Random initialization of the HillGene (inputIndexes need to be set already) */
void HillGene::randomStructureInitialization() {
	
	// initialize the number of inputs per module and the modules
	randomInitializationOfStructureAndModules();
	// initialize the gene activations
	randomInitializationOfAlpha();
}

// ----------------------------------------------------------------------------

/** Random initialization the number of inputs per module and the modules */
void HillGene::randomInitializationOfStructureAndModules() {

	GnwSettings *uni = GnwSettings::Instance();
	//Uniform uniform = uni.getUniformDistribution();
	int numInputsOfGene = inputGenes_.size();
	
	// RANDOM INITIALIZATION OF NUMBER OF MODULES AND NUMBER OF INPUTS PER MODULE

	// RANDOM INITIALIZATION OF THE MODULES

	// The inputGenes_ are ordered. Since the first m inputs go to the
	// first module, and the second n inputs to the second module, etc.,
	// we need to randomize the order of the inputs. We put the list in a
	// clone and will add one gene after the other back to the original list.
	std::vector<std::string> inputGenesClone = inputGenes_;
	//inputGenes_.clear();
	
	//regulatoryModules_ = new ArrayList<RegulatoryModule>(numInputsPerModule.size());
	std::vector<Edge::sign> inputSigns = getInputSigns();
	
    // count the input edge types of this module
    int numActivators = 0;
    int numInhibitors = 0;
    int numUnknown = 0;
	
    for (int i = 0; i < numInputsOfGene; ++i) {
		
        Edge::sign type = inputSigns.at(i);
        // for now, we don't model dual regulation and set the sign to unknown
        //if (type == Edge::UNKNOWN)
        //	type = Edge.UNKNOWN;
      
        if (type == Edge::ACTIVATOR)
            numActivators++;
        else if (type == Edge::INHIBITOR)
            numInhibitors++;
        else if (type == Edge::UNKNOWN) {
            numUnknown++;
            if (GnwSettings::Instance()->getRNG()->randInt(1) == 0) {
                inputSigns[i] = Edge::ACTIVATOR;
                numActivators++;
            } else {
                inputSigns[i] = Edge::INHIBITOR;
                numInhibitors++;
            }
        }
        else {
            ::logging::log::emit< ::logging::Error >() << "Invalid edge type!" 
                << ::logging::log::endl;
            throw Exception();			
        }

	}
    setInputSigns(inputSigns);
    setNumActivators(numActivators);
    setNumInhibitors(numInhibitors);
    
    // finally, initialize the numerical parameters of the module
    randomInitializationOfParameters();
		
    assert(numActivators + numInhibitors == numInputsOfGene);
    assert(inputGenes_.size() == numInputsOfGene);
}

// ----------------------------------------------------------------------------

/** Random initialization the alpha parameters for all possible states of the modules */
void HillGene::randomInitializationOfAlpha() {

    GnwSettings *uni = GnwSettings::Instance();
	double weakActivation = uni->getWeakActivation();
	
	int numInputs = getNumInputs();		
	int numStates = static_cast<int>(pow(2.0, static_cast<double>(numInputs)));
	Vec_DP dalpha(numInputs); // the effect on alpha_0 of each module individually
	//alpha_ = new double[numStates];
	
	// set the difference in gene activation due to any module alone
	for (int i=0; i<numInputs; i++) {
		dalpha[i] = uni->getRandomDeltaActivation();
		if (inputSigns_[i] == Edge::INHIBITOR)
			dalpha[i] *= -1;
	}
	
	// Compute the max positive and negative difference in gene activation, i.e.,
	// when all enhancers / repressors are active
	double maxDeltaPositive = 0;
	double maxDeltaNegative = 0;
	
	for (int i=0; i<numInputs; i++) {
		if (dalpha[i] > 0)
			maxDeltaPositive += dalpha[i];
		else
			maxDeltaNegative += dalpha[i];
	}

    std::stringstream ss;
    ss << numStates;
    //::logging::log::emit< ::logging::Debug >() << "\t\t" << ss.str().c_str() <<  ::logging::log::endl;
	alpha_ = Vec_DP(numStates);
	
	// Set alpha_0, the basal transcription rate
	if (numInputs == 0) 			// Case 1: No inputs
		alpha_[0] = 1;
	else if (maxDeltaPositive == 0) // Case 2: There are only repressors
		alpha_[0] = 1;
	else if (maxDeltaNegative == 0) // Case 3: There are only enhancers
		alpha_[0] = uni->getRandomLowBasalRate();
	else 							// Case 4: There are enhancers and repressors
		alpha_[0] = uni->getRandomMediumBasalRate();
	
	// make sure that the activation goes at least up to 1 in the maximally activated state
	// (if there is at least one activator)
	if (maxDeltaPositive > 0 && alpha_[0] + maxDeltaPositive < 1) {
		// find the smallest positive dalpha
		double minPos = 1;
		int indexMinPos = -1;
		
		for (int i=0; i<numInputs; i++) {
			if (dalpha[i] > 0 && dalpha[i] < minPos) {
				minPos = dalpha[i];
				indexMinPos = i;
			}
		}
		// increase the smallest dalpha so that: alpha_0 + maxDeltaPositive = 1
		dalpha[indexMinPos] += 1 - alpha_[0] - maxDeltaPositive;
	}
	
	// make sure that the activation falls within [0 weakActivation] in the maximally repressed state
	// (if there is a at least a repressor)
	if (maxDeltaNegative < 0 && alpha_[0] + maxDeltaNegative > weakActivation) {
		// find the weakest negative dalpha
		double minPos = -1;
		int indexMinPos = -1;
		
		for (int i=0; i<numInputs; i++) {
			if (dalpha[i] < 0 && dalpha[i] > minPos) {
				minPos = dalpha[i];
				indexMinPos = i;
			}
		}
		// increase the weakest dalpha so that: (alpha_[0] + maxDeltaNegative) in [0 weakActivation]
		dalpha[indexMinPos] += - alpha_[0] - maxDeltaNegative + uni->getRandomLowBasalRate();
			//weakActivation - alpha_[0] - maxDeltaNegative - uni.getRandomLowBasalRate();
	}
	
    ss.str("");
    ss << numStates;
    //::logging::log::emit< ::logging::Debug >() << "\t\t" << ss.str().c_str() 
    //    << ::logging::log::endl;
	
    // Set the alpha for all possible states
	// State s is interpreted as a binary number, bit k indicates whether module k
	// is active (1) or inactive (0) in this state. State 0 (alpha_0) has already
	// been set
	for (int i=1; i<numStates; i++) {
		alpha_[i] = alpha_[0];
		// note, leading 0's are not included
		char bits[numInputs];
		dec2bin(i, bits);
		std::stringstream ss;
		ss << bits;
		std::string s = ss.str();
		
		// if module j is active in this state, add its effect
		for (int j=0; j<numInputs; j++) {
			// if s.length()-j-1 smaller than zero, the corresponding module is off (it's one of the leading zeros)
			if (static_cast<int>(s.length())-j-1 >= 0 && s.at(static_cast<int>(s.length())-j-1) == '1')	
				alpha_[i] +=  dalpha[j];
		}
		
		// truncate to [0 1]
		if (alpha_[i] < 0)
			alpha_[i] = 0;
		else if (alpha_[i] > 1)
			alpha_[i] = 1;
	}
}

// ----------------------------------------------------------------------------

/**
 * Set the type (enhancing, inhibiting, dual, or unknown) of all input edges
 * of this gene according to the dynamically model.
 */
/*
void HillGene::setInputEdgeTypesAccordingToDynamicalModel() {
	
	// inputSigns is a vector of boolean indicating for each input the sign of its edge
	// (whether it has an enhancing or inhibitory effect on this gene)
	std::vector<bool> inputSigns;
	for (unsigned int i=0; i<regulatoryModules_.size(); i++) {
		std::vector<bool> signs = regulatoryModules_.at(i).getEdgeSigns();
		for (unsigned int j = 0; j < signs.size(); j++) 
			inputSigns.push_back(signs[j]);
	}
	
	// a sign must be specified for each input
	int numInputs = inputSigns.size();
	assert(numInputs == inputGenes_.size());
	
	// set the edges
	std::vector<Edge::sign> inputSigns = getInputSigns();
	for (int i=0; i<numInputs; i++) {
		if (inputSigns.at(i)) {
			grn_->getEdge(inputGenes_.at(i), *this, e);
			e.setSign(Edge::ACTIVATOR);
		}
		else {
			grn_->getEdge(inputGenes_.at(i), *this, e);
			e.setSign(Edge::INHIBITOR);
		}
	}
	
}
*/

// ----------------------------------------------------------------------------

void HillGene::pruneInput(HillGene& gene) {
	std::vector<std::string>::iterator it = inputGenes_.begin(); 
    for (; it != inputGenes_.end(); ++it) {
        if (*it == gene.getLabel())
            break;
    }
	assert(it != inputGenes_.end());
	
	inputGenes_.erase(it);
	inputSigns_.erase(inputSigns_.begin() + (it - inputGenes_.begin()));
	//regulatoryModules_.erase(regulatoryModules_.begin() + (it - inputGenes_.begin()));
	
	std::stringstream inputSize, moduleSize;
	inputSize << inputGenes_.size();
	//::logging::log::emit< ::logging::Debug >() << this->getLabel().c_str() << " input size " << inputSize.str().c_str() << " module size " << moduleSize.str().c_str() <<  ::logging::log::endl;
	
	// remove alpha
	std::vector<double> alpha_vec;
	int idx = it - inputGenes_.begin();
	//int counter = 0;
	for (int i = 0; i < alpha_.size(); ++i) {
		char bits[inputGenes_.size()];
		dec2bin(i, bits);
		std::stringstream ss;
		ss << bits;
		std::string s = ss.str();
		
		if (static_cast<int>(s.length())-idx-1 < 0 || s.at(static_cast<int>(s.length())-idx-1) == '0') {
			alpha_vec.push_back(alpha_[i]);
			//counter++;
		}
	}
	
	Vec_DP alpha_clone(alpha_vec.size());
	for (int i = 0; i < alpha_clone.size(); ++i)
		alpha_clone[i] = alpha_vec[i];
	alpha_ = alpha_clone;
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
void HillGene::randomInitializationOfParameters() {
	
	GnwSettings *uni = GnwSettings::Instance();
	int numInputs = numActivators_ + numInhibitors_;
	
	k_ = Vec_DP(numInputs);
	n_ = Vec_DP(numInputs);
	for (int i=0; i<numInputs; i++) {
		k_[i] = uni->getRandomK();
		n_[i] = uni->getRandomN();
	}
}

// ----------------------------------------------------------------------------

/** 
 * Compute the activation of this module as a function of the regulator concentrations x.
 * x must be ordered, first come the numActivators_ activators, then the deactivators.
 */
//inline
double HillGene::computeActivation(double x, int index) {

	// define xi_i := (x_i/k_i)^n_i
    assert(x >= 0.0);
    double xi = pow(x / k_[index], n_[index]);
    std::stringstream ssx,ssk,ssn;
    ssx << x;
    ssk << k_[index];
    ssn << n_[index];
    //::logging::log::emit< ::logging::Debug >() << "x=" << ssx.str().c_str() <<  
    //    " k=" << ssk.str().c_str() << " n=" << ssn.str().c_str() <<
    //    ::logging::log::endl;
	
	// compute the numerator
	double multiplyActivators = 1;
    if (inputSigns_[index] == Edge::ACTIVATOR)
		multiplyActivators *= xi;
	double numerator = multiplyActivators;
	
	// compute the partition function
	double denominator = 1;
	
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
    denominator *=  (xi + 1);
	
    double activation = numerator / denominator;
	
	std::stringstream num, denom;
	num << numerator;
	denom << denominator;
	//::logging::log::emit< ::logging::Debug >() << "numerator = " << 
    //    num.str().c_str() << " denominator = " << denom.str().c_str() << 
    //    ::logging::log::endl;
	
	assert(activation >= 0 && activation <= 1);
	return activation;
}

