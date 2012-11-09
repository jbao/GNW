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
*/

#include "GnwSettings.h"
#include "PerturbationMultifactorial.h"
#include "logging/logging.h"
using namespace ::logging;

/* function object to check the value of a map element
 */
template <class K, class V>
class value_equals {
private:
    V value;
    K key;
public:
    // constructor (initialize value to compare with)
    value_equals (const K& k) : key(k) {}
    // comparison
    bool operator() (pair<K, const V> elem) {
        return elem.first == key;
    }
};

// ============================================================================
// PUBLIC METHODS

PerturbationMultifactorial::PerturbationMultifactorial() : Perturbation() {
	perturbationProbability_ = GnwSettings::Instance()->getPerturbationProbability();
	stdev_ = GnwSettings::Instance()->getMultifactorialStdev();
	numGenes_ = 0;
}

PerturbationMultifactorial::PerturbationMultifactorial(GeneNetwork *grn) : Perturbation() {
	perturbationProbability_ = GnwSettings::Instance()->getPerturbationProbability();
	stdev_ = GnwSettings::Instance()->getMultifactorialStdev();
	grn_ = grn;
	numGenes_ = grn_->getSize();
}

PerturbationMultifactorial::~PerturbationMultifactorial() {

}

// ---------------------------------------------------------------------------

PerturbationMultifactorial::PerturbationMultifactorial(const PerturbationMultifactorial& multi) {
	numGenes_ = multi.numGenes_;
	perturbationProbability_ = multi.perturbationProbability_;
	stdev_ = multi.stdev_;
	grn_ = multi.grn_;
	wildType_ = multi.wildType_;
	numPerturbations_ = multi.numPerturbations_;
	perturbations_ = multi.perturbations_;
}

// ----------------------------------------------------------------------------

PerturbationMultifactorial& PerturbationMultifactorial::operator=(const PerturbationMultifactorial& rhs) {
	numGenes_ = rhs.numGenes_;
	perturbationProbability_ = rhs.perturbationProbability_;
	stdev_ = rhs.stdev_;
	grn_ = rhs.grn_;
	wildType_ = rhs.wildType_;
	numPerturbations_ = rhs.numPerturbations_;
	perturbations_ = rhs.perturbations_;
	
	return *this;
}

// ----------------------------------------------------------------------------

/**
 * The basal transcription rate alpha_0 of all genes is sampled from a normal
 * distribution with mean alpha_0 and standard deviation CV_.
 */
void PerturbationMultifactorial::multifactorialAllGenesWeak(int numPerturbations) {
	
	saveWildType();
	
	numPerturbations_ = numPerturbations;
	perturbations_ = Mat_DP(numPerturbations_, numGenes_);
	int numGenes_ = grn_->getSize();
	
	// generate perturbations
	for (int p=0; p<numPerturbations_; p++) {
		for (int g=0; g<numGenes_; g++) {
			perturbations_[p][g] = GnwSettings::Instance()->getNormalDistributionNumber(0, stdev_);
		}
	}
}

// ----------------------------------------------------------------------------

/**
 * The basal transcription rate alpha_0 of all genes is perturbed to the maximal
 * level 1.
 */
void PerturbationMultifactorial::multifactorialAllGenesMax(int numPerturbations) {
	
	saveWildType();
	
	numPerturbations_ = numPerturbations;
	perturbations_ = Mat_DP(numPerturbations_, numGenes_);
	int numGenes_ = grn_->getSize();
	
	// generate perturbations
	for (int p=0; p<numPerturbations_; p++) {
		for (int g=0; g<numGenes_; g++) {
            double alpha0 = grn_->getNodes().at(g).getAlpha()[0];
			perturbations_[p][g] = 1 - alpha0;
		}
	}
}

// ----------------------------------------------------------------------------

/**
 * The basal transcription rate alpha_0 of all genes is sampled from a normal
 * distribution with mean alpha_0 and standard deviation CV_.
 */
void PerturbationMultifactorial::multifactorialByRank(int numPerturbations) {
	
	saveWildType();
	
	numPerturbations_ = numPerturbations;
	perturbations_ = Mat_DP(numPerturbations_, numGenes_);
	int numGenes_ = grn_->getSize();
	
	//std::map<std::string,int> rank = grn_->getRank();
	
	// generate perturbations
	for (int p=0; p<numPerturbations_; p++) {
		for (int g=0; g<numGenes_; g++) {
			//perturbations_[p][g] = -0.01*rank[grn_->getNodes().at(g).getLabel()]+1+GnwSettings::Instance()->getNormalDistributionNumber(0, stdev_);
		}
	}
}

// ----------------------------------------------------------------------------

/**
 * A perturbation is applied to a gene with probability perturbationProbability_.
 * The basal activations alpha_0 are sampled from a uniform distribution in [0 1]
 * independently of their unperturbed value. 
 */
void PerturbationMultifactorial::multifactorialStrong(int numPerturbations) {
	
	saveWildType();
	
	numPerturbations_ = numPerturbations;
	perturbations_ = Mat_DP(numPerturbations_, numGenes_);

	//Uniform uniform = GnwSettings.getInstance().getUniformDistribution();
	int numGenes_ = grn_->getSize();

	// generate perturbations
	for (int p=0; p<numPerturbations_; p++) {
		for (int g=0; g<numGenes_; g++) {
			if (GnwSettings::Instance()->getUniformDistributionNumber(1) < perturbationProbability_) {
				// the new basal activation
				double delta = GnwSettings::Instance()->getUniformDistributionNumber(1);
				// but actually we need the difference to the wild-type
				delta = delta - grn_->getNodes().at(g).getBasalActivation();
				perturbations_[p][g] = delta;
			} else
				perturbations_[p][g] = 0;
		}
	}
}

// ----------------------------------------------------------------------------

/**
 * A maximal perturbation (up to 1) is applied to all genes from a certain
 * fraction with a sorted ascending degree 
 */
void PerturbationMultifactorial::multifactorialHub(int numPerturbations) {
	
	saveWildType();
	
	numPerturbations_ = numPerturbations;
	perturbations_ = Mat_DP(numPerturbations_, numGenes_);

	//Uniform uniform = GnwSettings.getInstance().getUniformDistribution();
	int numGenes_ = grn_->getSize();
    int low = GnwSettings::Instance()->getLowerPerturbationFraction() * numGenes_;
    int high = GnwSettings::Instance()->getUpperPerturbationFraction() * numGenes_;

	// generate perturbations
	for (int p=0; p<numPerturbations_; p++) {
		for (int g=0; g<numGenes_; g++) {
            std::vector< std::pair<std::string,int> >::iterator it = 
                find_if(grn_->degreesVec.begin()+low, grn_->degreesVec.begin()+high, 
                        value_equals<std::string,int>(grn_->getNodes().at(g).getLabel()));
            if (it != grn_->degreesVec.begin()+high) {
	            //::logging::log::emit<Debug>() << grn_->getNodes().at(g).getLabel().c_str()
                //    << ::logging::log::endl;
                double alpha0 = grn_->getNodes().at(g).getAlpha()[0];
			    perturbations_[p][g] = 0 - alpha0;
            }
            else
                perturbations_[p][g] = 0;
        }
	}
}

// ----------------------------------------------------------------------------

/**
 * A random perturbation is applied to a single gene
 */
void PerturbationMultifactorial::perturbSingleGene(int numPerturbations, std::string gene) {
	
	saveWildType();
	
	numPerturbations_ = numPerturbations;
	perturbations_ = Mat_DP(numPerturbations_, numGenes_);

	//Uniform uniform = GnwSettings.getInstance().getUniformDistribution();
	int numGenes_ = grn_->getSize();

	// generate perturbations
	for (int p=0; p<numPerturbations_; p++) {
	    int idx = grn_->getIndexOfNode(gene);
		for (int g=0; g<numGenes_; g++) {
            if (g == idx) {
                double delta = GnwSettings::Instance()->getUniformDistributionNumber(1);
                perturbations_[p][idx] = delta;
            }
            else
                perturbations_[p][g] = 0;
        }
    }
}

// ----------------------------------------------------------------------------

/** Apply the k'th perturbation to the grn_ */
void PerturbationMultifactorial::applyPerturbation(int k) {
	
	int numGenes_ = grn_->getSize();
	for (int i=0; i<numGenes_; i++)
		grn_->getNodes().at(i).perturbBasalActivation( perturbations_[k][i] );
}


// ----------------------------------------------------------------------------

/** Save the wild-type of the network grn_ in wildType_ */
void PerturbationMultifactorial::saveWildType() {
	
	int numGenes_ = grn_->getSize();
	wildType_ = Vec_DP(numGenes_);
	for (int i=0; i<numGenes_; i++)
		wildType_[i] = 0;
}


// ----------------------------------------------------------------------------

/** Restore the values before perturbations were applied */
void PerturbationMultifactorial::restoreWildType() {
	
	int numGenes_ = grn_->getSize();
	for (int i=0; i<numGenes_; i++)
		grn_->getNodes().at(i).restoreWildTypeBasalActivation();
}
	

	
	// ----------------------------------------------------------------------------

	/**
	 * The max transcription rates of all genes are sampled from a normal
	 * distribution with mean m_i and standard deviation m_i*CV_.
	 */
	/*
	@Deprecated
	public void multifactorialMax(int numPerturbations) {
		
		perturbBasalActivation_ = false;
		saveWildType();
		
		numPerturbations_ = numPerturbations;
		perturbations_ = new DenseDoubleMatrix2D(numPerturbations_, numGenes_);

		// Generates random numbers from a normal distribution, resamples if value < min=0
		// The max is set to Double.MAX_VALUE. The mean and stdev will have to be set each time 
		RandomParameterGaussian rand = new RandomParameterGaussian(0, Double.MAX_VALUE, 0, 0, false);
		
		// generate perturbations
		for (int p=0; p<numPerturbations_; p++) {
			for (int g=0; g<numGenes_; g++) {
				double m = wildType_.get(g);
				double stdev = CV_*m;
				rand.setMeanStdev(m, stdev);
				perturbations_.set(p, g, rand.getRandomValue());
			}
		}
	}
	*/

