#ifndef PERTURBATIONMULTIFACTORIAL_H
#define PERTURBATIONMULTIFACTORIAL_H

//#include "nr/nr.h"
#include "GeneNetwork.h"
#include "Perturbation.h"

class PerturbationMultifactorial : public Perturbation {

public:
	PerturbationMultifactorial();
	PerturbationMultifactorial(GeneNetwork *grn);
	~PerturbationMultifactorial();
	PerturbationMultifactorial(const PerturbationMultifactorial& multi);
	PerturbationMultifactorial& operator=(const PerturbationMultifactorial& rhs);
	
	void multifactorialAllGenesWeak(int numPerturbations);
	void multifactorialAllGenesMax(int numPerturbations);
	void multifactorialByRank(int numPerturbations);
	void multifactorialStrong(int numPerturbations);
	void multifactorialHub(int numPerturbations);
	void perturbSingleGene(int numPerturbations, std::string gene);
	void applyPerturbation(int k);
	void restoreWildType();
	
private:
	/** Number of genes */
	//int numGenes_;
	/** The probability that a gene is perturbed (for strong multifactorial perturbations) */
	double perturbationProbability_;
	/** The standard deviation for the weak multifactorial perturbations */
	double stdev_;
	/** The gene network to which the perturbations are being applied */
	//GeneNetwork *grn_;
	/** The wild-type (so that it can be restored after applying perturbations) */
	//Vec_DP wildType_;
	/** The number of different multifactorial perturbations */
	//int numPerturbations_;
	/** perturbations_(k, i) is the perturbed value of m_i in perturbation k */
	//Mat_DP perturbations_;
	
	void saveWildType();

};

#endif
