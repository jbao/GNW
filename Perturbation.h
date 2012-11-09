// $Id: Perturbation.h 26 2011-08-26 14:33:16Z jbao $

#ifndef PERTURBATION_H
#define PERTURBATION_H

#include "nr/nr.h"
#include "GeneNetwork.h"

/**
 * Class offers functionalities to apply perturbations to a gene network. 
 * The following types of perturbations are offered: single-gene, two-gene,
 * and multifactorial. Perturbations affect always only the maximum transcripion
 * rate m_i of the genes. In single-gene perturbation, only one m_i is perturbed at
 * a time. For two-gene perturbations, two m_i are perturbed. For multifactorial,
 * the max transcription rates of all genes are sampled from a normal
 * distribution with mean m_i and standard deviation m_i*CV_.
 * @author Daniel Marbach
 */
class Perturbation {
	
public:
	Perturbation();
	Perturbation(GeneNetwork& grn);
	virtual ~Perturbation();
	
	void printPerturbations(std::string postfix);
	void loadPerturbations(std::string label);
	
	int getNumPerturbations() { return numPerturbations_; }
	int getNumGenes() { return numGenes_; }

    Mat_DP& getPerturbations() { return perturbations_; }
    void setPerturbations(Mat_DP& perturbations) { perturbations_ = perturbations; }

	virtual void multifactorialAllGenesWeak(int numPerturbations) = 0;
	virtual void multifactorialAllGenesMax(int numPerturbations) = 0;
	virtual void multifactorialHub(int numPerturbations) = 0;
	virtual void perturbSingleGene(int numPerturbations, std::string gene) = 0;
	virtual void applyPerturbation(int k) = 0;
	
protected:
	/** The gene network to which the perturbations are being applied */
	GeneNetwork *grn_;
	/** The size of the network */
	int numGenes_;
	/** The wild-type (so that it can be restored after applying perturbations) */
	Vec_DP wildType_;
	/** The number of different multifactorial perturbations */
	int numPerturbations_;
	/** perturbations_(k, i) is the perturbed value of m_i in perturbation k */
	Mat_DP perturbations_;

};

#endif
