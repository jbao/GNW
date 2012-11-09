// $Id: HillGene.h 29 2012-01-04 17:06:55Z jbao $

#ifndef HILLGENE_H
#define HILLGENE_H

#include <vector>
#include <string>
#include <map>
#include "Node.h"
#include "Edge.h"

/**HillGene implements the gene regulation function of a gene.
 * 
 * The model is based on Hill type kinetics and principles of thermodynamics. It's quite tricky, 
 * especially the random initialization in order to achieve biologically plausible gene regulation
 * functions. This will be described in a paper as soon as possible.
 * 
 * @author Daniel Marbach (firstname.name@gmail.com)
 * @author Thomas Schaffter (firstname.name@gmail.com)
 * 
 */
 
 class HillGene : public Node {
 
 public:
 	HillGene();
 	HillGene(std::string label);
 	HillGene(const HillGene& hg);
 	HillGene& operator=(const HillGene& rhs);
 	~HillGene();
 	void initialization(std::map<std::string, double>& params, 
            const std::vector<std::string>& inputGenes);
 	void compileParameters(std::map<std::string, double>& params);
 	
 	std::vector<std::string>& getInputGenes() { return inputGenes_; }
 	void setInputGenes(const std::vector<std::string>& inputs) { inputGenes_ = inputs; }
 	
 	std::vector<std::string>& getOutputGenes() { return outputGenes_; }
 	void setOutputGenes(const std::vector<std::string>& outputs) { 
        outputGenes_ = outputs; 
    }
 	
 	std::vector<Edge::sign>& getInputSigns() { return inputSigns_; }
 	void setInputSigns(const std::vector<Edge::sign>& inputSigns) { inputSigns_ = inputSigns; }
 	
 	//std::vector<RegulatoryModule>& getRegulatoryModules() { return regulatoryModules_; }
 	//void setRegulatoryModules(const std::vector<RegulatoryModule>& modules) { regulatoryModules_ = modules; }
 	
 	void setMaxTranscription(double max) { max_ = max; }
	double getMaxTranscription() { return max_; }
 	
 	void setDelta(double d) { delta_ = d; }
	double getDelta() { return delta_; }
 	
 	void setMaxTranslation(double d) { maxTranslation_ = d; }
	double getMaxTranslation() { return maxTranslation_; }

	void setDeltaProtein(double d) { deltaProtein_ = d; }
	double getDeltaProtein() { return deltaProtein_; }
 	
 	Vec_DP& getAlpha() { return alpha_; }
 	int getNumInputs() { return inputGenes_.size(); }
 	int getNumOutputs() { return outputGenes_.size(); }

 	double getBasalActivation() { return alpha_[0]; }
 	void perturbBasalActivation(double deltaBasalActivation);
 	void restoreWildTypeBasalActivation();
 	
 	double computeMRnaProductionRate(Vec_INT& geneIndex, Vec_I_DP& c);
 	double computeMRnaDegradationRate(double x);
 	double computeProteinProductionRate(double x);
 	double computeProteinDegradationRate(double y);
 	
 	void randomKineticInitialization();
 	void randomStructureInitialization();
 	void randomInitializationOfAlpha();
 	
 	void setInputEdgeTypesAccordingToDynamicalModel();
 	
 	void pruneInput(HillGene& gene);
 	
 	void dec2bin(long decimal, char *binary);
 	
 	void randomInitializationOfParameters();
 	
	int getNumActivators() {return numActivators_; }
	void setNumActivators(int n) { numActivators_ = n; }

	int getNumInhibitors() { return numInhibitors_;	}
	void setNumInhibitors(int n) {	numInhibitors_ = n; }
	
	Vec_DP& getK() { return k_; }
	void setK(const Vec_DP& k) { k_ = k; }
	
	Vec_DP& getN() { return n_; }
	void setN(const Vec_DP& n) { n_ = n; }
 
 private:
 	/** Relative activations for all possible states of the regulatory modules **/
	Vec_DP alpha_;
	/** The wild-type alpha_ vector (used as backup when alpha_ is perturbed */
	Vec_DP alphaWildType_;
	/** The regulatory modules (they are activated independently from each other) */
	//std::vector<RegulatoryModule> regulatoryModules_;
	std::vector<std::string> inputGenes_;
	std::vector<std::string> outputGenes_;
	std::vector<Edge::sign> inputSigns_;
	//GeneNetwork *grn_;
	
	//std::vector<Node> inputGenes_;
	/** Maximum transcription rate */
	double max_;
	/** Degradation rate */
	double delta_;
	/** Maximum translation rate */
	double maxTranslation_;
	/** Protein degradation rate */
	double deltaProtein_;
	/** Dissociation constants k for each state */
	//double *k_;
	/** Hill coefficients for the regulators (activators and deactivators) */
	//double *n_;
	/** Perturbation of basal transcription rate */
	//protected double perturbationBasalActivation_;
	/** The activators of this module */
	int numActivators_;	
	/** The deactivators of this module */
	int numInhibitors_;
	/** Dissociation constants k for each state */
	Vec_DP k_;
	/** Hill coefficients for the regulators (activators and deactivators) */
	Vec_DP n_;
 
 	void initializeHalfLife(std::map<std::string, double>& params);
 	void initializeRegulation(std::map<std::string, double>& params);
 	
 	void randomInitializationOfStructureAndModules();
 
 	double computeActivation(double x, int idx);
 	
 };
 
 #endif
