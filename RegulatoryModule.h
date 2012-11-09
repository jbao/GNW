#ifndef REGULATORYMODULE_H
#define REGULATORYMODULE_H

/** A HillGene is independently regulated by one ore more regulatory modules.
 *  
 * A module is either an enhancer or a repressor (defined by isEnhancer_). The module is controlled
 * synergistically by its activators, it is active only if all of them
 * are bound together. Optionally, it can have deactivators, which synergistically
 * silence the module. Note that the activators of a repressor module are
 * actually repressors of the gene.
 * 
 * @author Daniel Marbach (firstname.name@gmail.com)
 * 
 */
 
 #include <string>
 #include <vector>
 #include "nr/nr.h"
 
 class RegulatoryModule {
 
 public:
 	RegulatoryModule();
 	RegulatoryModule(const std::vector<std::string>& labels);
 	RegulatoryModule(const RegulatoryModule& rm);
 	RegulatoryModule& operator=(const RegulatoryModule& rhs);
 	~RegulatoryModule() {};
 	
 	double computeActivation(Vec_I_DP& x);
 	
 	void randomInitializationOfParameters();
 	
 	std::vector<bool>& getEdgeSigns();
 	
 	// ============================================================================
	// SETTERS AND GETTERS

	bool isEnhancer() { return isEnhancer_; }
	void setIsEnhancer(bool b) { isEnhancer_ = b; }
	
	bool bindsAsComplex() { return bindsAsComplex_; }
	void setBindsAsComplex(bool b) { bindsAsComplex_ = b; }

	int getNumInputs() { return numActivators_ + numDeactivators_; }

	int getNumActivators() {return numActivators_; }
	void setNumActivators(int n) { numActivators_ = n; }

	int getNumDeactivators() { return numDeactivators_;	}
	void setNumDeactivators(int n) {	numDeactivators_ = n; }
	
	Vec_DP& getK() { return k_; }
	void setK(const Vec_DP& k) { k_ = k; }
	
	Vec_DP& getN() { return n_; }
	void setN(const Vec_DP& n) { n_ = n; }
 
 private:
 	/** The type of the module (true = enhancer, false = repressor) */
	bool isEnhancer_;
	/** True if the activators first form a hetero-oligomer, and this complex binds to the promoter */
	bool bindsAsComplex_;
	/** The activators of this module */
	int numActivators_;	
	/** The deactivators of this module */
	int numDeactivators_;
	/** Dissociation constants k for each state */
	Vec_DP k_;
	/** Hill coefficients for the regulators (activators and deactivators) */
	Vec_DP n_;
	/** Input labels of this module */
	std::vector<std::string> labels_;
 
 };
 
 #endif