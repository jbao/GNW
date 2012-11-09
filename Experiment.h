// $Id: Experiment.h 19 2011-06-17 08:20:59Z jbao $

#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <string>
#include "nr/nr.h"
#include "GeneNetwork.h"
#include "PerturbationMultifactorial.h"

/** Abstract class for an experiment type.
 * 
 * Subclasses are SteadyStateExperiment
 * and TimeSeriesExperiment. The subclasses define the data types and implement
 * the corresponding simulations, this class contains everything that is common
 * to the different experiment types.
 * 
 * @author Daniel Marbach (firstname.name@gmail.com)
 * 
 */
class Experiment {

public:
	Experiment(const std::string& label);
	Experiment(Perturbation *perturbation, std::string label);
	Experiment(const Experiment& e);
 	Experiment& operator=(const Experiment& rhs);
	~Experiment();

	std::string getLabel() { return label_; }

	//void setNumGenes(int numGenes) { numGenes_ = numGenes; }
	//public int getNumGenes() { return numGenes_; }
	
	Perturbation* getPerturbation() { return perturbation_; }
	int getNumExperiments() { return numExperiments_; }
	
	GeneNetwork& getGrn() { return grn_; }
	void setGrn(GeneNetwork& grn) {
		grn_ = grn;
		numGenes_ = grn.getSize();
	}

protected:
	/** Define which type of solver to use (ODE or SDE) */
	//Solver.type solverType_;
	
	/** The label of the experiment (appended to filenames when saving, e.g. wildtype, knockouts, ...) */
	std::string label_;
	
	/** Reference to GeneNetwork object */
	GeneNetwork grn_;
	/** Number of genes */
	int numGenes_;
	
	// for rk4 numerical integration method
	Vec_DP *xx_p;
	Mat_DP *y_p;
	
	/** Defines the perturbations to be applied (compute only wild-type if null) */
	Perturbation *perturbation_;
	/** Number of perturbations, if perturbation_ is set, otherwise 1 (only wild-type) */
	int numExperiments_;
	/** Initial conditions (mRNA and proteins in one concatenated vector) */
	Vec_DP xy0_;

	/** Set true if proteins are modeled */
	//protected boolean modelTranslation_;
	
	/** Normal distribution used to generate the different types of noise */
	//private Normal normalDistribution_;
		
	/** Flag, set true after noise has been added to the data */
	bool noiseHasBeenAdded_;
	
	double addNoise(double x);
	
	void rkdumb(Vec_I_DP &vstart, const DP x1, const DP x2,
		void (GeneNetwork::*derivs)(const DP, Vec_I_DP &, Vec_O_DP &));
	
private:
	/** Set true to add normal noise to the data */
	bool addNormalNoise_;
	/** Set true to add lognormal noise to the data */
	bool addLognormalNoise_;
	/** Set true to use a realistic model of microarray noise, similar to a mix of normal and lognormal */
	bool addMicroarrayNoise_;
	/** The standard deviation of the normal noise */
	double normalStdev_;
	/** The standard deviation of the lognormal noise */
	double lognormalStdev_;

	double addLogNormalNoise(double x);
	double addMicroarrayNoise(double x);
	
	void rk4(Vec_I_DP &y, Vec_I_DP &dydx, const DP x, const DP h,
		Vec_O_DP &yout, void (GeneNetwork::*derivs)(const DP, Vec_I_DP &, Vec_O_DP &));

};

#endif
