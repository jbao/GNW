// $Id: GnwSettings.h 20 2011-08-03 09:27:19Z jbao $

#ifndef GNWSETTINGS_H
#define GNWSETTINGS_H

#include <string>
#include "MersenneTwister.h"
#include "RandomParameter.h"
#include "RandomParameterGaussian.h"
#include "RandomParameterUniform.h"

/** 
 * Offers global parameters (settings) and functions used by all classes of the
 * gnw package.
 * 
 * GnwSettings makes use of the Singleton design pattern: There's at most one
 * instance present, which can only be accessed through getInstance().
 * 
 * @author Thomas Schaffter (firstname.name@gmail.com)
 * @author Daniel Marbach (firstname.name@gmail.com)
 */
class GnwSettings {

public:
	static GnwSettings *Instance();
	
	// ============================================================================
	// SETTERS AND GETTERS
	
	double getUniformDistributionNumber(const double& bound) {
		return r_->rand(bound); }
	double getNormalDistributionNumber(const double& mean, const double& std) { 
		return r_->randNormInv(mean, std); }
	
	int getRandomSeed() { return randomSeed_; }
	void setRandomSeed(int s) { 
		randomSeed_ = s; 
		r_->seed(randomSeed_);
		//mersenneTwister_ = new MersenneTwister(randomSeed_);
		//uniformDistribution_ = new Uniform(mersenneTwister_);
		//normalDistribution_ = new Normal(0, 1, mersenneTwister_); // mean 0, stdev 1
	}
	
	MTRand* getRNG() { return r_; }
	void setRNG(MTRand *rng) { r_ = rng; }
	
	int getPerturbationNumber() { return perturbationNumber_; }
	void setPerturbationNumber(int perturbationNumber) { perturbationNumber_ = perturbationNumber; }
	
	int getNumRegulators() { return numRegulators_; }
	void setNumRegulators(int numRegulators) { numRegulators_ = numRegulators; }
	
	double getTruncatedSelectionFraction() { return truncatedSelectionFraction_; }
	void setTruncatedSelectionFraction(double f) { truncatedSelectionFraction_ = f; }
	
	int getNumSeedsFromStronglyConnectedComponents() { return numSeedsFromStronglyConnectedComponents_; }
	void setNumSeedsFromStronglyConnectedComponents(int N) { numSeedsFromStronglyConnectedComponents_ = N; }
	
	void setAbsolutePrecision(double value) { absolutePrecision_ = value; }
	double getAbsolutePrecision() { return absolutePrecision_; }
	
	void setRelativePrecision(double value) { relativePrecision_ = value; }
	double getRelativePrecision() { return relativePrecision_; }
	
	void setAddNormalNoise(bool b) { addNormalNoise_ = b; }
	bool getAddNormalNoise() { return addNormalNoise_; }
	
	void setAddLognormalNoise(bool b) { addLognormalNoise_ = b; }
	bool getAddLognormalNoise() { return addLognormalNoise_; }
	
	void setAddMicroarrayNoise(bool b) { addMicroarrayNoise_ = b; }
	bool getAddMicroarrayNoise() { return addMicroarrayNoise_; }
	
	void setNormalStdev(double s) { normalStdev_ = s; }
	double getNormalStdev() { return normalStdev_; }
	
	void setLognormalStdev(double s) { lognormalStdev_ = s; }
	double getLognormalStdev() { return lognormalStdev_; }
	
	void setNormalizeAfterAddingNoise(bool b) { normalizeAfterAddingNoise_ = b; }
	bool getNormalizeAfterAddingNoise() { return normalizeAfterAddingNoise_; }
	
	void setNumTimeSeries(int n) { numTimeSeries_ = n; }
	int getNumTimeSeries() { return numTimeSeries_; }
	
	void setMaxtTimeSeries(double maxt) { maxtTimeSeries_ = maxt; }
	double getMaxtTimeSeries() { return maxtTimeSeries_; }
	
	void setMaxtSteadyStateODE(double maxt) { maxtSteadyStateODE_ = maxt; }
	double getMaxtSteadyStateODE() { return maxtSteadyStateODE_; }
	
	void setDt(double dt) { dt_ = dt; }
	double getDt() { return dt_; }
	
	double getMultifactorialStdev() { return multifactorialStdev_; }
	void setMultifactorialStdev(double cv) { multifactorialStdev_ = cv; }
	
	double getPerturbationProbability() { return perturbationProbability_; }
	void setPerturbationProbability(double p) { perturbationProbability_ = p; }
	
	void setModelTranslation(bool b) { modelTranslation_ = b; }
	bool getModelTranslation() { return modelTranslation_; }
	
	void setIgnoreAutoregulatoryInteractionsInEvaluation(bool b) { ignoreAutoregulatoryInteractionsInEvaluation_ = b; }
	bool getIgnoreAutoregulatoryInteractionsInEvaluation() { return ignoreAutoregulatoryInteractionsInEvaluation_; }
	
	//void setSimulateODE(boolean b) { simulateODE_ = b; }
	//public boolean getSimulateODE() { return simulateODE_; }
	
	double getRandomHalfLife() { return randomHalfLife_->getRandomValue();	}
	void setRandomHalfLife(RandomParameter *r) {randomHalfLife_ = r; }
	
	double getRandomK() { return randomK_->getRandomValue(); }
	void setRandomK(RandomParameter *r) { randomK_ = r; }
	
	double getRandomN() { return randomN_->getRandomValue(); }
	void setRandomN(RandomParameter *r) { randomN_ = r; }
	
	double getRandomDeltaActivation() { return randomDeltaActivation_->getRandomValue(); }
	void setRandomDeltaActivation(RandomParameter *r) { randomDeltaActivation_ = r; }
	
	double getWeakActivation() { return weakActivation_; }
	void setWeakActivation(double w) { weakActivation_ = w; }
	
	double getRandomLowBasalRate() { return randomLowBasalRate_->getRandomValue(); }
	void setRandomLowBasalRate(RandomParameter *r) { randomLowBasalRate_ = r; }
	
	double getRandomMediumBasalRate() { return randomMediumBasalRate_->getRandomValue(); }
	void setRandomMediumBasalRate(RandomParameter *r) { randomMediumBasalRate_ = r; }
	
	void generateSsKnockouts(bool b) { ssKnockouts_ = b; }
	bool generateSsKnockouts() { return ssKnockouts_; }
	
	void generateSsKnockdowns(bool b) { ssKnockdowns_ = b; }
	bool generateSsKnockdowns() { return ssKnockdowns_; }
	
	void generateSsMultifactorial(bool b) { ssMultifactorial_ = b; }
	bool generateSsMultifactorial() { return ssMultifactorial_; }
	
	void generateSsDualKnockouts(bool b) { ssDualKnockouts_ = b; }
	bool generateSsDualKnockouts() { return ssDualKnockouts_; }
	
	void generateTsKnockouts(bool b) { tsKnockouts_ = b; }
	bool generateTsKnockouts() { return tsKnockouts_; }
	
	void generateTsKnockdowns(bool b) { tsKnockdowns_ = b; }
	bool generateTsKnockdowns() { return tsKnockdowns_; }
	
	void generateTsMultifactorial(bool b) { tsMultifactorial_ = b; }
	bool generateTsMultifactorial() { return tsMultifactorial_; }
	
	void generateTsDualKnockouts(bool b) { tsDualKnockouts_ = b; }
	bool generateTsDualKnockouts() { return tsDualKnockouts_; }
	
	void generateTsConstantInput(bool b) { tsConstantInput_ = b; }
	bool generateTsConstantInput() { return tsConstantInput_; }
	
	//public void setNumMeasuredPoints(int num) { numMeasuredPoints_ = num; }
	//public int getNumMeasuredPoints() { return numMeasuredPoints_; }
	
	std::string getOutputDirectory() { return outputDirectory_; }
	void setOutputDirectory(std::string &absPath) {
		outputDirectory_ = absPath;
		std::string sep = "/";
		if (outputDirectory_.substr(outputDirectory_.length()-1,1) != sep)
			outputDirectory_ += sep;
	}
	
	void stopSubnetExtraction(bool b) { stopSubnetExtraction_ = b; }
	bool stopSubnetExtraction() { return stopSubnetExtraction_; }
	
	void stopBenchmarkGeneration(bool b) { stopBenchmarkGeneration_ = b; }
	bool stopBenchmarkGeneration() { return stopBenchmarkGeneration_; }

    void setHubThreshold(int hubThreshold) { hubThreshold_ = hubThreshold; }
    int getHubThreshold() { return hubThreshold_; }

    void setPerturbationFraction(double lower, double upper) {
        lowerPerturbationFraction_ = lower;
        upperPerturbationFraction_ = upper;
    }
    double getLowerPerturbationFraction() { return lowerPerturbationFraction_; }
    double getUpperPerturbationFraction() { return upperPerturbationFraction_; }

private:
	GnwSettings();  // Private so that it can  not be called
	static GnwSettings *m_pInstance;
	
	/** Mersenne Twister random engine (should be used by all other random number generators) */
	MTRand *r_;
	/** Uniform distribution random number generator */
	//private Uniform uniformDistribution_;
	/** Normal distribution random number generator */
	//private Normal normalDistribution_;

	// VARIOUS
	/** Seed for the random number generator. Set to -1 to use current time */
	int randomSeed_;
	/** Default output directory to save stuff */
	std::string outputDirectory_;
	/** Model proteins and translation */
	bool modelTranslation_;
	/** Set true to remove self-links (Gi->Gi) when generating kinetic models */
	//private boolean removeAutoregulatoryInteractionsFromGeneNetworks_;
	/** Set true to ignore self-links (Gi->Gi) when saving gold standards in DREAM format */
	bool ignoreAutoregulatoryInteractionsInEvaluation_;
	int perturbationNumber_;

	// SUBNETWORK EXTRACTION
	/** The number of regulators in the extracted networks, set to 0 to disable control of number of regulators */
	int numRegulators_;
	/** Vertices are added using truncated selection with the given fraction (0=greedy, 1=random selection) */
	double truncatedSelectionFraction_;
	/** Number of seeds to be sampled from strongly connected components */
	int numSeedsFromStronglyConnectedComponents_;
	
	// STEADY-STATE EXPERIMENTS
	/** Generate steady states for knockouts */
	bool ssKnockouts_;
	/** Generate steady states knockdowns */
	bool ssKnockdowns_;
	/** Generate steady states for multifactorial perturbations */
	bool ssMultifactorial_;
	/** Generate steady states for dual knockouts */
	bool ssDualKnockouts_;
	/**
	 * For deterministic simulations (ODEs), we return the steady-states as soon as convergence is reached.
	 * If there is no convergence until time maxtSteadyStateODE_, the values at this point are returned and a
	 * warning message is displayed.
	 */
	double maxtSteadyStateODE_;
	
	// TIME-SERIES EXPERIMENTS
	/** Generate time series for knockouts */
	bool tsKnockouts_;
	/** Generate time series knockdowns */
	bool tsKnockdowns_;
	/** Generate time series for multifactorial perturbations */
	bool tsMultifactorial_;
	/** Generate time series for dual knockouts */
	bool tsDualKnockouts_;
	/** Generate time series with constant input */
	bool tsConstantInput_;
	/** Number of time-series experiments from different initial conditions */
	int numTimeSeries_; 
	/** Number of measured points per time series (must be consistent with maxtTimeSeries_ and dt_, does *not* affect precision) */
	//private int numMeasuredPoints_ = 21;
	/** Default max duration time in time-series experiments (must be consistent with numTimePoints_ and dt_) */
	double maxtTimeSeries_;
	/** Time step for the time-series (must be consistent with numTimePoints_ and maxtTimeSeries_) */
	double dt_;
	
	// MULTIFACTORIAL PERTURBATIONS
	/** Standard deviation for multifactorial perturbations */
	double multifactorialStdev_;
	/** The probability that a gene is perturbed (for DREAM4 time series) */
	double perturbationProbability_;
	
	// DETERMINISTIC MODEL (ODE)
	/** If set true, a deterministic simulation of the experiments is done using the ODE model */
	//private boolean simulateODE_ = false;
	/** Absolute _or_ relative precision _per variable_ need to be satisfied for convergence */
	double absolutePrecision_;
	/** See absolutePrecision_, in addition, this is also the tolerance used for integration */ 
	double relativePrecision_;
	
	// EXPERIMENTAL NOISE
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
	/** Set true to normalize the datasets after adding the experimental noise */
	bool normalizeAfterAddingNoise_;
	
	// RANDOM PARAMETERS
	/** Half-lives in minutes, Dassow2000 use [5 100]. */
	RandomParameter *randomHalfLife_;
	/** Dissociation constants, Dassow2000 use [0.001 1] */
	RandomParameter *randomK_;// DREAM3: (0.01, 1, 0, 0.2, false); lognormal: (0.1, 1, 0.1, 3.1623, true);
	/** Hill coefficients, Dassow2000 use [1 10] */
	RandomParameter *randomN_;
	/** Threshold for setting weak activations in random initialization */
	double weakActivation_;
	/** The difference in gene activation due to a module */
	RandomParameter *randomDeltaActivation_;
	/** Initialization of low basal rates (leakage) */
	RandomParameter *randomLowBasalRate_;//(0.001, 0.25, true);
	/** Initialization of medium basal rates */
	RandomParameter *randomMediumBasalRate_;	
	
	// PROCESS STATE
	bool stopSubnetExtraction_;
	bool stopBenchmarkGeneration_;

    // perturb certain fraction of hubs
    int hubThreshold_;
    double lowerPerturbationFraction_;
    double upperPerturbationFraction_;
};

#endif
