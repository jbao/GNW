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

$Id: GnwSettings.cpp 20 2011-08-03 09:27:19Z jbao $
*/

#include "GnwSettings.h"

// Global static pointer used to ensure a single instance of the class.
GnwSettings *GnwSettings::m_pInstance = NULL; 
	 
/** This function is called to create an instance of the class.
	Calling the constructor publicly is not allowed. The constructor
	is private and is only called by this Instance function.
*/	   
GnwSettings *GnwSettings::Instance() {
	if (!m_pInstance)   // Only allow one instance of class to be generated.
		m_pInstance = new GnwSettings;
	 
	return m_pInstance;
}
	
/**
 * Default constructor
 */
GnwSettings::GnwSettings() {
	/*
	if (randomSeed_ == -1)
		mersenneTwister_ = new MersenneTwister(new java.util.Date());
	else
		mersenneTwister_ = new MersenneTwister(randomSeed_);
	uniformDistribution_ = new Uniform(mersenneTwister_);
	normalDistribution_ = new Normal(0, 1, mersenneTwister_); // mean 0, stdev 1
	*/
	/** The unique instance of Universal (Singleton design pattern) */
	//private static GnwSettings instance_ = null;
	
	/** Mersenne Twister random engine (should be used by all other random number generators) */
	//private MersenneTwister mersenneTwister_;
	/** Uniform distribution random number generator */
	//private Uniform uniformDistribution_;
	/** Normal distribution random number generator */
	//private Normal normalDistribution_;

	// VARIOUS
	/** Seed for the random number generator. Set to -1 to use current time */
	randomSeed_ = -1;
	r_ = new MTRand();
	r_->seed(randomSeed_);
	r_->randNormInit();
	/** Default output directory to save stuff */
	outputDirectory_ = "";
	/** Model proteins and translation */
	modelTranslation_ = true;
    perturbationNumber_ = 0;

	// SUBNETWORK EXTRACTION
	/** The number of regulators in the extracted networks, set to 0 to disable control of number of regulators */
	numRegulators_ = -1;
	/** Vertices are added using truncated selection with the given fraction (0=greedy, 1=random selection) */
	truncatedSelectionFraction_ = 0.1;
	/** Number of seeds to be sampled from strongly connected components */
	numSeedsFromStronglyConnectedComponents_ = 0;
	
	// STEADY-STATE EXPERIMENTS
	/** Generate steady states for knockouts */
	ssKnockouts_ = false;
	/** Generate steady states knockdowns */
	ssKnockdowns_ = false;
	/** Generate steady states for multifactorial perturbations */
	ssMultifactorial_ = false;
	/** Generate steady states for dual knockouts */
	ssDualKnockouts_ = false;
	/**
	 * For deterministic simulations (ODEs), we return the steady-states as soon as convergence is reached.
	 * If there is no convergence until time maxtSteadyStateODE_, the values at this point are returned and a
	 * warning message is displayed.
	 */
	maxtSteadyStateODE_ = 2000;
	
	// TIME-SERIES EXPERIMENTS
	/** Generate time series for knockouts */
	tsKnockouts_ = false;
	/** Generate time series knockdowns */
	tsKnockdowns_ = false;
	/** Generate time series for multifactorial perturbations */
	tsMultifactorial_ = false;
	/** Generate time series for dual knockouts */
	tsDualKnockouts_ = false;
	/** Number of time-series experiments from different initial conditions */
	numTimeSeries_ = 10; 
	/** Number of measured points per time series (must be consistent with maxtTimeSeries_ and dt_, does *not* affect precision) */
	//private int numMeasuredPoints_ = 21;
	/** Default max duration time in time-series experiments (must be consistent with numTimePoints_ and dt_) */
	maxtTimeSeries_ = 1000;
	/** Time step for the time-series (must be consistent with numTimePoints_ and maxtTimeSeries_) */
	dt_ = 50;
	
	// MULTIFACTORIAL PERTURBATIONS
	/** Standard deviation for multifactorial perturbations */
	multifactorialStdev_ = 0.33;
	/** The probability that a gene is perturbed (for DREAM4 time series) */
	perturbationProbability_ = 0.33;
	
	// DETERMINISTIC MODEL (ODE)
	/** If set true, a deterministic simulation of the experiments is done using the ODE model */
	//private boolean simulateODE_ = false;
	/** Absolute _or_ relative precision _per variable_ need to be satisfied for convergence */
	absolutePrecision_ = 0.00001;
	/** See absolutePrecision_, in addition, this is also the tolerance used for integration */ 
	relativePrecision_ = 0.001;
	
	// EXPERIMENTAL NOISE
	/** Set true to add normal noise to the data */
	addNormalNoise_ = false;
	/** Set true to add lognormal noise to the data */
	addLognormalNoise_ = false;
	/** Set true to use a realistic model of microarray noise, similar to a mix of normal and lognormal */
	addMicroarrayNoise_ = true;
	/** The standard deviation of the normal noise */
	normalStdev_ = 0.025;
	/** The standard deviation of the lognormal noise */
	lognormalStdev_ = 0.075;
	/** Set true to normalize the datasets after adding the experimental noise */
	normalizeAfterAddingNoise_ = false;
	
	// RANDOM PARAMETERS
	/** Half-lives in minutes, Dassow2000 use [5 100]. */
	randomHalfLife_ = new RandomParameterGaussian(5, 50, false);
	//randomHalfLife_ = new RandomParameterGaussian(250, 350, false);
	/** Dissociation constants, Dassow2000 use [0.001 1] */
	//randomK_ = new RandomParameterUniform(0.01, 1);// DREAM3: (0.01, 1, 0, 0.2, false); lognormal: (0.1, 1, 0.1, 3.1623, true);
	randomK_ = new RandomParameterGaussian(0.01, 1, 0, 0.2, false);
	/** Hill coefficients, Dassow2000 use [1 10] */
	randomN_ = new RandomParameterGaussian(1, 10, 2, 2, false);
	/** Threshold for setting weak activations in random initialization */
	double weakActivation_ = 0.25;
	/** The difference in gene activation due to a module */
	randomDeltaActivation_ = new RandomParameterGaussian(weakActivation_, 1, false);
	/** Initialization of low basal rates (leakage) */
	randomLowBasalRate_ = new RandomParameterGaussian(0, weakActivation_, 0, 0.05, false);//(0.001, 0.25, true);
	/** Initialization of medium basal rates */
	randomMediumBasalRate_ = new RandomParameterGaussian(weakActivation_, 1-weakActivation_, false);	
	
	// PROCESS STATE
	stopSubnetExtraction_ = false;
	stopBenchmarkGeneration_ = false;

    // perturb certain fraction of hubs
    hubThreshold_ = 100;
    lowerPerturbationFraction_ = 0;
    upperPerturbationFraction_ = 0.2;
}
	


