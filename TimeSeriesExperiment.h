// $Id: TimeSeriesExperiment.h 29 2012-01-04 17:06:55Z jbao $

#ifndef TIMESERIESEXPERIMENT_H
#define TIMESERIESEXPERIMENT_H

#include <vector>
#include "Experiment.h"

/** Time course experiments, see documentation for details.
 *
 * @author Thomas Schaffter (firstname.name@gmail.com)
 * @author Daniel Marbach (firstname.name@gmail.com)
 *
 */
 class TimeSeriesExperiment : public Experiment {
 
 public:
 	TimeSeriesExperiment(Perturbation *perturbation, bool restoreWildTypeAtHalftime, std::string label);
 	~TimeSeriesExperiment();
    TimeSeriesExperiment(const TimeSeriesExperiment& e);
    TimeSeriesExperiment& operator=(const TimeSeriesExperiment& rhs);
 	
 	void run(Vec_DP xy0);
 	void integrate(int i);
 	void printAll(std::string postfix);
 	
 	int getNumTimePoints() { return numTimePoints_; }
	std::vector<Mat_DP> getTimeSeries() { return timeSeries_; }
	std::vector<Mat_DP> getTimeSeriesProteins() { return timeSeriesProteins_; }
	bool getRestoreWildTypeAtHalftime() { return restoreWildTypeAtHalftime_; }
	
	double getMaximumConcentration();
	void addNoise();
	void normalize(double max);
 
 private:
 	/** Time series data */
	std::vector<Mat_DP> timeSeries_;
	/** Input function */
	std::vector<Mat_DP> inputs_;
	/** Output flux */
	std::vector<Mat_DP> outputs_;
	/** Protein data */
	std::vector<Mat_DP> timeSeriesProteins_;
	/** The duration of the experiment */
	double maxt_;
	/** Number of time points (maxt/dt + 1)*/
	int numTimePoints_;
	/** Set true to remove the perturbation after maxt/2 */
	bool restoreWildTypeAtHalftime_;
 
 	void setMaxtAndNumTimePoints();
 	
 	void printTrajectories(std::string postfix, std::vector<Mat_DP> timeSeries);
 	void printMatrix(std::string postfix, std::vector<Mat_DP> timeSeries);
 	
 	double getMaximumConcentration(const Mat_I_DP& ts);
 	void addNoise(Mat_DP& ts);
 	void normalize(Mat_DP& ts, double max);
    
 };
 
 #endif
