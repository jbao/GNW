#include <cmath>
#include "RandomParameterGaussian.h"
#include "GnwSettings.h"

// ============================================================================
// PUBLIC METHODS

/** Constructor, mean_ and stdev_ are set according to the specified interval (see class documentation above) */
RandomParameterGaussian::RandomParameterGaussian(double min, double max, bool logNormal) {
	logNormal_ = logNormal; 
	setMinMax(min, max); // also sets mean_ and stdev_ and transforms to log if necessary
}

/** Constructor, mean_ and stdev_ are set explicitly */
RandomParameterGaussian::RandomParameterGaussian(double min, double max, double mean, double stdev, bool logNormal) {
	logNormal_ = logNormal;
	setMinMax(min, max);
	setMeanStdev(mean, stdev);
}
	
// ----------------------------------------------------------------------------

/** 
 * Draw a new random number from the distribution. If the parameters were given on
 * log-scale and paramsOnLogScale_ is set, the result is transformed back to linear scale.
 */
double RandomParameterGaussian::getRandomValue() {
	
	//Normal normal = GnwSettings.getInstance().getNormalDistribution();
	double value;
	
	do {
		value = GnwSettings::Instance()->getNormalDistributionNumber(mean_, stdev_);
	} while (value < min_ || value > max_);

	if (logNormal_)
		value = pow(10.0, value);
	
	return value;
}

	
// ----------------------------------------------------------------------------

/** Set min_, max_, and also mean_ and stdev_ accordingly (transforms to log if logNormal_ is set) */
void RandomParameterGaussian::setMinMax(double min, double max) {
	
	if (logNormal_) {
		min_ = log10(min);
		max_ = log10(max);
	} else {
		min_ = min;
		max_ = max;	
	}
	mean_ = (min_+max_) / 2.0;
	stdev_ = (max_-min_) / 6.0;
}


// ----------------------------------------------------------------------------

/** Set mean_ and stdev_ (transforms to log if logNormal_ is set) */
void RandomParameterGaussian::setMeanStdev(double mean, double stdev) {
	
	if (logNormal_) {
		mean_ = log10(mean);
		stdev_ = log10(stdev);
	} else {
		mean_ = mean;
		stdev_ = stdev;	
	}
}


// ----------------------------------------------------------------------------

/** Set the mean, transform to log if logNormal_ is set */
void RandomParameterGaussian::setMean(double mean) {
	
	if (logNormal_)
		mean_ = log10(mean);
	else
		mean_ = mean;
}		


