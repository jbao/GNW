#include "RandomParameterUniform.h"
#include "GnwSettings.h"

/** Constructor */
RandomParameterUniform::RandomParameterUniform(double min, double max) {
	min_ = min;
	max_ = max;
}

/** Return a number in [min max] from a uniform distribution */
double RandomParameterUniform::getRandomValue() {
	return GnwSettings::Instance()->getUniformDistributionNumber(max_-min_)+min_;
}


