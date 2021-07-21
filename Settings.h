//
// Created by Edoardo Ramalli on 05/05/21.
// Modified by Timoteo Dinelli on 20/07/21.
//

#ifndef SPLINE_SETTINGS_H
#define SPLINE_SETTINGS_H
/* Order of the basis functions */
int m = 4;

/* Degree of the basis functions */
int g = 3;

/* Orders of magnitude of difference between the smallest and the largest
possible value of the smoothing parameter lambda */
double lambdaSearchInterval = 6;

/* Number of steps in the for cycle for minimizing the smoothing parameter
lambda */
int numberOfStepsLambda = 13;

/* Number of steps in the for cycle for minimizing the smoothing parameter
lambda */
double numberOfRatiolkForAICcUse = 40;

/* Fraction of the range of a spline on the y-axis for determining which
segments of the spline count as asymptotes. If the oscillations of the spline
at one of its extremities are contained within a horizontal area with size
determined by this value, the corresponding segment is identified as an
asymptote */
double fractionOfOrdinateRangeForAsymptoteIdentification = 0.005;

/* Fraction of the range of a spline on the y-axis for determining which points
count as well-defined maxima. In order to be considered a well-defined maximum,
a point in a spline must not only have first derivative equal to 0 and negative
second derivative, it must also be sufficiently distant from the two surrounding
minima. The minimum admissible distance is determined using this variable */
double fractionOfOrdinateRangeForMaximumIdentification = 0.025;

/* Specifies whether negative segments on the y-axis are admissible for the
splines or whether they should be replaced with straight lines with ordinate 0
*/
bool possibleNegativeOrdinates = false;

/* Number of trapezoids for the numerical calculation of the indexes */
int numberOfTrapezoids = 99;

/* Default value for the relative experimental error. Used in case relative
errors are not provided along with the experimental data */
double defaultRelativeError = 0.1;

/**/
bool removeAsymptotes = true;

/* Pascal's triangle */
vector<vector<double>> pascalsTriangle;

/* Number of points to be calculated for each spline when saving the spline to a
.R file or to a .txt for future plotting */
int graphPoints = 500;

/**/
string criterion = "AIC";

#endif //SPLINE_SETTINGS_H
