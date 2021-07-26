//
// Created by Edoardo Ramalli on 05/05/21.
// Modified by Timoteo Dinelli on 20/07/21.
//

#ifndef SPLINE_SETTINGS_H
#define SPLINE_SETTINGS_H
/* Order of the basis functions */
int m;

/* Degree of the basis functions */
int g;

/* Orders of magnitude of difference between the smallest and the largest
possible value of the smoothing parameter lambda */
int lambdaSearchInterval; // ATTENZIONE PRIMA ERA UN DOUBLE messi cast in Spline.h

/* Number of steps in the for cycle for minimizing the smoothing parameter
lambda */
int numberOfStepsLambda;

/* Number of steps in the for cycle for minimizing the smoothing parameter
lambda */
int numberOfRatiolkForAICcUse;

/* Fraction of the range of a spline on the y-axis for determining which
segments of the spline count as asymptotes. If the oscillations of the spline
at one of its extremities are contained within a horizontal area with size
determined by this value, the corresponding segment is identified as an
asymptote */
double fractionOfOrdinateRangeForAsymptoteIdentification;

/* Fraction of the range of a spline on the y-axis for determining which points
count as well-defined maxima. In order to be considered a well-defined maximum,
a point in a spline must not only have first derivative equal to 0 and negative
second derivative, it must also be sufficiently distant from the two surrounding
minima. The minimum admissible distance is determined using this variable */
double fractionOfOrdinateRangeForMaximumIdentification;

/* Specifies whether negative segments on the y-axis are admissible for the
splines or whether they should be replaced with straight lines with ordinate 0
*/
//bool possibleNegativeOrdinates;


/**/
//bool removeAsymptotes;

/* Pascal's triangle */
vector<vector<double>> pascalsTriangle;

/* Number of points to be calculated for each spline when saving the spline to a
.R file or to a .txt for future plotting */
int graphPoints;

/**/
string criterion;

#endif //SPLINE_SETTINGS_H
