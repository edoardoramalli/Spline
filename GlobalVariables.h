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

/* DA TOGLIERE Number of plausible versions of the experimental data generated during the
bootstrap procedure. numberOfBootstrapVariations = 1 : no bootstrap */
int numberOfBootstrapVariations = 20;

/* DA TOGLIERE Specifies whether the program should attempt to line up the maxima of
experimental data and models during the shift procedure */
bool lineUpMaxima = true;

/* DA TOGLIERE Specifies whether the program should use the sum of the four dissimilarity
indexes when calculating the shift, thus finding a single shift amount, or
whether the program should consider each of the four indexes separately, thus
obtaining four different values for the shift */
bool useSumOfIndexesForAlignment = true;

/* Number of trapezoids for the numerical calculation of the indexes */
int numberOfTrapezoids = 99;

/* Default value for the relative experimental error. Used in case relative
errors are not provided along with the experimental data */
double defaultRelativeError = 0.1;

/* Specifies whether to create the graphs comparing the spline for the
experimental data with the splines for the models */
bool graphsD0 = true;

/* Specifies whether to create the graphs comparing the first derivative of the
spline for the experimental data with the first derivatives of the splines for
the models */
bool graphsD1 = false;

/* Per stampare tutte le spline non solo la migliore */
bool allSplineExp = false;

/* Per stampare i nodi scelti per il calcolo della spline */
bool knotsToSave = false;

/* Per stampare i nodi scelti per il calcolo della spline */
bool saveExpData = false;

/* Stampa di un file txt con un recap sulle spline */
bool coseUtili = false;

/**/
bool criterionSSE = false;

/**/
bool criterionAIC = true;

/**/
bool criterionBIC = false;

/**/
bool removeAsymptotes = true;

/* Specifies whether at least one kind of graph is to be created */
bool graphs = true;

/* Number of points to be calculated for each spline when saving the spline to a
.R file or to a .txt for future plotting */
int graphPoints = 500;

/* Pascal's triangle */
vector<vector<double>> pascalsTriangle;