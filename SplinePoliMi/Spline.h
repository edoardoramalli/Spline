
#include "Settings.h"

class Spline {

public:

    /* Type of spline. 0: Experimental data;  1: Model;  2: Error spline */
    int splineType;

    /* Abscissae of the spline */
    vector<double> abscissae;

    /* Ordinates of the spline */
    vector<double> ordinates;

    /* Number of data points */
    int n;

    /* Specifies whether there are enough data points to calculate the spline,
    or whether the spline can be considered a flat line with ordinate = 0 when
    compared to the experimental data */
    bool possibleToCalculateSpline;

    /* Real knots of the spline */
    vector<double> knots;

    /* Number of real knots of the spline */
    int numberOfKnots;

    /* Number of polynomials of the spline */
    int numberOfPolynomials;

    /* Coefficients of the polynomials of the spline, excluding those
    corresponding to coincident knots at the end points. coeffD0[i][j] refers to
    polynomial i and the coefficient of x^j */
    vector<vector<double>> coeffD0;

    /* Coefficients of the first derivatives of the polynomials of the spline,
    excluding those corresponding to coincident knots at the end points.
    coeffD1[i][j] refers to polynomial i and the coefficient of x^j */
    vector<vector<double>> coeffD1;

    /* Coefficients of the second derivatives of the polynomials of the spline,
    excluding those corresponding to coincident knots at the end points.
    coeffD2[i][j] refers to polynomial i and the coefficient of x^j */
    vector<vector<double>> coeffD2;

    /* Abscissae of the spline, as initially obtained from the input file */
    vector<double> originalAbscissae;

    /* Ordinates of the spline, as initially obtained from the input file */
    vector<double> originalOrdinates;

    /* Distance between the biggest and the smallest abscissae of the spline */
    double xRange;

    /* Degrees of freedom of the spline */
    int K;

    ////////////////////////////////////////////////////////////////////////////

    /* Calculates the spline */
    void solve(const vector<double>& abscissae,
               const vector<double>& ordinates,
               int splineType,
               int numberOfAbscissaeSeparatingConsecutiveKnots);

    /* Calculates the ordinate of the spline at position x on the x-axis */
    double D0(double x);

    /* Calculates the ordinate of the first derivative of the spline at position
    x on the x-axis */
    double D1(double x);

    /* Calculates the ordinate of the second derivative of the spline at
    position x on the x-axis */
    double D2(double x);

    /* Calculates the ordinate of the spline at position powersOfX[1] on the
    x-axis, using the normalized spline coefficients */
    double D0(const vector<double>& powersOfX);

    /* Calculates the ordinate of the first derivative of the spline at position
    powersOfX[1] on the x-axis, using the normalized spline coefficients */
    double D1(const vector<double>& powersOfX);

    /* Calculates the ordinate of the spline at position powersOfX[1] on the
    x-axis, using the normalized and shifted spline coefficients */
    double D0Shift(const vector<double>& powersOfX);

    /* Calculates the ordinate of the first derivative of the spline at position
    powersOfX[1] on the x-axis, using the normalized and shifted spline
    coefficients */
    double D1Shift(const vector<double>& powersOfX);

    /* Calculates coeffD0_shift_normalized, coeffD1_shift_normalized and
    knots_shift */
    void calculateShift(double shift);

    /* Calculates the real different roots of the spline or of the first or of
    the second derivative of the spline. Returns a vector with the roots sorted
    from smallest to largest. Returns a size() = 0 vector if there are no real
    roots */
    vector<double> calculateRoots(double derivativeOrder);

    /* Powers of the abscissa being considered */
    vector<double> powers;

////////////////////////////////////////////////////////////////////////////////

private:

    /* Degrees of freedom of the spline minus 1 */
    int G;

    /* Smoothing parameter */
    double lambda;

    /* Base 10 logarithm of the smoothing parameter lambda */
    double log10lambda;

    /* Value of log10lambda for which the elements of the matrices FiTFi and R
    have the same order of magnitude, rounded to the nearest 0.5 */
    double log10lambdaForSameOrderOfMagnitude;

    /* Lowest value of the interval for the search of the minimum for
    log10lambda */
    double log10lambdaMin;

    /* Highest value of the interval for the search of the minimum for
    log10lambda */
    double log10lambdaMax;

    /* Spline coefficients for obtaining coeffD0, coeffD1 and coeffD2 */
    vector<double> splineCoefficients;

    /* Knots of the spline, including non-real ones at the end points */
    vector<double> knotsForCalculations;

    /* Number of polynomials which are asymptotes on the left of the spline */
    double numberOfAsymptotePolynomialsLeft;

    /* Number of polynomials which are asymptotes on the right of the spline */
    double numberOfAsymptotePolynomialsRight;

    /* coeffD0, normalized */
    vector<vector<double>> coeffD0_normalized;

    /* coeffD1, normalized */
    vector<vector<double>> coeffD1_normalized;

    /* Shift value applied to coeffD0_normalized and coeffD1_normalized to
    obtain coeffD0_shift_normalized and coeffD1_shift_normalized, and to 'knots'
    to obtain knots_shift */
    double shift;

    /* coeffD0, normalized and shifted */
    vector<vector<double>> coeffD0_shift_normalized;

    /* coeffD1, normalized and shifted */
    vector<vector<double>> coeffD1_shift_normalized;

    /* Real knots of the spline, shifted */
    vector<double> knots_shift;

    ////////////////////////////////////////////////////////////////////////////

    /* Chooses the knots for the spline */
    void chooseKnots(int numberOfAbscissaeSeparatingConsecutiveKnots);

    /* Calculates the coefficients of the polynomials of the spline and the
    coefficients of the first derivative of the polynomials of the spline, for
    the real knots of the spline */
    void calculateCoefficients();

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void Spline::solve(const vector<double>& Abscissae,
                   const vector<double>& Ordinates,
                   int SplineType,
                   int numberOfAbscissaeSeparatingConsecutiveKnots) {

    abscissae = Abscissae;
    ordinates = Ordinates;
    splineType = SplineType;

    n = abscissae.size();

    if (originalAbscissae.size() == 0) {
        originalAbscissae = abscissae;
        originalOrdinates = ordinates;
    }

    possibleToCalculateSpline = abscissae.size() > 1 ? true : false;

    if (!possibleToCalculateSpline)
        return;

    this->chooseKnots(numberOfAbscissaeSeparatingConsecutiveKnots);

    this->calculateCoefficients();


//    this->findMaximaBetweenExtremes();



}

double Spline::D0(double x) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (x > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates the powers of x
    for (int i=1; i<m; ++i)
        powers[i] = powers[i-1]*x;

    // Calculates D0(x)
    double y = 0;
    for (int i=0; i<m; ++i)
        y += coeffD0[indexOfPolynomial][i]*powers[i];

    return y;

}



double Spline::D1(double x) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (x > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates the powers of x
    for (int i=1; i<g; ++i)
        powers[i] = powers[i-1]*x;

    // Calculates D1(x)
    double y = 0;
    for (int i=0; i<g; ++i)
        y += coeffD1[indexOfPolynomial][i]*powers[i];

    return y;

}



double Spline::D2(double x) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (x > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates the powers of x
    for (int i=1; i<g-1; ++i)
        powers[i] = powers[i-1]*x;

    // Calculates D2(x)
    double y = 0;
    for (int i=0; i<g-1; ++i)
        y += coeffD2[indexOfPolynomial][i]*powers[i];

    return y;

}



double Spline::D0(const vector<double>& powersOfX) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (powersOfX[1] > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates D0(powersOfX[1])
    double y = 0;
    for (int i=0; i<m; ++i)
        y += coeffD0_normalized[indexOfPolynomial][i]*powersOfX[i];

    return y;

}



double Spline::D1(const vector<double>& powersOfX) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (powersOfX[1] > knots[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates D1(powersOfX[1])
    double y = 0;
    for (int i=0; i<g; ++i)
        y += coeffD1_normalized[indexOfPolynomial][i]*powersOfX[i];

    return y;

}



double Spline::D0Shift(const vector<double>& powersOfX) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (powersOfX[1] > knots_shift[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates D0(powersOfX[1])
    double y = 0;
    for (int i=0; i<m; ++i)
        y += coeffD0_shift_normalized[indexOfPolynomial][i]*powersOfX[i];

    return y;

}



double Spline::D1Shift(const vector<double>& powersOfX) {

    int indexOfPolynomial = 0;
    for (int i=0; i<numberOfKnots-1; ++i)
        if (powersOfX[1] > knots_shift[i])
            indexOfPolynomial = i;
        else
            break;

    // Calculates D1(powersOfX[1])
    double y = 0;
    for (int i=0; i<g; ++i)
        y += coeffD1_shift_normalized[indexOfPolynomial][i]*powersOfX[i];

    return y;

}



void Spline::calculateShift(double Shift) {

    shift = Shift;

    auto powersShifts = vector<double>(m,1);

    for (int a=1; a<m; ++a)
        powersShifts[a] = powersShifts[a-1]*shift;

    coeffD0_shift_normalized =
        vector<vector<double>>(numberOfPolynomials,vector<double>(m,0));

    coeffD1_shift_normalized =
        vector<vector<double>>(numberOfPolynomials,vector<double>(m,0));

    for (int a=0; a<numberOfPolynomials; ++a)
        for (int b=0; b<m; ++b)
            coeffD0_shift_normalized[a][b] = coeffD0_normalized[a][b];

    for (int a=0; a<numberOfPolynomials; ++a)
        for (int b=0; b<g; ++b)
            coeffD1_shift_normalized[a][b] = coeffD1_normalized[a][b];

    for (int a=0; a<numberOfPolynomials; ++a)
        for (int b=0; b<g; ++b)
            for (int c=1; c<=g-b; ++c) {
                if (c%2 != 0)
                    coeffD0_shift_normalized[a][b] -=
                    pascalsTriangle[b+c][b]*coeffD0_shift_normalized[a][b+c]*
                    powersShifts[c];
                else
                    coeffD0_shift_normalized[a][b] +=
                    pascalsTriangle[b+c][b]*coeffD0_shift_normalized[a][b+c]*
                    powersShifts[c];
            }

    for (int a=0; a<numberOfPolynomials; ++a)
        for (int b=0; b<g; ++b)
            for (int c=1; c<=g-b; ++c) {
                if (c%2 != 0)
                    coeffD1_shift_normalized[a][b] -=
                    pascalsTriangle[b+c][b]*coeffD1_shift_normalized[a][b+c]*
                    powersShifts[c];
                else
                    coeffD1_shift_normalized[a][b] +=
                    pascalsTriangle[b+c][b]*coeffD1_shift_normalized[a][b+c]*
                    powersShifts[c];
            }

    knots_shift = vector<double>(numberOfKnots);

    for (int a=0; a<numberOfKnots; ++a)
        knots_shift[a] = knots[a]+shift;

}



void Spline::chooseKnots(int numberOfAbscissaeSeparatingConsecutiveKnots) {

    int number = numberOfAbscissaeSeparatingConsecutiveKnots;

    double meanKnotDistance =
        (abscissae.back()-abscissae[0]) / (double)(abscissae.size()-1);

    if (splineType == 1 /*Model*/) {

        vector<double> newX;
        vector<double> newY;

        newX.push_back(abscissae[0]);
        newY.push_back(ordinates[0]);

        // If there are less than 30 points, adds enough points to the spline to
        // reach at least 30 points
        if (abscissae.size() < 30) {

            double abscissaeLength = (abscissae.back()-abscissae[0]);
            int minPointsToAdd = 30-abscissae.size();

            for (int a=1; a<(int)abscissae.size(); ++a) {
                double segmentLength = (abscissae[a]-abscissae[a-1]);
                int numberOfPointstoAdd =
                    segmentLength/abscissaeLength*(double)(minPointsToAdd+1);
                double distanceBetweenPoints =
                    segmentLength/(double)(numberOfPointstoAdd+1);
                double slope = (ordinates[a]-ordinates[a-1])/segmentLength;
                for (int b=0; b<numberOfPointstoAdd; ++b) {
                    newX.push_back(newX.back()+distanceBetweenPoints);
                    newY.push_back(
                        ordinates[a-1]+slope*(newX.back()-abscissae[a-1]));
                }
            newX.push_back(abscissae[a]);
            newY.push_back(ordinates[a]);
            }
        }

        // Adds extra points between consecutive data points with a distance on
        // the x-axis greater than 3.*meanKnotDistance
        if (abscissae.size() >= 30)
            for (int a=1; a<(int)abscissae.size(); ++a) {
                double segmentLength = (abscissae[a]-abscissae[a-1]);
                if (segmentLength > 3.*meanKnotDistance) {
                    int numberOfNewPoints =
                        (int)(segmentLength/meanKnotDistance);
                    double distanceBetweenPoints =
                        segmentLength / (double)(numberOfNewPoints+1);
                    double slope = (ordinates[a]-ordinates[a-1])/segmentLength;
                    for (int b=0; b<numberOfNewPoints; ++b) {
                        newX.push_back(newX.back()+distanceBetweenPoints);
                        newY.push_back(
                            ordinates[a-1]+slope*(newX.back()-abscissae[a-1]));
                    }
                }
                newX.push_back(abscissae[a]);
                newY.push_back(ordinates[a]);
            }

        abscissae = newX;
        ordinates = newY;

        n = abscissae.size();

		meanKnotDistance =
			(abscissae.back()-abscissae[0]) / (double)(abscissae.size()-1);

    }

	double maxOrdinate = ordinates[0];
	double minOrdinate = ordinates[0];
	for (int a=1; a<n; ++a) {
		if (ordinates[a] > maxOrdinate)
			maxOrdinate = ordinates[a];
		if (ordinates[a] < minOrdinate)
			minOrdinate = ordinates[a];
	}
	double height = maxOrdinate - minOrdinate;

    knots.push_back(abscissae[0]);

    if (abscissae.size() > 2) {

        double y = ordinates[0];
        int k = 0;
		int l = 0;
        for (int a=1; a<(int)abscissae.size()-1; ++a) {
            ++k;
            double difference = abscissae[a] - abscissae[a-1];
            if (k > number ||
				difference > 2.*meanKnotDistance ||
				fabs(ordinates[a]-ordinates[a-1]) > 0.1*height ||
                l > 0)
                if (difference > 0.2*meanKnotDistance)
                    if (ordinates[a] != y) {
                        knots.push_back(abscissae[a]);
                        y = ordinates[a];
                        k = 0;
                        --l;
						if (fabs(ordinates[a]-ordinates[a-1]) > 0.1*height)
							l = number+1;
                    }
        }

        // Improves the positioning of the knots if the data ends with a
        // horizontal asymptote. Some of the intervals chosen with this approach
        // might not contain any data points
        if (y == ordinates.back() && k > 10) { 
            double knot = knots.back();
            knots.push_back(knot+(abscissae.back()-knot)*4./12.);
            knots.push_back(knot+(abscissae.back()-knot)*8./12.);
            knots.push_back(knot+(abscissae.back()-knot)*10./12.);
            knots.push_back(knot+(abscissae.back()-knot)*11./12.);
        }

    }

    knots.push_back(abscissae.back());

    // Sets the values of numberOfKnots, numberOfPolynomials, K and G
    numberOfKnots = knots.size();
    numberOfPolynomials = numberOfKnots - 1;
    K = numberOfKnots - 2 + m;
    G = K-1;

    xRange = knots.back() - knots[0];

    // Fills knotsForCalculations with the current knots plus additional knots
    // on the left and on the right of the spline, each at a distance from the
    // nearest knot equal to the mean distance of the other knots

    double meanDistance = xRange / (double)numberOfPolynomials;

    knotsForCalculations = vector<double>(numberOfKnots+2*g,0);
    for (int i=0; i<g; ++i)
        knotsForCalculations[i] = knots[0] + (double)(i-g)*meanDistance;
    for (int i=0; i<numberOfKnots; ++i)
        knotsForCalculations[i+g] = knots[i];
    for (int i=1; i<m; ++i)
        knotsForCalculations[i+g+numberOfPolynomials] =
            knots.back()+(double)i*meanDistance;

}



void Spline::calculateCoefficients() {

    // Calculates the basis functions
    auto basisFunctions = vector<BasisFunction>(K);
    for (int j=0; j<K; ++j)
        basisFunctions[j].calculateCoefficients(j,knotsForCalculations);

    // Calculates the Fi matrix
    auto Fi = vector<vector<double>>(n,vector<double>(K,0));
    for (int i=0; i<n; ++i)
        for (int j=0; j<K; ++j)
            Fi[i][j] = basisFunctions[j].D0(abscissae[i]);

    // Finds the limits for the non-zero elements in Fi
    auto firstInFi = vector<int>(n,0);
    auto lastInFi = vector<int>(n,0);
    for (int i=0; i<n; ++i)
        for (int j=0; j<K; ++j)
            if (Fi[i][j] != 0) {
                firstInFi[i] = j;
                break;
            }
    for (int i=0; i<n; ++i)
        for (int j=G; j>-1; --j)
            if (Fi[i][j] != 0) {
                lastInFi[i] = j;
                break;
            }

    // Finds the limits for the non-zero elements in FiT
    auto firstInFiT = vector<int>(K,0);
    auto lastInFiT = vector<int>(K,0);
    for (int j=0; j<K; ++j)
        for (int i=0; i<n; ++i)
            if (Fi[i][j] != 0) {
                firstInFiT[j] = i;
                break;
            }
    for (int j=0; j<K; ++j)
        for (int i=n-1; i>-1; --i)
            if (Fi[i][j] != 0) {
                lastInFiT[j] = i;
                break;
            }

    // Finds the limits for the non-zero elements in M, R and FiTFi
    auto firstInBandMatrices = vector<int>(K,0);
    auto lastInBandMatrices = vector<int>(K,G);
    if (K > m)
        for (int i=m; i<K; ++i)
            firstInBandMatrices[i] = i-g;
    if (K > m)
        for (int i=0; i<K-m; ++i)
            lastInBandMatrices[i] = i+g;

    // Calculates the FiTFi matrix, equal to the product of FiT and Fi
    auto FiTFi = vector<vector<double>>(K,vector<double>(K,0));
    for (int i=0; i<K; ++i)
        for (int j=firstInBandMatrices[i]; j<=lastInBandMatrices[i]; ++j)
            for (int k=0; k<n; ++k)
                FiTFi[i][j] += Fi[k][i] * Fi[k][j];

    // Calculates the R matrix
    auto R = vector<vector<double>>(K,vector<double>(K,0));
    for (int i=0; i<K; ++i)
        for (int j=i; j<=lastInBandMatrices[i]; ++j)
            R[i][j] = basisFunctions[i].integralOfProductD2(basisFunctions[j]);
    for (int i=1; i<K; ++i)
        for (int j=firstInBandMatrices[i]; j<i; ++j)
            R[i][j] = R[j][i];

    // Calculates the FiTy vector
    auto FiTy = vector<double>(K,0);
    for (int i=0; i<K; ++i)
        for (int k=firstInFiT[i]; k<=lastInFiT[i]; ++k)
            FiTy[i] += Fi[k][i] * ordinates[k];

    // Estimates the first derivatives of the experimental data. Contains an
    // additional 0 at position 0
    auto estimatedD1 = vector<double>(n-1,0);
    for (int i=1; i<n-1; ++i)
        estimatedD1[i] =
            (ordinates[i+1]-ordinates[i-1]) / (abscissae[i+1]-abscissae[i-1]);

    // Calculates the square root of the sum of squares of the elements of FiTFi
    double indexFiTFi = 0;
    for (int i=0; i<K; ++i)
        for (int j=firstInBandMatrices[i]; j<=lastInBandMatrices[i]; ++j)
            indexFiTFi += FiTFi[i][j] * FiTFi[i][j];
    indexFiTFi = sqrt(indexFiTFi);

    // Calculates the square root of the sum of squares of the elements of R
    double indexR = 0;
    for (int i=0; i<K; ++i)
        for (int j=firstInBandMatrices[i]; j<=lastInBandMatrices[i]; ++j)
            indexR += R[i][j] * R[i][j];
    indexR = sqrt(indexR);

    // Calculates the value of log10lambda for which the elements of FiTFi and R
    // have the same order of magnitude, rounded to the nearest 0.5
    log10lambdaForSameOrderOfMagnitude =
        round(2.*(log10(indexFiTFi)-log10(indexR)))/2.;

    // Calculates the end points of the log10lambda minimization interval
    log10lambdaMin = log10lambdaForSameOrderOfMagnitude-(double)lambdaSearchInterval/2.;
    log10lambdaMax = log10lambdaForSameOrderOfMagnitude+(double)lambdaSearchInterval/2.;

    // Calculates the log10 of the distance between two consecutive steps in the
    // for cycle for minimizing log10lambda
    double log10lambdaStep =
        (double)lambdaSearchInterval/(double)(numberOfStepsLambda-1);

    // Initializes the elements necessary for the minimization
    auto M = vector<vector<double>>(K,vector<double>(K,0));
    auto zed = vector<double>(K,0);
    auto Z = vector<vector<double>>(K,vector<double>(K,0));
    auto Minv = vector<vector<double>>(K,vector<double>(K,0));
    auto MinvFiT = vector<vector<double>>(K,vector<double>(n,0));
    auto splineCoefficientsForVariousLambdas =
        vector<vector<double>>(numberOfStepsLambda,vector<double>(K,0));
    auto GCV1 = vector<double>(numberOfStepsLambda,0);

    // Calculates the spline coefficients and GCV1 for each lambda in the for
    // cycle
    for (int a=0; a<numberOfStepsLambda; ++a) {

        // Obtains the value of lambda from that of a
        lambda = pow(10., log10lambdaMin + (double)a * log10lambdaStep);

        // Calculates the M matrix, sum of FiTFi and the product of Lambda and R
        for (int i=0; i<K; ++i)
            for (int j=firstInBandMatrices[i];j<=lastInBandMatrices[i];++j)
                M[i][j] = FiTFi[i][j] + lambda * R[i][j];

        // Uses the Doolittle decomposition to decompose M, and obtains the
        // triangular matrices L and U (M = L*U). Saves the values of the L and
        // U matrices, except for the main diagonal of L (consisting entirely of
        // ones), in place of the respective values in matrix M
        for (int i=0; i<K; ++i) {
            for (int j=i; j<=lastInBandMatrices[i]; ++j)
                for (int k=0; k<i; ++k)
                    M[i][j] -= M[i][k] * M[k][j];
            for (int j=i+1; j<=lastInBandMatrices[i]; ++j) {
                for (int k=0; k<i; ++k)
                    M[j][i] -= M[j][k] * M[k][i];
                M[j][i] /= M[i][i];
            }
        }

        // Calculates the zed vector, from the expression L*zed = FiTy, using
        // the forward substitution technique
        zed[0] = FiTy[0];
        for (int i=1; i<K; ++i) {
            zed[i] = FiTy[i];
            for (int k=firstInBandMatrices[i]; k<i; ++k)
                zed[i] -= M[i][k] * zed[k];
        }

        // Calculates the spline coefficients from R*coefficients = zed
        splineCoefficientsForVariousLambdas[a][G] = zed[G] / M[G][G];
        for (int i=G-1; i>-1; --i) {
            splineCoefficientsForVariousLambdas[a][i] = zed[i];
            for (int k=lastInBandMatrices[i]; k>i; --k)
                splineCoefficientsForVariousLambdas[a][i] -=
                M[i][k] * splineCoefficientsForVariousLambdas[a][k];
            splineCoefficientsForVariousLambdas[a][i] /= M[i][i];
        }

        // Solves L*Z = I
        for (int j=0; j<K; ++j) {
            Z[j][j] = 1.;
            for (int i=j+1; i<K; ++i) {
                Z[i][j] = 0;
                for (int k=firstInBandMatrices[i]; k<i; ++k)
                    Z[i][j] -= M[i][k] * Z[k][j];
            }
        }

        // Solves R*Minv = Z
        for (int j=G; j>-1; --j) {
            Minv[G][j] = Z[G][j] / M[G][G];
            for (int i=G-1; i>-1; --i) {
                if (j < i)
                    Minv[i][j] = Z[i][j];
                else if (j > i)
                    Minv[i][j] = 0;
                else
                    Minv[i][j] = 1.;
                for (int k=lastInBandMatrices[i]; k>i; --k)
                    Minv[i][j] -= M[i][k] * Minv[k][j];
                Minv[i][j] /= M[i][i];
            }
        }

        // Calculates the MinvFiT matrix, equal to the product of Minv and FiT,
        // in the locations necessary for the calculation of the trace of S
        for (int j=0; j<n; ++j)
            for (int i=firstInFi[j]; i<=lastInFi[j]; ++i) {
                MinvFiT[i][j] = 0;
                for (int k=firstInFi[j]; k<=lastInFi[j]; ++k)
                    MinvFiT[i][j] += Minv[i][k] * Fi[j][k];
            }

        // Calculates the numerator of GCV1(Lambda)
        double SSE1 = 0; // Sum of squared errors between yi' and f'(xi)
        for (int i=1; i<n-1; ++i) {
            double difference = estimatedD1[i];
            for (int j=0; j<K; ++j)
                difference -=
                splineCoefficientsForVariousLambdas[a][j] *
                basisFunctions[j].D1(abscissae[i]);
            SSE1 += difference * difference;
        }
        GCV1[a] = (double)n * SSE1;

        // Calculates the trace of matrix S
        double traceS = 0;
        for (int i=0; i<n; ++i)
            for (int k=firstInFi[i]; k<=lastInFi[i]; ++k)
                traceS += Fi[i][k] * MinvFiT[k][i];

        // Calculates GCV1(Lambda)
        GCV1[a] /= ((double)n-traceS)*((double)n-traceS);

    } // End of the for cycle for each lambda

    // Finds the minimum value of GCV1(lambda) in the GCV1 vector, and saves the
    // corresponding spline coefficients to splineCoefficients and the
    // corresponding lambda and log10(lambda) to 'lambda' and 'log10lambda'
    double GCVOne = GCV1[0];
    int index = 0;
    log10lambda = log10lambdaMin;
    lambda = pow(10.,log10lambda);
    for (int i=1; i<numberOfStepsLambda; ++i)
        if (GCV1[i] < GCVOne) {
            GCVOne = GCV1[i];
            index = i;
            log10lambda = log10lambdaMin + (double)i * log10lambdaStep;
            lambda = pow(10.,log10lambda);
        }
    splineCoefficients = splineCoefficientsForVariousLambdas[index];

    // Calculates the coefficients of the polynomials of the spline
    coeffD0 = vector<vector<double>>(numberOfPolynomials,vector<double>(m,0));
    int firstBasis = 0;
    for (int a=g; a<g+numberOfPolynomials; ++a) {
        for (int b=firstBasis; b<firstBasis+m; ++b) {
            for (int c=0; c<m; ++c) {
                coeffD0[a-g][c] += basisFunctions[b].coeffD0[g+firstBasis-b][c]*
                                   splineCoefficients[b];
            }
        }
        ++firstBasis;
    }

    // Calculates the coefficients of the first derivative of the polynomials of
    // the spline
    coeffD1 = vector<vector<double>>(numberOfPolynomials,vector<double>(m,0));
    firstBasis = 0;
    for (int a=g; a<g+numberOfPolynomials; ++a) {
        for (int b=firstBasis; b<firstBasis+m; ++b) {
            for (int c=0; c<g; ++c) {
                coeffD1[a-g][c] += basisFunctions[b].coeffD1[g+firstBasis-b][c]*
                                   splineCoefficients[b];
            }
        }
        ++firstBasis;
    }

    // Calculates the coefficients of the second derivative of the polynomials
    // of the spline
    coeffD2 = vector<vector<double>>(numberOfPolynomials, vector<double>(m,0));
    firstBasis = 0;
    for (int a=g; a<g+numberOfPolynomials; ++a) {
        for (int b=firstBasis; b<firstBasis+m; ++b) {
            for (int c=0; c<g-1; ++c) {
                coeffD2[a-g][c] += basisFunctions[b].coeffD2[g+firstBasis-b][c]*
                                   splineCoefficients[b];
            }
        }
        ++firstBasis;
    }

    // Initializes the 'powers' vector
    powers = vector<double>(m,1);

}
