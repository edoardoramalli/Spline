#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <iomanip>

using namespace std;

#include "Settings.h"
#include "BasisFunction.h"
#include "Utilities.h"
#include "Spline.h"
#include "ComputeSpline.h"

/*
                                TODO LIST

1. abscissa e originalAbscissae: c'è differenza. Se non c'è elimina
    --> originalAbscissae accesso solo in lettura, abscissa viene modificata ed usata per calcolare knots
2. Ricorda che quando calcoli spline modello c'è da fare estensione e normalizzazione rispetto spline exp
3. Cosa succede ai dati in input? ordinati? media?
    --> Ordinati, media, media errori, rimozione asintoti (non solo roba a 0) Riga 670 Indexes.h originale
4. Pre-rimozione asintoti dai dati numerici --> Noi lo dobbiamo fare in python, teniamo un po di asintoto
    Ma se lo teniamo puoi Spline.removeAsymptotes() lo rimuove lo stesso?
5. FW (Future Work): Massimi? --> parabola + solo guardando dei pointi posso sapere se è un max o è un outlier
6. FW: Filtraggio
7. FW: Outlier detection. Calcolo spline con e senza outlier e magari utilizzo il modello

*/


/*
    x and y have to be sorted and without duplicates on the x-axis.
*/
extern "C"
int compute_spline_cpp(double* x, double* y, int length, int splineType,
            int* numberOfKnots, int* numberOfPolynomials,
            double* coeffDO, double* coeffD1, double* coeffD2, double* knots,
            bool verbose,
            int m_, int g_, int lambdaSearchInterval_, int numberOfStepsLambda_, int numberOfRatiolkForAICcUse_,
            double fractionOfOrdinateRangeForAsymptoteIdentification_, double fractionOfOrdinateRangeForMaximumIdentification_,
                         int graphPoints_, char* criterion_){
//            bool possibleNegativeOrdinates_, bool removeAsymptotes_,


    // ----------  SET VARIABLE  ----------

    vector<double> x_vector(x, x + length);
    vector<double> y_vector(y, y + length);


    // ----------  SET SETTINGS  ----------

    m = m_;
    g = g_;
    lambdaSearchInterval = lambdaSearchInterval_;
    numberOfStepsLambda = numberOfStepsLambda_;
    numberOfRatiolkForAICcUse = numberOfRatiolkForAICcUse_;
    fractionOfOrdinateRangeForAsymptoteIdentification = fractionOfOrdinateRangeForAsymptoteIdentification_;
    fractionOfOrdinateRangeForMaximumIdentification = fractionOfOrdinateRangeForMaximumIdentification_;
//    possibleNegativeOrdinates = possibleNegativeOrdinates_;
//    removeAsymptotes = removeAsymptotes_;
    graphPoints = graphPoints_;
    criterion = string(criterion_);


    // ----------  COMPUTE BEST SPLINE  ----------

    vector<Spline> possibleSplines = calculateSplines(x_vector, y_vector, splineType);

    int index_best = calculateBestSpline(possibleSplines, criterion);

    Spline best_spline = possibleSplines[index_best];


    // ---------- VERBOSE ----------
    if(verbose){
        vector<vector<double>> tmp;
        cout << "Spline Type: " << splineType << endl;

        cout << "Original X: ";
        printV_inLine(x_vector);
        cout << "Original Y: ";
        printV_inLine(y_vector);
        cout << endl;

        cout << "Spline X: ";
        printV_inLine(best_spline.abscissae);
        cout << "Spline Y: ";
        printV_inLine(best_spline.ordinates);
        cout << endl;

        cout << "KNOTS: ";
        printV_inLine(best_spline.knots);
        cout << endl;

        cout << "CoeffD0:" << endl;
        printM(best_spline.coeffD0);
        cout << "CoeffD1:" << endl;
        printM(best_spline.coeffD1);
        cout << "CoeffD2:" << endl;
        printM(best_spline.coeffD2);
        cout << "SETTINGS:" << endl;
        printSettings();

        cout << endl;

        cout << "D0:" << endl;
        tmp = evaluateSpline(best_spline, 0);
        cout << "\tx: ";
        printV_inLine(tmp[0]);
        cout << "\ty: ";
        printV_inLine(tmp[1]);

        cout << "D1:" << endl;
        tmp = evaluateSpline(best_spline, 1);
        cout << "\tx: ";
        printV_inLine(tmp[0]);
        cout << "\ty: ";
        printV_inLine(tmp[1]);

        cout << "D2:" << endl;
        tmp = evaluateSpline(best_spline, 2);
        cout << "\tx: ";
        printV_inLine(tmp[0]);
        cout << "\ty: ";
        printV_inLine(tmp[1]);

    }

    // ----------  PASS BACK THE RESULTS  ----------

    *numberOfKnots = best_spline.knots.size();
    *numberOfPolynomials = best_spline.numberOfPolynomials;

    for(int i = 0; i < (int)best_spline.coeffD0.size(); i++){
        for(int j = 0; j < (int)best_spline.coeffD0[i].size(); j++){
            coeffDO[i * best_spline.coeffD0[i].size() + j] = best_spline.coeffD0[i][j];
            coeffD1[i * best_spline.coeffD0[i].size() + j] = best_spline.coeffD1[i][j];
            coeffD2[i * best_spline.coeffD0[i].size() + j] = best_spline.coeffD2[i][j];
        }
    }

    for(int i = 0; i < (int)best_spline.knots.size(); i++){
        knots[i] = best_spline.knots[i];
    }

    possibleSplines.clear();
    possibleSplines.shrink_to_fit();

    return 0;
}

int main() {
//    vector<vector<double>> splineD0;
//    vector<vector<double>> splineD1;
//
//    cout << "\nRunning Spline Calculations\n";
//
//    // qui scrivo io
//    // DATI DI INPUT
//    // NB!!!!! devono essere ordinate rispetto x
//    vector<double> x = {0.0 ,0.1, 0.2, 0.4};
//    vector<double> y = {97,98, 99, 100};
//
//    vector<Spline> possibleSplines = calculateSplines(x,y,removeAsymptotes);
//    int index_best = calculateBestSpline(x,y,criterion, possibleSplines);
//    splineD0 = evaluateBestSplineD0(possibleSplines[index_best]);
//    splineD1 = evaluateBestSplineD1(possibleSplines[index_best]);

    return 0;

}

