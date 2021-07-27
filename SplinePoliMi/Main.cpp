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

5. FW (Future Work): Massimi? --> mettere punti su parabola + solo guardando dei pointi posso sapere se è un max o è un outlier
6. FW: Filtraggio --> chooseKnots, scegliere i nodi in base alla densita di punti nella zona
7. FW: Outlier detection. Calcolo spline con e senza outlier e magari utilizzo il modello
8. FW: Plateau: se riesco a riconoscere un plateau potrei mettere piatta la spline

*/


/*

1. calcutesRoots ha il problema che è limitato al 3° grado. E' stato rimosso e viene fatto in python.
Veniva usato per removeNegativeSegments e yAndAsymptoteAnalysis. In removeNegativeSegments veniva usato per trovare
gli zeri della spline e poi vedere se la spline è negativa. In yAndAsymptoteAnalysis per trovare gli zeri della derivata
prima per lo studio dei massimi e minimi

2. yAndAsymptoteAnalysis è stata rimossa. Il suo studio dei massimi venivano usati da normalizeCoefficients e removeAsymptotes
che non usiamo più vedi dopo.

3. normalizeCoefficients lo abbiamo rimosso perchè la normalizzazione veniva fatto ma non succedeva niente in quanto la spline
veniva normalizzata con se stessa. Aveva senso quando si faceva tutto insieme, i.e. normalizzo spline modello rispetto a
quella sperimentale

4. removeAsymptotes è stata rimossa perchè non in realtà non rimuove gli asintoti ma solo dei nodi e rispettivi coefficienti.
In realtà la spline può sempre essere calcolata anche in punti dove non si hanno nodi vicini. E' da fare poi in python
un ragionamento per intersezione domini e rimozione asintoti troppo grandi

5. extendSpline è stata tolta perchè serviva solo per quella del modello nel momento del curve matching per estendere il
dominio.

6. findMaximaBetweenExtremes è stata tolta perchè non veniva usato più da nessuno in questo momento.

7. updateVariables. Tolta perchè aggiornava numberOfKnots, K, G etc.. ma dal momento che non aggiungiamo più nodi
(venivano aggiunti da removeNegativeSegments per togliere parti negative, nodi messi a 0) o tolti da removeAsymptotes

8. removeNegativeSegments sul mac crashava quindi l'abbiamo messo in python... Aveva il problema ereditato da calculateRoots per
spline di grado maggiore di 3


*/


/*
    x and y have to be sorted and without duplicates on the x-axis.
*/
extern "C"
int compute_spline_cpp(double* x, double* y, int length, int splineType,
            int* numberOfKnots, int* numberOfPolynomials,
            double* coeffDO, double* coeffD1, double* coeffD2, double* knots,
            bool verbose,
            int g_, int lambdaSearchInterval_, int numberOfStepsLambda_, int numberOfRatiolkForAICcUse_,
            double fractionOfOrdinateRangeForAsymptoteIdentification_,
            double fractionOfOrdinateRangeForMaximumIdentification_,
            int graphPoints_, char* criterion_){


    // ----------  SET VARIABLE  ----------

    vector<double> x_vector(x, x + length);
    vector<double> y_vector(y, y + length);


    // ----------  SET SETTINGS  ----------

    g = g_;
    m = g + 1;
    lambdaSearchInterval = lambdaSearchInterval_;
    numberOfStepsLambda = numberOfStepsLambda_;
    numberOfRatiolkForAICcUse = numberOfRatiolkForAICcUse_;
    fractionOfOrdinateRangeForAsymptoteIdentification = fractionOfOrdinateRangeForAsymptoteIdentification_;
    fractionOfOrdinateRangeForMaximumIdentification = fractionOfOrdinateRangeForMaximumIdentification_;
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

