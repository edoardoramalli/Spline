#include "Settings.h"

/* Prints a vector */
auto printV = [](const auto& vector) {

    for (int i=0; i< (int)vector.size(); ++i)
        cout << "(" << i << ") " << vector[i] << endl;
    cout << endl;

};

auto printV_inLine = [](const auto& vector) {

    cout << "[";
    for (int i=0; i< (int)vector.size(); ++i)
        cout << vector[i] << ", ";
    cout << "]" << endl;

};

/* Print settings */
void printSettings(){
    cout << "m:  " << m << endl;
    cout << "g:  " << g << endl;
    cout << "lambdaSearchInterval:  " << lambdaSearchInterval << endl;
    cout << "numberOfStepsLambda:  " << numberOfStepsLambda << endl;
    cout << "numberOfRatiolkForAICcUse:  " << numberOfRatiolkForAICcUse << endl;
    cout << "fractionOfOrdinateRangeForAsymptoteIdentification:  " << fractionOfOrdinateRangeForAsymptoteIdentification << endl;
    cout << "fractionOfOrdinateRangeForMaximumIdentification:  " << fractionOfOrdinateRangeForMaximumIdentification << endl;
//    cout << "possibleNegativeOrdinates:  " << possibleNegativeOrdinates << endl;
//    cout << "removeAsymptotes:  " << removeAsymptotes << endl;
    cout << "graphPoints:  " << graphPoints << endl;
    cout << "criterion:  " << criterion << endl;

    // NB: pascalsTriangle is not printed

}



/* Prints a matrix */
auto printM = [](const auto& matrix) {

    for (int i=0; i< (int)matrix.size(); ++i) {
        for (int j=0; j< (int)matrix[i].size(); ++j)
            cout << "(" << i << "," << j << ") " << matrix[i][j] << "\t    ";
        cout << endl;
    }
    cout << endl;

};



/* Prints the matrix equal to the difference of the two input matrices */
void print(const vector<vector<double>>& A, const vector<vector<double>>& B) {

    for (int i = 0; i< (int)A.size(); ++i) {
        for (int j = 0; j< (int)A[i].size(); ++j)
            cout << "(" << i << "," << j << ") " << A[i][j]-B[i][j] << "\t    ";
        cout << endl;
    }
    cout << endl;

}



/* Prints the vector equal to the difference of the two input vectors */
void print(const vector<double>& a, const vector<double>& b) {

    for (int i=0; i< (int)a.size(); ++i)
        cout << "(" << i << ") " << a[i]-b[i] << endl;
    cout << endl;

}



/* Prints the GCV1 vector and the corresponding log10lambda values */
void printGCV1(const vector<double>& GCV1,
               double log10lambdaMin,
               double log10lambdaStep) {

    for (int i=0; i< (int)GCV1.size(); ++i) {
        cout << "(" << log10lambdaMin+i*log10lambdaStep << ")\t " << GCV1[i];
        if (i != (int)GCV1.size()-1) {
            if (GCV1[i] > GCV1[i+1]) cout << "   \t ++";
            else if (GCV1[i] < GCV1[i+1]) cout << "   \t--";
            else cout << "   \tequal!";
        }
        cout << endl;
    }
    cout << endl;

}



/* Finds the minimum and the maximum value of a matrix */
void minMax(const vector<vector<double>>& matrix) {

    double min = matrix[0][0];
    int rowMin = 0;
    int columnMin = 0;
    double max = matrix[0][0];
    int rowMax = 0;
    int columnMax = 0;
    for (int i=0; i< (int)matrix.size(); ++i) {
        for (int j=0; j< (int)matrix[i].size(); ++j) {
            if (matrix[i][j] < min) {
                min = matrix[i][j];
                rowMin = i;
                columnMin = j;
            }
            if (matrix[i][j] > max) {
                max = matrix[i][j];
                rowMax = i;
                columnMax = j;
            }
        }
    }
    cout << "Min: (" << rowMin << "," << columnMin << ") " << min << endl
         << "Max: (" << rowMax << "," << columnMax << ") " << max << endl
         << endl;

}



/* Finds the minimum and the maximum value of the difference of the two input
matrices */
void minMax(const vector<vector<double>>& A, const vector<vector<double>>& B) {

    auto matrix = vector<vector<double>>(A.size(),vector<double>(A[0].size()));
    for (int i=0; i< (int)A.size(); ++i)
        for (int j=0; j < (int)A[0].size(); ++j)
            matrix[i][j] = A[i][j] - B[i][j];

    double min = matrix[0][0];
    int rowMin = 0;
    int columnMin = 0;
    double max = matrix[0][0];
    int rowMax = 0;
    int columnMax = 0;
    for (int i=0; i < (int)matrix.size(); ++i) {
        for (int j=0; j < (int)matrix[i].size(); ++j) {
            if (matrix[i][j] < min) {
                min = matrix[i][j];
                rowMin = i;
                columnMin = j;
            }
            if (matrix[i][j] > max) {
                max = matrix[i][j];
                rowMax = i;
                columnMax = j;
            }
        }
    }
    cout << "Min: (" << rowMin << "," << columnMin << ") " << min << endl
         << "Max: (" << rowMax << "," << columnMax << ") " << max << endl
         << endl;

}



/* Calculates the inverse of matrix A using the Gauss-Jordan method. Modifies
matrix A while running */
void invertWithGaussJordan(vector<vector<double>>& A,
                           vector<vector<double>>& I) {

    int K = A.size();

    // Matrix I starts as an identity matrix with the same dimensions as A
    I = vector<vector<double>>(K,vector<double>(K));
    for (int i=0; i<K; ++i)
        I[i][i] = 1.;

    for (int i=0; i<K; ++i) {
        double alfa = A[i][i];
        for (int j=0; j<K; ++j) {
            A[i][j] /= alfa;
            I[i][j] /= alfa;
        }
        for (int k=0; k<K; ++k) {
            if (k != i && A[k][i] != 0) {
                double beta = A[k][i];
                for (int j=0; j<K; ++j) {
                    A[k][j] -= beta * A[i][j];
                    I[k][j] -= beta * I[i][j];
                }
            }
        }
    }

}

