from ctypes import *
import numpy as np
import sys
import os
import subprocess
from statistics import mean
from .CLibrary import CLibrary
from copy import deepcopy

c_float_p = POINTER(c_double)


def listToArray(ll):
    return (c_double * len(ll))(*ll)


class Spline:
    moduleVersion = '_0.0.0.10'
    binariesFileName = f'SplineGenerator{moduleVersion}.o'
    criterion_list = ["AIC", "BIC", "SSE"]
    possibleSplineType = [0, 1]

    @staticmethod
    def compileBinaries(module_path, compiler='g++'):
        print('Compiling binaries...')
        flags_compiler = '-std=c++17 -shared -fPIC -O3 -Wall -DNDEBUG'
        input_main = os.path.join(module_path, 'main.cpp')
        output_exec = os.path.join(module_path, Spline.binariesFileName)
        subprocess.check_call(f'{compiler} {flags_compiler} {input_main} -o {output_exec}', shell=True)

    def checkSettings(self):
        if self.splineType not in self.possibleSplineType:
            raise ValueError("The selected splineType doesn't exist")
        if self.criterion not in self.criterion_list:
            raise ValueError("The selected criterion doesn't exist")
        if not 0 <= self._g <= 6:
            raise ValueError("g must stay between 0 and 6")
        if self.lambdaSearchInterval <= 0:
            raise ValueError("lambdaSearchInterval cannot be less or equal than zero")
        if self.numberOfStepsLambda <= 0:
            raise ValueError("numberOfStepsLambda cannot be less or equal than zero")
        if self.numberOfRatiolkForAICcUse <= 0:
            raise ValueError("numberOfRatiolkForAICcUse cannot be less or equal than zero")
        if not 0 <= self.fractionOfOrdinateRangeForAsymptoteIdentification <= 5:
            raise ValueError("fractionOfOrdinateRangeForAsymptoteIdentification cannot be less or equal than zero")
        if self.fractionOfOrdinateRangeForMaximumIdentification <= 0:
            raise ValueError("fractionOfOrdinateRangeForMaximumIdentification cannot be less or equal than zero")
        if self.graphPoints <= 0:
            raise ValueError("graphPoints cannot be less or equal than zero")

    def filterInputData(self):
        """
        Filter the x, y in input.
        Set originalX and originalY s.t. originalX has unique values and its sorted
        """
        if len(self.originalX) != len(self.originalY):
            raise ValueError('X and Y have different lengths!')

        if len(self.originalX) <= 1:
            raise ValueError('X and Y need more points!')

        # TODO do this stuff in python or in c++?

        diz_data = {}

        for index, i in enumerate(self.originalX):
            diz_data[i] = diz_data.get(i, []) + [self.originalY[index]]

        # Dicts preserve insertion order in Python 3.7+
        # Ensure unique (do the average in case) and sorted x
        diz_data = {k: mean(v) for k, v in sorted(diz_data.items(), key=lambda item: item[0])}

        self.x, self.y = list(diz_data.keys()), list(diz_data.values())

    def __init__(self, x: list, y: list,
                 verbose: bool = False,
                 splineType: int = 0,
                 g: int = 3, lambdaSearchInterval: int = 6, numberOfStepsLambda: int = 13,
                 numberOfRatiolkForAICcUse: int = 40, fractionOfOrdinateRangeForAsymptoteIdentification: float = 0.005,
                 fractionOfOrdinateRangeForMaximumIdentification: float = 0.025,
                 possibleNegativeOrdinates: bool = False, removeAsymptotes: bool = False, graphPoints: int = 500,
                 criterion: str = 'AIC'
                 ):
        """

        :param x: input x-values
        :param y: input y-values
        :param verbose: default is False. If True, ask to the c++ library to print input data, settings, coefficients
        and evaluate D0, D1, D2 in #graphPoints equidistant x-values from the first knots to the last knots
        :param splineType: default 0. 0 means experimental data, 1 means model data. TODO Explain better
        :param g: Degree of the basis function by default 3
        :param lambdaSearchInterval: Lambda is the smoothing parameter for the spline, lambdaSearchInterval is the Orders of magnitude 
        of difference between the smallest and the largest possible value of the smoothing parameter lambda 
        :param numberOfStepsLambda:
        :param numberOfRatiolkForAICcUse:
        :param fractionOfOrdinateRangeForAsymptoteIdentification:
        :param fractionOfOrdinateRangeForMaximumIdentification:
        :param possibleNegativeOrdinates:
        :param removeAsymptotes:
        :param graphPoints:
        :param criterion:
        """
        self.module_path = os.path.dirname(sys.modules[self.__module__].__file__)

        if not os.path.isfile(os.path.join(self.module_path, self.binariesFileName)):
            self.compileBinaries(module_path=self.module_path)

        # Manage Input Data
        self.originalX = x
        self.originalY = y
        self.x = None
        self.y = None
        self.filterInputData()

        self.verbose = verbose

        # Settings
        self._g = g
        self._m = self._g + 1
        self.lambdaSearchInterval = lambdaSearchInterval
        self.numberOfStepsLambda = numberOfStepsLambda
        self.numberOfRatiolkForAICcUse = numberOfRatiolkForAICcUse
        self.fractionOfOrdinateRangeForAsymptoteIdentification = fractionOfOrdinateRangeForAsymptoteIdentification
        self.fractionOfOrdinateRangeForMaximumIdentification = fractionOfOrdinateRangeForMaximumIdentification
        self.graphPoints = graphPoints
        self.criterion = criterion
        self.splineType = splineType

        # Check Settings
        self.checkSettings()

        # Backwards
        self._numberOfPolynomials = None
        self._knots = None
        self._coeffD0 = None
        self._coeffD1 = None
        self._coeffD2 = None

        # Start
        self.computeSpline()

        if not possibleNegativeOrdinates:
            self.removeNegativeSegments()

    def computeSpline(self):
        try:
            c_library = CLibrary(os.path.join(self.module_path, self.binariesFileName))
        except OSError:
            raise OSError("Unable to load the system C library")

        c_library.compute_spline_cpp.argtypes = [c_float_p,  # x
                                                 c_float_p,  # y
                                                 c_int,  # length of x, y
                                                 c_int,  # splineType
                                                 c_float_p,  # numberOfKnots
                                                 c_float_p,  # numberOfPolynomials
                                                 c_float_p,  # coeffDO
                                                 c_float_p,  # coeffD1
                                                 c_float_p,  # coeffD2
                                                 c_float_p,  # knots
                                                 c_bool,  # verbose
                                                 c_int,  # g
                                                 c_int,  # lambdaSearchInterval
                                                 c_int,  # numberOfStepsLambda
                                                 c_int,  # numberOfRatiolkForAICcUse
                                                 c_double,  # fractionOfOrdinateRangeForAsymptoteIdentification
                                                 c_double,  # fractionOfOrdinateRangeForMaximumIdentification
                                                 c_int,  # graphPoints
                                                 c_char_p,  # criterion
                                                 ]

        c_library.compute_spline_cpp.restype = c_int

        x_c = listToArray(self.x)
        y_c = listToArray(self.y)
        # TODO se non cambiano i nodi (tolti o aggiunti) si può conoscere numberOfKnots e si può togliere!
        numberOfKnots_c = c_int()
        numberOfPolynomials_c = c_int()

        size_coeff_matrix = (len(self.x) * self._m)
        size_knots = len(self.x)

        coeffD0_c = (size_coeff_matrix * c_double)()
        coeffD1_c = (size_coeff_matrix * c_double)()
        coeffD2_c = (size_coeff_matrix * c_double)()
        knots_c = (size_knots * c_double)()

        c_library.compute_spline_cpp(x_c,  # x
                                     y_c,  # y
                                     c_int(len(self.x)),  # length of x, y
                                     c_int(self.splineType),  # splineType
                                     pointer(numberOfKnots_c),  # numberOfKnots
                                     pointer(numberOfPolynomials_c),  # numberOfPolynomials
                                     pointer(coeffD0_c),  # coeffDO
                                     pointer(coeffD1_c),  # coeffD1
                                     pointer(coeffD2_c),  # coeffD2
                                     pointer(knots_c),  # knots
                                     c_bool(self.verbose),  # verbose
                                     c_int(self._g),  # g
                                     c_int(self.lambdaSearchInterval),  # lambdaSearchInterval
                                     c_int(self.numberOfStepsLambda),  # numberOfStepsLambda
                                     c_int(self.numberOfRatiolkForAICcUse),  # numberOfRatiolkForAICcUse
                                     c_double(self.fractionOfOrdinateRangeForAsymptoteIdentification),
                                     c_double(self.fractionOfOrdinateRangeForMaximumIdentification),
                                     c_int(self.graphPoints),  # graphPoints
                                     c_char_p(self.criterion.encode('utf-8')),  # criterion
                                     )

        self._numberOfPolynomials = numberOfPolynomials_c.value
        self._coeffD0 = np.reshape(np.array(deepcopy(coeffD0_c[0: self._m * self._numberOfPolynomials])),
                                  (self._numberOfPolynomials, self._m))
        self._coeffD1 = np.reshape(np.array(deepcopy(coeffD1_c[0: self._m * self._numberOfPolynomials])),
                                  (self._numberOfPolynomials, self._m))
        self._coeffD2 = np.reshape(np.array(deepcopy(coeffD2_c[0: self._m * self._numberOfPolynomials])),
                                  (self._numberOfPolynomials, self._m))
        self._knots = np.array(deepcopy(knots_c[0: numberOfKnots_c.value]))

        # Free Memory
        del x_c
        del y_c
        del numberOfKnots_c
        del numberOfPolynomials_c
        del coeffD0_c
        del coeffD1_c
        del coeffD2_c
        del knots_c

        del c_library

    def compute(self, x, k, coeff):
        # TODO ctyhon immplementation?
        indexOfPolynomial = 0

        for i in range(0, len(self._knots) - 1):
            if x > self._knots[i]:
                indexOfPolynomial = i
            else:
                break

        powers = [1] * k

        for i in range(1, k):
            powers[i] = powers[i - 1] * x

        y_val = 0
        for i in range(0, k):
            y_val += coeff[indexOfPolynomial][i] * powers[i]

        return y_val

    def evaluate(self, x, der: int = 0):
        """
        TODO set warning if outside abscissae range
        :param x: x could be a float or int number or a list of them. evaluate the derivative of the spline in x.
        :param der: default 0. 0 means D0, 1 means D1 (first derivative), 2 means D2 (second derivative)
        :return: the evaluated derivative on the x-value(s) x
        """
        if not any([len(self._knots), self._m, self._coeffD0, self._coeffD1, self._coeffD2]):
            raise ValueError('Spline is not computed yet!')
        # numpy.float64 etc..?
        # if not any([(lambda element: type(element) == float or type(element) == int)(e) for e in x]):
        #     raise ValueError('Value of x is not valid!')
        if der == 0:
            k = self._m
            coeff = self._coeffD0
        elif der == 1:
            k = self._g
            coeff = self._coeffD1
        elif der == 2:
            k = self._g - 1
            coeff = self._coeffD2
        else:
            raise ValueError('Derivative does not exists!')
        return self.compute(x, k, coeff) if not hasattr(x, '__iter__') else [self.compute(e, k, coeff) for e in x]

    def removeAsymptotes(self):

        return

        # TODO serve davvero?

        # Evaluations of polynomials left and right

        x_range = np.linspace(min(self.originalX), max(self.originalX), 100)
        y_eval = self.evaluate(x_range, der=0)

        self.yD0range = max(y_eval) - min(y_eval)
        minVariation = self.yD0range * self.fractionOfOrdinateRangeForAsymptoteIdentification
        numberAsymptotesPolynomialsLeft = 0
        numberAsymptotesPolynomialsRight = 0
        yFront = self.evaluate(self.knots[0], der=0)
        yBack = self.evaluate(self.knots[-1], der=0)

        for i in range(1, len(self.knots)):
            y = self.evaluate(self.knots[i], der=0)
            if abs(yFront - y) < minVariation:
                numberAsymptotesPolynomialsLeft += 1
            else:
                break

        for i in range(len(self.knots) - 2, 0, -1):
            y = self.evaluate(self.knots[i], der=0)
            if abs(yBack - y) < minVariation:
                numberAsymptotesPolynomialsRight += 1
            else:
                break

        # REMOVING THE ASYMPTOTES

        # If the spline is completely flat, doesn't remove any segment
        if numberAsymptotesPolynomialsLeft == self.numberOfPolynomials:
            return
        # If there are no horizontal asymptotes, doesn't remove any segment
        elif numberAsymptotesPolynomialsLeft + numberAsymptotesPolynomialsRight == 0:
            return
        elif numberAsymptotesPolynomialsLeft < self.numberOfPolynomials - numberAsymptotesPolynomialsRight:
            newKnots = np.array([])
            newCoeffD0 = np.array([])
            newCoeffD1 = np.array([])
            newCoeffD2 = np.array([])
            print('dooo', numberAsymptotesPolynomialsLeft, self.numberOfPolynomials, numberAsymptotesPolynomialsRight)

            for i in range(numberAsymptotesPolynomialsLeft,
                           self.numberOfPolynomials - numberAsymptotesPolynomialsRight):
                newKnots = np.append(newKnots, self.knots[i])
                newCoeffD0 = np.append(newCoeffD0, self.coeffD0[i])
                newCoeffD1 = np.append(newCoeffD1, self.coeffD1[i])
                newCoeffD2 = np.append(newCoeffD2, self.coeffD2[i])

            newKnots = np.append(newKnots, self.knots[self.numberOfPolynomials - numberAsymptotesPolynomialsRight])

            self.knots = newKnots
            self.coeffD0 = newCoeffD0
            self.coeffD1 = newCoeffD1
            self.coeffD2 = newCoeffD2
            self.numberOfPolynomials = len(self.knots) - 1

            print(self.coeffD0)

    def removeNegativeSegments(self):

        roots = []

        for i in range(self._numberOfPolynomials):
            roots_coeff = np.array([])
            for m_val in range(self._g, -1, -1):
                roots_coeff = np.append(roots_coeff, self._coeffD0[i][m_val])
            roots += list(np.roots(roots_coeff))

        for root in roots:
            if not np.iscomplex(root):
                index = 0
                root = float(np.real(root))
                for index_knots in range(len(self._knots)):
                    if root < self._knots[index_knots]:
                        break
                    else:
                        index += 1
                if 1 <= index < len(self._knots):
                    self._knots = np.insert(self._knots, index, root)
                    self._coeffD0 = np.insert(self._coeffD0, index, self._coeffD0[index - 1], axis=0)
                    self._coeffD1 = np.insert(self._coeffD1, index, self._coeffD1[index - 1], axis=0)
                    self._coeffD2 = np.insert(self._coeffD2, index, self._coeffD2[index - 1], axis=0)

        if len(roots) > 0:
            for index_knot in range(len(self._knots) - 1):
                midpoint = (self._knots[index_knot] + self._knots[index_knot + 1]) / 2.0
                if self.evaluate(midpoint, der=0) <= 0:
                    zeros = np.zeros(self._m)
                    self._coeffD0[index_knot] = zeros
                    self._coeffD1[index_knot] = zeros
                    self._coeffD2[index_knot] = zeros

        self._numberOfPolynomials = len(self._knots)
