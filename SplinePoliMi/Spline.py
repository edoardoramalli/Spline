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
    moduleVersion = '_0.0.0.9'
    binariesFileName = f'SplineGenerator{moduleVersion}.o'
    criterion_list = ["AIC", "BIC", "SSE"]
    possibleSplineType = [0, 1]

    @staticmethod
    def compileBinaries(module_path, compiler='g++'):
        print('Compiling binaries...')
        flags_compiler = '-std=c++17 -shared -fPIC'
        input_main = os.path.join(module_path, 'main.cpp')
        output_exec = os.path.join(module_path, Spline.binariesFileName)
        subprocess.check_call(f'{compiler} {flags_compiler} {input_main} -o {output_exec}', shell=True)

    def checkSettings(self):
        if self.splineType not in self.possibleSplineType:
            raise ValueError("The selected splineType doesn't exist")
        if self.criterion not in self.criterion_list:
            raise ValueError("The selected criterion doesn't exist")
        if self.m <= 0:
            raise ValueError("m cannot be less or equal than zero")
        if self.g <= 0:
            raise ValueError("g cannot be less or equal than zero")
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
                 m: int = 4, g: int = 3, lambdaSearchInterval: int = 6, numberOfStepsLambda: int = 13,
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
        :param m:
        :param g:
        :param lambdaSearchInterval:
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
        self.m = m
        self.g = g
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
        self.numberOfPolynomials = None
        self.knots = None
        self.coeffD0 = None
        self.coeffD1 = None
        self.coeffD2 = None

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
                                                 c_int,  # m
                                                 c_int,  # g
                                                 c_int,  # lambdaSearchInterval
                                                 c_int,  # numberOfStepsLambda
                                                 c_int,  # numberOfRatiolkForAICcUse
                                                 c_double,  # fractionOfOrdinateRangeForAsymptoteIdentification
                                                 c_double,  # fractionOfOrdinateRangeForMaximumIdentification
                                                 # c_bool,  # possibleNegativeOrdinates
                                                 # c_bool,  # removeAsymptotes
                                                 c_int,  # graphPoints
                                                 c_char_p,  # criterion
                                                 ]

        c_library.compute_spline_cpp.restype = c_int

        x_c = listToArray(self.x)
        y_c = listToArray(self.y)
        # TODO se non cambiano i nodi (tolti o aggiunti) si può conoscere numberOfKnots e si può togliere!
        numberOfKnots_c = c_int()
        numberOfPolynomials_c = c_int()

        size_coeff_matrix = (len(self.x) * self.m)
        size_knots = len(self.x)

        coeffD0_c = [0] * size_coeff_matrix
        coeffD1_c = [0] * size_coeff_matrix
        coeffD2_c = [0] * size_coeff_matrix
        knots_c = [0] * size_knots
        coeffD0_c = (size_coeff_matrix * c_double)(*coeffD0_c)
        coeffD1_c = (size_coeff_matrix * c_double)(*coeffD1_c)
        coeffD2_c = (size_coeff_matrix * c_double)(*coeffD2_c)
        knots_c = (size_knots * c_double)(*knots_c)

        c_library.compute_spline_cpp(x_c,  # x
                                     y_c,  # y
                                     c_int(len(self.x)),  # length of x, y
                                     c_int(self.splineType),  # splineType
                                     pointer(numberOfKnots_c),  # numberOfKnots
                                     pointer(numberOfPolynomials_c),  # numberOfPolynomials
                                     byref(coeffD0_c),  # coeffDO
                                     byref(coeffD1_c),  # coeffD1
                                     byref(coeffD2_c),  # coeffD2
                                     byref(knots_c),  # knots
                                     # coeffD0_c,  # coeffDO
                                     # coeffD1_c,  # coeffD1
                                     # coeffD2_c,  # coeffD2
                                     # knots_c,  # knots
                                     c_bool(self.verbose),  # verbose
                                     c_int(self.m),  # m
                                     c_int(self.g),  # g
                                     c_int(self.lambdaSearchInterval),  # lambdaSearchInterval
                                     c_int(self.numberOfStepsLambda),  # numberOfStepsLambda
                                     c_int(self.numberOfRatiolkForAICcUse),  # numberOfRatiolkForAICcUse
                                     c_double(self.fractionOfOrdinateRangeForAsymptoteIdentification),
                                     # fractionOfOrdinateRangeForAsymptoteIdentification
                                     c_double(self.fractionOfOrdinateRangeForMaximumIdentification),
                                     # fractionOfOrdinateRangeForMaximumIdentification
                                     # c_bool(self.possibleNegativeOrdinates),  # possibleNegativeOrdinates
                                     # c_bool(self.removeAsymptotes),  # removeAsymptotes
                                     c_int(self.graphPoints),  # graphPoints
                                     c_char_p(self.criterion.encode('utf-8')),  # criterion
                                     )

        self.numberOfPolynomials = numberOfPolynomials_c.value
        self.coeffD0 = np.reshape(np.array(deepcopy(coeffD0_c[0: self.m * self.numberOfPolynomials])),
                                  (self.numberOfPolynomials, self.m))
        self.coeffD1 = np.reshape(np.array(deepcopy(coeffD1_c[0: self.m * self.numberOfPolynomials])),
                                  (self.numberOfPolynomials, self.m))
        self.coeffD2 = np.reshape(np.array(deepcopy(coeffD2_c[0: self.m * self.numberOfPolynomials])),
                                  (self.numberOfPolynomials, self.m))
        self.knots = np.array(deepcopy(knots_c[0: numberOfKnots_c.value]))

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

        for i in range(0, len(self.knots) - 1):
            if x > self.knots[i]:
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
        if not any([len(self.knots), self.m, self.coeffD0, self.coeffD1, self.coeffD2]):
            raise ValueError('Spline is not computed yet!')
        # numpy.float64 etc..?
        # if not any([(lambda element: type(element) == float or type(element) == int)(e) for e in x]):
        #     raise ValueError('Value of x is not valid!')
        if der == 0:
            k = self.m
            coeff = self.coeffD0
        elif der == 1:
            k = self.g
            coeff = self.coeffD1
        elif der == 2:
            k = self.g - 1
            coeff = self.coeffD2
        else:
            raise ValueError('Derivative does not exists!')
        return self.compute(x, k, coeff) if not hasattr(x, '__iter__') else [self.compute(e, k, coeff) for e in x]

    def removeAsymptotes(self):
        pass

    def removeNegativeSegments(self):

        roots = []

        for i in range(self.numberOfPolynomials):
            # TODO capire cosa accade se m > 3
            a0 = self.coeffD0[i][3]
            b0 = self.coeffD0[i][2]
            c0 = self.coeffD0[i][1]
            d0 = self.coeffD0[i][0]
            roots += list(np.roots([a0, b0, c0, d0]))

        for root in roots:
            if not np.iscomplex(root):
                index = 0
                root = float(np.real(root))
                for index_knots in range(len(self.knots)):
                    if root < self.knots[index_knots]:
                        break
                    else:
                        index += 1
                if 1 <= index < len(self.knots):
                    self.knots = np.insert(self.knots, index, root)
                    self.coeffD0 = np.insert(self.coeffD0, index, self.coeffD0[index - 1], axis=0)
                    self.coeffD1 = np.insert(self.coeffD1, index, self.coeffD1[index - 1], axis=0)
                    self.coeffD2 = np.insert(self.coeffD2, index, self.coeffD2[index - 1], axis=0)

        for index_knot in range(len(self.knots) - 1):
            midpoint = (self.knots[index_knot] + self.knots[index_knot + 1]) / 2.0
            if self.evaluate(midpoint, der=0) <= 0:
                zeros = np.zeros(self.m)
                self.coeffD0[index_knot] = zeros
                self.coeffD1[index_knot] = zeros
                self.coeffD2[index_knot] = zeros
