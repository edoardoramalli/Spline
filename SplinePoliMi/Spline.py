from ctypes import *
import numpy as np
import sys
import os
import subprocess


def listToArray(ll):
    return (c_double * len(ll))(*ll)


class Spline:
    binariesFileName = 'SplineGenerator.o'

    @staticmethod
    def compileBinaries(module_path, compiler='g++'):
        print('Compiling binaries...')
        flags_compiler = '-std=c++17 -O3 -Wall -DNDEBUG'
        input_main = os.path.join(module_path, 'main.cpp')
        output_exec = os.path.join(module_path, Spline.binariesFileName)
        subprocess.check_call(f'{compiler} {flags_compiler} {input_main} -o {output_exec}', shell=True)

    def checkSettings(self):
        # TODO inserire controlli settings in ingresso
        pass

    def __init__(self, x: list, y: list,
                 verbose: bool = False,
                 m: int = 4, g: int = 3, lambdaSearchInterval: int = 6, numberOfStepsLambda: int = 13,
                 numberOfRatiolkForAICcUse: int = 40, fractionOfOrdinateRangeForAsymptoteIdentification: float = 0.005,
                 fractionOfOrdinateRangeForMaximumIdentification: float = 0.025,
                 possibleNegativeOrdinates: bool = False, removeAsymptotes: bool = True, graphPoints: int = 500,
                 criterion: str = 'AIC'
                 ):

        self.module_path = os.path.dirname(sys.modules[self.__module__].__file__)

        if not os.path.isfile(os.path.join(self.module_path, self.binariesFileName)):
            self.compileBinaries(module_path=self.module_path)

        if len(x) != len(y):
            raise ValueError('X and Y have different lengths!')

        if len(x) <= 1:
            raise ValueError('X and Y need more points!')

        # Input x,y has to be sorted by x
        self.x, self.y = zip(*sorted(zip(x, y)))

        self.verbose = verbose

        # Settings
        self.m = m
        self.g = g
        self.lambdaSearchInterval = lambdaSearchInterval
        self.numberOfStepsLambda = numberOfStepsLambda
        self.numberOfRatiolkForAICcUse = numberOfRatiolkForAICcUse
        self.fractionOfOrdinateRangeForAsymptoteIdentification = fractionOfOrdinateRangeForAsymptoteIdentification
        self.fractionOfOrdinateRangeForMaximumIdentification = fractionOfOrdinateRangeForMaximumIdentification
        self.possibleNegativeOrdinates = possibleNegativeOrdinates
        self.removeAsymptotes = removeAsymptotes
        self.graphPoints = graphPoints
        self.criterion = criterion

        # Backwards
        self.numberOfKnots = None
        self.numberOfPolynomials = None
        self.knots = None
        self.coeffD0 = None
        self.coeffD1 = None
        self.coeffD2 = None

        # Start
        self.computeSpline()

    def computeSpline(self):
        try:
            c_library = cdll.LoadLibrary(os.path.join(self.module_path, self.binariesFileName))
        except OSError:
            raise OSError("Unable to load the system C library")

        c_library.compute_spline_cpp.argtypes = [c_void_p,  # x
                                                 c_void_p,  # y
                                                 c_int,  # length of x, y
                                                 c_void_p,  # numberOfKnots
                                                 c_void_p,  # numberOfPolynomials
                                                 c_void_p,  # coeffDO
                                                 c_void_p,  # coeffD1
                                                 c_void_p,  # coeffD2
                                                 c_void_p,  # knots
                                                 c_bool,  # verbose
                                                 c_int,  # m
                                                 c_int,  # g
                                                 c_int,  # lambdaSearchInterval
                                                 c_int,  # numberOfStepsLambda
                                                 c_int,  # numberOfRatiolkForAICcUse
                                                 c_double,  # fractionOfOrdinateRangeForAsymptoteIdentification
                                                 c_double,  # fractionOfOrdinateRangeForMaximumIdentification
                                                 c_bool,  # possibleNegativeOrdinates
                                                 c_bool,  # removeAsymptotes
                                                 c_int,  # graphPoints
                                                 c_char_p,  # criterion
                                                 ]

        x_c = listToArray(self.x)
        y_c = listToArray(self.y)
        numberOfKnots_c = c_int()
        numberOfPolynomials_c = c_int()
        coeffD0_c = (len(self.x) * self.m * c_double)()
        coeffD1_c = (len(self.x) * self.m * c_double)()
        coeffD2_c = (len(self.x) * self.m * c_double)()
        knots_c = (len(self.x) * c_double)()

        c_library.compute_spline_cpp(x_c,  # x
                                     y_c,  # y
                                     c_int(len(self.x)),  # length of x, y
                                     pointer(numberOfKnots_c),  # numberOfKnots
                                     pointer(numberOfPolynomials_c),  # numberOfPolynomials
                                     pointer(coeffD0_c),  # coeffDO
                                     pointer(coeffD1_c),  # coeffD1
                                     pointer(coeffD2_c),  # coeffD2
                                     pointer(knots_c),  # knots
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
                                     c_bool(self.possibleNegativeOrdinates),  # possibleNegativeOrdinates
                                     c_bool(self.removeAsymptotes),  # removeAsymptotes
                                     c_int(self.graphPoints),  # graphPoints
                                     c_char_p(self.criterion.encode('utf-8')),  # criterion
                                     )

        self.numberOfKnots = numberOfKnots_c.value
        self.numberOfPolynomials = numberOfPolynomials_c.value
        self.coeffD0 = np.reshape(np.array(coeffD0_c[0: self.m * self.numberOfPolynomials]),
                                  (self.numberOfPolynomials, self.m))
        self.coeffD1 = np.reshape(np.array(coeffD1_c[0: self.m * self.numberOfPolynomials]),
                                  (self.numberOfPolynomials, self.m))
        self.coeffD2 = np.reshape(np.array(coeffD2_c[0: self.m * self.numberOfPolynomials]),
                                  (self.numberOfPolynomials, self.m))
        self.knots = np.array(knots_c[0: self.numberOfKnots])

    def compute(self, x, k, coeff):
        indexOfPolynomial = 0

        for i in range(0, self.numberOfKnots - 1):
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

    def evaluate(self, x: int, der: int = 0):
        if not any([self.numberOfKnots, self.knots, self.m, self.coeffD0, self.coeffD1, self.coeffD2]):
            raise ValueError('Spline is not computed yet!')
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
        return self.compute(x=x, k=k, coeff=coeff)
