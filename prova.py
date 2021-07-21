from ctypes import c_int, c_double, c_int16, c_void_p, c_char_p, byref, cdll


def listToArray(ll):
    return (c_double * len(ll))(*ll)


try:
    c_library = cdll.LoadLibrary('/Users/edoardo/Documents/GitHub/Spline/edo.o')
except OSError:
    raise OSError("Unable to load the system C library")

c_library.edo.argtypes = [c_void_p,
                          c_void_p,
                          c_int]

x = [1, 2, 3, 4]
y = [1, 2, 3, 4]

x_c = listToArray(x)
print(x_c[3])
y_c = listToArray(y)

c_library.edo(x_c, y_c, c_int(len(x)))
