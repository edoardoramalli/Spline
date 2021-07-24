from sys import platform
import ctypes
import _ctypes


class CLibrary:

    def __init__(self, lib_path):
        self._ctypes_lib = ctypes.CDLL(lib_path)
        self._handle = self._ctypes_lib._handle

    def __getattr__(self, attr):
        return self._ctypes_lib.__getattr__(attr)

    def __del__(self):

        del self._ctypes_lib

        if platform == "linux" or platform == "linux2" or platform == "darwin":
            _ctypes.dlclose(self._handle)
        elif platform == "win32":  # Windows
            _ctypes.FreeLibrary(self._handle)
        else:
            raise NotImplementedError("Unknown platform '{}'.".format(platform))