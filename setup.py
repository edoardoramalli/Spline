import os
from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info
from setuptools.command.install import install
import subprocess
import sys



with open("README.md", "r") as fh:
    long_description = fh.read()


def custom_command():
    subprocess.check_call('g++ -std=c++17  ./SplinePoliMi/Main.cpp -o ./SplinePoliMi/SplineGenerator.o -O3 -Wall -DNDEBUG', shell=True)


class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        custom_command()


class CustomDevelopCommand(develop):
    def run(self):
        develop.run(self)
        custom_command()


class CustomEggInfoCommand(egg_info):
    def run(self):
        egg_info.run(self)
        custom_command()

def copy_dir():
    dir_path = 'C++_Library'
    base_dir = os.path.join('MODULE_DIR_HERE', dir_path)
    for (dirpath, dirnames, files) in os.walk(base_dir):
        for f in files:
            yield os.path.join(dirpath.split('/', 1)[1], f)

setup(
    name="SplinePoliMi",
    version="0.0.0.6",
    author="Edoardo Ramalli, Timoteo Dinelli",
    author_email="edoardo.ramalli@polimi.it",
    description="Spline generation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/edoardoramalli/Spline",
    packages=find_packages(),
    install_requires=['numpy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",

    ],
    include_package_data=True,
    package_data={
        'SplinePoliMi': ['*.h', '*.cpp'],
    },

    cmdclass={
        'install': CustomInstallCommand,
        'develop': CustomDevelopCommand,
        'egg_info': CustomEggInfoCommand,
    },
)

# python setup.py bdist_wheel
# python -m twine upload dist/*
