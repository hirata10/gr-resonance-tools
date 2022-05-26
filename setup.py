from setuptools import Extension, setup
from Cython.Build import cythonize
import os
files = list(filter(lambda x: x.endswith(".c") and not x == "c_to_python.c", os.listdir()))

sourcefiles = ["c_to_python.pyx"] + files

extensions = [Extension("c_to_python", sourcefiles)]

setup(ext_modules=cythonize(extensions))











