from setuptools import Extension, setup
from Cython.Build import cythonize
import os
files = list(filter(lambda x: x.endswith(".c") and not x == "gr_wrapper.c", os.listdir()))

sourcefiles = ["gr_wrapper.pyx"] + files

extensions = [Extension("gr_wrapper", sourcefiles)]

setup(ext_modules=cythonize(extensions))











