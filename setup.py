from setuptools import Extension, setup
from Cython.Build import cythonize
import os

files = ["kerrgwem.c", "kerrmode.c", "kerrtraj.c", "Gamma.c", "J_dot.c", "resonance_find.c"]
#files = list(filter(lambda x: x.endswith(".c") and not x == "gr_wrapper.c", os.listdir()))
print(files)
sourcefiles = ["gr_wrapper.pyx"] + files

extensions = [Extension("gr_wrapper", sourcefiles)]

setup(ext_modules=cythonize(extensions))











