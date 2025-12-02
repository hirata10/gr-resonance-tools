from setuptools import Extension, setup
from Cython.Build import cythonize
import os

# Directory containing setup.py
here = os.path.dirname(os.path.abspath(__file__))

# Path to C directory
c_dir = os.path.join(here, "..", "C")

files = [
    "kerrgwem.c", "kerrmode.c", "kerrtraj.c",
    "Gamma.c", "J_dot.c", "J2J_dot.c", "resonance_find.c"
]

# Use relative paths instead of absolute paths
c_sources = [os.path.join("..", "C", f) for f in files]

sourcefiles = ["gr_wrapper.pyx"] + c_sources

extensions = [
    Extension(
        "gr_wrapper",
        sourcefiles,
        include_dirs=[c_dir],  # headers location (this can be absolute)
    )
]

setup(ext_modules=cythonize(extensions))