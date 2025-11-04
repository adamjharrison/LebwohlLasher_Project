from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "LebwohlLasher_cython_mpi",
        ["LebwohlLasher_cython_mpi.pyx"],
        extra_compile_args=['-O3'],
        include_dirs=[numpy.get_include()]
    )
]

setup(name="LebwholLasher_cython_mpi",
      ext_modules=cythonize(ext_modules),install_requires=['numpy'],)
