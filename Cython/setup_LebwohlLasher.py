from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "LebwohlLasher_cython",
        ["LebwohlLasher_cython.pyx"],
        extra_compile_args=['-O3'],
        include_dirs=[numpy.get_include()]
    )
]

setup(name="LebwholLasher_cython",
      ext_modules=cythonize(ext_modules),install_requires=['numpy'],)
