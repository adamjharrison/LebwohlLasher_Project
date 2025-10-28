from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "LebwohlLasher_cython",
        ["LebwohlLasher_cython.pyx"],
        extra_compile_args=['-fopenmp',
            '-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/'],
        extra_link_args=['-lgomp',
            '-Wl,-rpath,/opt/homebrew/opt/gcc/lib/gcc/current/',
            '-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/'],
        include_dirs=[numpy.get_include()]
    )
]

setup(name="LebwholLasher_cython",
      ext_modules=cythonize(ext_modules),install_requires=['numpy'],)