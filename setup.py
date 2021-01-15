"""
    Setup file for the libsphrtoids package.
    Copyright Constantine Shmirko (c) 2021
"""
from numpy.distutils.core import Extension

extension = Extension(name='libspheroids',
                      sources=['src/mo_par_DLS.f90',
                               'src/mo_DLS.f90',
                               'src/mo_alloc1.f90',
                               'src/mo_alloc.f90',
                               # 'src/lognormal.f90',
                               'src/mo_intrpl_linear.f90',
                               'src/mo_intrpl_spline.f90',
                               'src/DLS_fixget.f90',
                               'src/DLS_intrpl.f90',
                               'src/DLS_optchr.f90',
                               'src/DLS_read_input.f90',
                               'src/phase_func.f90',
                               'src/sizedstr.f90'
                               ],
                      # extra_f90_compile_args=['-fdefault-real-8'],
                      # extra_f77_compile_args=['-fdefault-real-8'],
                      )

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name='libspheroids',
          description="Library that calculates optical properties" +
          " of spherical aerosols",
          author="Constantine Shmirko",
          author_email="kshmirko@gmail.com",
          ext_modules=[extension]
          )
