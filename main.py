#!/usr/bin/env python
"""
  main.py
  PythonFindSolution

  Created by Константин Шмирко on 2021-01-07.
  Copyright 2021 Константин Шмирко. All rights reserved.

"""

import sys
import importlib
import numpy as np
from matplotlib import pyplot as plt
import nlopt
import libspheroids
from dataclasses import dataclass

@dataclass
class Params:
    n_coefs: int
    extinction_coeff: np.ndarray
    absorbption_coeff: np.ndarray
    backscatter_coeff: np.ndarray
    lidar_depol_ratio: np.ndarray
    calc_vector: np.ndarray
    meas_vector: np.ndarray
    back_lidar_ratio: np.ndarray
    
    def __init__(self, N):
        self.n_coefs = N
        self.extinction_coeff = np.zeros(self.n_coefs, dtype='float64')
        self.absorbption_coeff = np.zeros_like(self.extinction_coeff)
        self.backscatter_coeff = np.zeros_like(self.extinction_coeff)
        self.lidar_depol_ratio = np.zeros_like(self.extinction_coeff)
        self.calc_vector = np.zeros_like(self.extinction_coeff)
        self.meas_vector = np.zeros_like(self.extinction_coeff)
        self.back_lidar_ratio = np.zeros_like(self.extinction_coeff)


def main():
    """
    Основная функцииция
    """
    print()
    print("==========================================================")
    print("= PythonFindSolution                                     =")
    print("= Created by Константин Шмирко on 2021-01-07.            =")
    print("= Copyright 2021 Константин Шмирко. All rights reserved. =")
    print("==========================================================")

    config_modname = sys.argv[1]
    c = load_config(config_modname)
    p = Params(5)
    print(p)
    #x = np.zeros(5, dtype='float64')
    # Начальное решение
    x = 0.5*(c.params_hi_boundary+c.params_lo_boundary)

    # настройки солвера
    #opt = nlopt.opt(nlopt.GN_CRS2_LM, 5)
    opt = nlopt.opt(nlopt.GN_ESCH, 5)
    #opt = nlopt.opt(nlopt.GN_ISRES, 5)
    opt.set_lower_bounds(c.params_lo_boundary)
    opt.set_upper_bounds(c.params_hi_boundary)
    opt.set_min_objective(lambda x, grad: objective_funct(x, grad, c, p))
    opt.set_maxeval(10000)
    opt.set_population(c.generations_count)

    # читаем файл с настройками
    libspheroids.dls_read_input(c.input_fname)
    # выделяем памят для хранения промежуточных массивов
    libspheroids.alloc_dls_array(libspheroids.mo_dls.key,
                                 libspheroids.mo_dls.keyel, 1)

    xopt = opt.optimize(x)
    fval = objective_funct(xopt, None, c, p)

    print()
    print(" ============================== ")
    print(" =   Результаты расчетов.     =")
    print(" ============================== ")
    print()

    print("Fmin  value : ", fval)
    print("Xoptimal    : ", xopt)
    print("meas_vect.  : ", c.meas_vector)
    print("calc_vect.  : ", p.calc_vector)

    print("%13s|%13s|%13s"%("A=meas_vect.","B=calc_vect.","(A-B)/A*100%" ))
    for i, _ in enumerate(c.meas_vector):
        print("%13.3e|%13.3e|%7.2f" % (c.meas_vector[i]*1e9, p.calc_vector[i]*1e9,
                                       (c.meas_vector[i]-p.calc_vector[i])/c.meas_vector[i]*100)
              )

    plt.figure()
    plt.semilogy(np.r_[c.wavelengths, c.wavelengths[:2]], p.calc_vector, 'r.')
    plt.semilogy(np.r_[c.wavelengths, c.wavelengths[:2]], c.meas_vector, 'bo')

    # освобождаем память
    libspheroids.alloc_dls_array(libspheroids.mo_dls.key,
                                 libspheroids.mo_dls.keyel, 2)
    plt.show()
    pass


def load_config(config_modname):
    try:
        config = importlib.import_module(config_modname)
    except Exception as e:
        print(e)
    print_module_content(config)
    return config


def print_module_content(config):
  # печатаем параметры запуска
    print()
    print("============================== ")
    print("= Параметры запуска расчетов = ")
    print("============================== ")
    print()

    for name in dir(config):
        # печатаем только переменные, явно указанные пользователем и
        # игнорируем приватные переменные
        if not name.startswith('__'):
            tmpval = getattr(config, name)
            if isinstance(tmpval, list):
                tmpval = np.array(tmpval)

            print("{0:30s}\t=\t{1}".format(name, tmpval))
            if name == 'meas_vector':
                tmpval = tmpval*1e-9
            setattr(config, name, tmpval)
    print()


def objective_funct(x, grad, c, p):
    """
    x     - вектор параметров
    grad  - Якобиан
    c     - параметры из файла
    p     - параметры для расчета оптических свойств
    """

    # Если выбрали логнормальное распределение, то
    # число параметров в векторе x должно быть равно 5
    if c.funct_type == 0 and len(x) == 5:

        rr, ar, ac = libspheroids.sizedis2(-c.knots_count, [x[0]],
                                           [x[1]], [x[2]], c.r_min, c.r_max)
        libspheroids.mo_dls.rn.flat[0] = x[3]
        libspheroids.mo_dls.rk.flat[0] = x[4]
        libspheroids.mo_dls.sd[:] = ar[:]
        libspheroids.mo_dls.rrr[:] = rr[:]
    else:
        raise Exception("Неверный funct_type или размер вектроа x")

    b_print = False

    # для каждой длины волны из нашего списка вычисляем коэффициент
    # обратого рассеяния и ослабления, а также деполяризаионное отношение
    for i, wl in enumerate(c.wavelengths):
        if libspheroids.mo_dls.ndp == 0:
            print("LOADING DATABASE OF SPHEROIDS...")
            b_print = True

        libspheroids.mo_dls.wl = wl

        libspheroids.optchar(libspheroids.mo_dls.ndp)

        if b_print == True:
            print("...DONE")
            b_print = False
        p.extinction_coeff[i] = libspheroids.mo_dls.xext
        p.absorbption_coeff[i] = libspheroids.mo_dls.xabs
        p.backscatter_coeff[i] = libspheroids.mo_dls.xext / libspheroids.mo_dls.xblr
        p.lidar_depol_ratio[i] = libspheroids.mo_dls.xldr
        p.calc_vector[i] = p.backscatter_coeff[i]
        p.back_lidar_ratio[i] = libspheroids.mo_dls.xblr

    A, B = c.wavelengths_count, c.extinction_count+c.wavelengths_count

    for i in range(A, B):
        p.calc_vector[i] = p.extinction_coeff[i-A]

    # инициализируем новую переменную
    func_val = 0.0
    if c.discrepancy_kind == 0:
        func_val = np.sum(((np.log(p.calc_vector)-np.log(c.meas_vector))/np.log(c.meas_vector))**2)
        func_val = np.sqrt(func_val/(c.wavelengths_count +
                                     c.extinction_count))
    func_val = func_val * 100
    # print(func_val)
    return func_val


if __name__ == '__main__':
  # Вызов главной функции
    main()
