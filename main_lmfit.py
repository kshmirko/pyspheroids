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

from lmfit import minimize, Parameters, fit_report

import libspheroids
from dataclasses import dataclass
import contextlib

import matplotlib.ticker as ticker
import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"})
    

@contextlib.contextmanager
def printoptions(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    try:
        yield
    finally: 
        np.set_printoptions(**original)


class Config:
    pass

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
        self.calc_vector = np.zeros(self.n_coefs*2, dtype='float64')
        self.meas_vector = np.zeros_like(self.extinction_coeff)
        self.back_lidar_ratio = np.zeros_like(self.extinction_coeff)


def ln_funct(x, lnN, sigma, rm):
    return np.exp(lnN)/(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5*((np.log(x)-np.log(rm))/sigma)**2)

def load_config_yaml(yaml_fname):
    """
    Читаем файл yaml и создаем объект,
    добвляем ему поля по ключам в читаемом файле
    """
    
    from yaml import load
    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader
    
    config=Config()
    try:
        with open(yaml_fname, 'rb') as fin:
            ctx = load(fin, Loader=Loader)
    except IOError as e:
        print(e)
    finally:
        for key, val in ctx.items():
            setattr(config, key, val)
    return config

def prepare_dataset(c):
    """
    Подготовка данных:
    1.  вычисление параметра ангстрема и выполнение инерполяции 
        коэффициентов обратного рассеяния на 6 точек
    2.  вычисление параметра ангстрема и выполнение инерполяции 
        коэффициентов ослабления на 6 точек
    
    Итого вместо 5 уравнений будем иметь 12 (по 6 на кадый коэффициент)
    """
    setattr(c,'c1',np.polyfit(np.log(c.wavelengths), np.log(c.meas_vector[:3]), deg=1))
    setattr(c,'c2',np.polyfit(np.log(c.wavelengths[:2]), np.log(c.meas_vector[-2:]), deg=1))
    y1 = np.exp(np.polyval(c.c1, np.log(c.interp_wavelengths)))
    y2 = np.exp(np.polyval(c.c2, np.log(c.interp_wavelengths)))
    setattr(c,'meas_vector_interp', np.r_[y1, y2])
    
    
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
                tmpval = tmpval
            setattr(config, name, tmpval)
    print()


def display_solution(c, p, xopt):
    """docstring for display_solution"""
    grid = plt.GridSpec(3,2)
    plt.subplot(grid[0,0])
    plt.semilogy(c.interp_wavelengths, p.calc_vector[:6], 'b')
    plt.semilogy(c.interp_wavelengths, c.meas_vector_interp[:6],'bo')
    plt.ylabel(r'Back. coeff. 1/km $\times$ sr', color='b')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    
    #print(c.wavelengths[:2], p.calc_vector[-2:], c.meas_vector[-2:])
    plt.subplot(grid[1,0])
    plt.semilogy(c.interp_wavelengths, p.calc_vector[6:],'g')
    plt.semilogy(c.interp_wavelengths, c.meas_vector_interp[6:], 'go')
    plt.ylabel(r'Ext. coeff. 1/km', color='g')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    
    
    plt.subplot(grid[2,0])
    plt.plot(c.interp_wavelengths, p.calc_vector[6:]/p.calc_vector[:6], 'r')
    plt.plot(c.interp_wavelengths, c.meas_vector_interp[6:]/c.meas_vector_interp[:6], 'ro')
    plt.ylabel('LR ratio sr', color='r')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    plt.xlabel('Wavelength, nm')
    
    plt.subplot(grid[0,1])
    plt.plot(c.interp_wavelengths, p.lidar_depol_ratio, 'kd-')
    plt.ylabel(r'lid. dep. rat., \%')
    plt.xlabel('Wavelength, nm')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    
    rr = np.logspace(np.log10(c.r_min), np.log10(c.r_max), 50)
    plt.subplot(grid[1:,1])
    plt.semilogx(rr, ln_funct(rr, xopt.params['lnN1'].value, xopt.params['sigma1'].value, xopt.params['rm1'].value), 'b-')
    plt.semilogx(rr, ln_funct(rr, xopt.params['lnN2'].value, xopt.params['sigma2'].value, xopt.params['rm2'].value), 'g-')
    plt.semilogx(libspheroids.mo_dls.rrr[:c.knots_count], libspheroids.mo_dls.sd[:c.knots_count], 'ro', mfc=None)
    plt.ylabel(r'$\frac{dV}{d\ln r}$, $\frac{\mu m^3}{\mu m^3}$')
    plt.xlabel(r'Radius, $\mu m$')
    plt.grid(which='major')
    plt.grid(which='minor', linestyle=':')
    
    plt.suptitle('Lidar measurements processing results', fontsize="x-large")
    plt.tight_layout()
    print(c.interp_wavelengths)
    print(p.lidar_depol_ratio)
    print(p.back_lidar_ratio)
    plt.show()


def objective_funct(params, c, p):
    """
    x     - вектор параметров
    c     - конфигурационные параметры из файла (config.py)
    p     - параметры для расчета оптических свойств (структура)
    """

    # Если выбрали логнормальное распределение, то
    # число параметров в векторе x должно быть равно 5
    lnN1 = params['lnN1']
    sigma1 = params['sigma1']
    rm1 = params['rm1']
    lnN2 = params['lnN2']
    sigma2 = params['sigma2']
    rm2 = params['rm2']
    rn = params['rn']
    rk = params['rk']
    
    if c.funct_type == 0:

        rr, ar, ac = libspheroids.sizedis2(-c.knots_count, [np.exp(lnN1), np.exp(lnN2)],
                                            [sigma1, sigma2], [rm1, rm2], c.r_min, c.r_max)
        # rr1, ar1, ac1 = libspheroids.sizedis2(-c.knots_count, [np.exp(lnN1) ],
        #                                    [sigma1], [rm1], c.r_min, c.r_max)
        # _, ar2, _ = libspheroids.sizedis2(-c.knots_count, [np.exp(lnN2) ],
        #                                    [sigma2], [rm2], c.r_min, c.r_max)
        libspheroids.mo_dls.rn.flat[0] = rn
        libspheroids.mo_dls.rk.flat[0] = np.exp(rk)
        libspheroids.mo_dls.sd[:] = ar[:]
        libspheroids.mo_dls.rrr[:] = rr[:]
        
    else:
        raise Exception("Неверный funct_type или размер вектроа x")

    b_print = False

    # для каждой длины волны из нашего списка вычисляем коэффициент
    # обратого рассеяния и ослабления, а также деполяризаионное отношение
    for i, wl in enumerate(c.interp_wavelengths):
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
        
        p.back_lidar_ratio[i] = libspheroids.mo_dls.xblr
     
    p.calc_vector[:] = np.r_[p.backscatter_coeff[:], p.extinction_coeff[:]]

    # инициализируем новую переменную
    #func_val = 0.0
    #if c.discrepancy_kind == 0:
    #    func_val = np.sum(((np.log(p.calc_vector)-np.log(c.meas_vector))/np.log(c.meas_vector))**2)
    #    func_val = np.sqrt(func_val/(c.wavelengths_count +
    #                                 c.extinction_count))
    #func_val = func_val * 100
    func_val = ((p.calc_vector-c.meas_vector_interp))/(c.meas_vector_interp)
    return func_val
    
    

def main():
    """
    Основная функция
    """
    print()
    print("==========================================================")
    print("= PythonFindSolution                                     =")
    print("= Created by Константин Шмирко on 2021-01-07.            =")
    print("= Copyright 2021 Константин Шмирко. All rights reserved. =")
    print("==========================================================")

    config_modname = sys.argv[1]
    c = load_config_yaml(config_modname)
    print_module_content(c)
    
    p = Params(6)
    
    
    # настройки солвера
    # Параметры
    params = Parameters()
    params.add('lnN1', value=0.5*(c.lnN1_lo+c.lnN1_hi), min = c.lnN1_lo, max=c.lnN1_hi,
                brute_step = 0.25)
    params.add('sigma1', value=0.5*(c.sigma1_lo+c.sigma1_hi), min = c.sigma1_lo, max=c.sigma1_hi,
                brute_step=0.065)
    params.add('rm1', value=0.5*(c.rm1_lo+c.rm1_hi), min = c.rm1_lo, max=c.rm1_hi, 
                brute_step=0.109)
    params.add('lnN2', value=0.5*(c.lnN2_lo+c.lnN2_hi), min = c.lnN2_lo, max=c.lnN2_hi,
                brute_step = 0.25)
    params.add('sigma2', value=0.5*(c.sigma2_lo+c.sigma2_hi), min = c.sigma2_lo, max=c.sigma2_hi,
                brute_step=0.065)
    params.add('rm2', value=0.5*(c.rm2_lo+c.rm2_hi), min = c.rm2_lo, max=c.rm2_hi, 
                brute_step=0.109)
    params.add('rn', value=0.5*(c.rn_lo+c.rn_hi), min = c.rn_lo, max=c.rn_hi,
                brute_step=0.5)
    params.add('rk', 0.5*(c.lnrk_lo+c.lnrk_hi), min = c.lnrk_lo, max=c.lnrk_hi,
                brute_step=2.6)
    
    #аргументы
    args = (c, p)
    

    # читаем файл с настройками
    libspheroids.dls_read_input(c.input_fname)
    # выделяем памят для хранения промежуточных массивов
    libspheroids.alloc_dls_array(libspheroids.mo_dls.key,
                                 libspheroids.mo_dls.keyel, 1)

    data_count = c.meas_data.shape[0]
    for i, v in enumerate(c.meas_data):
        print('Prosessing data {elem}/{tot}...'.format(elem=i, tot=data_count), file=sys.stderr)
        
        c.meas_vector = v
        # подготоавливаем данные к расчетам
        prepare_dataset(c)
    
        if c.solver['name']=='dual_annealing':
            xopt = minimize(objective_funct, params, args=args, method=c.solver['name'], 
                            no_local_search=c.solver['no_local_search'])
        else:
            xopt = minimize(objective_funct, params, args=args, method=c.solver['name'])
    
        print("000000000000000000000000000000000000000000000000000")
        print("Results of processing {elem} dataset".format(elem=i))
        
        print(fit_report(xopt))
        fval = objective_funct(xopt.params,  c, p)
    

        print()
        print(" ==============================")
        print(" =   Результаты расчетов.     =")
        print(" ==============================")
        print()
    
        with printoptions(formatter={'float': '{: 0.2e}'.format}):
            print("Fmin  value (R^2): ", (fval*fval).sum())
            print()
            print("Optimal parameters:")
            print("-------------------")
        
            xopt.params.pretty_print()
            print()
            print("{0:20s}|{1:20s}|{2:20s}".format("meas_vector","meas_vect. interp.","calc_vect."))
            for i,_ in enumerate(p.calc_vector):
                if i<len(c.meas_vector):
                    print("{0:20.3e} {1:20.3e} {2:20.3e}".format(c.meas_vector[i], 
                        c.meas_vector_interp[i], p.calc_vector[i]))
                else:
                    print("{0:>20s} {1:20.3e} {2:20.3e}".format("--", 
                        c.meas_vector_interp[i], p.calc_vector[i]))
            # print("meas_vect.           : ", c.meas_vector)
    #         print("meas_vect. interp.   : ", c.meas_vector_interp[:6])
    #         print("                     : ", c.meas_vector_interp[6:])
    #         print("calc_vect.           : ", p.calc_vector)
    
        print()
        print("%13s|%13s|%13s"%("A=meas_vect.","B=calc_vect.","(A-B)/A*100%" ))
        for i, _ in enumerate(c.meas_vector_interp):
            print("%13.3e|%13.3e|%7.2f" % (c.meas_vector_interp[i], p.calc_vector[i],
                                            (c.meas_vector_interp[i]-p.calc_vector[i])/c.meas_vector_interp[i]*100)
                                            )

    

    # освобождаем память
    libspheroids.alloc_dls_array(libspheroids.mo_dls.key,
                                 libspheroids.mo_dls.keyel, 2)
    display_solution(c, p, xopt)
    
    

if __name__ == '__main__':
  # Вызов главной функции
    main()
