#!/usr/bin/env python
#
#  main.py
#  PythonFindSolution
#
#  Created by Константин Шмирко on 2021-01-07.
#  Copyright 2021 Константин Шмирко. All rights reserved.
#

import sys
import importlib
import numpy as np
import libspheroids
from matplotlib import pyplot as plt
import nlopt



class Params:
  def __init__(self, N):
    self.Ext = np.zeros(N,dtype='float64')
    self.Absb= np.zeros_like(self.Ext)
    self.Bsc= np.zeros_like(self.Ext)
    self.Ldr= np.zeros_like(self.Ext)
    self.Yc= np.zeros_like(self.Ext)
    self.Ym= np.zeros_like(self.Ext)
    self.Blr= np.zeros_like(self.Ext)


def main():
  np.set_printoptions(precision=2, )
  print()
  print("==========================================================")
  print("= PythonFindSolution                                     =")
  print("= Created by Константин Шмирко on 2021-01-07.            =")
  print("= Copyright 2021 Константин Шмирко. All rights reserved. =")
  print("==========================================================")   

  config_modname = sys.argv[1]
  c = load_config(config_modname)
  p = Params(5)
  
  x = np.zeros(5, dtype='float64')
  x = 0.5*(c.params_hi_boundary+c.params_lo_boundary)
  
	# настройки солвера
  opt = nlopt.opt(nlopt.GN_ISRES, 5)
  opt.set_lower_bounds(c.params_lo_boundary)
  opt.set_upper_bounds(c.params_hi_boundary)
  opt.set_min_objective(lambda x, grad: objective_funct(x, grad, c, p))
  opt.set_maxeval(10000)
  opt.set_population(c.generations_count)
  
  # читаем файл с настройками
  libspheroids.dls_read_input(c.input_fname)
  # выделяем памят для хранения промежуточных массивов
  libspheroids.alloc_dls_array(libspheroids.mo_dls.key,
                               libspheroids.mo_dls.keyel,1)
  
  xopt = opt.optimize(x)
  fval = objective_funct(xopt, None, c, p)
  
  print()
  print(" ============================== ")
  print(" =   Результаты расчетов.     =")
  print(" ============================== ")
  print()
  
  print("Finit value : ", fval)
  print("Xopt.       : ", xopt)
  print("Ym          : ", c.Ym)
  print("Yc          : ", p.Yc)
  
  for i, _ in enumerate(c.Ym):
    print("%13.3e %13.3e %7.2f"%(c.Ym[i], p.Yc[i], (c.Ym[i]-p.Yc[i])/c.Ym[i]*100))
  

  
  # освобождаем память
  libspheroids.alloc_dls_array(libspheroids.mo_dls.key,
                               libspheroids.mo_dls.keyel,2)
  pass

def load_config(config_modname):
  try:
    config=importlib.import_module(config_modname)
    
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
    if not name.startswith('__') :
      tmpval = getattr(config, name)
      if isinstance(tmpval, list):
        tmpval = np.array(tmpval)
        
      #if name == 'Ym':
      #  tmpval = tmpval
      setattr(config, name, tmpval)  
      print("{0:30s}\t=\t{1}".format(name, tmpval))
     
    
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
      print("DONE")
      b_print = False
    p.Ext[i] = libspheroids.mo_dls.xext
    p.Absb[i] = libspheroids.mo_dls.xabs
    p.Bsc[i] = libspheroids.mo_dls.xext / libspheroids.mo_dls.xblr
    p.Ldr[i] = libspheroids.mo_dls.xldr
    p.Yc[i] = p.Bsc[i]
    p.Blr[i] = libspheroids.mo_dls.xblr
  
  
  A, B = c.wavelengths_count, c.extinction_count+c.wavelengths_count

  for i in range(A, B):
    p.Yc[i] = p.Ext[i-A]
  
  func_val = 0.0
  if c.discrepancy_kind == 0:
    
    func_val = np.sum((p.Yc-c.Ym)**2)
    func_val = np.sqrt(func_val/(c.wavelengths_count+
                                 c.extinction_count))
  func_val = func_val * 100
  #print(func_val)
  return func_val
  
  
if __name__ == '__main__':
  main()