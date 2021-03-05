# pyspheroids

обертка для программного кода для расчета оптических свойств сфероидов

я используем только ту часть кода, которая соответствует параметрам

key=2, key_org = 0, key_SD=0 или 1, ID=0 или 3

остальную часть кода я попросту удалил для улучшения читабельности 
программы исходных текстов программы.

В программе два основных компонента: модуль libspheroids, который устанавливается 
автоматически по команде `python setup.py install` и код программы восстановления 
main.py.



# Состав пакета libspheroids

Этот модуль активно использует глобальные переменные, которые f2py обернул в 
соответствующую область видимости

### Functions:
	
	alloc_dls_array(key,keyel,key_alloc) - эта функция пердазначена для выделения и удаления динамической памяти
	для хранения матриц с коэффициентами рассеяния.
  	
	matrix_fix(kn1,grid1,wavel,kre,kim,are,aim,ratio,ndp)
  
	dlnr1 = usmatrix(ndp)
  
	usuf_ls(kn1,grid1,wavel,kre,kim,are,aim)
  
	usu_ls(ratio,nratn,kn1,grid1,wavel,kre,kim,are,aim)
  
	usu_ls_rd(ratio,nratn,kn1,grid1,wavel,kre,kim,are,aim)
  
	optchar(ndp)
  
	dls_read_input(fullpath_input)
  
	rrr,ar,ac = sizedisdn(kn,ia,id,cm,sm,rmm,rmin,rmax,knpar,idbg,nmd=len(cm))
  
	ac = sdnorm(c,s,rm,rmin,rmax)
  
	slog = slog(c,s,rm,rr)
  
	sizedisdn1(kn,ia,id,cm,sm,rmm,rmin,rmax,rrr,ar,ac,idbg,nmd=len(cm),knpar=len(rrr))
  
	rrr,ar,ac = sizedis2(kn,cm,sm,rmm,rmin,rmax,nmd=len(cm))


### COMMON blocks:

	/us1/ us11(181,41)
	
	/us2/ us12(181,41)
  
	/us3/ us22(181,41)
  
	/us4/ us33(181,41)
  
	/us5/ us34(181,41)
  
	/us6/ us44(181,41)
  
	/us0/ usea(2,41)
	
### Fortran 90/95 modules:
	каждый из нижеприведенных модулей содержит определенное количество используемых глобальных переменных
	
  	mo_par_dls - kn1par,kr1par,km1par,krepar,kimpar,knpar,krpar,kmpar,kmd  
	
	mo_dls - key,key_rd,keyel,keysub,keyls,key_org,key_fx,key_grid1,key_rd1,kn,km,kr,nratn,ndp,wl,rn,rk,pomin,pomax,xext,xabs,xsca,albedo,r,grid,sd,rd,f11,f12,f22,f33,f34,f44,angle,xblr,xldr,distname_o,distname_f,distname_n,comm_name,key_sd,id,nmd,nsd,ksd,cm,sm,rmm,rrr,ar,xgrid,ac  
	
	alloc1 - u11,u12,u22,u33,u34,u44,uea  
	
	alloc - uf11,uf12,uf22,uf33,uf34,uf44,ufea  
	
	mo_intrpl_linear - linear()  
	
	mo_intrpl_spline - intrpl_spline(),e01baf(),e02baf(),e02bbf(),p01abf(),x04aaf(),x04baf(),p01abz()  
	
	phase_func - sint(),sint_spline(),asypar().


# Состав пакета main.py

main_lmfit.py - это уже программа для выполнения расчетов


Вызов программы осуществляется следующим образом:
```
❯ ./main_lmfit.py config.yaml

==========================================================
= PythonFindSolution                                     =
= Created by Константин Шмирко on 2021-01-07.            =
= Copyright 2021 Константин Шмирко. All rights reserved. =
==========================================================

============================== 
= Параметры запуска расчетов = 
============================== 

depolarization_count          	=	0
discrepancy_kind              	=	0
extinction_count              	=	2
funct_type                    	=	0
input_fname                   	=	input1.dat
interp_wavelengths            	=	[0.355 0.4   0.532 0.6   0.8   1.064]
knots_count                   	=	22
lnN1_hi                       	=	-20.0
lnN1_lo                       	=	-30
lnN2_hi                       	=	-20.0
lnN2_lo                       	=	-30
lnrk_hi                       	=	-3
lnrk_lo                       	=	-16
meas_data                     	=	[[1.313150e-03 1.055597e-03 1.295030e-04 5.153500e-02 9.548333e-03]
 [3.530000e-02 1.710000e-02 6.050000e-03 4.190000e-01 2.840000e-01]
 [6.900000e-03 5.640000e-03 3.200000e-03 1.530000e-01 1.220000e-01]
 [4.300000e-03 1.629406e-03 1.437189e-03 8.598000e-02 2.442000e-02]]
meas_vector                   	=	[0.0043     0.00162941 0.00143719 0.08598    0.02442   ]
r_max                         	=	5.0
r_min                         	=	0.05
rm1_hi                        	=	0.3
rm1_lo                        	=	0.05
rm2_hi                        	=	0.8
rm2_lo                        	=	0.31
rn_hi                         	=	1.69
rn_lo                         	=	1.3
sigma1_hi                     	=	0.8
sigma1_lo                     	=	0.15
sigma2_hi                     	=	0.8
sigma2_lo                     	=	0.15
solver                        	=	{'name': 'powel', 'no_local_search': False}
threshold                     	=	45.0
wavelengths                   	=	[0.355 0.532 1.064]
wavelengths_count             	=	3

LOADING DATABASE OF SPHEROIDS...
 Volume mixture of spheroids
...DONE
[[Fit Statistics]]
    # fitting method   = Powell
    # function evals   = 2579
    # data points      = 12
    # variables        = 8
    chi-square         = 1.32453814
    reduced chi-square = 0.33113453
    Akaike info crit   = -10.4461139
    Bayesian info crit = -6.56686072
##  Warning: uncertainties could not be estimated:
    sigma1:  at boundary
    sigma2:  at boundary
    rk:      at boundary
[[Variables]]
    lnN1:   -25.2691270 (init = -25)
    sigma1:  0.15000001 (init = 0.475)
    rm1:     0.06164973 (init = 0.175)
    lnN2:   -27.5815596 (init = -25)
    sigma2:  0.15000065 (init = 0.475)
    rm2:     0.79866649 (init = 0.555)
    rn:      1.67745524 (init = 1.495)
    rk:     -15.9999990 (init = -9.5)

 ==============================
 =   Результаты расчетов.     =
 ==============================

Fmin  value (R^2):  1.3245381712423643

Optimal parameters:
-------------------
Name       Value      Min      Max   Stderr     Vary     Expr Brute_Step
lnN1      -25.27      -30      -20     None     True     None     0.25
lnN2      -27.58      -30      -20     None     True     None     0.25
rk           -16      -16       -3     None     True     None      2.6
rm1      0.06165     0.05      0.3     None     True     None    0.109
rm2       0.7987     0.31      0.8     None     True     None    0.109
rn         1.677      1.3     1.69     None     True     None      0.5
sigma1      0.15     0.15      0.8     None     True     None    0.065
sigma2      0.15     0.15      0.8     None     True     None    0.065

meas_vector         |meas_vect. interp.  |calc_vect.          
           4.300e-12            3.406e-12            3.556e-12
           1.629e-12            3.055e-12            3.097e-12
           1.437e-12            2.357e-12            1.794e-12
           8.598e-11            2.113e-12            1.344e-12
           2.442e-11            1.626e-12            5.768e-13
                  --            1.254e-12            2.684e-13
                  --            8.598e-11            9.937e-11
                  --            5.931e-11            6.890e-11
                  --            2.442e-11            2.686e-11
                  --            1.680e-11            1.743e-11
                  --            6.862e-12            6.839e-12
                  --            2.825e-12            3.369e-12

 A=meas_vect.| B=calc_vect.| (A-B)/A*100%
    3.406e-03|    3.556e-03|  -4.42
    3.055e-03|    3.097e-03|  -1.35
    2.357e-03|    1.794e-03|  23.89
    2.113e-03|    1.344e-03|  36.36
    1.626e-03|    5.768e-04|  64.53
    1.254e-03|    2.684e-04|  78.61
    8.598e-02|    9.937e-02| -15.58
    5.931e-02|    6.890e-02| -16.17
    2.442e-02|    2.686e-02| -10.00
    1.680e-02|    1.743e-02|  -3.78
    6.862e-03|    6.839e-03|   0.33
    2.825e-03|    3.369e-03| -19.23
[0.355 0.4   0.532 0.6   0.8   1.064]
[2.59813571 2.36147404 2.20619869 2.35963774 3.78077698 7.28293514]
[27.94269371 22.24932861 14.97393227 12.96582603 11.8564806  12.55254745]

```
# Отрисовка
Отрисовка графиков предполагается с использованием пакета matplotlib.


# Основные настройки программы
все настройки программы определяются файлом config.yaml. 

Его структура следующая

```
input_fname: 			"input1.dat"
wavelengths_count:		3
extinction_count:		2
depolarization_count:	0
wavelengths:			[0.355, 0.532, 1.064]
discrepancy_kind:		0  
funct_type:				0  
knots_count:			22
r_min:					0.05
r_max:					5.0
threshold:				45.0  # %

solver:
  name: powel
  no_local_search: False

lnN1_lo:	-30
lnN1_hi: 2.0
sigma1_lo: 0.15
sigma1_hi: 0.8
rm1_lo: 0.05
rm1_hi: 0.3
lnN2_lo:	-40
lnN2_hi: 2.0
sigma2_lo: 0.15
sigma2_hi: 0.8
rm2_lo: 0.31
rm2_hi: 0.8
rn_lo: 1.3
rn_hi: 1.45
lnrk_lo: -16
lnrk_hi: -3
interp_wavelengths:		[0.355, 0.4, 0.532, 0.6, 0.8, 1.064]
meas_vector:			[3.53E-02,	1.71E-02,	6.05E-03,	4.19E-01,	2.84E-01]
plot_solution: False
meas_data:
  - [0.00131315,	0.001055597,	0.000129503,	0.051535,	0.009548333]
  - [3.53E-02,	1.71E-02,	6.05E-03,	4.19E-01,	2.84E-01]
  - [6.90E-03,	5.64E-03,	3.20E-03,	1.53E-01,	1.22E-01]
  - [0.0043,	0.001629406,	0.001437189,	0.08598,	0.02442]

```

`input_fname`  - определяет настройки самого алгоритма расчета микрофизических свойств, там частиц, их характеристики

`wavelengths_count` - количество длин волн

`extinction_count` - количество коэффициентов экстинкии используемых в решении задачи

`depolarization_count`  - количество коэффициентов деполяризации, учавствуюших в решении

`wavelengths` - список длин волн в порядке возрастания

`discrepancy_kind` - способ вычисления невязки 0|1

`func_type` - тип апроксимирующей функции 0|1

`knots_count` - число узловых точек для построения решения (22)






