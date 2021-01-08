input_fname = "input1.dat"
wavelengths_count = 3
extinction_count = 2
depolarization_count = 0
wavelengths = [0.355, 0.532, 1.064]
discrepancy_kind = 0 #0 = DiscrKindRMS, 1 = DiscrKindMAXABS
funct_type = 0 # 0 = FunctLogNormal, 1 = FunctPowerLaw
knots_count= 22
r_min, r_max = 0.1, 3.0
threshold = 45.0 #%
generations_count = 400#
params_lo_boundary = [1e-13, 0.15, 0.05, 1.30, 1e-7] 
params_hi_boundary = [1e-9, 0.8, 1.0, 1.69, 0.05] # Hi Params Boundary
Ym = [0.0018,6.38E-04,3.77E-04,0.22953,0.1935]


