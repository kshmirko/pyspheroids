import numpy as np
from dataclasses import dataclass


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