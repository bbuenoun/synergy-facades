# -*- coding: utf-8 -*-
from __future__ import division
from scipy.optimize import fsolve
import math
from math import exp
from numpy import *
import numpy as np
from scipy.optimize import *
import subprocess
# bombaso
# import matplotlib.pyplot as plt
# import pandas as pd
import os

SB = 5.6693e-8  # [SI] Stefan Boltzmann constant
T_ref = 273.15  # Conversion C --> K
g = 9.81  # [m/s²] Gravitationnal constant
R = 8.314462  #
sigma = 5.6693e-8


class Layer:
    def __init__(self, width, height, emi_inf_f, emi_inf_b, ref_inf_f, ref_inf_b, tra_inf_f, tra_inf_b, solar_tra, solar_ref, solar_abs, thickness, conductivity, gap):
        self.width = width
        self.height = height
        self.emi_inf_f = emi_inf_f
        self.emi_inf_b = emi_inf_b
        self.ref_inf_f = ref_inf_f
        self.ref_inf_b = ref_inf_b
        self.tra_inf_f = tra_inf_f
        self.tra_inf_b = tra_inf_b
        self.solar_tra = solar_tra
        self.solar_ref = solar_ref
        self.solar_abs = solar_abs
        self.thickness = thickness
        self.conductivity = conductivity
        self.gap = gap


def gas_mixture_properties(
    T_gas, air_ratio, argon_ratio, krypton_ratio, xenon_ratio, P=101325
):
    # Calculates the thermal properties of a gas mixture
    # From ISO15099 p. 18/19
    # Attention: T_gas in K
    # 0=Air, 1=Argon, 2=Krypton, 3=Xenon
    lambda_air, mu_air, Cp_air, rho_air, M_air = gas_properties(T_gas, 0)
    lambda_argon, mu_argon, Cp_argon, rho_argon, M_argon = gas_properties(
        T_gas, 1)
    lambda_krypton, mu_krypton, Cp_krypton, rho_krypton, M_krypton = gas_properties(
        T_gas, 2
    )
    lambda_xenon, mu_xenon, Cp_xenon, rho_xenon, M_xenon = gas_properties(
        T_gas, 3)

    # Molecular mass
    molecular_mass_mix = (
        (air_ratio * M_air)
        + (argon_ratio * M_argon)
        + (krypton_ratio * M_krypton)
        + (xenon_ratio * M_xenon)
    )

    # Density
    rho_mix = P * molecular_mass_mix / (R * T_gas * 1000)

    # Specifiec heat
    cpi_air = Cp_air * M_air
    cpi_argon = Cp_argon * M_argon
    cpi_krypton = Cp_krypton * M_krypton
    cpi_xenon = Cp_xenon * M_xenon
    cp_mix = (
        (air_ratio * cpi_air)
        + (argon_ratio * cpi_argon)
        + (krypton_ratio * cpi_krypton)
        + (xenon_ratio * cpi_xenon)
    )
    Cp_mix = cp_mix / molecular_mass_mix

    # Viscosity
    mu_mix = (
        (air_ratio * mu_air)
        + (argon_ratio * mu_argon)
        + (krypton_ratio * mu_krypton)
        + (xenon_ratio * mu_xenon)
    )

    # Thermal conductivity
    lambda_mix = (
        (air_ratio * lambda_air)
        + (argon_ratio * lambda_argon)
        + (krypton_ratio * lambda_krypton)
        + (xenon_ratio * lambda_xenon)
    )

    return lambda_mix, mu_mix, Cp_mix, rho_mix


def gas_properties(T_gas, gas=0, P=101325):
    # Calculates the thermal properties of different gases, based on Annex B of ISO15099
    # ThP 2015
    # Attention: T_gas in K
    # 0=Air, 1=Argon, 2=Krypton, 3=Xenon

    lambda_a = [2.873e-3, 2.285e-3, 9.443e-4, 4.538e-4]
    lambda_b = [7.760e-5, 5.149e-5, 2.826e-5, 1.723e-5]
    mu_a = [3.723e-6, 3.379e-6, 2.213e-6, 1.069e-6]
    mu_b = [4.94e-8, 6.451e-8, 7.777e-8, 7.414e-8]
    Cp_a = [1002.7370, 521.9285, 248.0907, 158.3397]
    Cp_b = [1.2324e-2, 0, 0, 0]
    molecular_masses = [
        28.97,
        39.948,
        83.80,
        131.30,
    ]  # [kg/kgmol] Molecular masses of the different gases

    lambda_gas = lambda_a[gas] + lambda_b[gas] * \
        T_gas  # [W/(m*K)] Thermal conductivity
    mu_gas = mu_a[gas] + mu_b[gas] * T_gas  # [Pa*s] Dynamic viscosity
    Cp_gas = (
        Cp_a[gas] + Cp_b[gas] * T_gas
    )  # [J/(kg*K)]  Specific heat capacity at constant pressure
    rho_gas = (
        P * molecular_masses[gas] / (R * T_gas * 1000)
    )  # [kg/m³] Density of the gas
    M_gas = molecular_masses[gas]
    #print ('gas properties')
    return lambda_gas, mu_gas, Cp_gas, rho_gas, M_gas


def h_cv_closed_cavity(T1, T2, height_glazed_area, d_gv, argon_ratio=0):
    # Calculates the convective heat transfer coefficient of vertical closed cavities
    # From ISO15099 p. 17 and 24
    # ThP 2015
    # d_gv [m] thickness of the gas layer
    # gas type: 0=Air, 1=Argon, 2=Krypton, 3=Xenon
    # aminextras
    #d_gv =0.07

    #height_glazed_area = 1.5
    #
    Tm = (T1 + T2) / 2  # Mean surface temperatures [K]

    air_ratio = 1 - argon_ratio
    # We get the gas properties:
    lambda_gas, mu_gas, Cp_gas, rho_gas = gas_mixture_properties(
        Tm,
        air_ratio,
        argon_ratio,
        krypton_ratio=0,
        xenon_ratio=0,
    )

    Ra = ((rho_gas ** 2) * ((d_gv) ** 3) * g * Cp_gas * abs(T1 - T2)) / \
        (mu_gas * lambda_gas *
         Tm)  # [-] Rayleigh number of the air in the cavity

    A_gv = (height_glazed_area) / d_gv  # [-] Aspect ratio of the cavity

    Nu_1 = 0

    if 0 < Ra <= 1e4:
        Nu_1 = 1 + ((1.7596678e-10) * (Ra ** 2.2984755))
    elif 1e4 < Ra <= 5e4:
        Nu_1 = 0.028154 * (Ra ** 0.4134)
    elif Ra > 5e4:
        Nu_1 = 0.0673838 * (Ra ** (1 / 3))

    Nu_2 = 0.242 * ((Ra / A_gv) ** 0.272)

    Nu = max(Nu_1, Nu_2)

    h_cv = Nu * (
        lambda_gas / d_gv

    )  # [W/(m²*K)] Convective heat transfer coefficient from one surface to the other

    return h_cv


def h_cv_vent_cavity(T1, T2, height_glazed_area, d_gv, V):
    # h_cv in [W/(m²*K)]P.41 V is the mean air velocity
    # This function calculates the convective heat transfer coefficient in ventilated gaps
    # It is the HTC between one surface and the average air temperature in the ventilated cavity

    h_cv_vent = 2 * \
        h_cv_closed_cavity(T1, T2, height_glazed_area, d_gv) + (4 * V)

    return h_cv_vent


def dP_driving_i_k(T_gap_i, T_gap_k, height_glazed_area):
    # This function calculates the driving pressure difference for natural convection through holes between two cavities
    T_0 = 0.5 * (T_gap_i + T_gap_k)  # [K] Reference temperature
    rho_0 = gas_mixture_properties(
        T_0,
        air_ratio=1,
        argon_ratio=0,
        krypton_ratio=0,
        xenon_ratio=0,)[
        -1
    ]  # [kg/m³] Density at the reference temperature

    dP_driving = (
        rho_0
        * T_0
        * g
        * height_glazed_area
        * abs(T_gap_i - T_gap_k)
        / (T_gap_i * T_gap_k)
    )
    # print("dP_driving="+str(dP_driving))
    return dP_driving


def ventilated_case(T1, T2, T_gap, T_gap_in, A_ho, height_glazed_area, width_glazed_area, d_gv, V):
    # ventilated_case(T_f_3_0, T_b_2_0, T_gap_1_2,
    # The following calculations of values are needed in the calculation of q_i and q_i_plus_1 P.46
    # The outputs are used in the calcultaions of q_i, q_i_plus_1 and T_gap[i]

    T_av_1_2 = 0.5 * (T2 + T1)
    # the average temp of the surfaces of layers i and i+1 (tb and tf,i+1)
    A_s_i = d_gv * (width_glazed_area)

    lambda_gas, mu_gas, Cp_gas, rho_gas = gas_mixture_properties(
        T_gap,
        air_ratio=1,
        argon_ratio=0,
        krypton_ratio=0,
        xenon_ratio=0,
    )

    if V < 0:
        A_equ_inl_i = (A_ho)
        A_equ_out_i = (A_ho)
        Z_inl_i = ((A_s_i / (0.6 * A_equ_inl_i)) - 1) ** 2
        Z_out_i = ((A_s_i / (0.6 * A_equ_out_i)) - 1) ** 2
        a11 = (rho_gas / 2) * (1 + Z_inl_i + Z_out_i)
        b11 = (12 * mu_gas * height_glazed_area) / (d_gv ** 2)
        c11 = -dP_driving_i_k(T_av_1_2, T_gap_in, height_glazed_area)
        delta = (b11**2) - (4*a11*c11)
        V1 = (-b11 - math.sqrt(delta))/(2*a11)
        V2 = (-b11 + math.sqrt(delta))/(2*a11)
        if V1 > V2:
            V = V1
        else:
            V = V2
    airflow = V*((d_gv * (width_glazed_area)))
    h_cv_1_2 = h_cv_vent_cavity(T1, T2, height_glazed_area, d_gv, V)
    # [m] The characteristic height of the temperatur profile
    H_0_1_2 = (rho_gas * Cp_gas * d_gv * V) / (2 * h_cv_1_2)
    alfa = 1 - exp(-height_glazed_area / H_0_1_2)
    beta = exp(-height_glazed_area / H_0_1_2)
    # [K] The gap outlet temperature
    T_gap_out = T_av_1_2 - ((T_av_1_2 - T_gap_in) *
                            exp((-height_glazed_area) / H_0_1_2))
    T_gap_middle = T_av_1_2 - ((T_av_1_2 - T_gap_in)
                               * exp((-height_glazed_area) / (2 * H_0_1_2)))
    T_gap = T_av_1_2 - ((H_0_1_2/height_glazed_area) * (T_gap_out - T_gap_in))
    q_v_g = (rho_gas * Cp_gas * airflow) / \
        (width_glazed_area * height_glazed_area)
    return h_cv_1_2, T_gap_out, T_gap, q_v_g  # pylint: disable=unbalanced-tuple-unpacking


def ventilated_case_stvb(T1, T2, T_gap, T_gap_in, A_ho, height_glazed_area, width_glazed_area, d_gv, V):
    # ventilated_case(T_f_3_0, T_b_2_0, T_gap_1_2,
    # The following calculations of values are needed in the calculation of q_i and q_i_plus_1 P.46
    # The outputs are used in the calcultaions of q_i, q_i_plus_1 and T_gap[i]

    T_av_1_2 = 0.5 * (T2 + T1)
    # the average temp of the surfaces of layers i and i+1 (tb and tf,i+1)
    A_s_i = d_gv * (width_glazed_area)

    lambda_gas, mu_gas, Cp_gas, rho_gas = gas_mixture_properties(
        T_gap,
        air_ratio=1,
        argon_ratio=0,
        krypton_ratio=0,
        xenon_ratio=0,
    )

    if V < 0:
        A_equ_inl_i = (A_ho)
        A_equ_out_i = (A_ho)
        Z_inl_i = ((A_s_i / (0.6 * A_equ_inl_i)) - 1) ** 2
        Z_out_i = ((A_s_i / (0.6 * A_equ_out_i)) - 1) ** 2
        a11 = (rho_gas / 2) * (1 + Z_inl_i + Z_out_i)
        b11 = (12 * mu_gas * height_glazed_area) / (d_gv ** 2)
        c11 = -dP_driving_i_k(T_av_1_2, T_gap_in, height_glazed_area)
        delta = (b11**2) - (4*a11*c11)
        solution1 = (-b11 - math.sqrt(delta))/(2*a11)
        solution2 = (-b11 + math.sqrt(delta))/(2*a11)
        V1 = solution1
        V2 = solution2
        if V1 > V2:
            V = V1
        else:
            V = V2
    airflow = V*((d_gv * (width_glazed_area)))
    h_cv_1_2 = h_cv_vent_cavity(T1, T2, height_glazed_area, d_gv, V)
    # [m] The characteristic height of the temperatur profile
    H_0_1_2 = (rho_gas * Cp_gas * d_gv * V) / (2 * h_cv_1_2)
    alfa = 1 - exp(-height_glazed_area / H_0_1_2)
    beta = exp(-height_glazed_area / H_0_1_2)
    # [K] The gap outlet temperature
    T_gap_out = T_av_1_2 - ((T_av_1_2 - T_gap_in) *
                            exp((-height_glazed_area) / H_0_1_2))
    T_gap_middle = T_av_1_2 - ((T_av_1_2 - T_gap_in)
                               * exp((-height_glazed_area) / (2 * H_0_1_2)))
    T_gap = T_av_1_2 - ((H_0_1_2/height_glazed_area) * (T_gap_out - T_gap_in))
    q_v_g = (rho_gas * Cp_gas * airflow) / \
        (width_glazed_area * height_glazed_area)
    return h_cv_1_2, T_gap_out, T_gap_middle, T_gap, q_v_g, T_av_1_2, H_0_1_2, alfa, beta  # pylint: disable=unbalanced-tuple-unpacking


def ventilated_case_esso(T1, T2, T_gap, T_gap_in, A_ho, height_glazed_area, width_glazed_area, d_gv, V, d_top, d_bot, d_left, d_right, d_su):
    # ventilated_case(T_f_3_0, T_b_2_0, T_gap_1_2,
    # The following calculations of values are needed in the calculation of q_i and q_i_plus_1 P.46
    # The outputs are used in the calcultaions of q_i, q_i_plus_1 and T_gap[i]

    T_av_1_2 = 0.5 * (T2 + T1)
    # the average temp of the surfaces of layers i and i+1 (tb and tf,i+1)
    A_s_i = d_gv * (width_glazed_area)

    lambda_gas, mu_gas, Cp_gas, rho_gas = gas_mixture_properties(
        T_gap,
        air_ratio=1,
        argon_ratio=0,
        krypton_ratio=0,
        xenon_ratio=0,
    )
    C_1 = 0.078
    C_2 = 1.20
    A_ho = C_1 * (d_su ** C_2) * (height_glazed_area * width_glazed_area)
    if d_top > 0 or d_bot > 0:
        A_tp = d_top * width_glazed_area
        A_bo = d_bot * width_glazed_area
        A_lf = d_left * height_glazed_area
        A_rt = d_right * height_glazed_area
        A_equ_inl_i = A_tp + (A_bo * (A_lf + A_rt + A_ho)
                              ) / (2 * (A_bo + A_tp))
        A_equ_out_i = A_bo + (A_tp * (A_lf + A_rt + A_ho)
                              ) / (2 * (A_bo + A_tp))
    # when there is no gap at the sides and theres only surface opennes
    if d_top == 0 and d_bot == 0 and d_left == 0 and d_right == 0:
        A_equ_inl_i = 0.25 * A_ho
        A_equ_out_i = 0.25 * A_ho
    # when there is no gap at the sides nor surface opennes and only bottom and top openning as big as the cavity
    if d_top == 0 and d_bot == 0 and d_left == 0 and d_right == 0 and d_su == 0:
        A_equ_inl_i = d_gv * width_glazed_area
        A_equ_out_i = d_gv * width_glazed_area
    if V < 0:
        Z_inl_i = ((A_s_i / (0.6 * A_equ_inl_i)) - 1) ** 2
        Z_out_i = ((A_s_i / (0.6 * A_equ_out_i)) - 1) ** 2
        a11 = (rho_gas / 2) * (1 + Z_inl_i + Z_out_i)
        b11 = (12 * mu_gas * height_glazed_area) / (d_gv ** 2)
        c11 = -dP_driving_i_k(T_av_1_2, T_gap_in, height_glazed_area)
        delta = (b11**2) - (4*a11*c11)
        solution1 = (-b11 - math.sqrt(delta))/(2*a11)
        solution2 = (-b11 + math.sqrt(delta))/(2*a11)
        V1 = solution1
        V2 = solution2
        if V1 > V2:
            V = V1
        else:
            V = V2

    airflow = V*((d_gv * (width_glazed_area)))
    h_cv_1_2 = h_cv_vent_cavity(T1, T2, height_glazed_area, d_gv, V)
    # [m] The characteristic height of the temperatur profile
    H_0_1_2 = (rho_gas * Cp_gas * d_gv * V) / (2 * h_cv_1_2)
    alfa = 1 - exp(-height_glazed_area / H_0_1_2)
    beta = exp(-height_glazed_area / H_0_1_2)
    # [K] The gap outlet temperature
    T_gap_out = T_av_1_2 - ((T_av_1_2 - T_gap_in) *
                            exp((-height_glazed_area) / H_0_1_2))
    T_gap_middle = T_av_1_2 - ((T_av_1_2 - T_gap_in)
                               * exp((-height_glazed_area) / (2 * H_0_1_2)))
    T_gap = T_av_1_2 - ((H_0_1_2/height_glazed_area) * (T_gap_out - T_gap_in))
    q_v_g = (rho_gas * Cp_gas * airflow) / \
        (width_glazed_area * height_glazed_area)
    return h_cv_1_2, T_gap_out, T_gap_middle, T_gap, q_v_g, T_av_1_2, H_0_1_2, alfa, beta  # pylint: disable=unbalanced-tuple-unpacking
