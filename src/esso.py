#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Sat Feb 21 12:20:21 2021

@author: Amin Kouti
"""

from scipy.optimize import fsolve
from numpy import ndarray, array, empty, zeros
import solar_func
import constant


def en673(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    conv_coeff_in,
    conv_coeff_out,
    layerVect,
):
    """ EN673 - 2 panes with closed cavity
    Glazing type = Double glazing (4 mm + 12 mm space + 4 mm) ,
    Takes the boundary conditions and returns the heat exchange with the interior,
    U_value and G_value """
    """ Layers properties """
    number_layer = 2
    layer = ndarray((number_layer, ), dtype=object)
    """ First layer """
    layer[0] = solar_func.Layer(
    layerVect[0, 0],
    layerVect[0, 1],
    layerVect[0, 2],
    layerVect[0, 3],
    layerVect[0, 4],
    layerVect[0, 5],
    layerVect[0, 6],
    layerVect[0, 7],
    layerVect[0, 8],
    layerVect[0, 9],
    layerVect[0, 10],
    layerVect[0, 11],
    layerVect[0, 12],
    layerVect[0, 13],)
    layer[1] = solar_func.Layer(
    layerVect[1, 0],
    layerVect[1, 1],
    layerVect[1, 2],
    layerVect[1, 3],
    layerVect[1, 4],
    layerVect[1, 5],
    layerVect[1, 6],
    layerVect[1, 7],
    layerVect[1, 8],
    layerVect[1, 9],
    layerVect[1, 10],
    layerVect[1, 11],
    layerVect[1, 12],
    layerVect[1, 13])
    """ Absorbed Irradiance in each layer """
    abs_irrad = zeros(number_layer)
    abs_irrad[0] = irradiation * layer[0].solar_abs
    abs_irrad[1] = irradiation * layer[0].solar_tra * layer[1].solar_abs
    """ Cavity heat transfer coefficients and initial values """
    temp_b_1_0 = temp_outdoor + 2
    temp_f_2_0 = temp_outdoor + 5
    hc_1 = solar_func.h_cv_closed_cavity(temp_b_1_0, temp_f_2_0,
                                         layer[0].height, layer[0].gap, 0.9)

    def iso15099(z):
        temp_f_1 = z[0]
        temp_b_1 = z[1]
        j_f_1 = z[2]
        j_b_1 = z[3]
        temp_f_2 = z[4]
        temp_b_2 = z[5]
        j_f_2 = z[6]
        j_b_2 = z[7]
        F = empty(8)
        """ First layer equations """
        F[0] = conv_coeff_out * (temp_f_1 - temp_outdoor) - (hc_1
                                                             * (temp_f_2 - temp_b_1) + j_f_2 - j_b_1)
        F[1] = j_f_1 - (layer[0].emi_inf_f * constant.SIGMA * temp_f_1
                        ** 4 + layer[0].tra_inf_f * j_f_2)
        F[2] = j_b_1 - (layer[0].emi_inf_b * constant.SIGMA * temp_b_1
                        ** 4 + layer[0].ref_inf_b * j_f_2)
        F[3] = temp_b_1 - temp_f_1 - layer[0].thickness / (2
                                                           * layer[0].conductivity) * (2 * (hc_1 * (temp_f_2
                                                                                                    - temp_b_1) + j_f_2 - j_b_1))
        """ Second layer equations """
        F[4] = hc_1 * (temp_f_2 - temp_b_1) + j_f_2 - j_b_1 \
            - (abs_irrad[1] + conv_coeff_in * (temp_indoor - temp_b_2))
        F[5] = j_f_2 - (layer[1].emi_inf_f * constant.SIGMA * temp_f_2
                        ** 4 + layer[1].ref_inf_f * j_b_1)
        F[6] = j_b_2 - (layer[1].emi_inf_b * constant.SIGMA * temp_b_2
                        ** 4 + layer[1].tra_inf_f * j_b_1)
        F[7] = temp_b_2 - temp_f_2 - layer[1].thickness / (2
                                                           * layer[1].conductivity) * (2 * (conv_coeff_in
                                                                                            * (temp_indoor - temp_b_2)))
        return F

    z_init_guess = array([
        temp_outdoor,
        temp_outdoor,
        300,
        300,
        temp_outdoor,
        temp_indoor,
        300,
        300,
    ])
    z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
    results = z
    """ Updating the results with corrected first guess """
    while abs(results[1] - temp_b_1_0) > 0.001 and abs(results[4]
                                                       - temp_f_2_0) > 0.001:
        temp_b_1_0 = results[1]
        temp_f_2_0 = results[4]
        hc_1 = solar_func.h_cv_closed_cavity(temp_b_1_0, temp_f_2_0,
                                             layer[0].height, layer[0].gap, 0.9)
        z_init_guess = results
        z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
        results = z

    energy_into_room = conv_coeff_in * (temp_indoor - z[5])
    if temp_indoor - temp_outdoor == 0:
        print(' U_value needs temperature difference greater than zero')
        u_value = -7777
    else:
        u_value = energy_into_room / \
            (temp_indoor - temp_outdoor)  # [-] U_value
    if irradiation == 0:
        print(' G_value needs irradiation greater than zero')
        g_value = -7777
    else:
        g_value = (energy_into_room + irradiation * layer[0].solar_tra
                   * layer[1].solar_tra) / irradiation  # [-] Gvalue
    return (energy_into_room, u_value, g_value)


def id4(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    conv_coeff_in,
    conv_coeff_out,
    layerVect,
):
    """ ID4 - internal roller blind inside reveal
    Glazing type = F (pane 1 / 16 mm space* / pane 2 (flipped) low e on position 3)
    Naturally ventilated gap
    Roller blind = Tempotest Star Screen 8024/400 
    Takes the boundary conditions and returns the heat exchange with the interior,
    U_value and G_value """
    """ Layers properties """
    number_layer = 2
    layer = ndarray((number_layer, ), dtype=object)
    """ First layer """
    layer[0] = solar_func.Layer(
    layerVect[0, 0],
    layerVect[0, 1],
    layerVect[0, 2],
    layerVect[0, 3],
    layerVect[0, 4],
    layerVect[0, 5],
    layerVect[0, 6],
    layerVect[0, 7],
    layerVect[0, 8],
    layerVect[0, 9],
    layerVect[0, 10],
    layerVect[0, 11],
    layerVect[0, 12],
    layerVect[0, 13],)
    layer[1] = solar_func.Layer(
    layerVect[1, 0],
    layerVect[1, 1],
    layerVect[1, 2],
    layerVect[1, 3],
    layerVect[1, 4],
    layerVect[1, 5],
    layerVect[1, 6],
    layerVect[1, 7],
    layerVect[1, 8],
    layerVect[1, 9],
    layerVect[1, 10],
    layerVect[1, 11],
    layerVect[1, 12],
    layerVect[1, 13])
    """ Absorbed Irradiance in each layer """
    abs_irrad = zeros(number_layer)
    abs_irrad[0] = irradiation * layer[0].solar_abs
    abs_irrad[1] = irradiation * layer[0].solar_tra * layer[1].solar_abs
    """ Cavity heat transfer coefficients and initial values """
    temp_b_1_0 = 300
    temp_f_2_0 = 295
    temp_gap_1_0 = 297
    temp_gap_in_1 = temp_indoor
    d_top = 0.01
    d_bot = 0.01
    d_left = 0.01
    d_right = 0.01
    d_su = 0.034
    a_ho_1 = layer[1].width * layer[1].height * d_su  # area * surface openness
    (
        hc_1,
        temp_gap_out_1,
        temp_gap_1,
        q_v_g_1,
    ) = solar_func.ventilated_case( 
        temp_b_1_0,
        temp_f_2_0,
        temp_gap_1_0,
        temp_gap_in_1,
        d_su,
        d_top,
        d_bot,
        d_left,
        d_right,
        layer[1].height,
        layer[1].width,
        layer[0].gap,
        -7777,
    )
    # [W/m²] Sky (ambient) long wave irradiation
    j_b_0 = constant.SIGMA * temp_sky ** 4
    # [W/m²] interior long wave irradiation
    j_f_3 = constant.SIGMA * temp_indoor ** 4

    def iso15099(z):
        temp_f_1 = z[0]
        temp_b_1 = z[1]
        j_f_1 = z[2]
        j_b_1 = z[3]
        temp_f_2 = z[4]
        temp_b_2 = z[5]
        j_f_2 = z[6]
        j_b_2 = z[7]
        F = empty(8)
        """ First layer equations """
        F[0] = conv_coeff_out * (temp_f_1 - temp_outdoor) + j_f_1 \
            - j_b_0 - (abs_irrad[0] + hc_1 * (temp_gap_1 - temp_b_1)
                       + j_f_2 - j_b_1 )
        F[1] = j_f_1 - (layer[0].emi_inf_f * constant.SIGMA * temp_f_1
                        ** 4 + layer[0].tra_inf_f * j_f_2
                        + layer[0].ref_inf_f * j_b_0)
        F[2] = j_b_1 - (layer[0].emi_inf_b * constant.SIGMA * temp_b_1
                        ** 4 + layer[0].tra_inf_f * j_b_0
                        + layer[0].ref_inf_b * j_f_2)
        F[3] = temp_b_1 - temp_f_1 - layer[0].thickness / (2
                        * layer[0].conductivity) * (2 * (hc_1 * (temp_gap_1
                        - temp_b_1) + j_f_2 - j_b_1 ) + abs_irrad[0])
        """ Second layer equations """
        F[4] = hc_1 * (temp_f_2 - temp_gap_1) + j_f_2 - j_b_1 \
            - (abs_irrad[1] + conv_coeff_in * (temp_indoor - temp_b_2)
               + j_f_3 - j_b_2)
        F[5] = j_f_2 - (layer[1].emi_inf_f * constant.SIGMA * temp_f_2
                        ** 4 + layer[1].tra_inf_f * j_f_3
                        + layer[1].ref_inf_f * j_b_1)
        F[6] = j_b_2 - (layer[1].emi_inf_b * constant.SIGMA * temp_b_1
                        ** 4 + layer[1].tra_inf_f * j_b_1
                        + layer[1].ref_inf_b * j_f_3)
        F[7] = temp_b_2 - temp_f_2 - layer[1].thickness / (2
                                                           * layer[1].conductivity) * (2 * (conv_coeff_in
                                                                                            * (temp_indoor - temp_b_2) + j_f_3 - j_b_2)
                                                                                       + abs_irrad[1])
        return F
    z_init_guess = array([
        300,
        300,
        300,
        300,
        300,
        300,
        300,
        300,
    ])
    z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
    results = z
    """ Updating the results with corrected first guess """
    while abs(results[1] - temp_b_1_0) > 0.1 and abs(results[4]
                                                     - temp_f_2_0) > 0.1:
        temp_b_1_0 = results[1]
        temp_f_2_0 = results[4]
        temp_gap_0_1 = (temp_b_1_0 + temp_f_2_0) * 0.5
        (
            hc_1,
            temp_gap_out_1,
            temp_gap_1,
            q_v_g_1,
        ) = solar_func.ventilated_case(
	        temp_b_1_0,
	        temp_f_2_0,
	        temp_gap_1_0,
	        temp_gap_in_1,
	        d_su,
	        d_top,
	        d_bot,
	        d_left,
	        d_right,
	        layer[1].height,
	        layer[1].width,
	        layer[0].gap,
	        -7777,
        )
        z_init_guess = results
        z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
        results = z
    energy_into_room = conv_coeff_in * (temp_indoor - z[5]) + j_f_3 \
        - z[7]   # [W/m²] Heat exchang with the interior
    if temp_indoor - temp_outdoor == 0:
        print(' U_value needs temperature difference greater than zero')
        u_value = -7777
    else:
        u_value = energy_into_room / \
            (temp_indoor - temp_outdoor)  # [-] U_value
    if irradiation == 0:
        print(' G_value needs irradiation greater than zero')
        g_value = -7777
    else:
        g_value = (energy_into_room + irradiation * layer[0].solar_tra
                   * layer[1].solar_tra) / irradiation  # [-] Gvalue
    return (energy_into_room, u_value, g_value)


def id38(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    conv_coeff_in,
    conv_coeff_out,
    layerVect,
):
    """ ID38 - external roller blind free hanging fitted outside reveal
    Roller blind = Tempotest Star Screen 8024/400
    Naturally ventilated gap
    Glazing type = F (pane 1 / 16 mm space* / pane 2 (flipped) low e on position 3) 
    Takes the boundary conditions and returns the heat exchange with the interior,
    U_value and G_value """
    """ Layers properties """
    number_layer = 2
    layer = ndarray((number_layer, ), dtype=object)
    """ First layer """
    layer[0] = solar_func.Layer(
    layerVect[0, 0],
    layerVect[0, 1],
    layerVect[0, 2],
    layerVect[0, 3],
    layerVect[0, 4],
    layerVect[0, 5],
    layerVect[0, 6],
    layerVect[0, 7],
    layerVect[0, 8],
    layerVect[0, 9],
    layerVect[0, 10],
    layerVect[0, 11],
    layerVect[0, 12],
    layerVect[0, 13],)
    layer[1] = solar_func.Layer(
    layerVect[1, 0],
    layerVect[1, 1],
    layerVect[1, 2],
    layerVect[1, 3],
    layerVect[1, 4],
    layerVect[1, 5],
    layerVect[1, 6],
    layerVect[1, 7],
    layerVect[1, 8],
    layerVect[1, 9],
    layerVect[1, 10],
    layerVect[1, 11],
    layerVect[1, 12],
    layerVect[1, 13])
    """ Absorbed Irradiance in each layer """
    abs_irrad = zeros(number_layer)
    abs_irrad[0] = irradiation * layer[0].solar_abs
    abs_irrad[1] = irradiation * layer[0].solar_tra * layer[1].solar_abs
    """ Cavity heat transfer coefficients and initial values """
    temp_b_1_0 = temp_outdoor + 2
    temp_f_2_0 = temp_outdoor + 5
    temp_gap_1_0 = temp_outdoor + 10
    temp_gap_in_1 = temp_outdoor
    d_top = 0.05
    d_bot = 0.05
    d_left = 0.05
    d_right = 0.05
    d_su = 0.034
    a_ho_1 = layer[0].width * layer[0].height * d_su  # area * surface openness
    (
        hc_1,
        temp_gap_out_1,
        temp_gap_1,
        q_v_g_1,
    ) = solar_func.ventilated_case(
        temp_b_1_0,
        temp_f_2_0,
        temp_gap_1_0,
        temp_gap_in_1,
        d_su,
        d_top,
        d_bot,
        d_left,
        d_right,
        layer[0].height,
        layer[0].width,
        layer[0].gap,
        -7777,
    )
    # [W/m²] Sky (ambient) long wave irradiation
    j_b_0 = constant.SIGMA * temp_sky ** 4
    # [W/m²] interior long wave irradiation
    j_f_3 = constant.SIGMA * temp_indoor ** 4

    def iso15099(z):
        temp_f_1 = z[0]
        temp_b_1 = z[1]
        j_f_1 = z[2]
        j_b_1 = z[3]
        temp_f_2 = z[4]
        temp_b_2 = z[5]
        j_f_2 = z[6]
        j_b_2 = z[7]
        F = empty(8)
        """ First layer equations """
        F[0] = conv_coeff_out * (temp_f_1 - temp_outdoor) + j_f_1 \
            - j_b_0 - (abs_irrad[0] + hc_1 * (temp_gap_1 - temp_b_1)
                       + j_f_2 - j_b_1 )
        F[1] = j_f_1 - (layer[0].emi_inf_f * constant.SIGMA * temp_f_1
                        ** 4 + layer[0].tra_inf_f * j_f_2
                        + layer[0].ref_inf_f * j_b_0)
        F[2] = j_b_1 - (layer[0].emi_inf_b * constant.SIGMA * temp_b_1
                        ** 4 + layer[0].tra_inf_f * j_b_0
                        + layer[0].ref_inf_b * j_f_2)
        F[3] = temp_b_1 - temp_f_1 - layer[0].thickness / (2
                        * layer[0].conductivity) * (2 * (hc_1 * (temp_gap_1
                        - temp_b_1) + j_f_2 - j_b_1) + abs_irrad[0])
        """ Second layer equations """
        F[4] = hc_1 * (temp_f_2 - temp_gap_1) + j_f_2 - j_b_1 \
            - (abs_irrad[1] + conv_coeff_in * (temp_indoor - temp_b_2)
               + j_f_3 - j_b_2)
        F[5] = j_f_2 - (layer[1].emi_inf_f * constant.SIGMA * temp_f_2
                        ** 4 + layer[1].tra_inf_f * j_f_3
                        + layer[1].ref_inf_f * j_b_1)
        F[6] = j_b_2 - (layer[1].emi_inf_b * constant.SIGMA * temp_b_2
                        ** 4 + layer[1].tra_inf_f * j_b_1
                        + layer[1].ref_inf_b * j_f_3)
        F[7] = temp_b_2 - temp_f_2 - layer[1].thickness / (2
                        * layer[1].conductivity) * (2 * (conv_coeff_in
                        * (temp_indoor - temp_b_2) + j_f_3 - j_b_2)
                        + abs_irrad[1])
        return F

    z_init_guess = array([
        temp_outdoor,
        temp_outdoor,
        300,
        300,
        temp_outdoor,
        temp_indoor,
        300,
        300,
    ])
    z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
    results = z
    """ Updating the results with corrected first guess """
    while abs(results[1] - temp_b_1_0) > 0.1 and abs(results[4]
                                                     - temp_f_2_0) > 0.1:
        temp_b_1_0 = results[1]
        temp_f_2_0 = results[4]
        temp_gap_0_1 = (temp_b_1_0 + temp_f_2_0) * 0.5
        (
            hc_1,
            temp_gap_out_1,
            temp_gap_1,
            q_v_g_1,
        ) = solar_func.ventilated_case(
	        temp_b_1_0,
	        temp_f_2_0,
	        temp_gap_1_0,
	        temp_gap_in_1,
	        d_su,
	        d_top,
	        d_bot,
	        d_left,
	        d_right,
	        layer[0].height,
	        layer[0].width,
	        layer[0].gap,
	        -7777,
        )
        z_init_guess = results
        z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
        results = z
    energy_into_room = conv_coeff_in * (temp_indoor - z[5]) + j_f_3 \
        - z[7]  # [W/m²] Heat exchang with the interior
    if temp_indoor - temp_outdoor == 0:
        print(' U_value needs temperature difference greater than zero')
        u_value = -7777
    else:
        u_value = energy_into_room / \
            (temp_indoor - temp_outdoor)  # [-] U_value
    if irradiation == 0:
        print(' G_value needs irradiation greater than zero')
        g_value = -7777
    else:
        g_value = (energy_into_room + irradiation * layer[0].solar_tra
                   * layer[1].solar_tra) / irradiation  # [-] Gvalue
    return (energy_into_room, u_value, g_value)


def id99(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    conv_coeff_in,
    conv_coeff_out,
    layerVect,
):
    """ ID99 - 
    Glazing type = (pane 1 / naturally ventilated gap from exterior / pane 2 
    Takes the boundary conditions and returns the heat exchange with the interior,
    U_value and G_value """
    """ Layers properties """
    number_layer = 2
    layer = ndarray((number_layer, ), dtype=object)
    """ First layer """
    layer[0] = solar_func.Layer(
    layerVect[0, 0],
    layerVect[0, 1],
    layerVect[0, 2],
    layerVect[0, 3],
    layerVect[0, 4],
    layerVect[0, 5],
    layerVect[0, 6],
    layerVect[0, 7],
    layerVect[0, 8],
    layerVect[0, 9],
    layerVect[0, 10],
    layerVect[0, 11],
    layerVect[0, 12],
    layerVect[0, 13],)
    layer[1] = solar_func.Layer(
    layerVect[1, 0],
    layerVect[1, 1],
    layerVect[1, 2],
    layerVect[1, 3],
    layerVect[1, 4],
    layerVect[1, 5],
    layerVect[1, 6],
    layerVect[1, 7],
    layerVect[1, 8],
    layerVect[1, 9],
    layerVect[1, 10],
    layerVect[1, 11],
    layerVect[1, 12],
    layerVect[1, 13])
    """ Absorbed Irradiance in each layer """
    abs_irrad = zeros(number_layer)
    abs_irrad[0] = irradiation * layer[0].solar_abs
    abs_irrad[1] = irradiation * layer[0].solar_tra * layer[1].solar_abs
    """ Cavity heat transfer coefficients and initial values """
    temp_b_1_0 = temp_outdoor + 2
    temp_f_2_0 = temp_outdoor + 5
    temp_gap_1_0 = temp_outdoor + 10
    temp_gap_in_1 = temp_outdoor
    d_top = 0.005
    d_bot = 0.005
    d_left = 0.005
    d_right = 0.005
    d_su = 0
    a_ho_1 = layer[0].width * layer[0].height * d_su  # area * surface openness
    (
        hc_1,
        temp_gap_out_1,
        temp_gap_1,
        q_v_g_1,
    ) = solar_func.ventilated_case(
        temp_b_1_0,
        temp_f_2_0,
        temp_gap_1_0,
        temp_gap_in_1,
        d_su,
        d_top,
        d_bot,
        d_left,
        d_right,
        layer[0].height,
        layer[0].width,
        layer[0].gap,
        -7777,
    )
    # [W/m²] Sky (ambient) long wave irradiation
    j_b_0 = constant.SIGMA * temp_sky ** 4
    # [W/m²] interior long wave irradiation
    j_f_3 = constant.SIGMA * temp_indoor ** 4

    def iso15099(z):
        temp_f_1 = z[0]
        temp_b_1 = z[1]
        j_f_1 = z[2]
        j_b_1 = z[3]
        temp_f_2 = z[4]
        temp_b_2 = z[5]
        j_f_2 = z[6]
        j_b_2 = z[7]
        F = empty(8)
        """ First layer equations """
        F[0] = conv_coeff_out * (temp_f_1 - temp_outdoor) + j_f_1 \
            - j_b_0 - (abs_irrad[0] + hc_1 * (temp_gap_1 - temp_b_1)
                       + j_f_2 - j_b_1 )
        F[1] = j_f_1 - (layer[0].emi_inf_f * constant.SIGMA * temp_f_1
                        ** 4 + layer[0].tra_inf_f * j_f_2
                        + layer[0].ref_inf_f * j_b_0)
        F[2] = j_b_1 - (layer[0].emi_inf_b * constant.SIGMA * temp_b_1
                        ** 4 + layer[0].tra_inf_f * j_b_0
                        + layer[0].ref_inf_b * j_f_2)
        F[3] = temp_b_1 - temp_f_1 - layer[0].thickness / (2
                        * layer[0].conductivity) * (2 * (hc_1 * (temp_gap_1
                        - temp_b_1) + j_f_2 - j_b_1) + abs_irrad[0])
        """ Second layer equations """
        F[4] = hc_1 * (temp_f_2 - temp_gap_1) + j_f_2 - j_b_1 \
            - (abs_irrad[1] + conv_coeff_in * (temp_indoor - temp_b_2)
               + j_f_3 - j_b_2)
        F[5] = j_f_2 - (layer[1].emi_inf_f * constant.SIGMA * temp_f_2
                        ** 4 + layer[1].tra_inf_f * j_f_3
                        + layer[1].ref_inf_f * j_b_1)
        F[6] = j_b_2 - (layer[1].emi_inf_b * constant.SIGMA * temp_b_2
                        ** 4 + layer[1].tra_inf_f * j_b_1
                        + layer[1].ref_inf_b * j_f_3)
        F[7] = temp_b_2 - temp_f_2 - layer[1].thickness / (2
                                                           * layer[1].conductivity) * (2 * (conv_coeff_in
                                                                                            * (temp_indoor - temp_b_2) + j_f_3 - j_b_2)
                                                                                       + abs_irrad[1])
        return F

    z_init_guess = array([
        temp_outdoor,
        temp_outdoor,
        300,
        300,
        temp_outdoor,
        temp_indoor,
        300,
        300,
    ])
    z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
    results = z
    """ Updating the results with corrected first guess """
    while abs(results[1] - temp_b_1_0) > 0.1 and abs(results[4]
                                                     - temp_f_2_0) > 0.1:
        temp_b_1_0 = results[1]
        temp_f_2_0 = results[4]
        temp_gap_0_1 = (temp_b_1_0 + temp_f_2_0) * 0.5
        (
            hc_1,
            temp_gap_out_1,
            temp_gap_1,
            q_v_g_1,
        ) = solar_func.ventilated_case(
	        temp_b_1_0,
	        temp_f_2_0,
	        temp_gap_1_0,
	        temp_gap_in_1,
	        d_su,
	        d_top,
	        d_bot,
	        d_left,
	        d_right,
	        layer[0].height,
	        layer[0].width,
	        layer[0].gap,
	        -7777,
        )
        z_init_guess = results
        z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
        results = z
    energy_into_room = conv_coeff_in * (temp_indoor - z[5]) + j_f_3 \
        - z[7]  # [W/m²] Heat exchang with the interior
    if temp_indoor - temp_outdoor == 0:
        print(' U_value needs temperature difference greater than zero')
        u_value = -7777
    else:
        u_value = energy_into_room / \
            (temp_indoor - temp_outdoor)  # [-] U_value
    if irradiation == 0:
        print(' G_value needs irradiation greater than zero')
        g_value = -7777
    else:
        g_value = (energy_into_room + irradiation * layer[0].solar_tra
                   * layer[1].solar_tra) / irradiation  # [-] Gvalue
    return (energy_into_room, u_value, g_value)
