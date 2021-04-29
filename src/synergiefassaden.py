#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on April 06 2021

@author: Bruno Bueno
"""

from scipy.optimize import fsolve
from numpy import ndarray, array, empty, zeros, ones, pi
import solar_func
import constant


def heat_transfer_through_textil(
    temp_air,
    tube_diameter,
    textil_thickness,
    textil_open_factor,
    open_area,
    textil_volume_flowrate
):
    # air properties
    air_thermal_conductivity, \
        air_viscosity, \
        air_specific_heat, \
        air_density, \
        M_gas = solar_func.gas_properties(temp_air)

    # !! Assumption
    tube_air_velocity = textil_volume_flowrate/open_area
    print("tube_air_velocity", tube_air_velocity)

    reynolds = air_density*tube_air_velocity*tube_diameter/air_viscosity
    prandtl = air_specific_heat*air_viscosity/air_thermal_conductivity
    print("reynolds", reynolds)
    if reynolds < 2300:
        print("LAMINAR")
    else:
        print("TURBULENT")
        exit()

    # Average Nusselt for a tube (entrance effects)
    # Eq 4.50, Mills 1999, p. 304
    entrance_ratio = tube_diameter/textil_thickness
    nusselt_tube = 3.66 + 0.065*entrance_ratio*reynolds*prandtl / \
        (1+0.04*(entrance_ratio*reynolds*prandtl)**(2./3.))

    tube_conv_coeff = nusselt_tube*air_thermal_conductivity/tube_diameter
    print("tube_conv_coeff", tube_conv_coeff)

    return tube_conv_coeff, air_density, air_specific_heat


def analytic(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    wind_speed,
    conv_coeff_in,
    conv_coeff_out,
    facade_width,
    facade_height,
    cavity_depth,
    volume_flowrate,
    cavity_flowrate_factor,
    textil_thickness,
    textil_thermal_conductivity,
    textil_solar_absortivity,
    textil_emissivity_front,
    textil_open_factor,
    temp_ini
):
    """ model inputs """
    # ! Assumption
    tube_diameter = textil_thickness/2
    number_layer = 4
    number_vertical_partitions = 5
    textil_emissivity_back = textil_emissivity_front
    # textil_air_permeability = textil_open_factor

    """ Dimensions and areas """
    facade_area = facade_width*facade_height  # m2
    partition_area = facade_area/number_vertical_partitions
    # ! Assumption
    open_area = textil_open_factor*partition_area
    hole_area = pi * (tube_diameter/2)**2
    number_holes = open_area / hole_area
    tube_area = pi * tube_diameter * textil_thickness * number_holes
    print('tube_area', tube_area)

    """ Flow rates and air velocities """
    textil_volume_flowrate = volume_flowrate * \
        (1-cavity_flowrate_factor)/number_vertical_partitions

    cavity_volume_flowrate = zeros(number_vertical_partitions)
    cavity_velocity = zeros(number_vertical_partitions)
    for i in range(number_vertical_partitions):
        cavity_volume_flowrate[i] = volume_flowrate * \
            cavity_flowrate_factor + textil_volume_flowrate*i
        cavity_velocity[i] = cavity_volume_flowrate[i] / \
            facade_width/cavity_depth

    """ Heat transfer through the textile """
    tube_conv_coeff, air_density, air_specific_heat = heat_transfer_through_textil(
        temp_ini,
        tube_diameter,
        textil_thickness,
        textil_open_factor,
        open_area,
        textil_volume_flowrate
    )

    """ Layers properties """
    layer = ndarray((number_layer, ), dtype=object)
    """ First layer """
    layer[0] = solar_func.Layer(
        facade_width,
        facade_height,
        0.89,
        0.08,
        0.11,
        0.92,
        0,
        0,
        0.75,
        0.12,
        0.13,
        0.004,
        1,
        0.014,
    )
    """ Second layer """
    layer[1] = solar_func.Layer(
        facade_width,
        facade_height,
        0.89,
        0.89,
        0.11,
        0.11,
        0,
        0,
        0.90,
        0.08,
        0.02,
        0.004,
        1,
        0.014,
    )
    """ Third layer """
    layer[2] = solar_func.Layer(
        facade_width,
        facade_height,
        0.08,
        0.89,
        0.92,
        0.11,
        0,
        0,
        0.75,
        0.14,
        0.11,
        0.004,
        1,
        0.210,
    )
    """ Textil layer """
    layer[3] = solar_func.Layer(
        facade_width,
        facade_height,
        textil_emissivity_front,
        textil_emissivity_back,
        1-textil_emissivity_front-textil_open_factor,
        1-textil_emissivity_back-textil_open_factor,
        textil_open_factor,
        textil_open_factor,
        1-textil_solar_absortivity,
        0,
        textil_solar_absortivity,
        textil_thickness,
        textil_thermal_conductivity,
        0,
    )

    """ Absorbed Irradiance in each layer """
    abs_irrad = zeros(number_layer)
    abs_irrad[0] = irradiation * layer[0].solar_abs
    abs_irrad[1] = irradiation * layer[0].solar_tra * layer[1].solar_abs
    abs_irrad[2] = irradiation * layer[0].solar_tra * layer[1].solar_tra\
        * layer[2].solar_abs
    abs_irrad[3] = irradiation * layer[0].solar_tra * layer[1].solar_tra\
        * layer[2].solar_tra * layer[3].solar_abs
    """ Cavity heat transfer coefficients and initial values """
    temp_b_1_0 = temp_ini + 2
    temp_f_2_0 = temp_b_1_0 + 2
    hc_1 = solar_func.h_cv_closed_cavity(temp_b_1_0, temp_f_2_0,
                                         layer[0].height, layer[0].gap, 0.9)
    temp_b_2_0 = temp_f_2_0 + 2
    temp_f_3_0 = temp_b_2_0 + 2
    hc_2 = solar_func.h_cv_closed_cavity(
        temp_b_2_0, temp_f_3_0, layer[1].height, layer[1].gap, 0.9)

    temp_b_3_0 = temp_f_3_0 + 2
    temp_f_4_0 = temp_b_3_0 + 2
    vent_cavity_conv_coeff = zeros(number_vertical_partitions)
    for i in range(number_vertical_partitions):
        vent_cavity_conv_coeff[i] = solar_func.h_cv_vent_cavity(
            temp_b_3_0, temp_f_4_0, facade_height, cavity_depth, cavity_velocity[i])

    def iso15099(z):
        temp_f_1 = z[0]
        temp_b_1 = z[1]
        j_f_1 = z[2]
        j_b_1 = z[3]
        temp_f_2 = z[4]
        temp_b_2 = z[5]
        j_f_2 = z[6]
        j_b_2 = z[7]
        temp_f_3 = z[8]
        temp_b_3 = z[9]
        j_f_3 = z[10]
        j_b_3 = z[11]
        temp_c = zeros(number_vertical_partitions)
        for i in range(number_vertical_partitions):
            temp_c[i] = z[12+i]
        index = 12+number_vertical_partitions
        temp_t = zeros(number_vertical_partitions)
        for i in range(number_vertical_partitions):
            temp_t[i] = z[index+i]
        index = index + number_vertical_partitions
        temp_f_4 = zeros(number_vertical_partitions)
        temp_b_4 = zeros(number_vertical_partitions)
        j_f_4 = zeros(number_vertical_partitions)
        j_b_4 = zeros(number_vertical_partitions)
        for i in range(number_vertical_partitions):
            temp_f_4[i] = z[index+i*4]
            temp_b_4[i] = z[index+i*4+1]
            j_f_4[i] = z[index+i*4+2]
            j_b_4[i] = z[index+i*4+3]
        F = empty(index+i*4+3+1)

        """ First glass layer equations """
        F[0] = conv_coeff_out * (temp_f_1 - temp_outdoor) \
            + layer[0].conductivity/layer[0].thickness * (temp_f_1 - temp_b_1) \
            - abs_irrad[0]/2
        F[1] = j_f_1 - layer[0].emi_inf_f * constant.SIGMA * \
            temp_f_1 ** 4 - layer[0].tra_inf_b * j_f_2
        F[2] = j_b_1 - layer[0].emi_inf_b * constant.SIGMA * \
            temp_b_1 ** 4 - layer[0].ref_inf_b * j_f_2
        F[3] = hc_1 * (temp_b_1 - temp_f_2) \
            + layer[0].conductivity/layer[0].thickness * (temp_b_1 - temp_f_1) \
            - abs_irrad[0]/2 \
            + j_b_1 - j_f_2
        """ Second glass layer equations """
        F[4] = hc_1 * (temp_f_2 - temp_b_1) \
            + layer[1].conductivity/layer[1].thickness * (temp_f_2 - temp_b_2) \
            - abs_irrad[1]/2 \
            + j_f_2 - j_b_1
        F[5] = j_f_2 - layer[1].emi_inf_f * constant.SIGMA * \
            temp_f_2 ** 4 \
            - layer[1].tra_inf_b * j_f_3 \
            - layer[1].ref_inf_f * j_b_1
        F[6] = j_b_2 - layer[1].emi_inf_b * constant.SIGMA * \
            temp_b_2 ** 4 \
            - layer[1].ref_inf_b * j_f_3 \
            - layer[1].tra_inf_f * j_b_1
        F[7] = hc_2 * (temp_b_2 - temp_f_3) \
            + layer[1].conductivity/layer[1].thickness * (temp_b_2 - temp_f_2) \
            - abs_irrad[1]/2 \
            + j_b_2 - j_f_3
        """ Third glass layer equations """
        F[8] = hc_2 * (temp_f_3 - temp_b_2) \
            + (temp_f_3 - temp_b_3)*layer[2].conductivity/layer[2].thickness \
            + j_f_3 - j_b_2 - abs_irrad[2]/2
        F[9] = j_f_3 - layer[2].emi_inf_f * constant.SIGMA * temp_f_3 ** 4 \
            - layer[2].ref_inf_f * j_b_2
        for i in range(number_vertical_partitions):
            F[9] = F[9] - layer[2].tra_inf_b * \
                j_f_4[i]/number_vertical_partitions
        F[10] = j_b_3 - layer[2].emi_inf_b * constant.SIGMA * temp_b_3 ** 4 \
            - layer[2].tra_inf_f * j_b_2
        for i in range(number_vertical_partitions):
            F[10] = F[10] - layer[2].ref_inf_b * \
                j_f_4[i]/number_vertical_partitions
        F[11] = (temp_b_3 - temp_f_3)*layer[2].conductivity/layer[2].thickness \
            - abs_irrad[2]/2 \
            + j_b_3
        for i in range(number_vertical_partitions):
            F[11] = F[11] - vent_cavity_conv_coeff[i]/number_vertical_partitions * (temp_c[i] - temp_b_3) \
                - j_f_4[i]/number_vertical_partitions

        """ Cavity equation """
        F[12] = vent_cavity_conv_coeff[0] * partition_area * (temp_c[0] - temp_b_3) \
            + vent_cavity_conv_coeff[0] * partition_area * (temp_c[0] - temp_f_4[0]) \
            + air_density * air_specific_heat * cavity_volume_flowrate[0] * (temp_c[0] - temp_indoor) \
            + air_density * air_specific_heat * cavity_volume_flowrate[0] * (temp_c[0] - temp_c[1]) \
            + air_density * air_specific_heat * \
            textil_volume_flowrate * (temp_c[0] - temp_t[0])
        for i in range(1, number_vertical_partitions-1):
            F[12+i] = vent_cavity_conv_coeff[i] * partition_area * (temp_c[i] - temp_b_3) \
                + vent_cavity_conv_coeff[i] * partition_area * (temp_c[i] - temp_f_4[i]) \
                + air_density * air_specific_heat * cavity_volume_flowrate[i] * (temp_c[i] - temp_c[i-1]) \
                + air_density * air_specific_heat * cavity_volume_flowrate[i] * (temp_c[i] - temp_c[i+1]) \
                + air_density * air_specific_heat * \
                textil_volume_flowrate * (temp_c[i] - temp_t[i])
        index = number_vertical_partitions-1
        # ! Check where the heat goes
        F[12+index] = vent_cavity_conv_coeff[index] * partition_area * (temp_c[index] - temp_b_3) \
            + vent_cavity_conv_coeff[i] * partition_area * (temp_c[index] - temp_f_4[index]) \
            + air_density * air_specific_heat * cavity_volume_flowrate[index] * (temp_c[index] - temp_c[index-1]) \
            + air_density * air_specific_heat * \
            textil_volume_flowrate * (temp_c[index] - temp_t[index])

        """ Tube equation """
        index = 12 + number_vertical_partitions
        for i in range(number_vertical_partitions):
            F[index+i] = air_density * air_specific_heat * textil_volume_flowrate * (temp_t[i] - temp_c[i]) \
                + tube_conv_coeff * tube_area * (temp_t[i] - temp_b_4[i])/2 \
                + tube_conv_coeff * tube_area * (temp_t[i] - temp_f_4[i])/2 \
                + air_density * air_specific_heat * \
                textil_volume_flowrate * (temp_t[i] - temp_indoor)

        """ textil equation """
        index = index + number_vertical_partitions
        for i in range(number_vertical_partitions):
            F[index+i*4] = vent_cavity_conv_coeff[i] * partition_area * (temp_f_4[i] - temp_c[i]) \
                + layer[3].conductivity/layer[3].thickness * partition_area * (temp_f_4[i] - temp_b_4[i]) \
                + tube_conv_coeff * tube_area * (temp_f_4[i] - temp_t[i])/2 \
                - abs_irrad[3] * partition_area \
                + j_f_4[i] * facade_area/number_vertical_partitions \
                - j_b_3 * facade_area/number_vertical_partitions
            F[index+i*4+1] = j_f_4[i] - layer[3].emi_inf_f * constant.SIGMA * temp_f_4[i] ** 4 \
                - layer[3].ref_inf_f * 1 * j_b_3
            F[index+i*4+2] = j_b_4[i] - layer[3].emi_inf_b * constant.SIGMA * temp_b_4[i] ** 4 \
                - layer[3].tra_inf_f * 1 * j_b_3
            F[index+i*4+3] = conv_coeff_in * partition_area * (temp_b_4[i] - temp_indoor) \
                + layer[3].conductivity/layer[3].thickness * partition_area * (temp_b_4[i] - temp_f_4[i]) \
                + tube_conv_coeff * tube_area * (temp_b_4[i] - temp_t[i])/2

        return F

    z_init_guess = ones(11+number_vertical_partitions*6+1)*temp_ini
    z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
    results = z
    """ Updating the results with corrected first guess """
    while abs(results[1] - temp_b_1_0) > 0.001 \
            or abs(results[4] - temp_f_2_0) > 0.001 \
            or abs(results[5] - temp_b_2_0) > 0.001 \
            or abs(results[8] - temp_f_3_0) > 0.001 \
            or abs(results[9] - temp_b_3_0) > 0.001:
        temp_b_1_0 = results[1]
        temp_f_2_0 = results[4]
        hc_1 = solar_func.h_cv_closed_cavity(temp_b_1_0, temp_f_2_0,
                                             layer[0].height, layer[0].gap, 0.9)
        temp_b_2_0 = results[5]
        temp_f_3_0 = results[8]
        hc_2 = solar_func.h_cv_closed_cavity(temp_b_2_0, temp_f_3_0,
                                             layer[1].height, layer[1].gap, 0.9)

        temp_b_3_0 = results[9]
        temp_f_4_0 = (results[18]+results[22]+results[26])/3.
        temp_f_4_0 = 0.
        for i in range(number_vertical_partitions):
            temp_f_4_0 = temp_f_4_0 + \
                results[12+number_vertical_partitions*2] / \
                number_vertical_partitions
        print('textil_temp', temp_f_4_0 - 273.15)

        vent_cavity_conv_coeff = zeros(number_vertical_partitions)
        for i in range(number_vertical_partitions):
            vent_cavity_conv_coeff[i] = solar_func.h_cv_vent_cavity(
                temp_b_3_0, temp_f_4_0, facade_height, cavity_depth, cavity_velocity[i])

        tube_conv_coeff, air_density, air_specific_heat = heat_transfer_through_textil(
            temp_f_4_0,
            tube_diameter,
            textil_thickness,
            textil_open_factor,
            open_area,
            textil_volume_flowrate
        )

        z_init_guess = results
        z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
        results = z

    temp_f_1 = z[0]
    temp_b_1 = z[1]
    j_f_1 = z[2]
    j_b_1 = z[3]
    temp_f_2 = z[4]
    temp_b_2 = z[5]
    j_f_2 = z[6]
    j_b_2 = z[7]
    temp_f_3 = z[8]
    temp_b_3 = z[9]
    j_f_3 = z[10]
    j_b_3 = z[11]
    temp_c = zeros(number_vertical_partitions)
    for i in range(number_vertical_partitions):
        temp_c[i] = z[12+i]
    index = 12+number_vertical_partitions
    temp_t = zeros(number_vertical_partitions)
    for i in range(number_vertical_partitions):
        temp_t[i] = z[index+i]
    index = index + number_vertical_partitions
    temp_f_4 = zeros(number_vertical_partitions)
    temp_b_4 = zeros(number_vertical_partitions)
    j_f_4 = zeros(number_vertical_partitions)
    j_b_4 = zeros(number_vertical_partitions)
    for i in range(number_vertical_partitions):
        temp_f_4[i] = z[index+i*4]
        temp_b_4[i] = z[index+i*4+1]
        j_f_4[i] = z[index+i*4+2]
        j_b_4[i] = z[index+i*4+3]
    print('glass_temp_out')
    print(temp_f_1-273.15)
    print('glass_temp_gap')
    print(temp_b_3-273.15)
    print('textil_temp_gap')
    print(temp_f_4-273.15)
    print('tube_temp')
    print(temp_t-273.15)
    print('textil_temp_in')
    print(temp_b_4-273.15)
    print('cavity_temp')
    print(temp_c-273.15)
    return (temp_b_4, temp_c)


def four_layer_model(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    wind_speed,
    conv_coeff_in,
    conv_coeff_out,
    facade_width,
    facade_height,
    cavity_depth,
    volume_flowrate,
    cavity_flowrate_factor,
    textil_thickness,
    textil_thermal_conductivity,
    textil_solar_absortivity,
    textil_emissivity_front,
    textil_open_factor,
    temp_ini
):
    """ tgu_hannover_sh - 3 panes with closed cavity according to email from M. Hiller
    Glazing type = Tripe glazing,
    Indoor ventilated roller blind
    Takes the boundary conditions and returns the heat exchange with the interior,
    U_value and G_value """
    """ Layers properties """
    textil_emissivity_back = textil_emissivity_front
    number_layer = 4
    layer = ndarray((number_layer, ), dtype=object)
    """ First layer """
    layer[0] = solar_func.Layer(
        facade_width,
        facade_height,
        0.86,
        0.01,
        0.14,
        0.99,
        0,
        0,
        0.275,
        1-0.275-0.4889,
        0.4889,
        0.016,
        1,
        0.016,
    )
    """ Second layer """
    layer[1] = solar_func.Layer(
        facade_width,
        facade_height,
        0.86,
        0.86,
        0.14,
        0.14,
        0,
        0,
        0.842,
        1-0.842-0.0136,
        0.0136,
        0.004,
        1,
        0.016,
    )
    """ Third layer """
    layer[2] = solar_func.Layer(
        facade_width,
        facade_height,
        0.03,
        0.86,
        0.97,
        0.14,
        0,
        0,
        0.586,
        1-0.586-0.0376,
        0.0376,
        0.008,
        1,
        0.210,
    )
    """ Textil layer """
    layer[3] = solar_func.Layer(
        facade_width,
        facade_height,
        textil_emissivity_front,
        textil_emissivity_back,
        1-textil_emissivity_front-textil_open_factor,
        1-textil_emissivity_back-textil_open_factor,
        textil_open_factor,
        textil_open_factor,
        1-0.14-textil_solar_absortivity,
        0.14,
        textil_solar_absortivity,
        textil_thickness,
        textil_thermal_conductivity,
        0,
    )
    """ Absorbed Irradiance in each layer !ADAPT for multiple reflections"""
    abs_irrad = zeros(number_layer)
    abs_irrad[0] = irradiation * layer[0].solar_abs
    abs_irrad[1] = irradiation * layer[1].solar_abs
    abs_irrad[2] = irradiation * layer[2].solar_abs
    abs_irrad[3] = irradiation * layer[3].solar_abs
    """ Cavity heat transfer coefficients and initial values """
    temp_b_1_0 = temp_ini + 20
    temp_f_2_0 = temp_b_1_0 + 2
    hc_1 = solar_func.h_cv_closed_cavity(temp_b_1_0, temp_f_2_0,
                                         layer[0].height, layer[0].gap, 0.9)
    temp_b_2_0 = temp_f_2_0 + 2
    temp_f_3_0 = temp_b_2_0 + 2
    hc_2 = solar_func.h_cv_closed_cavity(
        temp_b_2_0, temp_f_3_0, layer[1].height, layer[1].gap, 0.9)

    temp_b_3_0 = temp_f_3_0 + 2
    temp_f_4_0 = temp_b_3_0 + 2
    # hc_3 = solar_func.h_cv_closed_cavity(
    #     temp_b_3_0, temp_f_4_0, facade_height, cavity_depth)
    # cavity_velocity = volume_flowrate / facade_width / cavity_depth
    # hc_3 = solar_func.h_cv_vent_cavity(
    #     temp_b_3_0, temp_f_4_0, facade_height, cavity_depth, cavity_velocity)

    # """ Cavity heat transfer coefficients and initial values """
    temp_gap_3_0 = temp_b_3_0/2. + temp_f_4_0/2.
    temp_gap_in_3 = temp_indoor
    # d_top = 0.01
    # d_bot = 0.01
    # d_left = 0.01
    # d_right = 0.01
    d_su = 0.04
    a_ho_3 = layer[3].width * layer[3].height * d_su  # area * surface openness
    # (
    #     hc_3,
    #     temp_gap_out_3,
    #     temp_gap_middle_3,
    #     temp_gap_3,
    #     q_v_g_3,
    #     temp_av_3,
    #     h_0_3,
    #     alfa_3,
    #     beta_3,
    # ) = solar_func.ventilated_case_esso(
    #     temp_b_3_0,
    #     temp_f_4_0,
    #     temp_gap_3_0,
    #     temp_gap_in_3,
    #     a_ho_3,
    #     layer[2].height,
    #     layer[2].width,
    #     layer[2].gap,
    #     -7777,
    #     d_top,
    #     d_bot,
    #     d_left,
    #     d_right,
    #     d_su,
    # )
    (
        hc_3,
        temp_gap_out_3,
        temp_gap_3,
        q_v_g_3,
    ) = solar_func.ventilated_case(
        temp_b_3_0,
        temp_f_4_0,
        temp_gap_3_0,
        temp_gap_in_3,
        a_ho_3,
        layer[2].height,
        layer[2].width,
        layer[2].gap,
        -7777,
    )
    hc_4 = conv_coeff_in
    temp_f_5 = temp_indoor
    j_f_5 = constant.SIGMA * temp_f_5 ** 4
    hc_0 = conv_coeff_out
    temp_b_0 = temp_sky
    j_b_0 = constant.SIGMA * temp_b_0 ** 4

    def iso15099(z):
        temp_f_1 = z[0]
        temp_b_1 = z[1]
        j_f_1 = z[2]
        j_b_1 = z[3]
        temp_f_2 = z[4]
        temp_b_2 = z[5]
        j_f_2 = z[6]
        j_b_2 = z[7]
        temp_f_3 = z[8]
        temp_b_3 = z[9]
        j_f_3 = z[10]
        j_b_3 = z[11]
        temp_f_4 = z[12]
        temp_b_4 = z[13]
        j_f_4 = z[14]
        j_b_4 = z[15]
        F = empty(16)
        """ First glass layer equations """
        F[0] = hc_0 * (temp_f_1 - temp_b_0) \
            + layer[0].conductivity/layer[0].thickness * (temp_f_1 - temp_b_1) \
            - abs_irrad[0]/2 \
            + j_f_1 - j_b_0
        F[1] = j_f_1 - layer[0].emi_inf_f * constant.SIGMA * \
            temp_f_1 ** 4 \
            - layer[0].tra_inf_b * j_f_2 \
            - layer[0].ref_inf_f * j_b_0
        F[2] = j_b_1 - layer[0].emi_inf_b * constant.SIGMA * \
            temp_b_1 ** 4 \
            - layer[0].ref_inf_b * j_f_2 \
            - layer[0].tra_inf_f * j_b_0
        F[3] = hc_1 * (temp_b_1 - temp_f_2) \
            + layer[0].conductivity/layer[0].thickness * (temp_b_1 - temp_f_1) \
            - abs_irrad[0]/2 \
            + j_b_1 - j_f_2
        """ Second glass layer equations """
        F[4] = hc_1 * (temp_f_2 - temp_b_1) \
            + layer[1].conductivity/layer[1].thickness * (temp_f_2 - temp_b_2) \
            - abs_irrad[1]/2 \
            + j_f_2 - j_b_1
        F[5] = j_f_2 - layer[1].emi_inf_f * constant.SIGMA * \
            temp_f_2 ** 4 \
            - layer[1].tra_inf_b * j_f_3 \
            - layer[1].ref_inf_f * j_b_1
        F[6] = j_b_2 - layer[1].emi_inf_b * constant.SIGMA * \
            temp_b_2 ** 4 \
            - layer[1].ref_inf_b * j_f_3 \
            - layer[1].tra_inf_f * j_b_1
        F[7] = hc_2 * (temp_b_2 - temp_f_3) \
            + layer[1].conductivity/layer[1].thickness * (temp_b_2 - temp_f_2) \
            - abs_irrad[1]/2 \
            + j_b_2 - j_f_3
        """ Third glass layer equations """
        F[8] = hc_2 * (temp_f_3 - temp_b_2) \
            + layer[2].conductivity/layer[2].thickness * (temp_f_3 - temp_b_3) \
            - abs_irrad[2]/2 \
            + j_f_3 - j_b_2
        F[9] = j_f_3 - layer[2].emi_inf_f * constant.SIGMA * \
            temp_f_3 ** 4 \
            - layer[2].tra_inf_b * j_f_4 \
            - layer[2].ref_inf_f * j_b_2
        F[10] = j_b_3 - layer[2].emi_inf_b * constant.SIGMA * \
            temp_b_3 ** 4 \
            - layer[2].ref_inf_b * j_f_4 \
            - layer[2].tra_inf_f * j_b_2
        F[11] = hc_3 * (temp_b_3 - temp_f_4) \
            + layer[2].conductivity/layer[2].thickness * (temp_b_3 - temp_f_3) \
            - abs_irrad[2]/2 \
            + j_b_3 - j_f_4
        """ Four layer equations """
        F[12] = hc_3 * (temp_f_4 - temp_b_3) \
            + layer[3].conductivity/layer[3].thickness*(temp_f_4 - temp_b_4) \
            - abs_irrad[3]/2 \
            + j_f_4 - j_b_3
        F[13] = j_f_4 - layer[3].emi_inf_f * constant.SIGMA * temp_f_4**4 \
            - layer[3].tra_inf_b * j_f_5 \
            - layer[3].ref_inf_f * j_b_3
        F[14] = j_b_4 - layer[3].emi_inf_b * constant.SIGMA * temp_b_4**4 \
            - layer[3].ref_inf_b * j_f_5 \
            - layer[3].tra_inf_f * j_b_3
        F[15] = hc_4 * (temp_b_4 - temp_f_5) \
            + layer[3].conductivity/layer[3].thickness*(temp_b_4 - temp_f_4) \
            - abs_irrad[3]/2 \
            + j_b_4 - j_f_5
        return F
    z_init_guess = ones(16)*temp_ini
    z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
    results = z
    """ Updating the results with corrected first guess """
    while abs(results[1] - temp_b_1_0) > 0.001 and abs(results[4]
                                                       - temp_f_2_0) > 0.001:
        temp_b_1_0 = results[1]
        temp_f_2_0 = results[4]
        hc_1 = solar_func.h_cv_closed_cavity(temp_b_1_0, temp_f_2_0,
                                             layer[0].height, layer[0].gap, 0.9)
        temp_b_2_0 = results[5]
        temp_f_3_0 = results[8]
        hc_2 = solar_func.h_cv_closed_cavity(temp_b_2_0, temp_f_3_0,
                                             layer[1].height, layer[1].gap, 0.9)
        temp_b_3_0 = results[9]
        temp_f_4_0 = results[12]
        # hc_3 = solar_func.h_cv_closed_cavity(
        #     temp_b_3_0, temp_f_4_0, facade_height, cavity_depth)
        # hc_3 = solar_func.h_cv_vent_cavity(
        #     temp_b_3_0, temp_f_4_0, facade_height, cavity_depth, cavity_velocity)
        temp_gap_3_0 = temp_b_3_0/2. + temp_f_4_0/2.
        temp_gap_in_3 = temp_indoor
        # (
        #     hc_3,
        #     temp_gap_out_3,
        #     temp_gap_middle_3,
        #     temp_gap_3,
        #     q_v_g_3,
        #     temp_av_3,
        #     h_0_3,
        #     alfa_3,
        #     beta_3,
        # ) = solar_func.ventilated_case_esso(
        #     temp_b_3_0,
        #     temp_f_4_0,
        #     temp_gap_3_0,
        #     temp_gap_in_3,
        #     a_ho_3,
        #     layer[2].height,
        #     layer[2].width,
        #     layer[2].gap,
        #     -7777,
        #     d_top,
        #     d_bot,
        #     d_left,
        #     d_right,
        #     d_su,
        # )
        (
            hc_3,
            temp_gap_out_3,
            temp_gap_3,
            q_v_g_3,
        ) = solar_func.ventilated_case(
            temp_b_3_0,
            temp_f_4_0,
            temp_gap_3_0,
            temp_gap_in_3,
            a_ho_3,
            layer[2].height,
            layer[2].width,
            layer[2].gap,
            -7777,
        )
        print('hc_3', hc_3)
        z_init_guess = results
        z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
        results = z
    temp_f_1 = z[0]
    temp_b_1 = z[1]
    temp_f_2 = z[4]
    temp_b_2 = z[5]
    temp_f_3 = z[8]
    temp_b_3 = z[9]
    temp_f_4 = z[12]
    temp_b_4 = z[13]
    print('temp_f_1')
    print(temp_f_1-273.15)
    print('temp_b_1')
    print(temp_b_1-273.15)
    print('temp_f_2')
    print(temp_f_2-273.15)
    print('temp_b_2')
    print(temp_b_2-273.15)
    print('temp_f_3')
    print(temp_f_3-273.15)
    print('temp_b_3')
    print(temp_b_3-273.15)
    print('temp_f_4')
    print(temp_f_4-273.15)
    print('temp_b_4')
    print(temp_b_4-273.15)
    temp_t = [temp_b_4]
    temp_c = [temp_f_4/2.+temp_b_3/2.]
    return (temp_t, temp_c)


def tgu(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    wind_speed,
    conv_coeff_in,
    conv_coeff_out,
):
    """ tgu_hannover_sh - 3 panes with closed cavity according to email from M. Hiller
    Glazing type = Tripe glazing,
    Takes the boundary conditions and returns the heat exchange with the interior,
    U_value and G_value """
    """ Layers properties """
    number_layer = 3
    layer = ndarray((number_layer, ), dtype=object)
    """ First layer """
    layer[0] = solar_func.Layer(
        1.5,
        1,
        0.89,
        0.08,
        0.11,
        0.92,
        0,
        0,
        0.75,
        0.12,
        0.13,
        0.004,
        1,
        0.014,
    )
    """ Second layer """
    layer[1] = solar_func.Layer(
        1.5,
        1,
        0.89,
        0.89,
        0.11,
        0.11,
        0,
        0,
        0.90,
        0.08,
        0.02,
        0.004,
        1,
        0.014,
    )
    """ Third layer """
    layer[2] = solar_func.Layer(
        1.5,
        1,
        0.08,
        0.89,
        0.92,
        0.11,
        0,
        0,
        0.75,
        0.14,
        0.11,
        0.004,
        1,
        0,
    )
    """ Absorbed Irradiance in each layer """
    abs_irrad = zeros(number_layer)
    abs_irrad[0] = irradiation * layer[0].solar_abs
    abs_irrad[1] = irradiation * layer[0].solar_tra * layer[1].solar_abs
    abs_irrad[2] = irradiation * layer[0].solar_tra * layer[1].solar_tra\
        * layer[2].solar_abs
    """ Cavity heat transfer coefficients and initial values """
    temp_b_1_0 = temp_outdoor + 2
    temp_f_2_0 = temp_b_1_0 + 2
    hc_1 = solar_func.h_cv_closed_cavity(temp_b_1_0, temp_f_2_0,
                                         layer[0].height, layer[0].gap, 0.9)
    temp_b_2_0 = temp_f_2_0 + 2
    temp_f_3_0 = temp_b_2_0 + 3
    hc_2 = solar_func.h_cv_closed_cavity(
        temp_b_2_0, temp_f_3_0, layer[1].height, layer[1].gap, 0.9)

    def iso15099(z):
        temp_f_1 = z[0]
        temp_b_1 = z[1]
        j_f_1 = z[2]
        j_b_1 = z[3]
        temp_f_2 = z[4]
        temp_b_2 = z[5]
        j_f_2 = z[6]
        j_b_2 = z[7]
        temp_f_3 = z[8]
        temp_b_3 = z[9]
        j_f_3 = z[10]
        j_b_3 = z[11]
        F = empty(12)
        """ First layer equations """
        F[0] = conv_coeff_out * (temp_f_1 - temp_outdoor) - abs_irrad[0] \
            - hc_1 * (temp_f_2 - temp_b_1) - j_f_2 + j_b_1
        F[1] = j_f_1 - (layer[0].emi_inf_f * constant.SIGMA * temp_f_1
                        ** 4 + layer[0].tra_inf_f * j_f_2)
        F[2] = j_b_1 - (layer[0].emi_inf_b * constant.SIGMA * temp_b_1
                        ** 4 + layer[0].ref_inf_b * j_f_2)
        F[3] = layer[0].conductivity/layer[0].thickness * (temp_b_1 - temp_f_1) \
            - hc_1 * (temp_f_2 - temp_b_1) + j_f_2 - j_b_1
        """ Second layer equations """
        F[4] = hc_1 * (temp_f_2 - temp_b_1) + j_f_2 - j_b_1 \
            + hc_2 * (temp_b_2 - temp_f_3) + j_b_2 - j_f_3 \
            - abs_irrad[1]
        F[5] = j_f_2 - (layer[1].emi_inf_f * constant.SIGMA * temp_f_2 ** 4
                        + layer[1].ref_inf_f * j_b_1
                        + layer[1].tra_inf_b * j_f_3)
        F[6] = j_b_2 - (layer[1].emi_inf_b * constant.SIGMA * temp_b_2 ** 4
                        + layer[1].tra_inf_f * j_b_1
                        + layer[1].ref_inf_b * j_f_3)
        F[7] = temp_b_2 - temp_f_2 - layer[1].thickness / (2
                                                           * layer[1].conductivity) * (2 * (conv_coeff_in
                                                                                            * (temp_indoor - temp_b_2)))

        """ Third layer equations """
        F[8] = hc_2 * (temp_f_3 - temp_b_2) + j_f_3 - j_b_2  \
            - (abs_irrad[2] + conv_coeff_in * (temp_indoor - temp_b_3))
        F[9] = j_f_3 - (layer[2].emi_inf_f * constant.SIGMA * temp_f_3
                        ** 4 + layer[2].ref_inf_f * j_b_2)
        F[10] = j_b_3 - (layer[2].emi_inf_b * constant.SIGMA * temp_b_3
                         ** 4 + layer[2].tra_inf_f * j_b_2)
        F[11] = layer[2].conductivity/layer[2].thickness*(temp_b_3 - temp_f_3) \
            - conv_coeff_in * (temp_indoor - temp_b_3)
        return F

    z_init_guess = array([
        temp_outdoor,
        temp_outdoor,
        300,
        300,
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
        temp_b_2_0 = results[5]
        temp_f_3_0 = results[8]
        hc_2 = solar_func.h_cv_closed_cavity(temp_b_2_0, temp_f_3_0,
                                             layer[1].height, layer[1].gap, 0.9)
        z_init_guess = results
        z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
        results = z

    energy_into_room = conv_coeff_in * (z[9] - temp_indoor)
    if temp_indoor - temp_outdoor == 0:
        print(' U_value needs temperature difference greater than zero')
        u_value = -7777
    else:
        u_value = -energy_into_room / \
            (temp_indoor - temp_outdoor)  # [-] U_value
    if irradiation == 0:
        print(' G_value needs irradiation greater than zero')
        g_value = -7777
    else:
        g_value = (energy_into_room + irradiation * layer[0].solar_tra
                   * layer[1].solar_tra * layer[2].solar_tra) / irradiation  # [-] Gvalue
    return (energy_into_room, u_value, g_value)


def dgu(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    wind_speed,
    conv_coeff_in,
    conv_coeff_out,
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
        1.5,
        1,
        0.84,
        0.04,
        0.16,
        0.96,
        0,
        0,
        0.32,
        0.28,
        0.40,
        0.004,
        1,
        0.016,
    )
    """ Second layer """
    layer[1] = solar_func.Layer(
        1.5,
        1,
        0.84,
        0.84,
        0.16,
        0.16,
        0,
        0,
        0.83,
        0.08,
        0.09,
        0.004,
        1,
        0,
    )
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
        F[0] = conv_coeff_out * (temp_f_1 - temp_outdoor) - (abs_irrad[0] + hc_1
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
    while abs(results[1] - temp_b_1_0) > 0.001 or abs(results[4]
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


def tgu_hannover_sh(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    wind_speed,
    conv_coeff_in,
    conv_coeff_out,
):
    """ tgu_hannover_sh - 3 panes with closed cavity according to email from M. Hiller
    Glazing type = Tripe glazing,
    Takes the boundary conditions and returns the heat exchange with the interior,
    U_value and G_value """
    """ Layers properties """
    number_layer = 3
    layer = ndarray((number_layer, ), dtype=object)
    """ First layer """
    layer[0] = solar_func.Layer(
        1,
        3,
        0.86,
        0.01,
        0.14,
        0.99,
        0,
        0,
        0.275,
        0.279,
        0.4637,
        0.016,
        1,
        0.016,
    )
    """ Second layer """
    layer[1] = solar_func.Layer(
        1,
        3,
        0.86,
        0.86,
        0.14,
        0.14,
        0,
        0,
        0.842,
        0.076,
        0.0099,
        0.004,
        1,
        0.016,
    )
    """ Third layer """
    layer[2] = solar_func.Layer(
        1,
        3,
        0.03,
        0.86,
        0.97,
        0.14,
        0,
        0,
        0.586,
        0.210,
        0.0233,
        0.008,
        1,
        0,
    )
    """ Absorbed Irradiance in each layer """
    abs_irrad = zeros(number_layer)
    abs_irrad[0] = irradiation * layer[0].solar_abs
    abs_irrad[1] = irradiation * layer[1].solar_abs
    abs_irrad[2] = irradiation * layer[2].solar_abs
    """ Cavity heat transfer coefficients and initial values """
    temp_b_1_0 = temp_outdoor + 2
    temp_f_2_0 = temp_b_1_0 + 2
    hc_1 = solar_func.h_cv_closed_cavity(temp_b_1_0, temp_f_2_0,
                                         layer[0].height, layer[0].gap, 0.9)
    temp_b_2_0 = temp_f_2_0 + 2
    temp_f_3_0 = temp_b_2_0 + 3
    hc_2 = solar_func.h_cv_closed_cavity(
        temp_b_2_0, temp_f_3_0, layer[1].height, layer[1].gap, 0.9)

    hc_3 = conv_coeff_in
    temp_f_4 = temp_indoor
    j_f_4 = constant.SIGMA * temp_f_4 ** 4
    hc_0 = conv_coeff_out
    temp_b_0 = temp_sky
    j_b_0 = constant.SIGMA * temp_b_0 ** 4

    def iso15099(z):
        temp_f_1 = z[0]
        temp_b_1 = z[1]
        j_f_1 = z[2]
        j_b_1 = z[3]
        temp_f_2 = z[4]
        temp_b_2 = z[5]
        j_f_2 = z[6]
        j_b_2 = z[7]
        temp_f_3 = z[8]
        temp_b_3 = z[9]
        j_f_3 = z[10]
        j_b_3 = z[11]
        F = empty(12)
        """ First glass layer equations """
        F[0] = hc_0 * (temp_f_1 - temp_b_0) \
            + layer[0].conductivity/layer[0].thickness * (temp_f_1 - temp_b_1) \
            - abs_irrad[0]/2 \
            + j_f_1 - j_b_0
        F[1] = j_f_1 - layer[0].emi_inf_f * constant.SIGMA * \
            temp_f_1 ** 4 \
            - layer[0].tra_inf_b * j_f_2 \
            - layer[0].ref_inf_f * j_b_0
        F[2] = j_b_1 - layer[0].emi_inf_b * constant.SIGMA * \
            temp_b_1 ** 4 \
            - layer[0].ref_inf_b * j_f_2 \
            - layer[0].tra_inf_f * j_b_0
        F[3] = hc_1 * (temp_b_1 - temp_f_2) \
            + layer[0].conductivity/layer[0].thickness * (temp_b_1 - temp_f_1) \
            - abs_irrad[0]/2 \
            + j_b_1 - j_f_2
        """ Second glass layer equations """
        F[4] = hc_1 * (temp_f_2 - temp_b_1) \
            + layer[1].conductivity/layer[1].thickness * (temp_f_2 - temp_b_2) \
            - abs_irrad[1]/2 \
            + j_f_2 - j_b_1
        F[5] = j_f_2 - layer[1].emi_inf_f * constant.SIGMA * \
            temp_f_2 ** 4 \
            - layer[1].tra_inf_b * j_f_3 \
            - layer[1].ref_inf_f * j_b_1
        F[6] = j_b_2 - layer[1].emi_inf_b * constant.SIGMA * \
            temp_b_2 ** 4 \
            - layer[1].ref_inf_b * j_f_3 \
            - layer[1].tra_inf_f * j_b_1
        F[7] = hc_2 * (temp_b_2 - temp_f_3) \
            + layer[1].conductivity/layer[1].thickness * (temp_b_2 - temp_f_2) \
            - abs_irrad[1]/2 \
            + j_b_2 - j_f_3
        """ Third glass layer equations """
        F[8] = hc_2 * (temp_f_3 - temp_b_2) \
            + layer[2].conductivity/layer[2].thickness * (temp_f_3 - temp_b_3) \
            - abs_irrad[2]/2 \
            + j_f_3 - j_b_2
        F[9] = j_f_3 - layer[2].emi_inf_f * constant.SIGMA * \
            temp_f_3 ** 4 \
            - layer[2].tra_inf_b * j_f_4 \
            - layer[2].ref_inf_f * j_b_2
        F[10] = j_b_3 - layer[2].emi_inf_b * constant.SIGMA * \
            temp_b_3 ** 4 \
            - layer[2].ref_inf_b * j_f_4 \
            - layer[2].tra_inf_f * j_b_2
        F[11] = hc_3 * (temp_b_3 - temp_f_4) \
            + layer[2].conductivity/layer[2].thickness * (temp_b_3 - temp_f_3) \
            - abs_irrad[2]/2 \
            + j_b_3 - j_f_4

        return F

    z_init_guess = array([
        temp_outdoor,
        temp_outdoor,
        300,
        300,
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
        temp_b_2_0 = results[5]
        temp_f_3_0 = results[8]
        hc_2 = solar_func.h_cv_closed_cavity(temp_b_2_0, temp_f_3_0,
                                             layer[1].height, layer[1].gap, 0.9)
        z_init_guess = results
        z = fsolve(iso15099, z_init_guess, xtol=1.49012e-08, factor=0.1)
        results = z
    temp_f_1 = z[0]
    temp_b_1 = z[1]
    temp_f_2 = z[4]
    temp_b_2 = z[5]
    temp_f_3 = z[8]
    temp_b_3 = z[9]
    print('temp_f_1')
    print(temp_f_1-273.15)
    print('temp_b_1')
    print(temp_b_1-273.15)
    print('temp_f_2')
    print(temp_f_2-273.15)
    print('temp_b_2')
    print(temp_b_2-273.15)
    print('temp_f_3')
    print(temp_f_3-273.15)
    print('temp_b_3')
    print(temp_b_3-273.15)
    energy_into_room = conv_coeff_in * (z[9] - temp_indoor) + z[11] - j_f_4
    if temp_indoor - temp_outdoor == 0:
        print(' U_value needs temperature difference greater than zero')
        u_value = -7777
    else:
        u_value = -energy_into_room / \
            (temp_indoor - temp_outdoor)  # [-] U_value
    if irradiation == 0:
        print(' G_value needs irradiation greater than zero')
        g_value = -7777
    else:
        # g_value = (energy_into_room + irradiation * layer[0].solar_tra
        #            * layer[1].solar_tra * layer[2].solar_tra) / irradiation  # [-] Gvalue
        g_value = (energy_into_room + irradiation * layer[0].solar_tra
                   * layer[1].solar_tra * layer[2].solar_tra) / irradiation  # [-] Gvalue
    return (energy_into_room, u_value, g_value)


def en673(
    temp_outdoor,
    temp_indoor,
    temp_sky,
    irradiation,
    wind_speed,
    conv_coeff_in,
    conv_coeff_out,
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
        1,
        3,
        0.86,
        0.01,
        0.14,
        0.99,
        0,
        0,
        0.275,
        0.279,
        0.446,
        0.016,
        1,
        0.016,
    )
    """ Second layer """
    layer[1] = solar_func.Layer(
        1,
        3,
        0.86,
        0.86,
        0.14,
        0.14,
        0,
        0,
        0.842,
        0.076,
        0.082,
        0.004,
        1,
        0,
    )
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
        F[0] = conv_coeff_out * (temp_f_1 - temp_outdoor) - (abs_irrad[0] + hc_1
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
    while abs(results[1] - temp_b_1_0) > 0.001 or abs(results[4]
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
