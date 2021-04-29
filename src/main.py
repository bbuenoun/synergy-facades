"""Calculates performance curves of building integrated solar thermal collectors.

Description.

  Typical usage example:

  python3 ./src/main.py {directory_path}/config.scv -model {model_index}

"""

from subprocess import PIPE, run
import argparse
import sys
import os
import traceback
from configparser import ConfigParser
import numpy as np
import constant
import tabsolar
import stvb
import esso
import synergiefassaden


def parse_args():
    return _build_parser().parse_args()


def _build_parser():
    parser = argparse.ArgumentParser()
    _add_common_arguments(parser)
    return parser


def _add_common_arguments(parser):
    parser.add_argument("config", help="file of input parameters")
    parser.add_argument(
        "-model",
        action="store",
        type=int,
        default=1,
        help="model of solar thermal technology",
    )


class Config:
    pass


def parse_config(config_path, opts):
    parser = ConfigParser()
    parser.read(config_path)
    config = Config()
    _add_variables(parser, config, opts)
    _add_paths(parser, config)
    _add_inputs(parser, config, opts)
    return config


def _add_variables(parser, config, opts):
    """ Reads parameters from config file """
    if opts.model >= 13:
        config.temp_outdoor = parser.getfloat("VARIABLES", "temp_outdoor")
        config.temp_indoor = parser.getfloat("VARIABLES", "temp_indoor")
        config.irradiation = parser.getfloat("VARIABLES", "irradiation")
        config.conv_coeff_out = parser.getfloat("VARIABLES", "conv_coeff_out")
        config.conv_coeff_in = parser.getfloat("VARIABLES", "conv_coeff_in")
        config.facade_width = parser.getfloat("VARIABLES", "facade_width")
        config.facade_height = parser.getfloat("VARIABLES", "facade_height")
        config.cavity_depth = parser.getfloat("VARIABLES", "cavity_depth")
        config.volume_flowrate = parser.getfloat(
            "VARIABLES", "volume_flowrate")
        config.cavity_flowrate_factor = parser.getfloat(
            "VARIABLES", "cavity_flowrate_factor")
        config.textil_thickness = parser.getfloat(
            "VARIABLES", "textil_thickness")
        config.textil_thermal_conductivity = parser.getfloat(
            "VARIABLES", "textil_thermal_conductivity")
        config.textil_solar_absortivity = parser.getfloat(
            "VARIABLES", "textil_solar_absortivity")
        config.textil_emissivity_front = parser.getfloat(
            "VARIABLES", "textil_emissivity_front")
        config.textil_open_factor = parser.getfloat(
            "VARIABLES", "textil_open_factor")
        config.temp_ini = parser.getfloat("VARIABLES", "temp_ini")
    elif opts.model >= 7:
        config.temp_outdoor = parser.getfloat("VARIABLES", "temp_outdoor")
        config.temp_indoor = parser.getfloat("VARIABLES", "temp_indoor")
        config.irradiation = parser.getfloat("VARIABLES", "irradiation")
        config.conv_coeff_out = parser.getfloat("VARIABLES", "conv_coeff_out")
        config.conv_coeff_in = parser.getfloat("VARIABLES", "conv_coeff_in")
    else:
        config.temp_outdoor = parser.getfloat("VARIABLES", "temp_outdoor")
        config.temp_indoor = parser.getfloat("VARIABLES", "temp_indoor")
        config.irradiation = parser.getfloat("VARIABLES", "irradiation")
        config.fluidmass_flow = parser.getfloat("VARIABLES", "fluidmass_flow")
        config.wind_speed = parser.getfloat("VARIABLES", "wind_speed")


def _add_paths(parser, config):
    config.output_dir = parser.get("PATHS", "output_dir")


def _add_inputs(parser, config, opts):
    if opts.model < 7:
        config.temp_fluid_in = parser.get("PATHS", "temp_fluid_in")
        config.area = parser.get("PATHS", "area")


def _create_non_existing_directories(config):
    if hasattr(config, "output_dir") and not os.path.exists(config.output_dir):
        os.makedirs(config.output_dir)


def shell(cmd, output_flag=False):
    """ Runs shell commands """
    print(cmd)
    stdout = PIPE if output_flag else None
    _completed_process = run(
        cmd, stdout=stdout, stderr=PIPE, shell=True, check=True)
    if output_flag:
        return _completed_process.stdout.decode("utf-8")


class Parameter:
    """ Declares model parameters """

    def __init__(self, constant, config):
        # self.temp_sky = ((0.5+(0.5*(5.31*10**(-13))*(config.temp_outdoor)**6))/
        # (constant.SIGMA))**0.25 # [K]
        self.temp_sky = 0.0552*(config.temp_outdoor)**1.5
        self.conv_coeff_out = 5.7 + 3.8 * config.wind_speed  # [W/(m²*K)]
        self.conv_coeff_in = 7.7  # [W/(m²*K)]


def print_output(energy_yield_vect, efficiency_vect, temp_reduced_vect, front_loss_vect=0, back_loss_vect=0, edge_loss_vect=0, g_value_vect=0):
    """ Saves outputs as a text file """
    np.savetxt(
        "%s/energy_yield.out" % config.output_dir,
        np.transpose(energy_yield_vect),
        delimiter=",",
        fmt="%1.2f",
    )
    np.savetxt(
        "%s/efficiency.out" % config.output_dir,
        np.transpose(efficiency_vect),
        delimiter=",",
        fmt="%1.2f",
    )
    np.savetxt(
        "%s/temp_reduced.out" % config.output_dir,
        np.transpose(temp_reduced_vect),
        delimiter=",",
        fmt="%1.2f",
    )
    np.savetxt(
        "%s/g_value.out" % config.output_dir,
        np.transpose(g_value_vect),
        delimiter=",",
        fmt="%1.2f",
    )
    np.savetxt(
        "%s/front_loss.out" % config.output_dir,
        np.transpose(front_loss_vect),
        delimiter=",",
        fmt="%1.2f",
    )
    np.savetxt(
        "%s/back_loss.out" % config.output_dir,
        np.transpose(back_loss_vect),
        delimiter=",",
        fmt="%1.2f",
    )
    np.savetxt(
        "%s/edge_loss.out" % config.output_dir,
        np.transpose(edge_loss_vect),
        delimiter=",",
        fmt="%1.2f",
    )


def print_synergiefassaden(textil_temp_in, cavity_temp):
    np.savetxt(
        "%s/textil_temp_in.out" % config.output_dir,
        np.transpose(textil_temp_in),
        delimiter=",",
        fmt="%1.2f",
    )
    np.savetxt(
        "%s/cavity_temp.out" % config.output_dir,
        np.transpose(cavity_temp),
        delimiter=",",
        fmt="%1.2f",
    )


def curves(opts, config, constant):

    if opts.model < 7:
        temp_fluid_in = np.genfromtxt(
            config.temp_fluid_in,
            delimiter=",",
        )
        area = np.genfromtxt(
            config.area,
            delimiter=",",
        )
        curve_parameter = Parameter(constant, config)
        number_point = len(temp_fluid_in)
        energy_yield_vect = np.zeros(number_point)
        efficiency_vect = np.zeros(number_point)
        temp_reduced_vect = np.zeros(number_point)
        g_value_vect = np.zeros(number_point)
        front_loss_vect = np.zeros(number_point)
        back_loss_vect = np.zeros(number_point)
        edge_loss_vect = np.zeros(number_point)
    if opts.model == 1:
        print('*** TABSOLAR DESIGN ***')
        for i in range(number_point):
            energy_yield_tot, efficiency, temp_reduced, front_loss, back_loss, edge_loss = tabsolar.design(
                config.temp_outdoor,
                config.temp_indoor,
                curve_parameter.temp_sky,
                temp_fluid_in[i],
                config.fluidmass_flow,
                config.irradiation,
                config.wind_speed,
                curve_parameter.conv_coeff_in,
                curve_parameter.conv_coeff_out,
                area
            )
            print(energy_yield_tot, efficiency, temp_reduced)
            energy_yield_vect[i] = energy_yield_tot
            efficiency_vect[i] = efficiency
            temp_reduced_vect[i] = temp_reduced
            front_loss_vect[i] = front_loss
            back_loss_vect[i] = back_loss
            edge_loss_vect[i] = edge_loss
    if opts.model == 2:
        print('*** TABSOLAR PREMIUM ***')
        for i in range(number_point):
            energy_yield_tot, efficiency, temp_reduced, front_loss, back_loss, edge_loss = tabsolar.premium(
                config.temp_outdoor,
                config.temp_indoor,
                curve_parameter.temp_sky,
                temp_fluid_in[i],
                config.fluidmass_flow,
                config.irradiation,
                config.wind_speed,
                curve_parameter.conv_coeff_in,
                curve_parameter.conv_coeff_out,
                area
            )
            print(energy_yield_tot, efficiency, temp_reduced)
            energy_yield_vect[i] = energy_yield_tot
            efficiency_vect[i] = efficiency
            temp_reduced_vect[i] = temp_reduced
            front_loss_vect[i] = front_loss
            back_loss_vect[i] = back_loss
            edge_loss_vect[i] = edge_loss
    if opts.model == 3:
        print('*** TABSOLAR ECONOMY ***')
        for i in range(number_point):
            energy_yield_tot, efficiency, temp_reduced, front_loss, back_loss, edge_loss = tabsolar.economy(
                config.temp_outdoor,
                config.temp_indoor,
                curve_parameter.temp_sky,
                temp_fluid_in[i],
                config.fluidmass_flow,
                config.irradiation,
                config.wind_speed,
                curve_parameter.conv_coeff_in,
                curve_parameter.conv_coeff_out,
                area
            )
            print(energy_yield_tot, efficiency, temp_reduced)
            energy_yield_vect[i] = energy_yield_tot
            efficiency_vect[i] = efficiency
            temp_reduced_vect[i] = temp_reduced
            front_loss_vect[i] = front_loss
            back_loss_vect[i] = back_loss
            edge_loss_vect[i] = edge_loss
    if opts.model == 4:
        print('*** STVB 01A ***')
        for i in range(number_point):
            energy_yield, efficiency, temp_reduced, g_value = stvb.stvb_01a(
                config.temp_outdoor,
                config.temp_indoor,
                curve_parameter.temp_sky,
                temp_fluid_in[i],
                config.fluidmass_flow,
                config.irradiation,
                config.wind_speed,
                curve_parameter.conv_coeff_in,
                curve_parameter.conv_coeff_out
            )
            print(energy_yield, efficiency, temp_reduced, g_value)
            energy_yield_vect[i] = energy_yield
            efficiency_vect[i] = efficiency
            temp_reduced_vect[i] = temp_reduced
            g_value_vect[i] = g_value
    if opts.model == 5:
        print('*** STVB 01B ***')
        for i in range(number_point):
            energy_yield, efficiency, temp_reduced, g_value = stvb.stvb_01b(
                config.temp_outdoor,
                config.temp_indoor,
                curve_parameter.temp_sky,
                temp_fluid_in[i],
                config.fluidmass_flow,
                config.irradiation,
                config.wind_speed,
                curve_parameter.conv_coeff_in,
                curve_parameter.conv_coeff_out
            )
            print(energy_yield, efficiency, temp_reduced, g_value)
            energy_yield_vect[i] = energy_yield
            efficiency_vect[i] = efficiency
            temp_reduced_vect[i] = temp_reduced
            g_value_vect[i] = g_value
    if opts.model == 6:
        print('*** STVB 02 ***')
        for i in range(number_point):
            energy_yield, efficiency, temp_reduced, g_value = stvb.stvb_02(
                config.temp_outdoor,
                config.temp_indoor,
                curve_parameter.temp_sky,
                temp_fluid_in[i],
                config.fluidmass_flow,
                config.irradiation,
                config.wind_speed,
                curve_parameter.conv_coeff_in,
                curve_parameter.conv_coeff_out
            )
            print(energy_yield, efficiency, temp_reduced, g_value)
            energy_yield_vect[i] = energy_yield
            efficiency_vect[i] = efficiency
            temp_reduced_vect[i] = temp_reduced
            g_value_vect[i] = g_value
    if opts.model == 7:
        print('*** ESSO ID4 ***')
        for i in range(1):
            energy_into_room, u_value, g_value = esso.id4(
                config.temp_outdoor,
                config.temp_indoor,
                config.fluidmass_flow,
                config.irradiation,
                config.wind_speed,
                7.7,
                25,
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)
    if opts.model == 8:
        print('*** ESSO ID38 ***')
        for i in range(1):
            energy_into_room, u_value, g_value = esso.id38(
                config.temp_outdoor,
                config.temp_indoor,
                config.fluidmass_flow,
                config.irradiation,
                config.wind_speed,
                7.7,
                25,
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)
    if opts.model == 9:
        print('*** ESSO ID99 ***')
        for i in range(1):
            energy_into_room, u_value, g_value = esso.id99(
                config.temp_outdoor,
                config.temp_indoor,
                config.fluidmass_flow,
                config.irradiation,
                config.wind_speed,
                7.7,
                25,
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)
    if opts.model == 10:
        print('*** EN673 ***')
        for i in range(1):
            energy_into_room, u_value, g_value = esso.en673(
                273.15,
                288.15,
                273.15,
                config.irradiation,
                5.5,
                config.conv_coeff_in,
                config.conv_coeff_out,
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)
    if opts.model == 11:
        print('*** SYNERGIEFASSADEN TGU ***')
        for i in range(1):
            energy_into_room, u_value, g_value = synergiefassaden.tgu(
                config.temp_outdoor,
                config.temp_indoor,
                config.temp_outdoor,
                config.irradiation,
                5.5,
                config.conv_coeff_in,
                config.conv_coeff_out,
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)
    if opts.model == 12:
        print('*** SYNERGIEFASSADEN DGU ***')
        for i in range(1):
            energy_into_room, u_value, g_value = synergiefassaden.dgu(
                config.temp_outdoor,
                config.temp_indoor,
                config.temp_outdoor,
                config.irradiation,
                5.5,
                config.conv_coeff_in,
                config.conv_coeff_out,
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)
    if opts.model == 13:
        print('*** SYNERGIEFASSADEN Analytic ***')
        for i in range(1):
            textil_temp_in, cavity_temp = synergiefassaden.analytic(
                config.temp_outdoor,
                config.temp_indoor,
                config.temp_outdoor,
                config.irradiation,
                5.5,
                config.conv_coeff_in,
                config.conv_coeff_out,
                config.facade_width,
                config.facade_height,
                config.cavity_depth,
                config.volume_flowrate,
                config.cavity_flowrate_factor,
                config.textil_thickness,
                config.textil_thermal_conductivity,
                config.textil_solar_absortivity,
                config.textil_emissivity_front,
                config.textil_open_factor,
                config.temp_ini
            )
        print_synergiefassaden(textil_temp_in, cavity_temp)
    if opts.model == 14:
        print('*** SYNERGIEFASSADEN four_layer_model ***')
        for i in range(1):
            textil_temp_in, cavity_temp = synergiefassaden.four_layer_model(
                config.temp_outdoor,
                config.temp_indoor,
                config.temp_outdoor,
                config.irradiation,
                5.5,
                config.conv_coeff_in,
                config.conv_coeff_out,
                config.facade_width,
                config.facade_height,
                config.cavity_depth,
                config.volume_flowrate,
                config.cavity_flowrate_factor,
                config.textil_thickness,
                config.textil_thermal_conductivity,
                config.textil_solar_absortivity,
                config.textil_emissivity_front,
                config.textil_open_factor,
                config.temp_ini
            )
        print_synergiefassaden(textil_temp_in, cavity_temp)
    if opts.model == 15:
        print('*** SYNERGIEFASSADEN tgu_hannover_sh ***')
        for i in range(1):
            energy_into_room, u_value, g_value = synergiefassaden.tgu_hannover_sh(
                config.temp_outdoor,
                config.temp_indoor,
                config.temp_outdoor,
                config.irradiation,
                5.5,
                config.conv_coeff_in,
                config.conv_coeff_out,
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)
    if opts.model < 7:
        print_output(energy_yield_vect, efficiency_vect, temp_reduced_vect,
                     front_loss_vect, back_loss_vect, edge_loss_vect, g_value_vect)


if __name__ == "__main__":
    try:
        opts = parse_args()
        config = parse_config(opts.config, opts)
        _create_non_existing_directories(config)
        curves(opts, config, constant)
    except Exception as error:
        print(
            "".join(
                traceback.format_exception(
                    etype=type(error), value=error, tb=error.__traceback__
                )
            )
        )
        sys.exit("The following error occurred: %s" % error)
