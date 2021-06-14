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
        help="model of solar thermal technology",)
    parser.add_argument(
        "-model 7",
        action="store",
        type=int,
        help="ESSO ID4",)
    parser.add_argument(
        "-model 8",
        action="store",
        type=int,
        help="ESSO ID38",)
    parser.add_argument(
        "-model 9",
        action="store",
        type=int,
        help="ESSO ID99",)
    parser.add_argument(
        "-model 10",
        action="store",
        type=int,
        help="EN673",)
    parser.add_argument(
        "-model 11",
        action="store",
        type=int,
        help="SYNERGIEFASSADEN TGU",)
    parser.add_argument(
        "-model 12",
        action="store",
        type=int,
        help="SYNERGIEFASSADEN DGU",)
    parser.add_argument(
        "-model 13",
        action="store",
        type=int,
        help="SYNERGIEFASSADEN Analytic",)
    parser.add_argument(
        "-model 14",
        action="store",
        type=int,
        help="SYNERGIEFASSADEN four_layer_model",)
    parser.add_argument(
        "-model 15",
        action="store",
        type=int,
        help="SYNERGIEFASSADEN tgu_hannover_sh",)


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

def _add_inputs(parser, config, opts):
	config.layers = parser.get("PATHS", "layers")
    
        
def _add_variables(parser, config, opts):
    """ Reads parameters from config file """
    config.temp_outdoor = parser.getfloat("VARIABLES", "temp_outdoor")
    config.temp_indoor = parser.getfloat("VARIABLES", "temp_indoor")
    config.irradiation = parser.getfloat("VARIABLES", "irradiation")
    config.conv_coeff_out = parser.getfloat("VARIABLES", "conv_coeff_out")
    config.conv_coeff_in = parser.getfloat("VARIABLES", "conv_coeff_in")
    if opts.model >= 13:
        config.cavity_depth = parser.getfloat("VARIABLES", "cavity_depth")
        config.volume_flowrate = parser.getfloat(
            "VARIABLES", "volume_flowrate")
        config.cavity_flowrate_factor = parser.getfloat(
            "VARIABLES", "cavity_flowrate_factor")
        config.textil_open_factor = parser.getfloat(
            "VARIABLES", "textil_open_factor")
        config.textil_thickness = parser.getfloat(
            "VARIABLES", "textil_thickness")
        config.temp_ini = parser.getfloat("VARIABLES", "temp_ini")




def _add_paths(parser, config):
    config.output_dir = parser.get("PATHS", "output_dir")



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



def print_synergiefassaden(textil_temp_in, cavity_temp,layers_temp_radiosities =[0,0]):
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
    np.savetxt(
        "%s/layers_temp_radiosities.out" % config.output_dir,
        np.transpose(layers_temp_radiosities),
        delimiter=",",
        fmt="%1.2f",
    )

def curves(opts, config, constant):
    
    layerVect = np.atleast_2d(np.genfromtxt(config.layers,skip_header=1,delimiter=','))
    
    if opts.model == 7:
        print('*** ESSO ID4 ***')
        for i in range(1):
            energy_into_room, u_value, g_value = esso.id4(
                config.temp_outdoor,
                config.temp_indoor,
                273.15,
                config.irradiation,
                7.7,
                25,
                layerVect
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)
    if opts.model == 8:
        print('*** ESSO ID38 ***')
        for i in range(1):
            energy_into_room, u_value, g_value = esso.id38(
                config.temp_outdoor,
                config.temp_indoor,
                273.15,
                config.irradiation,
                7.7,
                25,
                layerVect
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)
    if opts.model == 9:
        print('*** ESSO ID99 ***')
        for i in range(1):
            energy_into_room, u_value, g_value = esso.id99(
                config.temp_outdoor,
                config.temp_indoor,
                273.15,
                config.irradiation,
                7.7,
                25,
                layerVect
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
                config.conv_coeff_in,
                config.conv_coeff_out,
                layerVect,
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
                layerVect,
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
                layerVect,
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
                config.cavity_depth,
                config.volume_flowrate,
                config.cavity_flowrate_factor,
                config.textil_open_factor,
                config.textil_thickness,
                config.temp_ini,
                layerVect,
            )
        print_synergiefassaden(textil_temp_in, cavity_temp)
    if opts.model == 14:
        print('*** SYNERGIEFASSADEN four_layer_model ***')
        for i in range(1):
            textil_temp_in, cavity_temp, layers_temp_radiosities = synergiefassaden.four_layer_model(
                config.temp_outdoor,
                config.temp_indoor,
                config.temp_outdoor,
                config.irradiation,
                5.5,
                config.conv_coeff_in,
                config.conv_coeff_out,
                config.cavity_depth,
                config.volume_flowrate,
                config.cavity_flowrate_factor,
                config.textil_open_factor,
                config.textil_thickness,
                config.temp_ini,
                layerVect,
            )
        print_synergiefassaden(textil_temp_in, cavity_temp,layers_temp_radiosities)
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
                layerVect,
            )
            print('Heat exchange with int, U_value, G_value')
            print(energy_into_room, u_value, g_value)


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
