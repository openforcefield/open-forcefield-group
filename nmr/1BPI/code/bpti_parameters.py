from simtk import unit as u
import numpy as np
import pandas as pd
import itertools

code = "1bpi"
padding = 0.95 * u.nanometers
cutoff = 0.9 * u.nanometers

temperature = 303 * u.kelvin
ionicStrength = 95 * u.millimolar

timestep = 2.0 * u.femtoseconds
equilibration_timestep = 1.0 * u.femtoseconds

barostat_frequency = 25
pressure = 1.0 * u.atmosphere
friction = 0.25 / u.picoseconds

n_steps = 30000000
n_equil_steps = 100000
output_frequency = 1000
output_frequency_allatoms = output_frequency * 50
equilibrate_output_frequency = n_equil_steps - 1

forcefields = ["amber99sbildn", "amber96", "amber99sbnmr"]
water_models = ["tip3p", "tip4pew", "tip3p-fb", "tip4p-fb"]

products = itertools.product(forcefields, water_models, sequences)


base_waters = {"tip3p":"tip3p", "tip4pew":"tip4pew", "tip3p-fb":"tip3p", "tip4p-fb":"tip4pew"}

