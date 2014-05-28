from simtk import unit as u
import numpy as np
import pandas as pd
import itertools

padding = 0.95 * u.nanometers
cutoff = 0.9 * u.nanometers

temperature = 298 * u.kelvin

timestep = 2.0 * u.femtoseconds
equilibration_timestep = 1.0 * u.femtoseconds

barostat_frequency = 25
pressure = 1.0 * u.atmosphere
friction = 0.25 / u.picoseconds

n_steps = 30000000
n_equil_steps = 1000000
output_frequency = 1000
output_frequency_allatoms = output_frequency * 50
equilibrate_output_frequency = n_equil_steps - 1

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]

sequences = ["ACE_%sG_NH2" % (aa) for aa in amino_acids]
sequences.extend(["ACE_G%s_NH2" % (aa) for aa in amino_acids])

#forcefields = ["amber99sbildn", "amber96", "amber99sbnmr"]
forcefields = ["amber99sbildn"]  # For now let's look at a single FF
water_models = ["tip3p", "tip4pew", "tip3p-fb", "tip4p-fb"]

products = itertools.product(forcefields, water_models, sequences)


base_waters = {"tip3p":"tip3p", "tip4pew":"tip4pew", "tip3p-fb":"tip3p", "tip4p-fb":"tip4pew"}
