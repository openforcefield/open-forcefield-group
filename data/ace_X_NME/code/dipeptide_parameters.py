from simtk import unit as u
import itertools

padding = 1.0 * u.nanometers
cutoff = 0.9 * u.nanometers

temperature = 303 * u.kelvin

timestep = 2.0 * u.femtoseconds
equilibration_timestep = 1.0 * u.femtoseconds

barostat_frequency = 25
pressure = 1.0 * u.atmosphere
friction = 0.25 / u.picoseconds

n_steps = 100000000
n_equil_steps = 5000000
output_frequency = 1000
equilibrate_output_frequency = n_equil_steps - 1

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]

sequences = ["ACE_%s_NME" % aa for aa in amino_acids]

forcefields = ["amber99sbildn", "amber96", "amber99sbnmr"]
water_models = ["tip3p", "tip4pew", "tip3p-fb", "tip4p-fb"]

products = itertools.product(forcefields, water_models, sequences)


base_waters = {"tip3p":"tip3p", "tip4pew":"tip4pew", "tip3p-fb":"tip3p", "tip4p-fb":"tip4pew"}

