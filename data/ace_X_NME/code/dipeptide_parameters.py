from simtk import unit as u
import itertools

cutoff = 0.9 * u.nanometers
temperature = 298 * u.kelvin
timestep = 2.0 * u.femtoseconds
pressure = 1.0 * u.atmosphere
friction = 0.25 / u.picoseconds
n_steps = 100000000
output_frequency = 1000

ff_name = "amber99sbnmr"
water_name = "tip3p"

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]

capped_string = "capped"
sequences = ["ACE_%s_NME" % aa for aa in amino_acids]

products = itertools.product(forcefields, water_models, sequences)
