import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.modeling.physical_models import BlackBody



Sv = 891 * u.Jy  
d = 388 * u.pc  
d_cm = d.to(u.cm)  
kv = 1.84 * (u.cm**2 / u.g)  
Tdust = 36.1 * u.K  

Bv = BlackBody(Tdust)  

wavelength_micron = 850 * u.micron  
wavelength_cm = wavelength_micron.to(u.cm)
frequency = (const.c / wavelength_cm).to(u.Hz) 

bbRadiation = Bv(frequency).to(u.erg / (u.cm**2 * u.s * u.Hz * u.sr))  
bbRadiation = bbRadiation * u.sr  # Remove per-steradian factor
print(f"Blackbody Radiation: {bbRadiation}")

# Mdust = Sνd^2 / κνBν(Tdust)
mass_dust = (Sv * d_cm**2 / (kv * bbRadiation))
mass_dust = mass_dust.to(u.g)

print(f"Mass of dust: {mass_dust:.3e}")
print(f"Mass of gas: {mass_dust * 100:.3e}")
print(f"Mass of gas in solar masses: {(mass_dust).value * 100 / 1.989e+33:.3}")

volume = 3.32100273247897e+53 * u.cm**3
gas_density = mass_dust*100 / volume
Q = .5 # unitless
σv = 1.32e5 * (u.cm / u.s)
σθ = .0745 # radians but acts unitless

# B = Q * sqrt(2πp) * (σv / σθ)
magnetic_field_pos = (Q * np.sqrt(4 * np.pi * gas_density) * (σv/σθ)).value
print(f"Plane of sky magnetic field strength: {magnetic_field_pos * 1000:.5} mG")

# Bpos = (pi/4) |B| --> |B| = (4/pi) Bpos
magnetic_field = magnetic_field_pos * (4 / np.pi)
print(f"Magnetic field strength: {magnetic_field * 1000:.5} mG")

# (M / I) / (.12 / sqrt(G))     I = B * A
area = 5.797894231541227e+35
upperRatio = mass_dust.value*100 / (magnetic_field * area)
print(f"Ratio: {upperRatio / (.12 / np.sqrt(const.G.value * 1000)):.3}")