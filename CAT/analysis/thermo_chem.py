""" A module related to calculating thermochemical properties. """

__all__ = ['get_thermo', 'get_entropy']

import numpy as np

from scm.plams.tools.units import Units


def get_entropy(mol, freqs, T=298.15):
    """
    Calculate the translational, vibrational and rotational entropy.
    All units and constants are in SI units.

    mol <plams.Molecule>: A PLAMS molecule.
    freqs <np.ndarray>: An iterable consisting of vibrational frequencies in units of s**-1.
    T <float>: The temperature in Kelvin
    Return <np.ndarray>: A numpy array containing the translational, rotational and
        vibrational contributions to the entropy
    """
    if not isinstance(freqs, np.ndarray):
        freqs = np.array(freqs)

    # Define constants
    kT = 1.380648 * 10**-23 * T  # Boltzmann constant * temperature
    h = 6.6260701 * 10**-34  # Planck constant
    hv_kT = (h * freqs) / kT  # (Planck constant * frequencies) / (Boltzmann * temperature)
    R = 8.31445  # Gas constant
    V_Na = ((R * T) / 10**5) / Units.constants['NA']  # Volume(1 mol ideal gas) / Avogadro's number
    pi = np.pi

    # Extract atomic masses and coordinates
    m = np.array([at.mass for at in mol]) * 1.6605390 * 10**-27
    x, y, z = mol.as_array().T * 10**-10

    # Calculate the rotational partition function
    inertia = np.array([sum(m*(y**2 + z**2)), -sum(m*x*y), -sum(m*x*z),
                        -sum(m*x*y), sum(m*(x**2 + z**2)), -sum(m*y*z),
                        -sum(m*x*z), -sum(m*y*z), sum(m*(x**2 + y**2))]).reshape(3, 3)
    inertia = np.product(np.linalg.eig(inertia)[0])
    q_rot = pi**0.5 * ((8 * pi**2 * kT) / h**2)**1.5 * inertia**0.5

    # Calculate the translational, rotational and vibrational entropy (divided by R)
    S_trans = 1.5 + np.log(V_Na * ((2 * pi * sum(m) * kT) / h**2)**1.5)
    S_rot = 1.5 + np.log(q_rot)
    S_vib = sum(hv_kT / np.expm1(hv_kT) - np.log(1 - np.exp(-hv_kT)))

    return R * np.array([S_trans, S_rot, S_vib])


def get_thermo(mol, freqs, E, T=298.15, export=['E', 'H', 'S', 'G'], unit='kcal/mol'):
    """
    Extract and return Gibbs free energies, entropies and/or enthalpies from an AMS KF file.
    All vibrational frequencies smaller than 100 cm**-1 are set to 100 cm**-1.

    mol <plams.Molecule>: A PLAMS molecule.
    freqs <np.array>: An iterable consisting of vibrational frequencies in units of cm**-1.
    E <float>: The eletronic energy in kcal/mol.
    T <float>: The temperature in Kelvin
    export <list>[<str>]: An iterable containing strings of the to be exported energies:
        'E': Electronic energy
        'U': Interal energy (E + U_nuc)
        'H': Enthalpy (U + pV)
        'S': Entropy
        'G': Gibbs free energy (H - T*S)
    unit <str>: The unit of the to be returned energies.
    Return <float> or <dict>[<float>]: An energy or dictionary of energies
    """
    # Get frequencies; set all frequencies smaller than 100 cm**-1 to 100 cm**-1
    freqs = np.array(freqs)
    freqs[freqs < 100] = 100
    freqs *= 100 * Units.constants['c']

    # hv_kT = (Planck constant * frequencies) / (Boltzmann constant * temperature)
    hv_kT = (6.6260701 * 10**-34 * freqs) / (1.380648 * 10**-23 * T)
    RT = 8.31445 * T  # Gas constant * temperature

    # Extract and/or calculate the various energies
    E = E * Units.conversion_ratio('kcal/mol', 'kj/mol') * 1000
    U = E + RT * (3.0 + sum(0.5 * hv_kT + hv_kT / np.expm1(hv_kT)))
    H = U + RT
    S = sum(get_entropy(mol, freqs, T=T))
    G = H - T * S

    ret = {'E': E, 'U': U, 'H': H, 'S': S, 'G': G}

    if len(export) == 1:
        return Units.convert(ret[export[0]], 'kj/mol', unit) / 1000
    return {i: Units.convert(ret[i], 'kj/mol', unit) / 1000 for i in ret if i in export}
