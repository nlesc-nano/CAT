""" A module designed for running COSMO-RS Jobs. """

__all__ = ['CRSResults', 'CRSJob']

import numpy as np

try:
    import pandas as pd
except ImportError:
    pass

from scm.plams.core.basejob import SingleJob
from scm.plams.tools.units import Units
from scm.plams.interfaces.adfsuite.scmjob import (SCMJob, SCMResults)


class CRSResults(SCMResults):
    """
    A class for accessing results of COSMO-RS jobs.
    """
    _kfext = '.crskf'
    _rename_map = {'CRSKF': '$JN.crskf'}

    def get_energy(self, unit='kcal/mol'):
        """ Returns the solvation energy from a Activity Coefficients calculation. """
        E = self.readkf('ACTIVITYCOEF', 'deltag')[0]
        return Units.convert(E, 'kcal/mol', unit)

    def get_sigma(self, section, unit='kcal/mol'):
        """ Grab all values of sigma and the sigmapotential/profile;
        combine them into a dictionary or pandas dataframe. """
        sigma = self._sigma_y(section, unit)
        if self.readkf('PURE' + section) is not None:
            sigma['mixture'] = self._sigma_y(section, unit)
        sigma['sigma'] = self._sigma_x(section)
        try:
            return sigma.set_index('sigma')
        except AttributeError:
            return sigma

    def get_sigma_profile(self, unit='kcal/mol'):
        """ Returns all sigma profiles, expressed in *unit*.
        Returns a dictionary of numpy arrays or, if available, a pandas dataframe. """
        return self.get_sigma('SIGMAPOTENTIAL', unit)

    def get_sigma_potential(self):
        """ Returns all sigma profiles, expressed in *unit*.
        Returns a dictionary of numpy arrays or, if available, a pandas dataframe. """
        return self.get_sigma('SIGMAPROFILE')

    def _sigma_x(self, section):
        """ Construct all values of sigma. """
        min_max = self.readkf(section, 'sigmax')
        nitems = self.readkf(section, 'nitems')
        step = int((1 + 2 * min_max) / nitems)
        return np.arange(-min_max, min_max, step)

    def _sigma_y(self, section, unit='kcal/mol'):
        """ Get all values of . """
        values = np.array(self.readkf(section, 'profil'))
        values *= Units.conversion_ratio('kcal/mol', unit)
        if 'PURE' in section:
            ncomp = self.readkf(section, 'ncomp')
            values.shape = len(values) // ncomp, ncomp
        keys = self.readkf(section, 'filename')
        ret = dict(zip(keys, values))
        try:
            return pd.DataFrame(ret).set_index('sigma')
        except NameError:
            return ret


class CRSJob(SCMJob):
    """ A class for running COSMO-RS jobs. """
    _command = 'crs'
    _result_type = CRSResults

    def __init__(self, **kwargs):
        SingleJob.__init__(self, **kwargs)
        self.ignore_molecule = True
