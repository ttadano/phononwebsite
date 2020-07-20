# Copyright (c) 2020, Terumasa Tadano
# All rights reserved.
#
#
""" Read phonon dispersion from alamode output """
from math import pi
import numpy as np
from .pw import *
from .phononweb import Phonon, bohr_angstroem, atomic_numbers
import h5py

class AlamodePhonon(Phonon):
    """
    Class to read phonons from ALAMODE

    Input:
        prefix: <prefix>.evec.hdf5 file containing polarization vectors
    """

    def __init__(self, filename, name, reps=(3, 3, 3), reorder=False, highsym_qpts=None, folder='.'):
        self.reps = reps
        self.name = name
        self.folder = folder
        self.filename = "%s/%s" % (folder, filename)
        self.highsym_qpts = highsym_qpts
        self.atomic_masses = None

        # if the file already exists then we read it
        if os.path.isfile(self.filename):
            self.read_data_from_hdf5()
        else:
            raise ValueError('File {} does not exist.'.format(self.filename))

        # reorder eigenvales
        if reorder:
            self.reorder_eigenvalues()
        self.get_distances_qpts()
        self.labels_qpts = None

    def read_data_from_hdf5(self):
        """
        Function to read the eigenvalues and eigenvectors from ALAMODE evec.HDF5 format
        """

        hdf5file = h5py.File(self.filename, 'r')

        self.eigenvalues = hdf5file['Eigenvalues']['frequencies'][()]
        vectors = hdf5file['Eigenvalues']['polarization_vectors'][()]
        self.qpoints = hdf5file['Kpoints']['kpoint_coordinates'][()] # reduced coordinates
        self.pos = hdf5file['PrimitiveCell']['fractional_coordinate'][()]
        self.cell = hdf5file['PrimitiveCell']['lattice_vector'][()] * bohr_angstroem
        self.nqpoints = len(self.qpoints)
        self.chemical_symbols = np.array(hdf5file['PrimitiveCell']['elements'][()], dtype=str)
        kinds =  hdf5file['PrimitiveCell']['atomic_kinds'][()]
        masskd = hdf5file['PrimitiveCell']['masses'][()]
        self.atomic_masses = [masskd[i] for i in kinds]
        self.atom_types = [self.chemical_symbols[i] for i in kinds]
        self.atom_numbers = [atomic_numbers[x] for x in self.atom_types]
        self.atomic_numbers = np.unique(self.atom_numbers)
        self.natoms = len(self.pos)
        _, self.nphons = np.shape(self.eigenvalues)
        hdf5file.close()

        # the abinit eigenvectors are scaled with the atomic masses but the phonopy ones are not
        # so we always scale the eigenvectors with the atomic masses in the javascript of the website
        # here we scale then with sqrt(m) so that we recover the correct scalling on the website
        vectors = np.reshape(vectors, (self.nqpoints, self.nphons, self.natoms, 3, 2))

        for na in range(self.natoms):
           vectors[:,:,na,:,:] /= sqrt(self.atomic_masses[na])

        # normalize the eigenvectors
        self.eigenvectors = vectors / np.linalg.norm(vectors[0, 0])
        self.chemical_formula = self.get_chemical_formula()
