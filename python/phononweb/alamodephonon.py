# Copyright (c) 2020, Terumasa Tadano
# All rights reserved.
#
#
""" Read phonon dispersion from alamode output """
import os
from math import sqrt
import numpy as np
from phononweb.phonon import Phonon
from phononweb.units import bohr_angstroem, atomic_numbers

class AlamodePhonon(Phonon):
    """
    Class to read phonons from ALAMODE

    Input:
        filename: HDF5 file containing polarization vectors written by anphon with PRINTEVEC = 1,
                  e.g. <prefix>.band.evec.h5 (ALAMODE >= 2.0) or <prefix>.band.evec.hdf5 (ALAMODE < 2.0)
    """

    def __init__(self,
                 filename,
                 name,
                 reps=(3, 3, 3),
                 reorder=False,
                 highsym_qpts=None,
                 band_range=None,
                 folder='.'):

        self.reps = reps
        self.name = name
        self.folder = folder
        self.filename = os.path.join(folder, filename)
        self.highsym_qpts = highsym_qpts
        self.atomic_masses = None
        self.band_range = band_range

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
        Function to read the eigenvalues and eigenvectors from ALAMODE evec.h5 (evec.hdf5) format
        """
        import h5py

        with h5py.File(self.filename, 'r') as hdf5file:
            self.eigenvalues = hdf5file['Eigenvalues']['frequencies'][()]
            vectors = hdf5file['Eigenvalues']['polarization_vectors'][()]
            self.qpoints = hdf5file['Kpoints']['kpoint_coordinates'][()] # reduced coordinates
            self.pos = hdf5file['PrimitiveCell']['fractional_coordinate'][()]
            self.cell = hdf5file['PrimitiveCell']['lattice_vector'][()] * bohr_angstroem
            self.chemical_symbols = np.array(hdf5file['PrimitiveCell']['elements'][()], dtype=str)
            kinds = hdf5file['PrimitiveCell']['atomic_kinds'][()]
            masskd = hdf5file['PrimitiveCell']['masses'][()]

        self.nqpoints = len(self.qpoints)
        self.atomic_masses = [masskd[i] for i in kinds]
        self.atom_types = [self.chemical_symbols[i] for i in kinds]
        self.atom_numbers = [atomic_numbers[x] for x in self.atom_types]
        self.atomic_numbers = [atomic_numbers[x] for x in self.chemical_symbols]
        self.natoms = len(self.pos)
        _, self.nphons = np.shape(self.eigenvalues)

        # the abinit eigenvectors are scaled with the atomic masses but the phonopy ones are not
        # so we always scale the eigenvectors with the atomic masses in the javascript of the website
        # here we scale then with sqrt(m) so that we recover the correct scalling on the website
        vectors = np.reshape(vectors, (self.nqpoints, self.nphons, self.natoms, 3, 2))

        if self.band_range:
            iband_s, iband_e = [int(i) for i in self.band_range.strip().split()]
            iband_s -= 1

            if iband_s < 0:
                iband_s = 0

            if iband_e > self.nphons:
                iband_e = self.nphons

            vectors = vectors[:,iband_s:iband_e,:,:,:]
            self.nphons = iband_e - iband_s
            self.eigenvalues = self.eigenvalues[:,iband_s:iband_e]

        for na in range(self.natoms):
           vectors[:,:,na,:,:] /= sqrt(self.atomic_masses[na])

        # normalize the eigenvectors
        self.eigenvectors = vectors / np.linalg.norm(vectors[0, 0])
        self.chemical_formula = self.get_chemical_formula()
