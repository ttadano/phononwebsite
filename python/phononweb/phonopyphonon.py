# Copyright (c) 2018, Henrique Miranda
# All rights reserved.
#
# This file is part of the phononwebsite project
#
""" Helpper function to create band-structures with phonopy """

import os
import json
import re
import copy
import numpy as np

from phonopy import Phonopy
from phonopy.units import Hartree, Bohr
from phonopy.interface.phonopy_yaml import *
import phonopy.file_IO as file_IO

class PhonopyPhonon():
    """
    Calculate the phonon dispersion from phononpy
    """
    def __init__(self,phonon):
        self.phonon = phonon

    @classmethod
    def from_files(self,phonon_yaml_filename,force_sets_filename,nac_filename=None):
        """initialize the PhonopyPhonon"""
        #get phonon_yaml
        ph_yaml = PhonopyYaml()
        ph_yaml.read(phonon_yaml_filename)
        # Phonopy API changed over time; support both old and current attributes.
        if hasattr(ph_yaml, "get_unitcell"):
            atoms = ph_yaml.get_unitcell()
        else:
            atoms = ph_yaml.unitcell

        if hasattr(ph_yaml, "supercell_matrix"):
            supercell_matrix = ph_yaml.supercell_matrix
        else:
            supercell_matrix = ph_yaml._data['supercell_matrix']

        #get force_sets
        force_sets = file_IO.parse_FORCE_SETS(filename=force_sets_filename)

        phonon = Phonopy(atoms,supercell_matrix, primitive_matrix='auto')
        if hasattr(phonon, "set_displacement_dataset"):
            phonon.set_displacement_dataset(force_sets)
        else:
            phonon.dataset = force_sets
        phonon.produce_force_constants()
        phonon.symmetrize_force_constants_by_space_group()

        #get NAC
        if nac_filename:
            primitive = phonon.get_primitive() if hasattr(phonon, "get_primitive") else phonon.primitive
            nac_params = file_IO.parse_BORN(primitive, filename=nac_filename)
            nac_factor = Hartree * Bohr
            if nac_params.get('factor') is None:
                nac_params['factor'] = nac_factor
            if hasattr(phonon, "set_nac_params"):
                phonon.set_nac_params(nac_params)
            else:
                phonon.nac_params = nac_params

        return PhonopyPhonon(phonon)

    def set_bandstructure_mp(self,mp_id,mp_api_key=None,band_points=5,verbose=False):
        """
        get bandstructure from the materials project
        """
        from pymatgen.ext.matproj import MPRester
        
        #start mprester
        self.mprester = MPRester(mp_api_key)

        #get bandstruccture
        bs = self.mprester.get_bandstructure_by_material_id(mp_id)

        #get high symmetry k-points
        if verbose: print("nkpoints:", len(bs.kpoints))
        branches = bs.as_dict()['branches']

        self.bands = []
        self.labels = [bs.kpoints[0].label]
        for path in branches:
            start = path['start_index'] 
            end   = path['end_index']

            start_kpoint = bs.kpoints[start].frac_coords
            end_kpoint   = bs.kpoints[end].frac_coords
            step_kpoint  = end_kpoint-start_kpoint

            self.labels.append(bs.kpoints[end].label)

            branch = []
            for i in range(band_points+1):
                branch.append(start_kpoint + float(i)/band_points*step_kpoint )
            self.bands.append(np.array(branch))

    def set_bandstructure_seekpath(self,reference_distance=0.1):
        """Get the bandstructure using seekpath"""
        import seekpath

        unitcell = self.phonon.get_unitcell() if hasattr(self.phonon, "get_unitcell") else self.phonon.unitcell
        cell = unitcell.get_cell() if hasattr(unitcell, "get_cell") else unitcell.cell
        atoms = unitcell.get_atomic_numbers() if hasattr(unitcell, "get_atomic_numbers") else unitcell.numbers
        pos = unitcell.get_scaled_positions() if hasattr(unitcell, "get_scaled_positions") else unitcell.scaled_positions
        
        path = seekpath.get_explicit_k_path((cell,pos,atoms),reference_distance=reference_distance)
        kpath  = path['explicit_kpoints_rel']
        explicit_labels = path['explicit_kpoints_labels']
        segments = path['segments'] if 'segments' in path else path['explicit_segments']
        self.bands = []
        self.labels = []
        for segment in segments:
            start_k, end_k = segment
            self.labels.append(explicit_labels[start_k])
            self.bands.append(kpath[start_k:end_k])
        self.labels.append(explicit_labels[-1])

    def set_bandstructure_seekpath_points(self, band_points=21):
        """Get the bandstructure using seekpath with a fixed number of points per segment."""
        import seekpath

        unitcell = self.phonon.get_unitcell() if hasattr(self.phonon, "get_unitcell") else self.phonon.unitcell
        cell = unitcell.get_cell() if hasattr(unitcell, "get_cell") else unitcell.cell
        atoms = unitcell.get_atomic_numbers() if hasattr(unitcell, "get_atomic_numbers") else unitcell.numbers
        pos = unitcell.get_scaled_positions() if hasattr(unitcell, "get_scaled_positions") else unitcell.scaled_positions

        path = seekpath.get_path((cell, pos, atoms))
        point_coords = path['point_coords']
        path_segments = path['path']

        segment_points = max(2, int(band_points))
        self.bands = []
        self.labels = []
        for start_label, end_label in path_segments:
            start_k = np.array(point_coords[start_label], dtype=float)
            end_k = np.array(point_coords[end_label], dtype=float)
            branch = np.linspace(start_k, end_k, num=segment_points, endpoint=True)
            self.bands.append(branch)
            self.labels.append(start_label)
        self.labels.append(path_segments[-1][1])

    def get_frequencies_with_eigenvectors(self,qpoint=(0,0,0)):
        """calculate the eigenvalues and eigenvectors at a specific qpoint"""
        frequencies, eigenvectors = self.phonon.get_frequencies_with_eigenvectors(qpoint)
        return frequencies, eigenvectors

    def get_bandstructure(self, is_eigenvectors=True, is_band_connection=True):
        """calculate the bandstructure"""
        if hasattr(self.phonon, "set_band_structure"):
            self.phonon.set_band_structure(
                self.bands,
                is_eigenvectors=is_eigenvectors,
                is_band_connection=is_band_connection,
            )
            return self.phonon.get_band_structure()

        self.phonon.run_band_structure(
            self.bands,
            with_eigenvectors=is_eigenvectors,
            is_band_connection=is_band_connection,
        )
        if hasattr(self.phonon, "get_band_structure_dict"):
            return self.phonon.get_band_structure_dict()
        return self.phonon.band_structure
 
    def write_disp_yaml(self,filename='disp.yaml'):
        """write disp yaml file"""
        displacements = self.phonon.get_displacements()
        directions = self.phonon.get_displacement_directions()
        supercell = self.phonon.get_supercell()
        file_IO.write_disp_yaml(displacements, supercell, directions=directions, filename=filename)

    def write_band_yaml(self,eigenvectors=True,filename='band.yaml'):
        """export a yaml file with the band-structure data"""
        if eigenvectors:
            phonon = self.phonon
        else:
            phonon = copy.deepcopy(self.phonon)
            phonon._band_structure._eigenvectors = None
        phonon.write_yaml_band_structure(filename=filename)
