import assert from 'node:assert/strict';
import { afterEach, beforeEach, describe, it } from 'node:test';

import {
  setupLegacyTestEnv,
  teardownLegacyTestEnv,
} from './helpers/legacy-test-env.mjs';

import { PhononJson } from '../src/phononjson.js';

describe('PhononJson runtime eigenvectors', () => {
  let dom;

  beforeEach(() => {
    ({ dom } = setupLegacyTestEnv());
  });

  afterEach(() => {
    teardownLegacyTestEnv(dom);
  });

  it('computes missing q-point eigenvectors from a runtime dynamical-matrix payload', async () => {
    const phonon = new PhononJson();

    await new Promise((resolve) => {
      phonon.getFromInternalJson({
        name: 'Runtime',
        natoms: 1,
        atom_types: ['C'],
        atom_numbers: [6],
        atom_pos_car: [[0, 0, 0]],
        atom_pos_red: [[0, 0, 0]],
        lattice: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        qpoints: [[0, 0, 0]],
        distances: [0],
        formula: 'C',
        eigenvalues: [[10, 20, 30]],
        repetitions: [1, 1, 1],
        masses: [12],
        highsym_qpts: [[0, 'GAMMA']],
        line_breaks: [[0, 1]],
        dynamical_matrix: {
          format: 'phonopy-dynamical-matrix-v1',
          acoustic_sum_rule: 'off',
          primitive_natoms: 1,
          supercell_natoms: 1,
          primitive_lattice: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
          masses: [1],
          force_constants_compact: [[[
            [1, 0, 0],
            [0, 2, 0],
            [0, 0, 3],
          ]]],
          shortest_vectors: [[0, 0, 0]],
          multiplicity: [[[1, 0]]],
          s2pp_map: [0],
        },
      }, resolve);
    });

    assert.ok(Array.isArray(phonon.vec[0]));
    assert.equal(phonon.ensureQpointEigenvectors(0), true);
    assert.ok(Array.isArray(phonon.vec[0]));
    assert.equal(phonon.vec[0].length, 3);

    const firstMode = phonon.vec[0][0][0];
    assert.ok(Math.abs(firstMode[0][0]) > 0.99);
    assert.ok(Math.abs(firstMode[1][0]) < 1e-6);
    assert.ok(Math.abs(firstMode[2][0]) < 1e-6);
    assert.ok(Math.abs(phonon.eigenvalues[0][0] - 33.35641) < 1e-4);
  });

  it('uses the same band-connection mapping as the pre-generated path', async () => {
    const phonon = new PhononJson();

    const prevEigenvectors = [
      [[1, 0], [0, 0], [0, 0]],
      [[0, 0], [1, 0], [0, 0]],
      [[0, 0], [0, 0], [1, 0]],
    ];
    const currentEigenvectors = [
      [[0, 0], [1, 0], [0, 0]],
      [[1, 0], [0, 0], [0, 0]],
      [[0, 0], [0, 0], [1, 0]],
    ];
    const prevBandOrder = [0, 1, 2];

    assert.deepEqual(
      phonon.estimateRuntimeBandConnection(prevEigenvectors, currentEigenvectors, prevBandOrder),
      [1, 0, 2]
    );
  });

  it('keeps the Python-style continuous band connection across the full path', async () => {
    const phonon = new PhononJson();

    let calls = [];
    phonon.estimateRuntimeBandConnection = (prevEigenvectors, eigenvectors, prevBandOrder) => {
      calls.push(prevBandOrder.slice());
      return [1, 0, 2];
    };
    phonon.convertComplexEigenvectorsToInternal = (vectors) => vectors;
    phonon.applyModeAmplitudeConventionToEigenvectors = (vectors) => vectors;
    phonon.normalizeEigenvectorSet = () => {};
    phonon.invalidateEigenvectorCaches = () => {};
    phonon.canComputeEigenvectorsOnDemand = () => true;
    phonon.getNacQDirection = () => null;

    phonon.kpoints = [[0, 0, 0], [0.1, 0, 0], [0.1, 0, 0], [0, 0.1, 0]];
    phonon.line_breaks = [[0, 2], [2, 4]];
    phonon.vec = new Array(4).fill(null);
    phonon.eigenvalues = [
      [1, 2, 3],
      [1, 2, 3],
      [1, 2, 3],
      [1, 2, 3],
    ];
    phonon.dynamical_matrix = {};

    phonon.solveRuntimeEigenSystem = async () => ({
      eigenvectors: [
        [[1, 0], [0, 0], [0, 0]],
        [[0, 0], [1, 0], [0, 0]],
        [[0, 0], [0, 0], [1, 0]],
      ],
      eigenvaluesCm1: [1, 2, 3],
    });

    await phonon.computeAllRuntimeEigenvectorsWithBandConnectionAsync();

    assert.equal(calls.length, 3);
    assert.deepEqual(calls[0], [0, 1, 2]);
    assert.deepEqual(calls[1], [1, 0, 2]);
    assert.deepEqual(calls[2], [1, 0, 2]);
  });

  it('supports Gonze NAC payloads through the isolated correction path', async () => {
    const phonon = new PhononJson();

    await new Promise((resolve) => {
      phonon.getFromInternalJson({
        name: 'Runtime Gonze',
        natoms: 1,
        atom_types: ['C'],
        atom_numbers: [6],
        atom_pos_car: [[0, 0, 0]],
        atom_pos_red: [[0, 0, 0]],
        lattice: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        qpoints: [[0, 0, 0], [0.2, 0, 0]],
        distances: [0, 1],
        formula: 'C',
        eigenvalues: [[10, 20, 30], [10, 20, 30]],
        repetitions: [1, 1, 1],
        masses: [12],
        highsym_qpts: [[0, 'GAMMA']],
        line_breaks: [[0, 2]],
        dynamical_matrix: {
          format: 'phonopy-dynamical-matrix-v1',
          acoustic_sum_rule: 'off',
          primitive_natoms: 1,
          supercell_natoms: 1,
          primitive_lattice: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
          masses: [1],
          force_constants_compact: [[[
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
          ]]],
          shortest_vectors: [[0, 0, 0]],
          multiplicity: [[[1, 0]]],
          s2pp_map: [0],
          nac: {
            method: 'gonze',
            born: [[
              [1, 0, 0],
              [0, 0, 0],
              [0, 0, 0],
            ]],
            dielectric: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            nac_factor: 1,
            positions_car: [[0, 0, 0]],
            g_list: [[0, 0, 0]],
            dd_q0: {
              real: [[[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
              imag: [[[0, 0, 0], [0, 0, 0], [0, 0, 0]]],
            },
            lambda: 1,
            q_direction_tolerance: 1e-5,
          },
        },
      }, resolve);
    });

    assert.equal(phonon.ensureQpointEigenvectors(0), true);
    const hasLongitudinalMode = phonon.vec[0].some((mode) => {
      const atom = mode[0];
      return Math.abs(atom[0][0]) > 0.99
        && Math.abs(atom[1][0]) < 1e-6
        && Math.abs(atom[2][0]) < 1e-6;
    });
    assert.ok(hasLongitudinalMode);
  });
});
