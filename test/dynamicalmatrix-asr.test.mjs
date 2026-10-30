import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

import {
  buildDynamicalMatrixBlocks,
  getAsrCorrectedCompactForceConstants,
} from '../src/dynamicalmatrix.js';

describe('dynamicalmatrix acoustic sum rule', () => {
  it('imposes translational invariance on compact force constants', () => {
    const payload = {
      format: 'phonopy-dynamical-matrix-v1',
      primitive_natoms: 1,
      supercell_natoms: 2,
      masses: [1],
      force_constants_compact: [[
        [
          [2, 0, 0],
          [0, 0, 0],
          [0, 0, 0],
        ],
        [
          [0, 0, 0],
          [0, 0, 0],
          [0, 0, 0],
        ],
      ]],
      shortest_vectors: [[0, 0, 0], [0, 0, 0]],
      multiplicity: [
        [[1, 0]],
        [[1, 1]],
      ],
      s2pp_map: [0, 0],
    };

    const corrected = getAsrCorrectedCompactForceConstants(payload);
    assert.equal(corrected[0][0][0][0], 1);
    assert.equal(corrected[0][1][0][0], -1);

    const rowSum = corrected[0][0][0][0] + corrected[0][1][0][0];
    assert.ok(Math.abs(rowSum) < 1e-12);

    const matrix = buildDynamicalMatrixBlocks(payload, [0, 0, 0]);
    assert.ok(Math.abs(matrix.real[0][0]) < 1e-12);
  });

  it('allows ASR to be disabled explicitly', () => {
    const payload = {
      format: 'phonopy-dynamical-matrix-v1',
      acoustic_sum_rule: 'off',
      primitive_natoms: 1,
      supercell_natoms: 2,
      masses: [1],
      force_constants_compact: [[
        [
          [2, 0, 0],
          [0, 0, 0],
          [0, 0, 0],
        ],
        [
          [0, 0, 0],
          [0, 0, 0],
          [0, 0, 0],
        ],
      ]],
      shortest_vectors: [[0, 0, 0], [0, 0, 0]],
      multiplicity: [
        [[1, 0]],
        [[1, 1]],
      ],
      s2pp_map: [0, 0],
    };

    const corrected = getAsrCorrectedCompactForceConstants(payload);
    assert.equal(corrected, payload.force_constants_compact);

    const matrix = buildDynamicalMatrixBlocks(payload, [0, 0, 0]);
    assert.equal(matrix.real[0][0], 2);
  });
});
