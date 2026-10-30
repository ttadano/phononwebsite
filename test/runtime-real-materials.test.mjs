import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import { afterEach, beforeEach, describe, it } from 'node:test';

import {
  setupLegacyTestEnv,
  teardownLegacyTestEnv,
} from './helpers/legacy-test-env.mjs';

import { PhononJson } from '../src/phononjson.js';

function loadGzJson(filePath) {
  return JSON.parse(zlib.gunzipSync(fs.readFileSync(filePath)).toString('utf8'));
}

function loadInternalJson(filePath) {
  return new Promise((resolve) => {
    const phonon = new PhononJson();
    phonon.getFromInternalJson(loadGzJson(filePath), () => resolve(phonon));
  });
}

function getMaxMismatch(summary) {
  let maxMismatch = 0;
  for (const item of summary) {
    if (item && item.maxDeltaCm1 !== null && Number.isFinite(item.maxDeltaCm1)) {
      maxMismatch = Math.max(maxMismatch, item.maxDeltaCm1);
    }
  }
  return maxMismatch;
}

describe('Runtime real-material diagnostics', () => {
  let dom;

  beforeEach(() => {
    ({ dom } = setupLegacyTestEnv());
  });

  afterEach(() => {
    teardownLegacyTestEnv(dom);
  });

  it('keeps a testable comparison harness for real runtime materials', async () => {
    const phonon = new PhononJson();
    phonon.kpoints = [[0, 0, 0], [0.1, 0, 0]];
    phonon.eigenvalues = [[1, 2], [3, 4]];
    phonon.canComputeEigenvectorsOnDemand = () => true;
    phonon.solveRuntimeEigenSystem = (qIndex) => ({
      eigenvectors: qIndex === 0
        ? [
          [[1, 0], [0, 0]],
          [[0, 0], [1, 0]],
        ]
        : [
          [[0, 0], [1, 0]],
          [[1, 0], [0, 0]],
        ],
      eigenvaluesCm1: qIndex === 0 ? [1, 2] : [4, 3],
    });

    const solution = await phonon.computeRuntimeBandConnectedEigenSolutionAsync();

    assert.ok(solution);
    assert.deepEqual(solution.eigenvalues, [[1, 2], [3, 4]]);
    assert.equal(getMaxMismatch(solution.mismatchSummary), 0);
  });

  it('compares computed runtime eigenvalues against stored values for local materials when enabled', async (t) => {
    if (process.env.PHONON_RUNTIME_REAL_MATERIALS !== '1') {
      t.skip('Set PHONON_RUNTIME_REAL_MATERIALS=1 to run local real-material diagnostics.');
      return;
    }

    const materials = ['mp-1000', 'mp-1175', 'mp-2918'];
    const results = [];

    for (const material of materials) {
      const runtimePath = path.join(process.cwd(), 'data/runtime-test', `${material}.json.gz`);
      const referencePath = path.join(process.cwd(), 'data/runtime-reference', `${material}.json.gz`);
      if (!fs.existsSync(runtimePath) || !fs.existsSync(referencePath)) {
        t.skip(`Missing local fixture for ${material}.`);
        return;
      }

      const runtime = await loadInternalJson(runtimePath);
      const reference = await loadInternalJson(referencePath);
      const solution = await runtime.computeRuntimeBandConnectedEigenSolutionAsync();

      assert.ok(solution, `${material} should produce a runtime band-connected solution`);
      assert.equal(
        solution.eigenvalues.length,
        reference.eigenvalues.length,
        `${material} should have the same q-point count`
      );

      results.push({
        material,
        maxMismatchCm1: getMaxMismatch(solution.mismatchSummary),
      });
    }

    for (const result of results) {
      assert.ok(
        result.maxMismatchCm1 < 0.5,
        `${result.material} runtime eigenvalues drift by ${result.maxMismatchCm1.toFixed(3)} cm^-1`
      );
    }
  });
});
