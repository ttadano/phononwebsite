import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

import { getReasonableRepetitions } from '../src/utils.js';

describe('utils.getReasonableRepetitions', () => {
  it('returns [3,3,3] for tiny systems', () => {
    assert.deepEqual(getReasonableRepetitions(2), [3, 3, 3]);
    assert.deepEqual(getReasonableRepetitions(4), [3, 3, 3]);
  });

  it('returns [2,2,2] for small systems', () => {
    assert.deepEqual(getReasonableRepetitions(5), [2, 2, 2]);
    assert.deepEqual(getReasonableRepetitions(15), [2, 2, 2]);
  });

  it('returns [2,2,1] for medium systems', () => {
    assert.deepEqual(getReasonableRepetitions(16), [2, 2, 1]);
    assert.deepEqual(getReasonableRepetitions(50), [2, 2, 1]);
  });

  it('returns [1,1,1] for large systems', () => {
    assert.deepEqual(getReasonableRepetitions(51), [1, 1, 1]);
    assert.deepEqual(getReasonableRepetitions(500), [1, 1, 1]);
  });
});
