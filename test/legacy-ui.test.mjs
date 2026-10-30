import assert from 'node:assert/strict';
import { afterEach, beforeEach, describe, it } from 'node:test';

import {
  loadPhononClasses,
  setupLegacyTestEnv,
  teardownLegacyTestEnv,
} from './helpers/legacy-test-env.mjs';

describe('Legacy compatibility: UI wiring', () => {
  let dom;
  let PhononWebpage;

  beforeEach(async () => {
    ({ dom } = setupLegacyTestEnv());
    ({ PhononWebpage } = await loadPhononClasses());
  });

  afterEach(() => {
    teardownLegacyTestEnv(dom);
  });

  it('updateMenu populates materials and references', () => {
    const visualizer = { update() {} };
    const dispersion = { setClickEvent() {}, update() {} };
    const p = new PhononWebpage(visualizer, dispersion);
    p.setMaterialsList(global.$('#mat'));
    p.setReferencesList(global.$('#ref'));

    p.updateMenu();

    assert.ok(global.$('#mat').children().length > 0, 'materials list should not be empty');
    assert.ok(global.$('#ref').children().length > 0, 'references list should not be empty');
  });
});
