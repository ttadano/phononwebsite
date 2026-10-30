import assert from 'node:assert/strict';
import { afterEach, beforeEach, describe, it } from 'node:test';

import {
  loadPhononClasses,
  setupLegacyTestEnv,
  teardownLegacyTestEnv,
} from './helpers/legacy-test-env.mjs';

describe('Legacy compatibility: construction', () => {
  let dom;
  let VibCrystal;
  let PhononHighcharts;
  let PhononWebpage;

  beforeEach(async () => {
    ({ dom } = setupLegacyTestEnv());
    ({ VibCrystal, PhononHighcharts, PhononWebpage } = await loadPhononClasses());
  });

  afterEach(() => {
    teardownLegacyTestEnv(dom);
  });

  it('initializes VibCrystal, PhononHighcharts and PhononWebpage', () => {
    const fakeContainer = {
      width: () => 640,
      height: () => 480,
      get: () => ({ getBoundingClientRect: () => ({ width: 640, height: 480 }) }),
    };

    const v = new VibCrystal(fakeContainer);
    assert.ok(v instanceof VibCrystal);

    const d = new PhononHighcharts(global.$('#highcharts'));
    assert.ok(d instanceof PhononHighcharts);

    const p = new PhononWebpage(v, d);
    assert.ok(p instanceof PhononWebpage);
  });
});
