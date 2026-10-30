import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';

function read(relPath) {
  const absPath = path.join(process.cwd(), relPath);
  return fs.readFileSync(absPath, 'utf8');
}

describe('Build smoke tests', () => {
  it('produces minified bundle artifacts', () => {
    const requiredFiles = [
      'build/main.js',
      'build/main.js.map',
      'build/main.min.js',
      'build/main.min.js.map',
      'build/exciton.js',
      'build/exciton.js.map',
      'build/exciton.min.js',
      'build/exciton.min.js.map',
    ];

    for (const file of requiredFiles) {
      assert.ok(fs.existsSync(path.join(process.cwd(), file)), `Missing ${file}`);
    }
  });

  it('does not leak Node fs imports into browser bundle', () => {
    const mainJs = read('build/main.js');
    const mainMinJs = read('build/main.min.js');
    assert.ok(!mainJs.includes("from 'fs'"), 'build/main.js must not import fs');
    assert.ok(!mainMinJs.includes("from 'fs'"), 'build/main.min.js must not import fs');
  });

  it('uses minified production entrypoint', () => {
    const html = read('phonon.html');
    assert.ok(
      html.includes('<script src="main.min.js" type="module"></script>'),
      'phonon.html should load main.min.js'
    );
  });

  it('uses minified exciton entrypoint', () => {
    const html = read('exciton.html');
    assert.ok(
      html.includes('<script src="exciton.min.js" type="module"></script>'),
      'exciton.html should load exciton.min.js'
    );
  });

  it('keeps exciton datasets in the shared standard format', () => {
    const datasets = [
      'data/excitondb/bn/absorptionspectra.json',
      'data/excitondb/mos2/absorptionspectra.json',
      'data/excitondb/mote2/absorptionspectra.json',
    ];

    for (const dataset of datasets) {
      const data = JSON.parse(read(dataset));
      assert.ok(Array.isArray(data['E/ev[1]']), `${dataset} should provide E/ev[1]`);
      assert.ok(Array.isArray(data['EPS-Im[2]']), `${dataset} should provide EPS-Im[2]`);
      assert.ok(Array.isArray(data['EPS-Re[3]']), `${dataset} should provide EPS-Re[3]`);
      assert.ok(Array.isArray(data['EPSo-Im[4]']), `${dataset} should provide EPSo-Im[4]`);
      assert.ok(Array.isArray(data['EPSo-Re[5]']), `${dataset} should provide EPSo-Re[5]`);
      assert.ok(Array.isArray(data.lattice), `${dataset} should provide a primitive lattice`);
      assert.ok(Array.isArray(data.supercell_lattice), `${dataset} should provide a supercell lattice for the grid`);
      assert.ok(!('eps' in data), `${dataset} should not use the legacy eps matrix format`);
    }
  });
});
