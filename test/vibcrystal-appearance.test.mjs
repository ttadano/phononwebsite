import assert from 'node:assert/strict';
import { describe, it } from 'node:test';
import { JSDOM } from 'jsdom';
import jqueryFactory from 'jquery';

import { VibCrystal } from '../src/vibcrystal.js';
import { Complex } from '../src/legacycomplex.js';

function makeContainer() {
  return {
    width: () => 640,
    height: () => 480,
    get: () => ({
      getBoundingClientRect: () => ({ width: 640, height: 480 }),
      clientHeight: 480,
      style: {},
    }),
  };
}

function makeAppearanceDom() {
  return new JSDOM(`<!doctype html><html><body>
    <div id="atom-list"></div>
    <select id="display"><option value="jmol">Jmol</option><option value="vesta">Vesta</option></select>
    <input id="atom-color" type="color">
    <input id="arrow-color" type="color">
    <input id="bond-color" type="color">
    <input id="bond-color-by-atom" type="checkbox">
    <input id="atom-radius" type="number">
    <input id="bond-radius" type="number">
    <input id="arrow-radius" type="number">
    <div id="bond-rules-list"></div>
    <select id="bond-add-a"></select>
    <select id="bond-add-b"></select>
    <input id="bond-add-cutoff" type="number">
    <button id="reset-atom" type="button"></button>
    <button id="reset-bonds" type="button"></button>
    <button id="reset-vectors" type="button"></button>
  </body></html>`);
}

describe('VibCrystal advanced appearance', () => {
  it('supports atom color overrides and reset', () => {
    const v = new VibCrystal(makeContainer());
    v.display = 'jmol';

    const defaultColor = v.getDefaultAtomColor(6);
    assert.equal(v.getAtomColorHex(6), defaultColor);

    v.setAtomColorOverride(6, '#123456');
    assert.equal(v.getAtomColorHex(6), 0x123456);

    v.clearAtomColorOverride(6);
    assert.equal(v.getAtomColorHex(6), defaultColor);
  });

  it('preserves selected atom when rebuilding selectors', () => {
    const dom = makeAppearanceDom();
    const $ = jqueryFactory(dom.window);

    const v = new VibCrystal(makeContainer());
    v.atom_numbers = [6, 8, 6];
    v.updatelocal = () => {};
    v.setAdvancedAppearanceControls(
      $('#atom-list'),
      $('#display'),
      $('#atom-color'),
      $('#arrow-color'),
      $('#bond-color'),
      $('#bond-color-by-atom'),
      $('#atom-radius'),
      $('#bond-radius'),
      $('#arrow-radius'),
      $('#bond-rules-list'),
      $('#bond-add-a'),
      $('#bond-add-b'),
      $('#bond-add-cutoff'),
      $('#reset-atom'),
      $('#reset-bonds'),
      $('#reset-vectors')
    );

    v.adjustCovalentRadiiSelect();
    $('#atom-list').find('button[data-atom-number="8"]').trigger('click');
    v.adjustCovalentRadiiSelect();

    assert.equal(v.getSelectedAppearanceAtomNumber(), 8);
    dom.window.close();
  });

  it('resets atom, bonds and vectors sections independently', () => {
    const dom = makeAppearanceDom();
    const $ = jqueryFactory(dom.window);

    const v = new VibCrystal(makeContainer());
    v.atom_numbers = [6];
    v.updatelocal = () => {};
    v.setAdvancedAppearanceControls(
      $('#atom-list'),
      $('#display'),
      $('#atom-color'),
      $('#arrow-color'),
      $('#bond-color'),
      $('#bond-color-by-atom'),
      $('#atom-radius'),
      $('#bond-radius'),
      $('#arrow-radius'),
      $('#bond-rules-list'),
      $('#bond-add-a'),
      $('#bond-add-b'),
      $('#bond-add-cutoff'),
      $('#reset-atom'),
      $('#reset-bonds'),
      $('#reset-vectors')
    );
    v.adjustCovalentRadiiSelect();

    const defaultAtom = v.getDefaultAtomColor(6);
    v.setAtomColorOverride(6, '#123456');
    v.setAtomRadiusScaleOverride(6, 2.0);
    v.arrowcolor = 0x111111;
    v.bondscolor = 0x222222;
    v.bondColorByAtom = true;
    v.bondRadius = 0.3;
    v.arrowRadius = 0.4;
    v.arrowScale = 3.5;
    v.arrows = true;

    $('#reset-atom').trigger('click');
    assert.equal(v.getAtomColorHex(6), defaultAtom);
    assert.equal(v.getAtomRadiusScale(6), v.defaultAtomRadiusScale);
    assert.equal(v.bondscolor, 0x222222);

    $('#reset-bonds').trigger('click');
    assert.equal(v.bondscolor, v.defaultBondsColor);
    assert.equal(v.bondColorByAtom, v.defaultBondColorByAtom);
    assert.equal(v.bondRadius, v.defaultBondRadius);
    assert.equal(v.arrowcolor, 0x111111);

    $('#reset-vectors').trigger('click');
    assert.equal(v.arrowcolor, v.defaultArrowColor);
    assert.equal(v.arrowRadius, v.defaultArrowRadius);
    assert.equal(v.arrowScale, v.defaultArrowScale);
    assert.equal(v.arrows, false);
    dom.window.close();
  });

  it('adds and removes bond rules from the UI controls', () => {
    const dom = makeAppearanceDom();
    const $ = jqueryFactory(dom.window);

    const v = new VibCrystal(makeContainer());
    v.atom_numbers = [6, 8];
    v.updatelocal = () => {};
    v.setAdvancedAppearanceControls(
      $('#atom-list'),
      $('#display'),
      $('#atom-color'),
      $('#arrow-color'),
      $('#bond-color'),
      $('#bond-color-by-atom'),
      $('#atom-radius'),
      $('#bond-radius'),
      $('#arrow-radius'),
      $('#bond-rules-list'),
      $('#bond-add-a'),
      $('#bond-add-b'),
      $('#bond-add-cutoff'),
      $('#reset-atom'),
      $('#reset-bonds'),
      $('#reset-vectors')
    );
    v.adjustCovalentRadiiSelect();

    $('#bond-add-a').val('6');
    $('#bond-add-b').val('8');
    $('#bond-add-cutoff').val('1.23');
    $('#bond-add-cutoff').trigger($.Event('keydown', { key: 'Enter' }));
    assert.equal(v.hasBondRule(6, 8), true);
    const key = v.getBondRuleKey(6, 8);
    assert.equal(v.bondRules[key].cutoff, 1.23);

    $('#bond-rules-list').find(`button[data-remove-key="${key}"]`).trigger('click');
    assert.equal(v.hasBondRule(6, 8), false);
    dom.window.close();
  });

  it('updates split bond colors mode from the bonds controls', () => {
    const dom = makeAppearanceDom();
    const $ = jqueryFactory(dom.window);

    const v = new VibCrystal(makeContainer());
    v.atom_numbers = [6, 8];
    v.updatelocal = () => {};
    v.setAdvancedAppearanceControls(
      $('#atom-list'),
      $('#display'),
      $('#atom-color'),
      $('#arrow-color'),
      $('#bond-color'),
      $('#bond-color-by-atom'),
      $('#atom-radius'),
      $('#bond-radius'),
      $('#arrow-radius'),
      $('#bond-rules-list'),
      $('#bond-add-a'),
      $('#bond-add-b'),
      $('#bond-add-cutoff'),
      $('#reset-atom'),
      $('#reset-bonds'),
      $('#reset-vectors')
    );

    assert.equal(v.bondColorByAtom, true);
    $('#bond-color-by-atom').prop('checked', false).trigger('change');
    assert.equal(v.bondColorByAtom, false);

    $('#reset-bonds').trigger('click');
    assert.equal(v.bondColorByAtom, true);
    dom.window.close();
  });

  it('disables bond color input while split bond colors mode is enabled', () => {
    const dom = makeAppearanceDom();
    const $ = jqueryFactory(dom.window);

    const v = new VibCrystal(makeContainer());
    v.atom_numbers = [6, 8];
    v.updatelocal = () => {};
    v.setAdvancedAppearanceControls(
      $('#atom-list'),
      $('#display'),
      $('#atom-color'),
      $('#arrow-color'),
      $('#bond-color'),
      $('#bond-color-by-atom'),
      $('#atom-radius'),
      $('#bond-radius'),
      $('#arrow-radius'),
      $('#bond-rules-list'),
      $('#bond-add-a'),
      $('#bond-add-b'),
      $('#bond-add-cutoff'),
      $('#reset-atom'),
      $('#reset-bonds'),
      $('#reset-vectors')
    );

    assert.equal($('#bond-color').prop('disabled'), true);
    $('#bond-color-by-atom').prop('checked', false).trigger('change');
    assert.equal($('#bond-color').prop('disabled'), false);
    $('#bond-color').val('#123456').trigger('change');
    assert.equal(v.bondscolor, 0x123456);

    $('#bond-color-by-atom').prop('checked', true).trigger('change');
    assert.equal($('#bond-color').prop('disabled'), true);
    dom.window.close();
  });

  it('keeps default atom colors dynamic across display modes', () => {
    const dom = makeAppearanceDom();
    const $ = jqueryFactory(dom.window);

    const v = new VibCrystal(makeContainer());
    v.atom_numbers = [6];
    v.updatelocal = () => {};
    v.setAdvancedAppearanceControls(
      $('#atom-list'),
      $('#display'),
      $('#atom-color'),
      $('#arrow-color'),
      $('#bond-color'),
      $('#bond-color-by-atom'),
      $('#atom-radius'),
      $('#bond-radius'),
      $('#arrow-radius'),
      $('#bond-rules-list'),
      $('#bond-add-a'),
      $('#bond-add-b'),
      $('#bond-add-cutoff'),
      $('#reset-atom'),
      $('#reset-bonds'),
      $('#reset-vectors')
    );
    v.adjustCovalentRadiiSelect();

    const jmolDefault = v.getDefaultAtomColor(6);
    $('#atom-color').val(v.colorToInputHex(jmolDefault));
    $('#atom-radius').trigger($.Event('keydown', { key: 'Enter' }));

    v.display = 'vesta';
    const vestaDefault = v.getDefaultAtomColor(6);
    assert.equal(v.getAtomColorHex(6), vestaDefault);
    dom.window.close();
  });

  it('keeps manually added bond rules across phonon mode updates', () => {
    const v = new VibCrystal(makeContainer());
    v.initialized = true;
    v.updatelocal = () => {};

    const phonon = { atom_numbers: [6, 8] };
    const phononweb = {
      phonon,
      vibrations: [[
        { real: () => 0, imag: () => 0 },
        { real: () => 0, imag: () => 0 },
        { real: () => 0, imag: () => 0 },
      ]],
      atoms: [
        [0, 0, 0, 0],
        [1, 1.6, 0, 0],
      ],
      k: 0,
      n: 0,
    };

    v.update(phononweb);
    v.setBondRule(6, 8, 2.25);
    v.update(phononweb);

    assert.equal(v.bondRules[v.getBondRuleKey(6, 8)].cutoff, 2.25);
  });

  it('uses the closest observed pair distance as the default bond cutoff', () => {
    const v = new VibCrystal(makeContainer());
    v.phonon = { atom_numbers: [6, 8] };
    v.atoms = [
      [0, 0, 0, 0],
      [1, 2.05, 0, 0],
      [1, 4.40, 0, 0],
    ];

    const cutoff = v.getDefaultBondCutoff(6, 8);
    assert.equal(Number(cutoff.toFixed(2)), 2.11);
  });

  it('preserves real and imaginary vibration components from the local complex helper', () => {
    const v = new VibCrystal(makeContainer());
    v.init = () => {};
    v.updatelocal = () => {};
    v.initialized = true;

    v.update({
      phonon: { atom_numbers: [6] },
      atoms: [[0, 0, 0, 0]],
      vibrations: [[
        Complex(0.25, -0.5),
        Complex(-0.75, 0.125),
        Complex(0.5, 0.625),
      ]],
      k: 0,
      n: 0,
    });

    assert.deepEqual(v.vibrationComponents[0], [
      [0.25, -0.5],
      [-0.75, 0.125],
      [0.5, 0.625],
    ]);
  });
});
