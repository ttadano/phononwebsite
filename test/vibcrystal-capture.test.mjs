import assert from 'node:assert/strict';
import { afterEach, beforeEach, describe, it } from 'node:test';
import { JSDOM } from 'jsdom';

import { VibCrystal } from '../src/vibcrystal.js';

describe('VibCrystal capture lifecycle', () => {
  let dom;
  let fakeContainer;
  let captureInstances;
  let lastCapture;
  let originalCreateElement;
  let createdUrls;

  beforeEach(() => {
    dom = new JSDOM(`<!doctype html><html><body><div id="progress"></div></body></html>`);
    global.window = dom.window;
    global.document = dom.window.document;
    originalCreateElement = document.createElement.bind(document);

    fakeContainer = {
      width: () => 640,
      height: () => 480,
      get: () => ({
        getBoundingClientRect: () => ({ width: 640, height: 480 }),
        clientHeight: 480,
        style: {},
      }),
    };

    captureInstances = 0;
    lastCapture = null;
    createdUrls = [];
    global.GIF = function FakeGIF() {};
    global.MediaRecorder = class FakeMediaRecorder {
      constructor(stream, options = {}) {
        this.stream = stream;
        this.options = options;
        this.mimeType = options.mimeType || 'video/webm';
        this.state = 'inactive';
        this.ondataavailable = null;
        this.onstop = null;
      }
      static isTypeSupported(type) {
        return type.indexOf('video/webm') === 0;
      }
      start() {
        this.state = 'recording';
      }
      stop() {
        this.state = 'inactive';
        if (this.ondataavailable) {
          this.ondataavailable({ data: new Blob(['test'], { type: this.mimeType }), size: 4 });
        }
        if (this.onstop) {
          this.onstop();
        }
      }
    };
    global.window.URL.createObjectURL = function (blob) {
      createdUrls.push(blob);
      return 'blob:mock-capture';
    };
    global.CCapture = class FakeCCapture {
      constructor(options) {
        captureInstances += 1;
        this.options = options;
        this.started = false;
        this.stopped = false;
        this.saved = false;
        lastCapture = this;
      }
      start() {
        this.started = true;
      }
      stop() {
        this.stopped = true;
      }
      save(callback) {
        this.saved = true;
        callback('blob:mock-capture');
      }
      capture() {}
    };
  });

  afterEach(() => {
    document.createElement = originalCreateElement;
    delete global.window.URL.createObjectURL;
    dom.window.close();
    delete global.window;
    delete global.document;
    delete global.GIF;
    delete global.MediaRecorder;
    delete global.CCapture;
  });

  it('shows the webm button when the browser reports webm support', () => {
    const v = new VibCrystal(fakeContainer);
    const button = {
      hidden: false,
      hide() {
        this.hidden = true;
      },
      click(handler) {
        this.handler = handler;
      },
    };

    document.createElement = function(tagName) {
      if (tagName === 'canvas') {
        return {
          captureStream() {
            return {
              getTracks() {
                return [{ stop() {} }];
              },
            };
          },
          toDataURL(type) {
            return type === 'image/webp' ? 'data:image/webp;base64,AAAA' : '';
          },
        };
      }
      if (tagName === 'video') {
        return {
          canPlayType(type) {
            return type.indexOf('video/webm') !== -1 ? 'probably' : '';
          },
        };
      }
      return originalCreateElement(tagName);
    };

    v.setWebmButton(button);

    assert.equal(button.hidden, false);
    assert.equal(typeof button.handler, 'function');
  });

  it('uses MediaRecorder for webm capture when available', () => {
    const v = new VibCrystal(fakeContainer);
    v.speed = 1;
    v.time = 1.25;
    v.canvas = {
      captureStream() {
        return {
          getTracks() {
            return [{ stop() {} }];
          },
        };
      },
    };
    v.phonon = { name: 'Graphene demo' };

    v.capturestart('webm');

    assert.equal(v.captureState, 'capturing');
    assert.ok(v.captureRecorder);
    assert.equal(v.capturer, null);
    assert.equal(v.captureFrameTarget, v.fps);
    assert.equal(v.captureDurationMs, 1000);
    assert.equal(v.capturePhaseStart, 0.25);

    v.captureend('webm');

    assert.equal(v.captureState, 'idle');
    assert.equal(v.captureRecorder, null);
    assert.equal(createdUrls.length, 1);
  });

  it('ignores duplicate capture start while already capturing', () => {
    const v = new VibCrystal(fakeContainer);

    v.capturestart('gif');
    v.capturestart('gif');

    assert.equal(captureInstances, 1);
    assert.equal(v.captureState, 'capturing');
    assert.ok(v.capturer);
  });

  it('resets state/progress cleanly when capture ends', () => {
    const v = new VibCrystal(fakeContainer);
    v.phonon = { name: 'Graphene demo' };
    v.speed = 1;

    v.capturestart('gif');
    assert.ok(lastCapture);
    assert.equal(v.captureFrameCount, 0);
    assert.equal(v.captureFrameTarget, v.fps);

    lastCapture.options.onProgress(0.5);
    assert.equal(document.getElementById('progress').style.width, '50%');

    v.captureend('gif');

    assert.equal(v.captureState, 'idle');
    assert.equal(v.capturer, null);
    assert.equal(v.captureFrameCount, 0);
    assert.equal(v.captureFrameTarget, 0);
    assert.equal(document.getElementById('progress').style.width, '0%');
    assert.ok(lastCapture.stopped);
    assert.ok(lastCapture.saved);

    v.captureend('gif');
    assert.equal(v.captureState, 'idle');
  });

  it('includes k-point and mode indices in capture filename when available', () => {
    const v = new VibCrystal(fakeContainer);
    v.phonon = { name: 'Graphene demo' };
    v.captureK = 7;
    v.captureN = 2;

    assert.equal(v.getCaptureFilename('gif'), 'Graphene_demo_k7_n2.gif');
  });

  it('captures one full oscillation cycle based on speed', () => {
    const v = new VibCrystal(fakeContainer);
    v.speed = 2.0;

    assert.equal(v.getCaptureFrameTarget(), 30);
    assert.equal(v.getCaptureDurationMs(30), 500);
  });

  it('freezes animation time during capture and samples exact loop phases', () => {
    const v = new VibCrystal(fakeContainer);
    v.speed = 1;
    v.time = 2.4;
    v.captureState = 'capturing';
    v.captureFrameTarget = 60;
    v.captureFrameCount = 15;
    v.capturePhaseStart = 0.4;

    const timestamp = 1000;
    v.lastFrameTime = timestamp - 16.6666666667;
    v.controls = { update() {} };
    v.camera = {};
    v.scene = {};
    v.renderer = { render() {} };
    v.stats = null;
    v.atomobjects = [];
    v.atompos = [];
    v.vibrationComponents = [];
    v.atommeshes = [];
    v.bondmeshes = [];
    v.bonds = [];
    v.bondmesh = null;
    v.canvas = {};
    v.needsRender = true;
    v.paused = false;

    v.animate(timestamp);

    assert.equal(v.time, 2.4);
    assert.equal(v.capturePhaseStart + (v.captureFrameCount / v.captureFrameTarget), 0.65);
  });
});
