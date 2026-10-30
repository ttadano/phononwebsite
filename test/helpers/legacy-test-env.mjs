import fs from 'node:fs';
import path from 'node:path';
import { JSDOM } from 'jsdom';
import jqueryFactory from 'jquery';
import jsyaml from 'js-yaml';

import { Complex } from '../../src/legacycomplex.js';

function resolveAsset(url) {
  const clean = decodeURIComponent(url.split('?')[0]);
  return path.join(process.cwd(), clean);
}

function installAjaxStubs($) {
  $.getJSON = function (url, callback) {
    const text = fs.readFileSync(resolveAsset(url), 'utf8');
    callback(JSON.parse(text));
  };

  $.get = function (url, callback) {
    const filePath = resolveAsset(url);
    const text = fs.readFileSync(filePath, 'utf8');
    if (filePath.endsWith('.json')) {
      callback(JSON.parse(text));
    } else {
      callback(text);
    }
  };
}

export function setupLegacyTestEnv() {
  const dom = new JSDOM(
    `
    <!doctype html>
    <html>
    <body>
      <div id="vibcrystal"></div>
      <div id="highcharts"></div>
      <ul id="mat"></ul>
      <ul id="ref"></ul>
    </body>
    </html>
    `,
    { url: 'http://localhost/phonon.html' }
  );

  global.window = dom.window;
  global.document = dom.window.document;
  global.location = dom.window.location;
  Object.defineProperty(global, 'navigator', {
    configurable: true,
    writable: true,
    value: dom.window.navigator,
  });
  global.alert = () => {};
  global.jsyaml = jsyaml;
  global.Complex = Complex;

  const $ = jqueryFactory(dom.window);
  installAjaxStubs($);
  global.$ = $;
  global.jQuery = $;

  return { dom, $ };
}

export function teardownLegacyTestEnv(dom) {
  if (dom && dom.window) {
    dom.window.close();
  }
  delete global.window;
  delete global.document;
  delete global.location;
  Reflect.deleteProperty(global, 'navigator');
  delete global.alert;
  delete global.jsyaml;
  delete global.Complex;
  delete global.$;
  delete global.jQuery;
}

export async function loadPhononClasses() {
  const moduleUrl = new URL('../../src/phononwebsite.js', import.meta.url);
  moduleUrl.searchParams.set('t', String(Date.now()));
  return import(moduleUrl.href);
}
