import assert from 'node:assert/strict';
import { after, before, describe, it } from 'node:test';
import { spawn } from 'node:child_process';

let playwright = null;
try {
  playwright = await import('playwright');
} catch {
  // Playwright is optional for now.
}

const runBrowserSmoke = process.env.PLAYWRIGHT_SMOKE === '1' && !!playwright;

describe('Browser smoke (optional)', { concurrency: false }, () => {
  const port = 8123;
  const baseUrl = `http://127.0.0.1:${port}`;
  let server;

  before({ skip: !runBrowserSmoke }, async () => {
    server = spawn('python3', ['-m', 'http.server', String(port)], {
      cwd: process.cwd(),
      stdio: 'ignore',
    });

    await new Promise((resolve) => setTimeout(resolve, 800));
  });

  after(() => {
    if (server && !server.killed) {
      server.kill('SIGTERM');
    }
  });

  it('loads phonon.html with no page errors', { skip: !runBrowserSmoke, timeout: 30000 }, async () => {
    const browser = await playwright.chromium.launch({ headless: true });
    const page = await browser.newPage();
    const pageErrors = [];

    page.on('pageerror', (error) => {
      pageErrors.push(error.message || String(error));
    });

    await page.goto(`${baseUrl}/phonon.html`, { waitUntil: 'networkidle' });
    await page.waitForSelector('#highcharts');

    await browser.close();

    assert.equal(pageErrors.length, 0, `Browser page errors:\n${pageErrors.join('\n')}`);
  });

  it('loads exciton.html with no page errors', { skip: !runBrowserSmoke, timeout: 30000 }, async () => {
    const browser = await playwright.chromium.launch({ headless: true });
    const page = await browser.newPage();
    const pageErrors = [];

    page.on('pageerror', (error) => {
      pageErrors.push(error.message || String(error));
    });

    await page.goto(`${baseUrl}/exciton.html`, { waitUntil: 'networkidle' });
    await page.waitForSelector('#highcharts');

    await browser.close();

    assert.equal(pageErrors.length, 0, `Browser page errors:\n${pageErrors.join('\n')}`);
  });

  it('set PLAYWRIGHT_SMOKE=1 and install playwright to run', { skip: runBrowserSmoke }, () => {});
});
