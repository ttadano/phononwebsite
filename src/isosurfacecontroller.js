import { MarchingCubesIsosurfaceBackend, RaymarchIsosurfaceBackend } from './isosurfacebackends.js';

export class IsosurfaceController {
    constructor(host) {
        this.host = host;
        this.mode = 'marching-cubes';
        this.backends = {
            'marching-cubes': new MarchingCubesIsosurfaceBackend(this),
            'raymarch-trilinear': new RaymarchIsosurfaceBackend(this, 'trilinear'),
            'raymarch-tricubic': new RaymarchIsosurfaceBackend(this, 'tricubic'),
        };
        this.backend = this.backends[this.mode];
    }

    clearPreviewCache() {
        this.backend.clearCache();
    }

    setMode(mode) {
        const normalizedMode = mode === 'raymarch' ? 'raymarch-trilinear' : mode;
        const nextMode = this.backends[normalizedMode] ? normalizedMode : 'marching-cubes';
        if (nextMode === this.mode) {
            return;
        }
        if (this.backend) {
            this.backend.dispose();
        }
        this.mode = nextMode;
        this.backend = this.backends[this.mode];
        this.host.removeNamedSceneObjects('isosurface');
    }

    getMode() {
        return this.mode;
    }

    getAvailableModes() {
        return Object.keys(this.backends);
    }

    requestUpdate() { this.backend.requestUpdate(); }
    updateSync() { this.backend.updateSync(); }
    updatePreview() { this.backend.updatePreview(); }
}
