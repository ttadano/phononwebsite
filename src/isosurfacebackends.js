import * as THREE from 'three';

import { buildMarchingCubesGeometry } from './marchingcubesgeometry.js';
import { createVolumeTexture, disposeVolumeTexture } from './volumetexture.js';
import { createRaymarchedIsosurface, supportsRaymarchedIsosurface, updateRaymarchedIsosurface } from './raymarchedisosurface.js';

class BaseIsosurfaceBackend {
    constructor(controller) {
        this.controller = controller;
    }

    get host() {
        return this.controller.host;
    }

    clearCache() {}

    dispose() {}

    applyGeometry(geometry) {
        this.host.removeNamedSceneObjects('isosurface');
        this.host.addMarchingCubesGeometry(geometry);
        this.host.render();
    }
}

export class MarchingCubesIsosurfaceBackend extends BaseIsosurfaceBackend {
    constructor(controller) {
        super(controller);
        this.worker = null;
        this.workerFailed = false;
        this.requestId = 0;
        this.workerBusy = false;
        this.previewCache = null;
    }

    clearCache() {
        this.previewCache = null;
    }

    dispose() {
        this.resetWorker();
        this.previewCache = null;
    }

    getPreviewStride() {
        const { sizex = 1, sizey = 1, sizez = 1 } = this.host;
        const voxelCount = sizex * sizey * sizez;
        const minSize = Math.min(sizex, sizey, sizez);
        if (minSize < 12) return 1;
        if (voxelCount > 800000 && minSize >= 24) return 4;
        if (voxelCount > 250000 && minSize >= 18) return 3;
        if (voxelCount > 80000 && minSize >= 12) return 2;
        return 1;
    }

    buildDownsampledField(values, sizex, sizey, sizez, stride, periodic) {
        if (stride <= 1) {
            return { values, sizex, sizey, sizez };
        }

        const reducedX = periodic ? Math.max(2, Math.floor(sizex / stride)) : Math.max(2, Math.floor((sizex - 1) / stride) + 1);
        const reducedY = periodic ? Math.max(2, Math.floor(sizey / stride)) : Math.max(2, Math.floor((sizey - 1) / stride) + 1);
        const reducedZ = periodic ? Math.max(2, Math.floor(sizez / stride)) : Math.max(2, Math.floor((sizez - 1) / stride) + 1);
        const reducedValues = new Float32Array(reducedX * reducedY * reducedZ);
        let writeIndex = 0;

        for (let z = 0; z < reducedZ; z++) {
            const sourceZ = periodic ? (z * stride) % sizez : Math.min(z * stride, sizez - 1);
            for (let y = 0; y < reducedY; y++) {
                const sourceY = periodic ? (y * stride) % sizey : Math.min(y * stride, sizey - 1);
                for (let x = 0; x < reducedX; x++) {
                    const sourceX = periodic ? (x * stride) % sizex : Math.min(x * stride, sizex - 1);
                    reducedValues[writeIndex] = values[sourceX + sizex * sourceY + sizex * sizey * sourceZ];
                    writeIndex += 1;
                }
            }
        }

        return { values: reducedValues, sizex: reducedX, sizey: reducedY, sizez: reducedZ };
    }

    getInteractiveInput() {
        const stride = this.getPreviewStride();
        const { values, sizex, sizey, sizez } = this.host;
        if (stride <= 1) {
            return { values, sizex, sizey, sizez };
        }

        const options = this.host.getMarchingCubesOptions();
        const periodic = !!options.periodic;
        const cacheKey = `${stride}:${sizex}:${sizey}:${sizez}:${periodic ? 'p' : 'n'}`;
        if (!this.previewCache || this.previewCache.key !== cacheKey) {
            this.previewCache = {
                key: cacheKey,
                data: this.buildDownsampledField(values, sizex, sizey, sizez, stride, periodic),
            };
        }

        return this.previewCache.data;
    }

    getWorker() {
        if (this.workerFailed || typeof Worker === 'undefined') {
            return null;
        }
        if (this.worker) {
            return this.worker;
        }

        try {
            this.worker = new Worker(new URL('./marchingcubesworker.js', import.meta.url), { type: 'module' });
            this.worker.onmessage = (event) => {
                this.workerBusy = false;
                const { requestId, positions, normals } = event.data;
                this.applyBuffers(requestId, new Float32Array(positions), new Float32Array(normals));
            };
            this.worker.onerror = () => {
                this.workerBusy = false;
                this.workerFailed = true;
                if (this.worker) {
                    this.worker.terminate();
                    this.worker = null;
                }
            };
        } catch (error) {
            this.workerFailed = true;
            this.worker = null;
        }

        return this.worker;
    }

    resetWorker() {
        this.workerBusy = false;
        if (this.worker) {
            this.worker.terminate();
            this.worker = null;
        }
    }

    createGeometryFromBuffers(positions, normals) {
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        geometry.setAttribute('normal', new THREE.BufferAttribute(normals, 3));
        return geometry;
    }

    applyBuffers(requestId, positions, normals) {
        if (requestId !== this.requestId || !this.host.scene) {
            return;
        }
        this.applyGeometry(this.createGeometryFromBuffers(positions, normals));
    }

    buildGeometry(values, sizex, sizey, sizez) {
        return buildMarchingCubesGeometry(
            values,
            sizex,
            sizey,
            sizez,
            this.host.gridCell,
            this.host.isolevel,
            this.host.getMarchingCubesOptions(),
        );
    }

    requestUpdate() {
        if (!Array.isArray(this.host.values) || !this.host.values.length || !this.host.scene) {
            return;
        }

        const requestId = ++this.requestId;
        const payload = {
            requestId,
            values: Float32Array.from(this.host.values),
            sizex: this.host.sizex,
            sizey: this.host.sizey,
            sizez: this.host.sizez,
            gridCell: this.host.gridCell,
            isolevel: this.host.isolevel,
            options: this.host.getMarchingCubesOptions(),
        };

        const worker = this.getWorker();
        if (worker) {
            if (this.workerBusy) {
                this.resetWorker();
            }
            const activeWorker = this.getWorker();
            if (activeWorker) {
                this.workerBusy = true;
                activeWorker.postMessage(payload, [payload.values.buffer]);
                return;
            }
        }

        this.applyGeometry(this.buildGeometry(
            this.host.values,
            this.host.sizex,
            this.host.sizey,
            this.host.sizez,
        ));
    }

    updateSync() {
        if (!this.host.scene || !Array.isArray(this.host.values) || !this.host.values.length) {
            return;
        }
        this.requestId += 1;
        if (this.workerBusy) {
            this.resetWorker();
        }
        this.applyGeometry(this.buildGeometry(
            this.host.values,
            this.host.sizex,
            this.host.sizey,
            this.host.sizez,
        ));
    }

    updatePreview() {
        if (!this.host.scene || !Array.isArray(this.host.values) || !this.host.values.length) {
            return;
        }
        this.requestId += 1;
        if (this.workerBusy) {
            this.resetWorker();
        }
        const preview = this.getInteractiveInput();
        this.applyGeometry(this.buildGeometry(
            preview.values,
            preview.sizex,
            preview.sizey,
            preview.sizez,
        ));
    }
}

export class RaymarchIsosurfaceBackend extends BaseIsosurfaceBackend {
    constructor(controller, interpolation = 'trilinear') {
        super(controller);
        this.interpolation = interpolation;
        this.volumeTexture = null;
        this.fallbackBackend = new MarchingCubesIsosurfaceBackend(controller);
        this.warned = false;
        this.volumeTextureKey = null;
    }

    clearCache() {
        this.fallbackBackend.clearCache();
    }

    dispose() {
        disposeVolumeTexture(this.volumeTexture);
        this.volumeTexture = null;
        this.volumeTextureKey = null;
        this.fallbackBackend.dispose();
    }

    ensureVolumeTexture() {
        const nextKey = [
            this.host.values,
            this.host.sizex,
            this.host.sizey,
            this.host.sizez,
        ];
        if (
            this.volumeTexture &&
            this.volumeTextureKey &&
            this.volumeTextureKey[0] === nextKey[0] &&
            this.volumeTextureKey[1] === nextKey[1] &&
            this.volumeTextureKey[2] === nextKey[2] &&
            this.volumeTextureKey[3] === nextKey[3]
        ) {
            return;
        }
        disposeVolumeTexture(this.volumeTexture);
        this.volumeTexture = createVolumeTexture(
            this.host.values,
            this.host.sizex,
            this.host.sizey,
            this.host.sizez,
        );
        this.volumeTextureKey = nextKey;
    }

    useFallback() {
        if (!this.warned) {
            console.warn('Raymarched isosurface backend is unavailable here; falling back to marching cubes.');
            this.warned = true;
        }
        return true;
    }

    hasIsosurfaceObjects() {
        if (!this.host.scene) {
            return false;
        }
        for (let i = 0; i < this.host.scene.children.length; i++) {
            if (this.host.scene.children[i] && this.host.scene.children[i].name === 'isosurface') {
                return true;
            }
        }
        return false;
    }

    getRenderConfig() {
        if (typeof this.host.getRaymarchRenderConfig === 'function') {
            return this.host.getRaymarchRenderConfig(this.interpolation);
        }
        return {
            interpolation: this.interpolation,
            stepCount: Math.max(96, Math.min(384, Math.round(Math.max(this.host.sizex, this.host.sizey, this.host.sizez) * 2.0))),
        };
    }

    updateExistingObjects() {
        if (!this.hasIsosurfaceObjects()) {
            return false;
        }

        const renderConfig = this.getRenderConfig();
        let updated = false;
        this.host.forEachNamedSceneObject('isosurface', (object) => {
            const textureRepeat = object.userData && object.userData.textureRepeat ? object.userData.textureRepeat : [1, 1, 1];
            const repeatScale = Math.max(textureRepeat[0] || 1, textureRepeat[1] || 1, textureRepeat[2] || 1);
            updated = updateRaymarchedIsosurface(object, {
                texture: this.volumeTexture,
                isolevel: this.host.isolevel,
                opacity: this.host.isosurfaceOpacity,
                color: 0xffff00,
                periodic: !!this.host.getMarchingCubesOptions().periodic,
                gridSize: [this.host.sizex, this.host.sizey, this.host.sizez],
                interpolation: renderConfig.interpolation,
                textureRepeat: textureRepeat,
                stepCount: renderConfig.stepCount * repeatScale,
                activeRayHits: renderConfig.activeRayHits,
            }) || updated;
        });

        if (updated) {
            this.host.render();
        }
        return updated;
    }

    requestUpdate() {
        if (!Array.isArray(this.host.values) || !this.host.values.length || !this.host.scene) {
            return;
        }
        this.ensureVolumeTexture();
        if (!this.volumeTexture || !supportsRaymarchedIsosurface(this.host.renderer)) {
            this.useFallback();
            this.fallbackBackend.requestUpdate();
            return;
        }
        if (this.updateExistingObjects()) {
            return;
        }
        const renderConfig = this.getRenderConfig();
        const object = createRaymarchedIsosurface({
            texture: this.volumeTexture,
            gridCell: this.host.gridCell,
            isolevel: this.host.isolevel,
            opacity: this.host.isosurfaceOpacity,
            color: 0xffff00,
            periodic: !!this.host.getMarchingCubesOptions().periodic,
            gridSize: [this.host.sizex, this.host.sizey, this.host.sizez],
            interpolation: renderConfig.interpolation,
            textureRepeat: [1, 1, 1],
            stepCount: renderConfig.stepCount,
            activeRayHits: renderConfig.activeRayHits,
        });
        if (!object) {
            this.useFallback();
            this.fallbackBackend.requestUpdate();
            return;
        }
        this.host.removeNamedSceneObjects('isosurface');
        this.host.addIsosurfaceObject(object);
        this.host.render();
    }

    updateSync() {
        this.requestUpdate();
    }

    updatePreview() {
        this.requestUpdate();
    }
}
