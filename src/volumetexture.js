import * as THREE from 'three';

export function createVolumeTexture(values, sizex, sizey, sizez) {
    const textureValues = values instanceof Float32Array ? values : Float32Array.from(values);
    const DataTexture3DCtor = Reflect.get(THREE, 'Data3DTexture') || Reflect.get(THREE, 'DataTexture3D');
    const texture = new DataTexture3DCtor(textureValues, sizex, sizey, sizez);
    texture.format = THREE.RedFormat;
    texture.type = THREE.FloatType;
    texture.minFilter = THREE.LinearFilter;
    texture.magFilter = THREE.LinearFilter;
    texture.unpackAlignment = 1;
    texture.needsUpdate = true;
    return texture;
}

export function disposeVolumeTexture(texture) {
    if (texture && typeof texture.dispose === 'function') {
        texture.dispose();
    }
}
