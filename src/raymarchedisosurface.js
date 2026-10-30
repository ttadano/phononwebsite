import * as THREE from 'three';
import { RAYMARCH_PERFORMANCE_SETTINGS } from './raymarchperformance.js';

const MAX_STEPS = 384;
const MAX_REFINEMENT_STEPS = 2;
const MAX_RAY_HITS = 2;
const RAYMARCH_OPACITY_COMPENSATION = RAYMARCH_PERFORMANCE_SETTINGS.opacityCompensation;
const BoxGeometryCtor = Reflect.get(THREE, 'BoxGeometry') || Reflect.get(THREE, 'BoxBufferGeometry');

export function clampRaymarchStepCount(stepCount, gridSize = [1, 1, 1]) {
    return Math.max(24, Math.min(MAX_STEPS, Math.round(stepCount || (Math.max(...gridSize) * 2.0))));
}

function createCellMatrix(gridCell) {
    const a = gridCell[0];
    const b = gridCell[1];
    const c = gridCell[2];

    return new THREE.Matrix4().set(
        a[0], b[0], c[0], 0,
        a[1], b[1], c[1], 0,
        a[2], b[2], c[2], 0,
        0, 0, 0, 1,
    );
}

function updateCellMatrix(matrix, gridCell) {
    const a = gridCell[0];
    const b = gridCell[1];
    const c = gridCell[2];

    return matrix.set(
        a[0], b[0], c[0], 0,
        a[1], b[1], c[1], 0,
        a[2], b[2], c[2], 0,
        0, 0, 0, 1,
    );
}

function createRaymarchGeometry(gridCell) {
    const geometry = new BoxGeometryCtor(1, 1, 1);
    geometry.translate(0.5, 0.5, 0.5);
    geometry.setAttribute('texturePosition', geometry.attributes.position.clone());
    geometry.applyMatrix4(createCellMatrix(gridCell));
    return geometry;
}

function getInterpolationMode(interpolation) {
    return interpolation === 'tricubic' ? 1 : 0;
}

function createMaterial({ texture, gridCell, isolevel, opacity, color, periodic, gridSize, interpolation, textureRepeat, stepCount, activeRayHits }) {
    const cellMatrix4 = createCellMatrix(gridCell);
    const objectToTextureMatrix = new THREE.Matrix4().copy(cellMatrix4).invert();
    const gradientMatrix = new THREE.Matrix3().setFromMatrix4(objectToTextureMatrix).transpose();
    const colorValue = new THREE.Color(color);

    return new THREE.ShaderMaterial({
        glslVersion: THREE.GLSL3,
        side: THREE.DoubleSide,
        transparent: true,
        depthWrite: opacity >= 0.999,
        uniforms: {
            uVolume: { value: texture },
            uIsolevel: { value: isolevel },
            uOpacity: { value: opacity },
            uColor: { value: colorValue },
            uGridSize: { value: new THREE.Vector3(gridSize[0], gridSize[1], gridSize[2]) },
            uPeriodic: { value: periodic ? 1 : 0 },
            uGradientMatrix: { value: gradientMatrix },
            uWorldToTexture: { value: new THREE.Matrix4() },
            uTextureToWorld: { value: new THREE.Matrix4() },
            uStepCount: { value: clampRaymarchStepCount(stepCount, gridSize) },
            uActiveRayHits: { value: Math.max(1, Math.min(MAX_RAY_HITS, Math.round(activeRayHits || MAX_RAY_HITS))) },
            uInterpolationMode: { value: getInterpolationMode(interpolation) },
            uTextureRepeat: { value: new THREE.Vector3(textureRepeat[0], textureRepeat[1], textureRepeat[2]) },
        },
        vertexShader: `
            in vec3 texturePosition;

            out vec3 vTexturePosition;
            void main() {
                vTexturePosition = texturePosition;
                gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
            }
        `,
        fragmentShader: `
            precision highp float;
            precision highp sampler3D;

            uniform sampler3D uVolume;
            uniform float uIsolevel;
            uniform float uOpacity;
            uniform vec3 uColor;
            uniform vec3 uGridSize;
            uniform int uPeriodic;
            uniform mat3 uGradientMatrix;
            uniform mat4 uWorldToTexture;
            uniform mat4 uTextureToWorld;
            uniform int uStepCount;
            uniform int uActiveRayHits;
            uniform int uInterpolationMode;
            uniform vec3 uTextureRepeat;

            in vec3 vTexturePosition;
            out vec4 outColor;

            vec3 wrapSamplePosition(vec3 coord) {
                coord *= uTextureRepeat;
                if (uPeriodic == 1) {
                    return fract(coord);
                }
                return clamp(coord, vec3(0.0), vec3(1.0));
            }

            int wrapIndex(int index, int size) {
                if (uPeriodic == 1) {
                    int wrapped = index % size;
                    return wrapped < 0 ? wrapped + size : wrapped;
                }
                return clamp(index, 0, size - 1);
            }

            float fetchVoxel(ivec3 index) {
                ivec3 dims = ivec3(max(uGridSize, vec3(1.0)));
                ivec3 wrapped = ivec3(
                    wrapIndex(index.x, dims.x),
                    wrapIndex(index.y, dims.y),
                    wrapIndex(index.z, dims.z)
                );
                return texelFetch(uVolume, wrapped, 0).r;
            }

            float cubicCatmullRom(float p0, float p1, float p2, float p3, float t) {
                float a = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
                float b = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3;
                float c = -0.5 * p0 + 0.5 * p2;
                float d = p1;
                float value = ((a * t + b) * t + c) * t + d;
                float lower = min(min(p0, p1), min(p2, p3));
                float upper = max(max(p0, p1), max(p2, p3));
                return clamp(value, lower, upper);
            }

            float sampleVolumeTrilinear(vec3 coord) {
                vec3 wrappedCoord = wrapSamplePosition(coord);
                vec3 dims = max(uGridSize, vec3(2.0));
                vec3 gridCoord = (uPeriodic == 1)
                    ? wrappedCoord * dims
                    : wrappedCoord * (dims - vec3(1.0));

                vec3 baseCoord = floor(gridCoord);
                vec3 fraction = fract(gridCoord);

                if (uPeriodic != 1) {
                    vec3 maxBase = dims - vec3(2.0);
                    baseCoord = min(baseCoord, maxBase);
                    fraction = gridCoord - baseCoord;
                }

                ivec3 base = ivec3(baseCoord);

                float c000 = fetchVoxel(base + ivec3(0, 0, 0));
                float c100 = fetchVoxel(base + ivec3(1, 0, 0));
                float c010 = fetchVoxel(base + ivec3(0, 1, 0));
                float c110 = fetchVoxel(base + ivec3(1, 1, 0));
                float c001 = fetchVoxel(base + ivec3(0, 0, 1));
                float c101 = fetchVoxel(base + ivec3(1, 0, 1));
                float c011 = fetchVoxel(base + ivec3(0, 1, 1));
                float c111 = fetchVoxel(base + ivec3(1, 1, 1));

                float c00 = mix(c000, c100, fraction.x);
                float c10 = mix(c010, c110, fraction.x);
                float c01 = mix(c001, c101, fraction.x);
                float c11 = mix(c011, c111, fraction.x);
                float c0 = mix(c00, c10, fraction.y);
                float c1 = mix(c01, c11, fraction.y);
                return mix(c0, c1, fraction.z);
            }

            float sampleVolumeTricubic(vec3 coord) {
                vec3 wrappedCoord = wrapSamplePosition(coord);
                vec3 dims = max(uGridSize, vec3(2.0));
                vec3 gridCoord = (uPeriodic == 1)
                    ? wrappedCoord * dims
                    : wrappedCoord * (dims - vec3(1.0));
                vec3 baseCoord = floor(gridCoord);
                vec3 fraction = fract(gridCoord);

                if (uPeriodic != 1) {
                    vec3 maxBase = dims - vec3(2.0);
                    baseCoord = min(baseCoord, maxBase);
                    fraction = gridCoord - baseCoord;
                }

                ivec3 base = ivec3(baseCoord);

                float yzPlane[4];
                for (int kz = 0; kz < 4; kz++) {
                    float yLine[4];
                    for (int ky = 0; ky < 4; ky++) {
                        float xLine[4];
                        for (int kx = 0; kx < 4; kx++) {
                            ivec3 sampleIndex = base + ivec3(kx - 1, ky - 1, kz - 1);
                            xLine[kx] = fetchVoxel(sampleIndex);
                        }
                        yLine[ky] = cubicCatmullRom(xLine[0], xLine[1], xLine[2], xLine[3], fraction.x);
                    }
                    yzPlane[kz] = cubicCatmullRom(yLine[0], yLine[1], yLine[2], yLine[3], fraction.y);
                }

                return cubicCatmullRom(yzPlane[0], yzPlane[1], yzPlane[2], yzPlane[3], fraction.z);
            }

            float sampleVolume(vec3 coord) {
                if (uInterpolationMode == 1) {
                    return sampleVolumeTricubic(coord);
                }
                return sampleVolumeTrilinear(coord);
            }

            bool intersectBox(vec3 rayOrigin, vec3 rayDirection, out float tMin, out float tMax) {
                vec3 invDir = 1.0 / max(abs(rayDirection), vec3(1e-6)) * sign(rayDirection);
                vec3 t0 = (vec3(0.0) - rayOrigin) * invDir;
                vec3 t1 = (vec3(1.0) - rayOrigin) * invDir;
                vec3 tSmall = min(t0, t1);
                vec3 tLarge = max(t0, t1);
                tMin = max(max(tSmall.x, tSmall.y), tSmall.z);
                tMax = min(min(tLarge.x, tLarge.y), tLarge.z);
                return tMax > max(tMin, 0.0);
            }

            vec3 computeTextureGradient3(vec3 coord) {
                vec3 eps = 1.0 / max(uGridSize, vec3(2.0));
                vec3 e1 = vec3( eps.x, -eps.y, -eps.z);
                vec3 e2 = vec3(-eps.x, -eps.y,  eps.z);
                vec3 e3 = vec3(-eps.x,  eps.y, -eps.z);
                vec3 e4 = vec3( eps.x,  eps.y,  eps.z);
                float s1 = sampleVolume(coord + e1);
                float s2 = sampleVolume(coord + e2);
                float s3 = sampleVolume(coord + e3);
                float s4 = sampleVolume(coord + e4);
                return e1 * s1 + e2 * s2 + e3 * s3 + e4 * s4;
            }

            vec3 estimateNormal(vec3 coord) {
                vec3 gradient = computeTextureGradient3(coord);
                return normalize(uGradientMatrix * gradient);
            }

            vec3 refineZeroCrossing(vec3 lowPos, vec3 highPos, float lowValue, float highValue) {
                vec3 bestPos = mix(lowPos, highPos, 0.5);

                for (int refine = 0; refine < ${MAX_REFINEMENT_STEPS}; refine++) {
                    float denom = highValue - lowValue;
                    float secantWeight = 0.5;
                    if (abs(denom) > 1e-7) {
                        secantWeight = clamp(-lowValue / denom, 0.0, 1.0);
                    }

                    vec3 secantPos = mix(lowPos, highPos, secantWeight);
                    float secantValue = sampleVolume(secantPos) - uIsolevel;
                    bestPos = secantPos;

                    if (abs(secantValue) < 1e-5) {
                        return secantPos;
                    }

                    vec3 gradient = computeTextureGradient3(secantPos);
                    float derivative = dot(gradient, normalize(highPos - lowPos));
                    if (abs(derivative) > 1e-6) {
                        float newtonStep = clamp(secantValue / derivative, -0.5, 0.5);
                        vec3 newtonPos = secantPos - normalize(highPos - lowPos) * newtonStep * length(highPos - lowPos);
                        newtonPos = clamp(newtonPos, min(lowPos, highPos), max(lowPos, highPos));
                        float newtonValue = sampleVolume(newtonPos) - uIsolevel;
                        if (abs(newtonValue) < abs(secantValue)) {
                            bestPos = newtonPos;
                            secantPos = newtonPos;
                            secantValue = newtonValue;
                        }
                    }

                    if ((lowValue <= 0.0 && secantValue <= 0.0) || (lowValue >= 0.0 && secantValue >= 0.0)) {
                        lowPos = secantPos;
                        lowValue = secantValue;
                    } else {
                        highPos = secantPos;
                        highValue = secantValue;
                    }
                }

                return bestPos;
            }

            vec4 shadeHit(vec3 hitPos, vec3 rayDirection) {
                vec3 normal = estimateNormal(hitPos);
                if (dot(normal, rayDirection) > 0.0) {
                    normal = -normal;
                }
                vec3 lightDir = normalize(vec3(0.45, 0.55, 1.0));
                vec3 hitWorld = (uTextureToWorld * vec4(hitPos, 1.0)).xyz;
                vec3 viewDir = normalize(cameraPosition - hitWorld);
                float diffuse = 0.28 + 0.72 * max(dot(normal, lightDir), 0.0);
                vec3 halfVector = normalize(lightDir + viewDir);
                float specular = pow(max(dot(normal, halfVector), 0.0), 24.0);
                vec3 shadedColor = uColor * diffuse;
                shadedColor += uColor * (0.035 * specular) + vec3(0.012 * specular);
                shadedColor = min(shadedColor, vec3(1.0));
                float clampedOpacity = clamp(uOpacity, 0.0, 0.9999);
                float targetOpacity = min(0.9999, clampedOpacity * ${RAYMARCH_OPACITY_COMPENSATION.toFixed(2)});
                float hitCount = float(max(uActiveRayHits, 1));
                float effectiveOpacity = 1.0 - pow(max(1.0 - targetOpacity, 1e-4), 1.0 / hitCount);
                return vec4(shadedColor, effectiveOpacity);
            }

            void main() {
                vec3 cameraTex = (uWorldToTexture * vec4(cameraPosition, 1.0)).xyz;
                vec3 rayDirection = normalize(vTexturePosition - cameraTex);

                float tEnter;
                float tExit;
                if (!intersectBox(cameraTex, rayDirection, tEnter, tExit)) {
                    discard;
                }

                float startT = max(tEnter, 0.0);
                float endT = tExit;
                float span = endT - startT;
                if (span <= 0.0) {
                    discard;
                }

                float stepSize = span / float(max(uStepCount, 1));
                float voxelAdvance = 1.5 / max(max(uGridSize.x, uGridSize.y), uGridSize.z);
                float minHitAdvance = max(stepSize * 2.5, voxelAdvance);
                float currentStartT = startT;
                vec4 accumulated = vec4(0.0);

                for (int hitIndex = 0; hitIndex < ${MAX_RAY_HITS}; hitIndex++) {
                    if (hitIndex >= uActiveRayHits) {
                        break;
                    }
                    vec3 previousPos = cameraTex + rayDirection * currentStartT;
                    float previousValue = sampleVolume(previousPos) - uIsolevel;
                    bool foundHit = false;
                    float nextStartT = endT;

                    for (int stepIndex = 1; stepIndex <= ${MAX_STEPS}; stepIndex++) {
                        if (stepIndex > uStepCount) {
                            break;
                        }

                        float t = currentStartT + float(stepIndex) * stepSize;
                        if (t > endT) {
                            break;
                        }

                        vec3 currentPos = cameraTex + rayDirection * t;
                        float currentValue = sampleVolume(currentPos) - uIsolevel;

                        if ((previousValue <= 0.0 && currentValue >= 0.0) || (previousValue >= 0.0 && currentValue <= 0.0)) {
                            vec3 hitPos = refineZeroCrossing(previousPos, currentPos, previousValue, currentValue);
                            float hitT = dot(hitPos - cameraTex, rayDirection);
                            vec4 hitColor = shadeHit(hitPos, rayDirection);
                            accumulated.rgb += (1.0 - accumulated.a) * hitColor.a * hitColor.rgb;
                            accumulated.a += (1.0 - accumulated.a) * hitColor.a;
                            nextStartT = hitT + minHitAdvance;
                            foundHit = true;
                            break;
                        }

                        previousPos = currentPos;
                        previousValue = currentValue;
                    }

                    if (!foundHit || accumulated.a >= 0.94 || nextStartT >= endT) {
                        break;
                    }

                    currentStartT = nextStartT;
                }

                if (accumulated.a <= 0.0) {
                    discard;
                }
                vec3 finalColor = accumulated.rgb / max(accumulated.a, 1e-5);
                outColor = vec4(clamp(finalColor, vec3(0.0), vec3(1.0)), accumulated.a);
            }
        `,
    });
}

function attachRaymarchUniformUpdater(mesh, gridCell) {
    const state = mesh.userData.raymarchState || {
        cellMatrix: createCellMatrix(gridCell),
        textureToWorld: new THREE.Matrix4(),
        worldToTexture: new THREE.Matrix4(),
    };
    mesh.userData.raymarchState = state;
    mesh.onBeforeRender = function() {
        refreshRaymarchedIsosurfaceUniforms(this, gridCell);
    };
}

export function refreshRaymarchedIsosurfaceUniforms(mesh, fallbackGridCell = null) {
    if (!mesh || !mesh.material || !mesh.material.uniforms) {
        return;
    }
    const activeGridCell = (mesh.userData && mesh.userData.gridCell) ? mesh.userData.gridCell : fallbackGridCell;
    if (!Array.isArray(activeGridCell) || activeGridCell.length !== 3) {
        return;
    }
    const textureRepeat = (mesh.userData && mesh.userData.textureRepeat) ? mesh.userData.textureRepeat : [1, 1, 1];
    const state = mesh.userData && mesh.userData.raymarchState ? mesh.userData.raymarchState : null;
    const cellMatrix = state ? state.cellMatrix : createCellMatrix(activeGridCell);
    const textureToWorld = state ? state.textureToWorld : new THREE.Matrix4();
    const worldToTexture = state ? state.worldToTexture : new THREE.Matrix4();
    updateCellMatrix(cellMatrix, activeGridCell);
    mesh.updateMatrixWorld(true);
    textureToWorld.multiplyMatrices(mesh.matrixWorld, cellMatrix);
    worldToTexture.copy(textureToWorld).invert();
    mesh.material.uniforms.uTextureToWorld.value.copy(textureToWorld);
    mesh.material.uniforms.uWorldToTexture.value.copy(worldToTexture);
    if (mesh.material.uniforms.uTextureRepeat) {
        mesh.material.uniforms.uTextureRepeat.value.set(textureRepeat[0], textureRepeat[1], textureRepeat[2]);
    }
}

export function createRaymarchedIsosurface({
    texture,
    gridCell,
    isolevel,
    opacity,
    color = 0xffff00,
    periodic = false,
    gridSize = [1, 1, 1],
    interpolation = 'trilinear',
    textureRepeat = [1, 1, 1],
    stepCount = null,
    activeRayHits = MAX_RAY_HITS,
}) {
    if (!texture || !Array.isArray(gridCell) || gridCell.length !== 3) {
        return null;
    }

    const geometry = createRaymarchGeometry(gridCell);
    const material = createMaterial({
        texture,
        gridCell,
        isolevel,
        opacity,
        color,
        periodic,
        gridSize,
        interpolation,
        textureRepeat,
        stepCount,
        activeRayHits,
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.name = 'isosurface';
    mesh.frustumCulled = false;
    mesh.userData.isRaymarchedIsosurface = true;
    mesh.userData.gridCell = gridCell.map((vector) => vector.slice());
    mesh.userData.textureRepeat = textureRepeat.slice();
    attachRaymarchUniformUpdater(mesh, mesh.userData.gridCell);
    refreshRaymarchedIsosurfaceUniforms(mesh, mesh.userData.gridCell);
    return mesh;
}

export function updateRaymarchedIsosurface(object, {
    texture,
    isolevel,
    opacity,
    color = 0xffff00,
    periodic = false,
    gridSize = [1, 1, 1],
    interpolation = 'trilinear',
    textureRepeat = [1, 1, 1],
    stepCount = null,
    activeRayHits = MAX_RAY_HITS,
}) {
    if (!object || !object.material || !object.material.uniforms) {
        return false;
    }

    const material = object.material;
    const uniforms = object.material.uniforms;
    if (uniforms.uVolume) {
        uniforms.uVolume.value = texture;
    }
    if (uniforms.uIsolevel) {
        uniforms.uIsolevel.value = isolevel;
    }
    if (uniforms.uOpacity) {
        uniforms.uOpacity.value = opacity;
    }
    if (uniforms.uColor) {
        uniforms.uColor.value.set(color);
    }
    if (uniforms.uGridSize) {
        uniforms.uGridSize.value.set(gridSize[0], gridSize[1], gridSize[2]);
    }
    if (uniforms.uPeriodic) {
        uniforms.uPeriodic.value = periodic ? 1 : 0;
    }
    if (uniforms.uStepCount) {
        uniforms.uStepCount.value = clampRaymarchStepCount(stepCount, gridSize);
    }
    if (uniforms.uActiveRayHits) {
        uniforms.uActiveRayHits.value = Math.max(1, Math.min(MAX_RAY_HITS, Math.round(activeRayHits || MAX_RAY_HITS)));
    }
    if (uniforms.uInterpolationMode) {
        uniforms.uInterpolationMode.value = getInterpolationMode(interpolation);
    }
    if (uniforms.uTextureRepeat) {
        uniforms.uTextureRepeat.value.set(textureRepeat[0], textureRepeat[1], textureRepeat[2]);
    }
    if (object.userData) {
        object.userData.textureRepeat = textureRepeat.slice();
    }
    material.opacity = opacity;
    if (material.transparent !== true) {
        material.transparent = true;
        material.needsUpdate = true;
    }
    const nextDepthWrite = opacity >= 0.999;
    if (material.depthWrite !== nextDepthWrite) {
        material.depthWrite = nextDepthWrite;
        material.needsUpdate = true;
    }
    refreshRaymarchedIsosurfaceUniforms(object);
    return true;
}

export function supportsRaymarchedIsosurface(renderer) {
    return !!(renderer && renderer.capabilities && renderer.capabilities.isWebGL2);
}
