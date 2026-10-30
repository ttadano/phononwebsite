import * as THREE from 'three';

import { edgeTable, triTable } from './static_libs/MarchingCubesData.js';

const triangleCountTable = new Uint8Array(256);
for (let cubeindex = 0; cubeindex < 256; cubeindex++) {
    let count = 0;
    const triOffset = cubeindex << 4;
    while (triTable[triOffset + count * 3] !== -1) {
        count += 1;
    }
    triangleCountTable[cubeindex] = count;
}

function getWrappedIndexAndShift(index, size) {
    if (size <= 0) {
        return { index: 0, shift: 0 };
    }
    const shift = Math.floor(index / size);
    const wrappedIndex = ((index % size) + size) % size;
    return { index: wrappedIndex, shift };
}

function getGridValue(values, x, y, z, sizex, sizey, sizez, size2, periodic) {
    if (!periodic) {
        return values[x + sizex * y + size2 * z];
    }

    const wx = getWrappedIndexAndShift(x, sizex);
    const wy = getWrappedIndexAndShift(y, sizey);
    const wz = getWrappedIndexAndShift(z, sizez);
    return values[wx.index + sizex * wy.index + size2 * wz.index];
}

function getGridPointCoords(target, targetOffset, x, y, z, sizex, sizey, sizez, periodic, gridCell) {
    if (!periodic) {
        const fx = sizex > 1 ? x / (sizex - 1) : 0;
        const fy = sizey > 1 ? y / (sizey - 1) : 0;
        const fz = sizez > 1 ? z / (sizez - 1) : 0;
        target[targetOffset] = fx * gridCell[0][0] + fy * gridCell[1][0] + fz * gridCell[2][0];
        target[targetOffset + 1] = fx * gridCell[0][1] + fy * gridCell[1][1] + fz * gridCell[2][1];
        target[targetOffset + 2] = fx * gridCell[0][2] + fy * gridCell[1][2] + fz * gridCell[2][2];
        return;
    }

    const wx = getWrappedIndexAndShift(x, sizex);
    const wy = getWrappedIndexAndShift(y, sizey);
    const wz = getWrappedIndexAndShift(z, sizez);
    const fx = wx.index / sizex + wx.shift;
    const fy = wy.index / sizey + wy.shift;
    const fz = wz.index / sizez + wz.shift;
    target[targetOffset] = fx * gridCell[0][0] + fy * gridCell[1][0] + fz * gridCell[2][0];
    target[targetOffset + 1] = fx * gridCell[0][1] + fy * gridCell[1][1] + fz * gridCell[2][1];
    target[targetOffset + 2] = fx * gridCell[0][2] + fy * gridCell[1][2] + fz * gridCell[2][2];
}

function interpolateVertex(target, targetOffset, pointA, pointAOffset, pointB, pointBOffset, valueA, valueB, isolevel) {
    const delta = valueB - valueA;
    const mu = delta === 0 ? 0.5 : (isolevel - valueA) / delta;
    target[targetOffset] = pointA[pointAOffset] + (pointB[pointBOffset] - pointA[pointAOffset]) * mu;
    target[targetOffset + 1] = pointA[pointAOffset + 1] + (pointB[pointBOffset + 1] - pointA[pointAOffset + 1]) * mu;
    target[targetOffset + 2] = pointA[pointAOffset + 2] + (pointB[pointBOffset + 2] - pointA[pointAOffset + 2]) * mu;
}

function computeCubeIndex(value0, value1, value2, value3, value4, value5, value6, value7, isolevel, insideIsAbove) {
    let cubeindex = 0;
    if (insideIsAbove) {
        if (value0 > isolevel) cubeindex |= 1;
        if (value1 > isolevel) cubeindex |= 2;
        if (value2 > isolevel) cubeindex |= 8;
        if (value3 > isolevel) cubeindex |= 4;
        if (value4 > isolevel) cubeindex |= 16;
        if (value5 > isolevel) cubeindex |= 32;
        if (value6 > isolevel) cubeindex |= 128;
        if (value7 > isolevel) cubeindex |= 64;
        return cubeindex;
    }

    if (value0 < isolevel) cubeindex |= 1;
    if (value1 < isolevel) cubeindex |= 2;
    if (value2 < isolevel) cubeindex |= 8;
    if (value3 < isolevel) cubeindex |= 4;
    if (value4 < isolevel) cubeindex |= 16;
    if (value5 < isolevel) cubeindex |= 32;
    if (value6 < isolevel) cubeindex |= 128;
    if (value7 < isolevel) cubeindex |= 64;
    return cubeindex;
}

function writeTriangle(positions, normals, writeOffset, vlist, vertexIndex1, vertexIndex2, vertexIndex3) {
    const v1 = vertexIndex1 * 3;
    const v2 = vertexIndex2 * 3;
    const v3 = vertexIndex3 * 3;

    positions[writeOffset] = vlist[v1];
    positions[writeOffset + 1] = vlist[v1 + 1];
    positions[writeOffset + 2] = vlist[v1 + 2];
    positions[writeOffset + 3] = vlist[v2];
    positions[writeOffset + 4] = vlist[v2 + 1];
    positions[writeOffset + 5] = vlist[v2 + 2];
    positions[writeOffset + 6] = vlist[v3];
    positions[writeOffset + 7] = vlist[v3 + 1];
    positions[writeOffset + 8] = vlist[v3 + 2];

    const abx = vlist[v2] - vlist[v1];
    const aby = vlist[v2 + 1] - vlist[v1 + 1];
    const abz = vlist[v2 + 2] - vlist[v1 + 2];
    const acx = vlist[v3] - vlist[v1];
    const acy = vlist[v3 + 1] - vlist[v1 + 1];
    const acz = vlist[v3 + 2] - vlist[v1 + 2];
    let nx = aby * acz - abz * acy;
    let ny = abz * acx - abx * acz;
    let nz = abx * acy - aby * acx;
    const length = Math.hypot(nx, ny, nz) || 1;
    nx /= length;
    ny /= length;
    nz /= length;

    for (let i = 0; i < 9; i += 3) {
        normals[writeOffset + i] = nx;
        normals[writeOffset + i + 1] = ny;
        normals[writeOffset + i + 2] = nz;
    }
}

export function buildMarchingCubesBuffers(values, sizex, sizey, sizez, gridCell, isolevel, options = {}) {
    const periodic = !!options.periodic;
    const insideIsAbove = !!options.insideIsAbove;
    const size2 = sizex * sizey;
    const vlist = new Float32Array(36);
    const points = new Float32Array(24);

    const xLimit = periodic ? sizex : sizex - 1;
    const yLimit = periodic ? sizey : sizey - 1;
    const zLimit = periodic ? sizez : sizez - 1;
    let triangleCount = 0;

    for (let z = 0; z < zLimit; z++) {
        for (let y = 0; y < yLimit; y++) {
            for (let x = 0; x < xLimit; x++) {
                const value0 = getGridValue(values, x, y, z, sizex, sizey, sizez, size2, periodic);
                const value1 = getGridValue(values, x + 1, y, z, sizex, sizey, sizez, size2, periodic);
                const value2 = getGridValue(values, x, y + 1, z, sizex, sizey, sizez, size2, periodic);
                const value3 = getGridValue(values, x + 1, y + 1, z, sizex, sizey, sizez, size2, periodic);
                const value4 = getGridValue(values, x, y, z + 1, sizex, sizey, sizez, size2, periodic);
                const value5 = getGridValue(values, x + 1, y, z + 1, sizex, sizey, sizez, size2, periodic);
                const value6 = getGridValue(values, x, y + 1, z + 1, sizex, sizey, sizez, size2, periodic);
                const value7 = getGridValue(values, x + 1, y + 1, z + 1, sizex, sizey, sizez, size2, periodic);
                const cubeindex = computeCubeIndex(value0, value1, value2, value3, value4, value5, value6, value7, isolevel, insideIsAbove);
                triangleCount += triangleCountTable[cubeindex];
            }
        }
    }

    const positions = new Float32Array(triangleCount * 9);
    const normals = new Float32Array(triangleCount * 9);
    let writeOffset = 0;

    for (let z = 0; z < zLimit; z++) {
        for (let y = 0; y < yLimit; y++) {
            for (let x = 0; x < xLimit; x++) {
                const value0 = getGridValue(values, x, y, z, sizex, sizey, sizez, size2, periodic);
                const value1 = getGridValue(values, x + 1, y, z, sizex, sizey, sizez, size2, periodic);
                const value2 = getGridValue(values, x, y + 1, z, sizex, sizey, sizez, size2, periodic);
                const value3 = getGridValue(values, x + 1, y + 1, z, sizex, sizey, sizez, size2, periodic);
                const value4 = getGridValue(values, x, y, z + 1, sizex, sizey, sizez, size2, periodic);
                const value5 = getGridValue(values, x + 1, y, z + 1, sizex, sizey, sizez, size2, periodic);
                const value6 = getGridValue(values, x, y + 1, z + 1, sizex, sizey, sizez, size2, periodic);
                const value7 = getGridValue(values, x + 1, y + 1, z + 1, sizex, sizey, sizez, size2, periodic);
                const cubeindex = computeCubeIndex(value0, value1, value2, value3, value4, value5, value6, value7, isolevel, insideIsAbove);
                const bits = edgeTable[cubeindex];
                if (bits === 0) {
                    continue;
                }

                getGridPointCoords(points, 0, x, y, z, sizex, sizey, sizez, periodic, gridCell);
                getGridPointCoords(points, 3, x + 1, y, z, sizex, sizey, sizez, periodic, gridCell);
                getGridPointCoords(points, 6, x, y + 1, z, sizex, sizey, sizez, periodic, gridCell);
                getGridPointCoords(points, 9, x + 1, y + 1, z, sizex, sizey, sizez, periodic, gridCell);
                getGridPointCoords(points, 12, x, y, z + 1, sizex, sizey, sizez, periodic, gridCell);
                getGridPointCoords(points, 15, x + 1, y, z + 1, sizex, sizey, sizez, periodic, gridCell);
                getGridPointCoords(points, 18, x, y + 1, z + 1, sizex, sizey, sizez, periodic, gridCell);
                getGridPointCoords(points, 21, x + 1, y + 1, z + 1, sizex, sizey, sizez, periodic, gridCell);

                if (bits & 1) interpolateVertex(vlist, 0, points, 0, points, 3, value0, value1, isolevel);
                if (bits & 2) interpolateVertex(vlist, 3, points, 3, points, 9, value1, value3, isolevel);
                if (bits & 4) interpolateVertex(vlist, 6, points, 6, points, 9, value2, value3, isolevel);
                if (bits & 8) interpolateVertex(vlist, 9, points, 0, points, 6, value0, value2, isolevel);
                if (bits & 16) interpolateVertex(vlist, 12, points, 12, points, 15, value4, value5, isolevel);
                if (bits & 32) interpolateVertex(vlist, 15, points, 15, points, 21, value5, value7, isolevel);
                if (bits & 64) interpolateVertex(vlist, 18, points, 18, points, 21, value6, value7, isolevel);
                if (bits & 128) interpolateVertex(vlist, 21, points, 12, points, 18, value4, value6, isolevel);
                if (bits & 256) interpolateVertex(vlist, 24, points, 0, points, 12, value0, value4, isolevel);
                if (bits & 512) interpolateVertex(vlist, 27, points, 3, points, 15, value1, value5, isolevel);
                if (bits & 1024) interpolateVertex(vlist, 30, points, 9, points, 21, value3, value7, isolevel);
                if (bits & 2048) interpolateVertex(vlist, 33, points, 6, points, 18, value2, value6, isolevel);

                const triOffset = cubeindex << 4;
                for (let i = 0; triTable[triOffset + i] !== -1; i += 3) {
                    writeTriangle(
                        positions,
                        normals,
                        writeOffset,
                        vlist,
                        triTable[triOffset + i],
                        triTable[triOffset + i + 1],
                        triTable[triOffset + i + 2],
                    );
                    writeOffset += 9;
                }
            }
        }
    }

    return {
        positions,
        normals,
    };
}

export function buildMarchingCubesGeometry(values, sizex, sizey, sizez, gridCell, isolevel, options = {}) {
    const buffers = buildMarchingCubesBuffers(values, sizex, sizey, sizez, gridCell, isolevel, options);
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(buffers.positions, 3));
    geometry.setAttribute('normal', new THREE.BufferAttribute(buffers.normals, 3));
    return geometry;
}
