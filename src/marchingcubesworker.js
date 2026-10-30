import { buildMarchingCubesBuffers } from './marchingcubesgeometry.js';

self.onmessage = (event) => {
    const {
        requestId,
        values,
        sizex,
        sizey,
        sizez,
        gridCell,
        isolevel,
        options,
    } = event.data;

    const buffers = buildMarchingCubesBuffers(values, sizex, sizey, sizez, gridCell, isolevel, options);
    self.postMessage(
        {
            requestId,
            positions: buffers.positions.buffer,
            normals: buffers.normals.buffer,
        },
        [buffers.positions.buffer, buffers.normals.buffer],
    );
};
