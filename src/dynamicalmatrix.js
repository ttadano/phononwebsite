import { getGonzeReciprocalCorrection } from './gonze.js';
import { solveComplexHermitianWithEigenWasm } from './eigenwasm.js';

function dotShortestVectorWithQ(shortestVector, qpoint) {
    return shortestVector[0] * qpoint[0] + shortestVector[1] * qpoint[1] + shortestVector[2] * qpoint[2];
}

function determinant3x3(matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

function invert3x3(matrix) {
    let det = determinant3x3(matrix);
    if (Math.abs(det) < 1e-16) {
        return null;
    }
    let invDet = 1 / det;
    return [
        [
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) * invDet,
            (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) * invDet,
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * invDet,
        ],
        [
            (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) * invDet,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * invDet,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) * invDet,
        ],
        [
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) * invDet,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) * invDet,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * invDet,
        ],
    ];
}

function cloneCompactForceConstants(forceConstants) {
    let cloned = new Array(forceConstants.length);
    for (let i = 0; i < forceConstants.length; i++) {
        cloned[i] = new Array(forceConstants[i].length);
        for (let j = 0; j < forceConstants[i].length; j++) {
            cloned[i][j] = [
                forceConstants[i][j][0].slice(),
                forceConstants[i][j][1].slice(),
                forceConstants[i][j][2].slice(),
            ];
        }
    }
    return cloned;
}

function imposeTranslationalInvariancePerIndex(forceConstants, index) {
    let primaryLength = index === 0 ? forceConstants[0].length : forceConstants.length;
    let correctionLength = index === 0 ? forceConstants.length : forceConstants[0].length;

    for (let primaryIndex = 0; primaryIndex < primaryLength; primaryIndex++) {
        for (let axisRow = 0; axisRow < 3; axisRow++) {
            for (let axisCol = 0; axisCol < 3; axisCol++) {
                let drift = 0;
                for (let correctionIndex = 0; correctionIndex < correctionLength; correctionIndex++) {
                    if (index === 0) {
                        drift += forceConstants[correctionIndex][primaryIndex][axisRow][axisCol];
                    } else {
                        drift += forceConstants[primaryIndex][correctionIndex][axisRow][axisCol];
                    }
                }
                drift /= correctionLength;
                for (let correctionIndex = 0; correctionIndex < correctionLength; correctionIndex++) {
                    if (index === 0) {
                        forceConstants[correctionIndex][primaryIndex][axisRow][axisCol] -= drift;
                    } else {
                        forceConstants[primaryIndex][correctionIndex][axisRow][axisCol] -= drift;
                    }
                }
            }
        }
    }
}

function getAcousticSumRuleMode(payload) {
    if (!payload || payload.acoustic_sum_rule === false || payload.acoustic_sum_rule === 'off') {
        return 'off';
    }
    if (typeof payload.acoustic_sum_rule === 'string') {
        return payload.acoustic_sum_rule;
    }
    return payload && payload.format === 'phonopy-dynamical-matrix-v1'
        ? 'translational'
        : 'off';
}

export function getAsrCorrectedCompactForceConstants(payload) {
    if (!payload || !payload.force_constants_compact) {
        return null;
    }
    if (getAcousticSumRuleMode(payload) !== 'translational') {
        return payload.force_constants_compact;
    }
    if (!payload._cached_asr_force_constants_compact) {
        let corrected = cloneCompactForceConstants(payload.force_constants_compact);
        // Compact force constants are stored as (primitive atom, supercell atom).
        // The dominant translational drift relevant for matching the pre-generated
        // path is along the supercell-atom index, which is the one directly used
        // in the dynamical-matrix sum.
        imposeTranslationalInvariancePerIndex(corrected, 1);
        payload._cached_asr_force_constants_compact = corrected;
    }
    return payload._cached_asr_force_constants_compact;
}

export function buildDynamicalMatrixBlocks(payload, qpoint, qDirection = null) {
    let masses = payload.masses;
    let forceConstants = getAsrCorrectedCompactForceConstants(payload);
    let shortestVectors = payload.shortest_vectors;
    let multiplicity = payload.multiplicity;
    let s2ppMap = payload.s2pp_map;
    let natoms = masses.length;
    let primitiveToSupercellRatio = payload.supercell_natoms / payload.primitive_natoms;
    let nacCorrection = null;

    if (payload.nac) {
        if (payload.nac.method !== 'gonze') {
            throw new Error(`Unsupported NAC method: ${payload.nac.method}`);
        }
        nacCorrection = getGonzeReciprocalCorrection(payload, qpoint, qDirection);
    }

    let matrixReal = [];
    let matrixImag = [];
    for (let row = 0; row < natoms * 3; row++) {
        matrixReal.push(new Array(natoms * 3).fill(0));
        matrixImag.push(new Array(natoms * 3).fill(0));
    }

    for (let i = 0; i < natoms; i++) {
        for (let j = 0; j < natoms; j++) {
            let sqrtMass = Math.sqrt(masses[i] * masses[j]);
            let blockReal = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
            let blockImag = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

            if (nacCorrection) {
                for (let axisRow = 0; axisRow < 3; axisRow++) {
                    for (let axisCol = 0; axisCol < 3; axisCol++) {
                        blockReal[axisRow][axisCol] += nacCorrection[i][j][axisRow][axisCol][0];
                        blockImag[axisRow][axisCol] += nacCorrection[i][j][axisRow][axisCol][1];
                    }
                }
            }

            for (let superIndex = 0; superIndex < s2ppMap.length; superIndex++) {
                if (s2ppMap[superIndex] !== j) {
                    continue;
                }

                let count = multiplicity[superIndex][i][0];
                let address = multiplicity[superIndex][i][1];
                let phaseReal = 0;
                let phaseImag = 0;
                for (let vectorIndex = 0; vectorIndex < count; vectorIndex++) {
                    let phase = 2 * Math.PI * dotShortestVectorWithQ(shortestVectors[address + vectorIndex], qpoint);
                    phaseReal += Math.cos(phase);
                    phaseImag += Math.sin(phase);
                }

                let factor = 1 / (sqrtMass * count);
                let fcBlock = forceConstants[i][superIndex];
                for (let axisRow = 0; axisRow < 3; axisRow++) {
                    for (let axisCol = 0; axisCol < 3; axisCol++) {
                        let value = fcBlock[axisRow][axisCol];
                        value *= factor;
                        blockReal[axisRow][axisCol] += value * phaseReal;
                        blockImag[axisRow][axisCol] += value * phaseImag;
                    }
                }
            }

            for (let axisRow = 0; axisRow < 3; axisRow++) {
                for (let axisCol = 0; axisCol < 3; axisCol++) {
                    let row = i * 3 + axisRow;
                    let col = j * 3 + axisCol;
                    matrixReal[row][col] = blockReal[axisRow][axisCol];
                    matrixImag[row][col] = blockImag[axisRow][axisCol];
                }
            }
        }
    }

    for (let i = 0; i < natoms * 3; i++) {
        for (let j = i + 1; j < natoms * 3; j++) {
            let symReal = 0.5 * (matrixReal[i][j] + matrixReal[j][i]);
            let symImag = 0.5 * (matrixImag[i][j] - matrixImag[j][i]);
            matrixReal[i][j] = symReal;
            matrixReal[j][i] = symReal;
            matrixImag[i][j] = symImag;
            matrixImag[j][i] = -symImag;
        }
        matrixImag[i][i] = 0;
    }

    return { real: matrixReal, imag: matrixImag };
}

function eigenvalueToFrequencyCm1(value, payload = null) {
    if (!Number.isFinite(value)) {
        return NaN;
    }
    let conversionFactor = payload && Number.isFinite(Number(payload.frequency_conversion_factor))
        ? Number(payload.frequency_conversion_factor)
        : 1.0;
    let magnitude = Math.sqrt(Math.abs(value));
    let signed = value < 0 ? -magnitude : magnitude;
    return signed * conversionFactor * 33.35641;
}

export async function solveHermitianEigenSystem(payload, qpoint, qDirection = null) {
    let blocks = buildDynamicalMatrixBlocks(payload, qpoint, qDirection);
    let eigenSystem = await solveComplexHermitianWithEigenWasm(blocks.real, blocks.imag);
    return {
        eigenvectors: eigenSystem.vectors,
        eigenvaluesRaw: eigenSystem.values,
        eigenvaluesCm1: eigenSystem.values.map((value) => eigenvalueToFrequencyCm1(value, payload)),
    };
}
