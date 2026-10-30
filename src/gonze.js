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

function multiplyRowByMatrix(row, matrix) {
    return [
        row[0] * matrix[0][0] + row[1] * matrix[1][0] + row[2] * matrix[2][0],
        row[0] * matrix[0][1] + row[1] * matrix[1][1] + row[2] * matrix[2][1],
        row[0] * matrix[0][2] + row[1] * matrix[1][2] + row[2] * matrix[2][2],
    ];
}

function multiplyMatrixByVector(matrix, vector) {
    return [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ];
}

function dot3(left, right) {
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

function complexAddAssign(target, sourceReal, sourceImag) {
    target[0] += sourceReal;
    target[1] += sourceImag;
}

function complexScaleAssign(target, scale) {
    target[0] *= scale;
    target[1] *= scale;
}

function createComplexBlockTensor(natoms) {
    let tensor = [];
    for (let i = 0; i < natoms; i++) {
        let row = [];
        for (let j = 0; j < natoms; j++) {
            let block = [];
            for (let alpha = 0; alpha < 3; alpha++) {
                let components = [];
                for (let beta = 0; beta < 3; beta++) {
                    components.push([0, 0]);
                }
                block.push(components);
            }
            row.push(block);
        }
        tensor.push(row);
    }
    return tensor;
}

function getDielectricPart(q, dielectric) {
    let total = 0;
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            total += q[i] * dielectric[i][j] * q[j];
        }
    }
    return total;
}

function getKKTensor(gVector, qCart, qDirectionCart, dielectric, lambda, tolerance) {
    let qK = [
        gVector[0] + qCart[0],
        gVector[1] + qCart[1],
        gVector[2] + qCart[2],
    ];
    let norm = Math.sqrt(dot3(qK, qK));
    let kk = [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
    ];

    if (norm < tolerance) {
        if (!qDirectionCart) {
            return kk;
        }
        let dielectricPart = getDielectricPart(qDirectionCart, dielectric);
        if (!(Math.abs(dielectricPart) > 1e-16)) {
            return kk;
        }
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                kk[i][j] = qDirectionCart[i] * qDirectionCart[j] / dielectricPart;
            }
        }
        return kk;
    }

    let dielectricPart = getDielectricPart(qK, dielectric);
    if (!(Math.abs(dielectricPart) > 1e-16)) {
        return kk;
    }
    let l2 = 4 * lambda * lambda;
    let prefactor = Math.exp(-dielectricPart / l2) / dielectricPart;
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            kk[i][j] = qK[i] * qK[j] * prefactor;
        }
    }
    return kk;
}

function getDdTensorPart(gList, qCart, qDirectionCart, dielectric, positionsCar, lambda, tolerance) {
    let natoms = positionsCar.length;
    let ddPart = createComplexBlockTensor(natoms);

    for (let gIndex = 0; gIndex < gList.length; gIndex++) {
        let gVector = gList[gIndex];
        let kk = getKKTensor(gVector, qCart, qDirectionCart, dielectric, lambda, tolerance);

        for (let i = 0; i < natoms; i++) {
            for (let j = 0; j < natoms; j++) {
                let phase = 0;
                for (let axis = 0; axis < 3; axis++) {
                    phase += (positionsCar[i][axis] - positionsCar[j][axis]) * gVector[axis];
                }
                phase *= 2 * Math.PI;
                let cosPhase = Math.cos(phase);
                let sinPhase = Math.sin(phase);

                for (let alpha = 0; alpha < 3; alpha++) {
                    for (let beta = 0; beta < 3; beta++) {
                        let value = kk[alpha][beta];
                        complexAddAssign(ddPart[i][j][alpha][beta], value * cosPhase, value * sinPhase);
                    }
                }
            }
        }
    }

    return ddPart;
}

function multiplyBorns(ddIn, born) {
    let natoms = born.length;
    let dd = createComplexBlockTensor(natoms);

    for (let i = 0; i < natoms; i++) {
        for (let j = 0; j < natoms; j++) {
            for (let alpha = 0; alpha < 3; alpha++) {
                for (let beta = 0; beta < 3; beta++) {
                    let target = dd[i][j][alpha][beta];
                    for (let alphaPrime = 0; alphaPrime < 3; alphaPrime++) {
                        for (let betaPrime = 0; betaPrime < 3; betaPrime++) {
                            let zz = born[i][alphaPrime][alpha] * born[j][betaPrime][beta];
                            target[0] += ddIn[i][j][alphaPrime][betaPrime][0] * zz;
                            target[1] += ddIn[i][j][alphaPrime][betaPrime][1] * zz;
                        }
                    }
                }
            }
        }
    }

    return dd;
}

function getQCart(direction, primitiveLattice) {
    let inverseLattice = invert3x3(primitiveLattice);
    if (!inverseLattice) {
        return null;
    }
    return multiplyMatrixByVector(inverseLattice, direction);
}

export function getGonzeReciprocalCorrection(payload, qpoint, qDirection = null) {
    if (!payload.nac || payload.nac.method !== 'gonze') {
        return null;
    }

    let qCart = getQCart(qpoint, payload.primitive_lattice);
    if (!qCart) {
        return null;
    }
    let qDirectionCart = qDirection ? getQCart(qDirection, payload.primitive_lattice) : null;
    let ddPart = getDdTensorPart(
        payload.nac.g_list,
        qCart,
        qDirectionCart,
        payload.nac.dielectric,
        payload.nac.positions_car,
        payload.nac.lambda,
        payload.nac.q_direction_tolerance || 1e-5
    );
    let dd = multiplyBorns(ddPart, payload.nac.born);
    let natoms = payload.masses.length;

    for (let i = 0; i < natoms; i++) {
        for (let alpha = 0; alpha < 3; alpha++) {
            for (let beta = 0; beta < 3; beta++) {
                dd[i][i][alpha][beta][0] -= payload.nac.dd_q0.real[i][alpha][beta];
                dd[i][i][alpha][beta][1] -= payload.nac.dd_q0.imag[i][alpha][beta];
            }
        }
    }

    for (let i = 0; i < natoms; i++) {
        for (let j = 0; j < natoms; j++) {
            let scale = payload.nac.nac_factor / Math.sqrt(payload.masses[i] * payload.masses[j]);
            for (let alpha = 0; alpha < 3; alpha++) {
                for (let beta = 0; beta < 3; beta++) {
                    complexScaleAssign(dd[i][j][alpha][beta], scale);
                }
            }
        }
    }

    return dd;
}
