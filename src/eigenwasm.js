import createEigenModule from './eigen_solver.generated.mjs';

let modulePromise = null;
let moduleInstance = null;
let solveHermitianEigen = null;

export async function initEigenWasmSolver() {
    if (moduleInstance) {
        return moduleInstance;
    }
    if (!modulePromise) {
        modulePromise = Promise.resolve(createEigenModule())
            .then((module) => {
                moduleInstance = module;
                solveHermitianEigen = module.cwrap(
                    'solve_hermitian_eigen',
                    'number',
                    ['number', 'number', 'number', 'number']
                );
                return moduleInstance;
            })
            .catch((error) => {
                modulePromise = null;
                throw error;
            });
    }
    return modulePromise;
}

export async function solveComplexHermitianWithEigenWasm(matrixReal, matrixImag) {
    let module = await initEigenWasmSolver();
    let size = matrixReal.length;
    if (!size) {
        return {
            values: [],
            vectors: [],
        };
    }

    let matrixPtr = module._malloc(size * size * 16);
    let valuesPtr = module._malloc(size * 8);
    let vectorsPtr = module._malloc(size * size * 16);

    try {
        for (let col = 0; col < size; col++) {
            for (let row = 0; row < size; row++) {
                let offset = 16 * (col * size + row);
                module.setValue(matrixPtr + offset, matrixReal[row][col], 'double');
                module.setValue(matrixPtr + offset + 8, matrixImag[row][col], 'double');
            }
        }

        let status = solveHermitianEigen(size, matrixPtr, valuesPtr, vectorsPtr);
        if (status !== 0) {
            throw new Error(`Eigen wasm solver failed with status ${status}`);
        }

        let values = [];
        for (let modeIndex = 0; modeIndex < size; modeIndex++) {
            values.push(module.getValue(valuesPtr + modeIndex * 8, 'double'));
        }
        let vectors = [];
        for (let modeIndex = 0; modeIndex < size; modeIndex++) {
            let vector = [];
            for (let row = 0; row < size; row++) {
                let offset = 16 * (modeIndex * size + row);
                vector.push([
                    module.getValue(vectorsPtr + offset, 'double'),
                    module.getValue(vectorsPtr + offset + 8, 'double'),
                ]);
            }
            vectors.push(vector);
        }

        return {
            values: values,
            vectors: vectors,
        };
    } finally {
        module._free(vectorsPtr);
        module._free(valuesPtr);
        module._free(matrixPtr);
    }
}
