import * as mat from './mat.js';
import * as utils from './utils.js';
import { atomic_symbol } from './atomic_data.js';

const elementToAtomicNumber = new Map();
for (let i = 1; i < atomic_symbol.length; i++) {
    elementToAtomicNumber.set(String(atomic_symbol[i]).toLowerCase(), i);
}

function tokenize(line) {
    return String(line || '').trim().split(/\s+/).filter(Boolean);
}

function parseTriple(line, label) {
    const values = tokenize(line).slice(0, 3).map(Number);
    if (values.length !== 3 || values.some((value) => !Number.isFinite(value))) {
        throw new Error(`Invalid ${label} in POSCAR`);
    }
    return values;
}

function isIntegerTokenList(tokens) {
    return tokens.length > 0 && tokens.every((token) => /^[+-]?\d+$/.test(token));
}

function getScaleFactor(rawScale, lattice) {
    if (!Number.isFinite(rawScale) || rawScale === 0) {
        throw new Error('Invalid POSCAR scale factor');
    }
    if (rawScale > 0) {
        return rawScale;
    }
    const volume = Math.abs(mat.matrix_determinant(lattice));
    if (!Number.isFinite(volume) || volume <= 0) {
        throw new Error('Invalid POSCAR lattice');
    }
    return Math.cbrt(Math.abs(rawScale) / volume);
}

function normalizeElementSymbol(token, index) {
    const text = String(token || '').trim();
    const match = text.match(/[A-Z][a-z]?/) || text.match(/[a-z]+/);
    if (!match) {
        return `Type${index + 1}`;
    }
    const value = match[0];
    return value.charAt(0).toUpperCase() + value.slice(1).toLowerCase();
}

function expandTypeLabels(speciesLabels, counts) {
    const expanded = [];
    for (let i = 0; i < counts.length; i++) {
        const label = speciesLabels[i] || `Type${i + 1}`;
        for (let j = 0; j < counts[i]; j++) {
            expanded.push(label);
        }
    }
    return expanded;
}

function scaleLattice(lattice, factor) {
    return lattice.map((vector) => vector.map((value) => value * factor));
}

function cartesianToReduced(position, inverseLattice) {
    return [
        position[0] * inverseLattice[0][0] + position[1] * inverseLattice[1][0] + position[2] * inverseLattice[2][0],
        position[0] * inverseLattice[0][1] + position[1] * inverseLattice[1][1] + position[2] * inverseLattice[2][1],
        position[0] * inverseLattice[0][2] + position[1] * inverseLattice[1][2] + position[2] * inverseLattice[2][2],
    ];
}

export function parsePoscarStructureLines(lines) {
    if (!Array.isArray(lines) || lines.length < 8) {
        throw new Error('POSCAR is too short');
    }

    let cursor = 0;
    const title = String(lines[cursor++] || '').trim() || 'Structure';
    const rawScale = Number(String(lines[cursor++] || '').trim());
    const latticeRaw = [
        parseTriple(lines[cursor++], 'lattice vector a'),
        parseTriple(lines[cursor++], 'lattice vector b'),
        parseTriple(lines[cursor++], 'lattice vector c'),
    ];

    let speciesLabels = [];
    let counts = [];
    const speciesTokens = tokenize(lines[cursor]);
    if (isIntegerTokenList(speciesTokens)) {
        counts = speciesTokens.map((token) => parseInt(token, 10));
        speciesLabels = counts.map((value, index) => `Type${index + 1}`);
        cursor += 1;
    } else {
        speciesLabels = speciesTokens.map((token, index) => normalizeElementSymbol(token, index));
        cursor += 1;
        const countTokens = tokenize(lines[cursor]);
        if (!isIntegerTokenList(countTokens)) {
            throw new Error('Missing POSCAR atom counts');
        }
        counts = countTokens.map((token) => parseInt(token, 10));
        cursor += 1;
    }

    if (counts.some((value) => !Number.isFinite(value) || value <= 0)) {
        throw new Error('Invalid POSCAR atom counts');
    }

    let selectiveDynamics = false;
    const coordinateLine = String(lines[cursor] || '').trim();
    if (/^[sS]/.test(coordinateLine)) {
        selectiveDynamics = true;
        cursor += 1;
    }

    const coordinateToken = String(lines[cursor] || '').trim();
    if (!coordinateToken) {
        throw new Error('Missing POSCAR coordinate system');
    }
    const isDirect = /^[dD]/.test(coordinateToken);
    const isCartesian = /^[cCkK]/.test(coordinateToken);
    if (!isDirect && !isCartesian) {
        throw new Error('Unsupported POSCAR coordinate system');
    }
    cursor += 1;

    const natoms = counts.reduce((sum, value) => sum + value, 0);
    const coordinates = [];
    for (let i = 0; i < natoms; i++) {
        const tokens = tokenize(lines[cursor + i]);
        if (tokens.length < 3) {
            throw new Error(`Invalid POSCAR atomic position on line ${cursor + i + 1}`);
        }
        coordinates.push(tokens.slice(0, 3).map(Number));
    }
    cursor += natoms;

    const scaleFactor = getScaleFactor(rawScale, latticeRaw);
    const lat = scaleLattice(latticeRaw, scaleFactor);
    let atom_pos_red;
    let atom_pos_car;

    if (isDirect) {
        atom_pos_red = coordinates;
        atom_pos_car = utils.red_car_list(atom_pos_red, lat);
    } else {
        atom_pos_car = coordinates.map((position) => position.map((value) => value * scaleFactor));
        const inverseLattice = mat.matrix_inverse(lat);
        atom_pos_red = atom_pos_car.map((position) => cartesianToReduced(position, inverseLattice));
    }

    const atom_types = expandTypeLabels(speciesLabels, counts);
    const uniqueAtomNumbers = speciesLabels.map((label) => elementToAtomicNumber.get(String(label).toLowerCase()) || 0);
    const atomTypeIndexPerSite = [];
    const atoms = [];
    let expandedIndex = 0;
    for (let typeIndex = 0; typeIndex < counts.length; typeIndex++) {
        for (let countIndex = 0; countIndex < counts[typeIndex]; countIndex++) {
            atomTypeIndexPerSite.push(typeIndex);
            const position = atom_pos_car[expandedIndex];
            atoms.push([typeIndex, position[0], position[1], position[2]]);
            expandedIndex += 1;
        }
    }

    const formula = utils.get_formula(atom_types);
    const repetitions = utils.getReasonableRepetitions(natoms, lat);

    return {
        title,
        scaleFactor,
        lat,
        natoms,
        counts,
        speciesLabels,
        atomTypeIndexPerSite,
        atom_types,
        atom_numbers: uniqueAtomNumbers,
        atom_pos_red,
        atom_pos_car,
        atoms,
        formula,
        name: formula || title,
        repetitions,
        selectiveDynamics,
        nextLineIndex: cursor,
    };
}

export function parsePoscar(text) {
    const normalized = String(text || '').replace(/^\uFEFF/, '');
    const lines = normalized.split(/\r?\n/);
    const structure = parsePoscarStructureLines(lines);

    return {
        type: 'poscar',
        title: structure.title,
        name: structure.name,
        formula: structure.formula,
        lat: structure.lat,
        natoms: structure.natoms,
        atom_types: structure.atom_types,
        atom_numbers: structure.atom_numbers,
        atom_pos_red: structure.atom_pos_red,
        atom_pos_car: structure.atom_pos_car,
        atoms: structure.atoms,
        repetitions: structure.repetitions,
    };
}
