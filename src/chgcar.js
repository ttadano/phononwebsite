import * as mat from './mat.js';
import { parsePoscarStructureLines } from './poscarformat.js';

function tokenize(line) {
    return String(line || '').trim().split(/\s+/).filter(Boolean);
}

function nextNonEmptyLine(lines, startIndex) {
    let index = startIndex;
    while (index < lines.length && !String(lines[index] || '').trim()) {
        index += 1;
    }
    return index;
}

export function parseChgcar(text) {
    const normalized = String(text || '').replace(/^\uFEFF/, '');
    const lines = normalized.split(/\r?\n/);
    const structure = parsePoscarStructureLines(lines);

    let cursor = nextNonEmptyLine(lines, structure.nextLineIndex);
    const dimensionTokens = tokenize(lines[cursor]);
    if (dimensionTokens.length < 3 || !dimensionTokens.slice(0, 3).every((token) => /^[+-]?\d+$/.test(token))) {
        throw new Error('Missing CHGCAR grid dimensions');
    }

    const sizex = parseInt(dimensionTokens[0], 10);
    const sizey = parseInt(dimensionTokens[1], 10);
    const sizez = parseInt(dimensionTokens[2], 10);
    if (![sizex, sizey, sizez].every((value) => Number.isFinite(value) && value > 0)) {
        throw new Error('Invalid CHGCAR grid dimensions');
    }

    cursor += 1;
    const totalValues = sizex * sizey * sizez;
    const values = [];

    while (cursor < lines.length && values.length < totalValues) {
        const tokens = tokenize(lines[cursor]);
        for (let i = 0; i < tokens.length && values.length < totalValues; i++) {
            const value = Number(tokens[i]);
            if (!Number.isFinite(value)) {
                throw new Error('Invalid CHGCAR charge density value');
            }
            values.push(value);
        }
        cursor += 1;
    }

    if (values.length !== totalValues) {
        throw new Error('CHGCAR ended before the full charge density grid was read');
    }

    const cellVolume = Math.abs(mat.matrix_determinant(structure.lat));
    const normalization = cellVolume;
    const normalizedValues = normalization > 0
        ? values.map((value) => value / normalization)
        : values;

    return {
        type: 'chgcar',
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
        values: normalizedValues,
        sizex,
        sizey,
        sizez,
        gridCell: structure.lat,
        isolevel: 0.02,
    };
}
