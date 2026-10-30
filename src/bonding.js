import * as THREE from 'three';

import { getCombinations } from './utils.js';

export function getCovalentBondLength(atomNumberA, atomNumberB, covalentRadii) {
    return (covalentRadii[atomNumberA] || 0) + (covalentRadii[atomNumberB] || 0);
}

export function getChemicalBondLimit(atomNumberA, atomNumberB, covalentRadii) {
    const covalent = getCovalentBondLength(atomNumberA, atomNumberB, covalentRadii);
    return Math.max(0.4, covalent + 0.45);
}

export function getBondSearchLimit(atomNumberA, atomNumberB, covalentRadii) {
    const chemical = getChemicalBondLimit(atomNumberA, atomNumberB, covalentRadii);
    return chemical + Math.max(0.5, chemical * 0.2);
}

export function getEffectiveCoordinationCandidates(siteNeighbors) {
    if (!siteNeighbors || !siteNeighbors.length) {
        return [];
    }

    const sorted = siteNeighbors
        .filter((neighbor) => Number.isFinite(neighbor.distance) && neighbor.distance >= 0.4)
        .sort((a, b) => a.distance - b.distance);
    if (!sorted.length) {
        return [];
    }

    const shortest = sorted[0].distance;
    const accepted = [];

    for (let i = 0; i < sorted.length; i++) {
        const neighbor = sorted[i];
        if (neighbor.distance > shortest + Math.max(1.0, shortest * 0.6)) {
            break;
        }

        const relative = neighbor.distance / shortest;
        const weight = Math.exp(1.0 - Math.pow(relative, 6));
        if (weight < 0.35) {
            continue;
        }

        accepted.push({
            index: neighbor.index,
            atom_number: neighbor.atom_number,
            distance: neighbor.distance,
            weight,
        });
    }

    return accepted;
}

export function buildCrystalBondRules(atoms, atomNumbers, covalentRadii, getBondRuleKey) {
    const rules = {};
    if (!atoms || !atoms.length || !atomNumbers || !atomNumbers.length) {
        return rules;
    }

    const probeAtoms = atoms.map((atom, index) => ({
        index,
        atom_number: atomNumbers[atom[0]],
        position: new THREE.Vector3(atom[1], atom[2], atom[3]),
    }));

    const pairDistances = {};
    const siteNeighbors = probeAtoms.map(() => []);
    const combinations = getCombinations(probeAtoms);

    for (let i = 0; i < combinations.length; i++) {
        const a = combinations[i][0];
        const b = combinations[i][1];
        const length = a.position.distanceTo(b.position);
        const searchLimit = getBondSearchLimit(a.atom_number, b.atom_number, covalentRadii);
        if (length <= searchLimit && length >= 0.4) {
            siteNeighbors[a.index].push({
                index: b.index,
                atom_number: b.atom_number,
                distance: length,
            });
            siteNeighbors[b.index].push({
                index: a.index,
                atom_number: a.atom_number,
                distance: length,
            });
        }
    }

    const acceptedNeighbors = siteNeighbors.map((neighbors) => getEffectiveCoordinationCandidates(neighbors));
    const acceptedNeighborMaps = acceptedNeighbors.map((neighbors) => {
        const lookup = {};
        for (let i = 0; i < neighbors.length; i++) {
            lookup[neighbors[i].index] = neighbors[i];
        }
        return lookup;
    });

    for (let i = 0; i < combinations.length; i++) {
        const a = combinations[i][0];
        const b = combinations[i][1];
        const length = a.position.distanceTo(b.position);
        const chemicalLimit = getChemicalBondLimit(a.atom_number, b.atom_number, covalentRadii);
        const acceptedByA = acceptedNeighborMaps[a.index][b.index];
        const acceptedByB = acceptedNeighborMaps[b.index][a.index];

        if (!acceptedByA || !acceptedByB || length > chemicalLimit) {
            continue;
        }

        const key = getBondRuleKey(a.atom_number, b.atom_number);
        if (!pairDistances[key]) {
            pairDistances[key] = {
                a: Math.min(a.atom_number, b.atom_number),
                b: Math.max(a.atom_number, b.atom_number),
                distances: [],
            };
        }
        pairDistances[key].distances.push(length);
    }

    const keys = Object.keys(pairDistances);
    for (let i = 0; i < keys.length; i++) {
        const pair = pairDistances[keys[i]];
        if (!pair.distances.length) {
            continue;
        }

        let cutoff = Math.max.apply(null, pair.distances) + 0.04;
        cutoff = Math.min(cutoff, getChemicalBondLimit(pair.a, pair.b, covalentRadii));
        if (Number.isFinite(cutoff) && cutoff >= 0.4) {
            rules[keys[i]] = {
                a: pair.a,
                b: pair.b,
                cutoff,
            };
        }
    }

    return rules;
}
