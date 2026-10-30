import { atomic_symbol } from './atomic_data.js';
import { applyAtomBadgeStyle } from './atomcolors.js';
import * as mat from './mat.js';

function clearDom(domNode) {
    if (domNode && domNode.length) {
        domNode.empty();
    }
}

export function renderLatticeTable(domLattice, lattice) {
    if (!domLattice || !domLattice.length || !lattice) {
        return;
    }

    clearDom(domLattice);
    for (let i = 0; i < 3; i++) {
        const tr = document.createElement('tr');
        for (let j = 0; j < 3; j++) {
            const td = document.createElement('td');
            td.appendChild(document.createTextNode(Number(lattice[i][j]).toPrecision(4)));
            tr.append(td);
        }
        domLattice.append(tr);
    }
}

export function renderAtomPositionsTable(domAtompos, positionsReduced, siteLabels, siteAtomNumbers, getAtomColorHex) {
    if (!domAtompos || !domAtompos.length || !positionsReduced || !siteLabels || !siteAtomNumbers) {
        return;
    }

    clearDom(domAtompos);
    for (let i = 0; i < positionsReduced.length; i++) {
        const tr = document.createElement('tr');

        const typeCell = document.createElement('td');
        const badge = document.createElement('span');
        badge.className = 'atom-type-badge';
        badge.textContent = siteLabels[i];
        if (Number.isFinite(siteAtomNumbers[i]) && typeof getAtomColorHex === 'function') {
            applyAtomBadgeStyle(badge, siteAtomNumbers[i], getAtomColorHex);
        }
        typeCell.className = 'ap atom-type-cell';
        typeCell.appendChild(badge);
        tr.append(typeCell);

        for (let j = 0; j < 3; j++) {
            const td = document.createElement('td');
            td.appendChild(document.createTextNode(Number(positionsReduced[i][j]).toFixed(4)));
            tr.append(td);
        }

        domAtompos.append(tr);
    }
}

export function cartesianPositionsToReduced(cartesianPositions, lattice) {
    if (!Array.isArray(cartesianPositions) || !Array.isArray(lattice)) {
        return [];
    }

    const inverseLattice = mat.matrix_inverse(lattice);
    return cartesianPositions.map((position) => [
        position[0] * inverseLattice[0][0] + position[1] * inverseLattice[1][0] + position[2] * inverseLattice[2][0],
        position[0] * inverseLattice[0][1] + position[1] * inverseLattice[1][1] + position[2] * inverseLattice[2][1],
        position[0] * inverseLattice[0][2] + position[1] * inverseLattice[1][2] + position[2] * inverseLattice[2][2],
    ]);
}

export function buildSiteLabelsFromAtoms(atoms, atomNumbers) {
    if (!Array.isArray(atoms) || !Array.isArray(atomNumbers)) {
        return { labels: [], numbers: [] };
    }

    const labels = [];
    const numbers = [];
    for (let i = 0; i < atoms.length; i++) {
        const atomTypeIndex = atoms[i][0];
        const atomNumber = atomNumbers[atomTypeIndex];
        numbers.push(atomNumber);
        labels.push(atomic_symbol[atomNumber] || `Type${atomTypeIndex + 1}`);
    }
    return { labels, numbers };
}
