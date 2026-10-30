export function getViewerAtomRadius(display, atomNumber, atomScale, sphereRadius, covalentRadii) {
    if (display === 'vesta') {
        return ((covalentRadii[atomNumber] || 0) / 2.3) * atomScale;
    }
    return sphereRadius * atomScale;
}

export function buildBondList(atomobjects, bondRules, getBondRuleKey, getDefaultBondCutoff, options = {}) {
    const bonds = [];
    const requireRule = !!options.requireRule;

    for (let i = 0; i < atomobjects.length; i++) {
        const atomA = atomobjects[i];
        for (let j = i + 1; j < atomobjects.length; j++) {
            const atomB = atomobjects[j];
            const distance = atomA.position.distanceTo(atomB.position);
            const key = getBondRuleKey(atomA.atom_number, atomB.atom_number);
            const rule = bondRules[key];
            if (requireRule && !rule) {
                continue;
            }
            const cutoff = rule ? rule.cutoff : getDefaultBondCutoff(atomA.atom_number, atomB.atom_number);
            if (distance < cutoff) {
                bonds.push({
                    a: atomA.position,
                    b: atomB.position,
                    a_atom_number: atomA.atom_number,
                    b_atom_number: atomB.atom_number,
                    baseLength: distance,
                });
            }
        }
    }

    return bonds;
}
