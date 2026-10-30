export function atomColorHexToCss(colorHex) {
    return '#' + Number(colorHex).toString(16).padStart(6, '0');
}

export function getAtomBadgeTextColor(colorHex) {
    let color = Number(colorHex);
    let red = (color >> 16) & 255;
    let green = (color >> 8) & 255;
    let blue = color & 255;
    let luminance = (0.2126 * red + 0.7152 * green + 0.0722 * blue) / 255;
    return luminance > 0.6 ? '#111827' : '#ffffff';
}

export function applyAtomBadgeStyle(element, atomNumber, getAtomColorHex) {
    if (!element || !Number.isFinite(atomNumber) || typeof getAtomColorHex !== 'function') {
        return;
    }
    let colorHex = getAtomColorHex(atomNumber);
    element.style.setProperty('--atom-badge-bg', atomColorHexToCss(colorHex));
    element.style.setProperty('--atom-badge-fg', getAtomBadgeTextColor(colorHex));
}

export function createAtomBadgeHtml(label, atomNumber, getAtomColorHex) {
    let style = '';
    if (Number.isFinite(atomNumber) && typeof getAtomColorHex === 'function') {
        let colorHex = getAtomColorHex(atomNumber);
        style =
            ' style="--atom-badge-bg: ' + atomColorHexToCss(colorHex) +
            '; --atom-badge-fg: ' + getAtomBadgeTextColor(colorHex) + ';"';
    }
    return '<span class="atom-type-badge"' + style + '>' + label + '</span>';
}
