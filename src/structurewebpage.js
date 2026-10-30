import { parsePoscar } from './poscarformat.js';
import { parseChgcar } from './chgcar.js';
import { format_formula_html } from './utils.js';
import { renderLatticeTable, renderAtomPositionsTable } from './structureinfo.js';

function getFilenameStem(filename) {
    return String(filename || '').replace(/\.[^.]+$/u, '');
}

function parseStructureFile(filename, text) {
    const lowerName = String(filename || '').toLowerCase();
    if (lowerName.includes('chgcar')) {
        return parseChgcar(text);
    }
    if (lowerName.includes('poscar') || lowerName.includes('contcar') || lowerName.endsWith('.vasp')) {
        return parsePoscar(text);
    }

    try {
        return parseChgcar(text);
    } catch (chgcarError) {
        return parsePoscar(text);
    }
}

function quantile(sortedValues, fraction) {
    if (!sortedValues.length) {
        return 0;
    }
    if (sortedValues.length === 1) {
        return sortedValues[0];
    }
    const clamped = Math.max(0, Math.min(1, fraction));
    const position = clamped * (sortedValues.length - 1);
    const lower = Math.floor(position);
    const upper = Math.ceil(position);
    const weight = position - lower;
    if (lower === upper) {
        return sortedValues[lower];
    }
    return sortedValues[lower] * (1 - weight) + sortedValues[upper] * weight;
}

function formatIsolevelValue(value) {
    const numeric = Number(value);
    if (!Number.isFinite(numeric)) {
        return '';
    }
    if (numeric === 0) {
        return '0';
    }
    if (Math.abs(numeric) >= 1) {
        return numeric.toFixed(3).replace(/\.?0+$/, '');
    }
    return numeric.toPrecision(4).replace(/\.?0+$/, '');
}

export class StructureWebpage {

    constructor(viewer) {
        this.viewer = viewer;
        this.currentData = null;
    }

    setTitle(domTitle) {
        this.domTitle = domTitle;
    }

    setFileInput(domInput) {
        this.domFileInput = domInput;
        domInput.on('change', this.loadCustomFile.bind(this));
        domInput.on('click', function() { this.value = ''; });
    }

    setLattice(domLattice) {
        this.domLattice = domLattice;
    }

    setAtomPositions(domAtompos) {
        this.domAtompos = domAtompos;
    }

    setIsosurfaceSection(domSection) {
        this.domIsosurfaceSection = domSection;
    }

    setIsosurfaceModeInput(domInput) {
        this.domIsosurfaceModeInput = domInput;
        if (!domInput || !domInput.length) {
            return;
        }
        domInput.val(this.viewer.getIsosurfaceRenderMode());
        domInput.on('change', () => {
            this.viewer.setIsosurfaceRenderMode(domInput.val());
        });
    }

    setIsolevelInput(domInput, domValue = null) {
        this.domIsolevelInput = domInput;
        this.domIsolevelValue = domValue;

        const updateLabel = () => {
            const value = Number(domInput.val());
            if (this.domIsolevelValue) {
                this.domIsolevelValue.text(formatIsolevelValue(value));
            }
            return value;
        };

        const previewHandler = () => {
            const value = updateLabel();
            this.viewer.previewIsolevel(value);
        };

        const finalHandler = () => {
            const value = updateLabel();
            this.viewer.changeIsolevel(value);
        };

        domInput.on('input', previewHandler);
        domInput.on('change', finalHandler);
        updateLabel();
    }

    setIsosurfaceOpacityInput(domInput, domValue = null) {
        this.domOpacityInput = domInput;
        this.domOpacityValue = domValue;

        const handler = () => {
            const value = Number(domInput.val());
            if (this.domOpacityValue) {
                this.domOpacityValue.text(formatIsolevelValue(value));
            }
            this.viewer.changeIsosurfaceOpacity(value);
        };

        domInput.on('input change', handler);
        handler();
    }

    setRepetitionControls(domNx, domNy, domNz, domButton) {
        this.domNx = domNx;
        this.domNy = domNy;
        this.domNz = domNz;
        this.domRepeatButton = domButton;

        const apply = () => {
            this.viewer.setRepetitions(domNx.val(), domNy.val(), domNz.val());
        };

        domButton.on('click', apply);
        [domNx, domNy, domNz].forEach((domInput) => {
            domInput.on('keydown', (event) => {
                if (event.key === 'Enter') {
                    event.preventDefault();
                    apply();
                }
            });
        });
    }

    setCameraDirectionButton(domButton, direction) {
        domButton.on('click', () => this.viewer.setCameraDirection(direction));
    }

    init() {
        this.setTitleText('Structure');
        this.setChargeDensityVisibility(false);
    }

    loadCustomFile(event) {
        const file = event.target.files[0];
        if (!file) {
            return;
        }

        const reader = new FileReader();
        reader.readAsText(file);
        reader.onloadend = () => {
            try {
                const data = parseStructureFile(file.name, reader.result);
                if (data.values) {
                    this.configureChargeDensityRange(data);
                }
                this.currentData = data;
                this.viewer.setData(data);
                this.applyDefaultRepetitions(data.repetitions || [1, 1, 1]);
                this.setChargeDensityVisibility(!!data.values);
                if (this.domIsolevelInput && data.values) {
                    this.domIsolevelInput.val(this.viewer.getIsolevel());
                    if (this.domIsolevelValue) {
                        this.domIsolevelValue.text(formatIsolevelValue(this.viewer.getIsolevel()));
                    }
                }
                this.setTitleText(data.formula || getFilenameStem(file.name));
                this.updateStructureInfo(data);
            } catch (error) {
                const message = error && error.message ? error.message : 'Unable to read this structure file.';
                if (typeof alert === 'function') {
                    alert(message);
                } else {
                    console.error(message);
                }
            }
        };
    }

    applyDefaultRepetitions(repetitions) {
        const nx = repetitions[0] || 1;
        const ny = repetitions[1] || 1;
        const nz = repetitions[2] || 1;

        if (this.domNx) { this.domNx.val(nx); }
        if (this.domNy) { this.domNy.val(ny); }
        if (this.domNz) { this.domNz.val(nz); }

        this.viewer.setRepetitions(nx, ny, nz);
    }

    setChargeDensityVisibility(isVisible) {
        if (this.domIsosurfaceSection) {
            this.domIsosurfaceSection.toggle(!!isVisible);
        }
    }

    configureChargeDensityRange(data) {
        if (!this.domIsolevelInput || !this.domIsolevelInput.length || !Array.isArray(data.values) || !data.values.length) {
            return;
        }

        const positiveValues = data.values.filter((value) => Number.isFinite(value) && value > 0).sort((a, b) => a - b);
        if (!positiveValues.length) {
            return;
        }

        const minValue = 0;
        const maxValue = Math.max(
            quantile(positiveValues, 0.999),
            positiveValues[positiveValues.length - 1] * 0.02
        );
        const defaultValue = Math.min(
            Math.max(quantile(positiveValues, 0.95), maxValue * 0.02),
            maxValue * 0.9
        );
        const stepValue = Math.max(maxValue / 500, 1e-4);

        data.isolevel = defaultValue;
        this.domIsolevelInput.attr('min', minValue);
        this.domIsolevelInput.attr('max', maxValue);
        this.domIsolevelInput.attr('step', stepValue);
        this.domIsolevelInput.attr('value', defaultValue);
        if (this.domIsolevelValue) {
            this.domIsolevelValue.text(formatIsolevelValue(defaultValue));
        }
    }

    setTitleText(title) {
        if (!this.domTitle) {
            return;
        }
        this.domTitle.html(format_formula_html(title));
    }

    updateStructureInfo(data) {
        if (!data) {
            return;
        }

        renderLatticeTable(this.domLattice, data.lat);
        renderAtomPositionsTable(
            this.domAtompos,
            data.atom_pos_red,
            data.atom_types,
            data.atoms.map((atom) => data.atom_numbers[atom[0]]),
            this.viewer.getAtomColorHex.bind(this.viewer)
        );
    }
}
