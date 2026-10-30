import $ from 'jquery';

import { format_formula_html } from './utils.js';
import { renderLatticeTable } from './structureinfo.js';

export class ExcitonWebpage {

    constructor(viewer, spectra) {
        this.viewer = viewer;
        this.spectra = spectra;
        this.materials = {
            bn: {
                name: 'Boron Nitride',
                file: 'data/excitondb/bn/absorptionspectra.json',
                visible: true,
            },
            mos2: {
                name: 'MoS2',
                file: 'data/excitondb/mos2/absorptionspectra.json',
                visible: false,
            },
            mote2: {
                name: 'MoTe2',
                file: 'data/excitondb/mote2/absorptionspectra.json',
                visible: false,
            },
        };
        this.currentMaterial = 'bn';

        this.spectra.setSelectionHandler((index) => {
            this.selectExciton(index);
        });
        this.viewer.onIsolevelRangeChanged = this.updateIsolevelControls.bind(this);
    }

    setTitle(domTitle) {
        this.domTitle = domTitle;
    }

    setMaterialsList(domList) {
        this.domMaterials = domList;
    }

    setLattice(domLattice) {
        this.domLattice = domLattice;
    }

    setFileInput(domInput) {
        this.domFileInput = domInput;
        domInput.change(this.loadCustomFile.bind(this));
        domInput.click(function () { this.value = ''; });
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
                this.domIsolevelValue.text(value.toFixed(3));
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
                this.domOpacityValue.text(value.toFixed(2).replace(/\.?0+$/, ''));
            }
            this.viewer.changeIsosurfaceOpacity(value);
        };

        domInput.on('input change', handler);
        handler();
    }

    updateIsolevelControls(range) {
        if (!this.domIsolevelInput || !this.domIsolevelInput.length || !range) {
            return;
        }

        this.domIsolevelInput.attr('min', range.min);
        this.domIsolevelInput.attr('max', range.max);
        this.domIsolevelInput.attr('step', range.step);
        this.domIsolevelInput.val(range.value);
        if (this.domIsolevelValue) {
            this.domIsolevelValue.text(Number(range.value).toPrecision(4).replace(/\.?0+$/, ''));
        }
    }

    setCameraDirectionButton(domButton, direction) {
        domButton.click(() => this.viewer.setCameraDirection(direction));
    }

    init() {
        this.updateMaterialsMenu();
        return this.loadMaterial(this.currentMaterial);
    }

    updateMaterialsMenu() {
        if (!this.domMaterials) {
            return;
        }

        this.domMaterials.empty();

        Object.keys(this.materials).forEach((key) => {
            const material = this.materials[key];
            if (material.visible === false) {
                return;
            }
            const link = $('<a href="#"></a>');
            link.html(format_formula_html(material.name));
            link.on('click', (event) => {
                event.preventDefault();
                this.loadMaterial(key);
            });

            const item = $('<li></li>');
            if (key === this.currentMaterial) {
                item.addClass('active');
            }
            item.append(link);
            this.domMaterials.append(item);
        });
    }

    loadMaterial(key) {
        const material = this.materials[key];
        if (!material) {
            return Promise.resolve();
        }

        this.currentMaterial = key;
        this.updateMaterialsMenu();
        this.setTitleText(material.name);

        return this.spectra.getDataFilename(material.file).then(() => {
            this.syncViewerFromSpectra();
        });
    }

    loadCustomFile(event) {
        const file = event.target.files[0];
        if (!file) {
            return;
        }

        const reader = new FileReader();
        reader.readAsText(file);
        reader.onloadend = () => {
            const data = JSON.parse(reader.result);
            this.currentMaterial = '';
            this.updateMaterialsMenu();
            this.setTitleText(file.name.replace(/\.json$/i, ''));
            this.spectra.getDataObject(data);
            this.spectra.render();
            this.syncViewerFromSpectra();
        };
    }

    selectExciton(index) {
        this.spectra.setExcitonIndex(index, false);
        this.viewer.setExcitonIndex(index);
        this.viewer.updateStructure();
    }

    setTitleText(title) {
        if (this.domTitle) {
            this.domTitle.html(format_formula_html(title));
        }
    }

    updateStructureInfo() {
        if (!this.spectra) {
            return;
        }

        const lattice = this.spectra.gridCell || this.spectra.cell;
        renderLatticeTable(this.domLattice, lattice);
    }

    syncViewerFromSpectra() {
        if (!this.viewer || !this.spectra) {
            return;
        }

        this.viewer.setData(this.spectra);
        this.updateStructureInfo();
        this.viewer.updateStructure();
    }
}
