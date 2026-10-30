import * as THREE from 'three';

import { StructureViewerBase } from './structureviewerbase.js';

export class ExcitonWf extends StructureViewerBase {

    constructor() {
        super();
        this.excitonIndex = 0;
        this.excitons = null;
    }

    setData(absorption) {
        this.excitons = absorption.excitons;
        this.excitonIndex = absorption.excitonIndex;
        this.values = absorption.excitons[this.excitonIndex].datagrid;
        this.sizex = absorption.sizex;
        this.sizey = absorption.sizey;
        this.sizez = absorption.sizez;
        this.cell = absorption.cell;
        this.gridCell = absorption.gridCell || absorption.cell;
        this.atoms = absorption.atoms;
        this.natoms = absorption.atoms.length;
        this.atom_numbers = absorption.atom_numbers;
        this.clearIsosurfacePreviewCache();
        this.updateRecommendedIsolevel();

        this.geometricCenter = new THREE.Vector3(0, 0, 0);
        for (let i = 0; i < this.atoms.length; i++) {
            const pos = new THREE.Vector3(this.atoms[i][1], this.atoms[i][2], this.atoms[i][3]);
            this.geometricCenter.add(pos);
        }
        this.geometricCenter.multiplyScalar(1.0 / this.atoms.length);

        this.updateViewParameters();
        this.initializeBondRulesFromAtoms();
        this.refreshAppearanceControls();
    }

    setExcitonIndex(index) {
        this.excitonIndex = index;
        if (this.excitons && this.excitons[index]) {
            this.values = this.excitons[index].datagrid;
            this.clearIsosurfacePreviewCache();
            this.updateRecommendedIsolevel();
        }
    }

    updateRecommendedIsolevel(force = false) {
        if (!Array.isArray(this.values) || !this.values.length) {
            return;
        }

        const positiveValues = this.values.filter((value) => Number.isFinite(value) && value > 0).sort((a, b) => a - b);
        if (!positiveValues.length) {
            return;
        }

        const quantile = (fraction) => {
            const position = Math.max(0, Math.min(1, fraction)) * (positiveValues.length - 1);
            const lower = Math.floor(position);
            const upper = Math.ceil(position);
            const weight = position - lower;
            if (lower === upper) {
                return positiveValues[lower];
            }
            return positiveValues[lower] * (1 - weight) + positiveValues[upper] * weight;
        };

        const recommended = Math.min(Math.max(quantile(0.99), positiveValues[positiveValues.length - 1] * 0.02), 0.9);
        if (
            force ||
            !this.hasManualIsolevel() ||
            !Number.isFinite(this.getIsolevel()) ||
            this.getIsolevel() <= 0
        ) {
            this.resetIsolevelState(recommended);
        }

        if (this.onIsolevelRangeChanged) {
            this.onIsolevelRangeChanged({
                min: 0,
                max: Math.max(quantile(0.999), positiveValues[positiveValues.length - 1] * 0.1),
                step: Math.max(1e-4, positiveValues[positiveValues.length - 1] / 500),
                value: this.getIsolevel(),
            });
        }
    }

    updateStructure() {
        if (!this.scene || !this.excitons || !this.excitons.length) {
            return;
        }

        this.values = this.excitons[this.excitonIndex].datagrid;
        this.removeStructure();

        this.addLights();
        this.addMarchingCubes();
        this.getAtomMaterials();
        this.addStructure();
        this.render();
    }
}
