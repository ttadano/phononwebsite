import * as THREE from 'three';

import { StructureViewerBase } from './structureviewerbase.js';
import { covalent_radii } from './atomic_data.js';
import { buildCrystalBondRules, getChemicalBondLimit } from './bonding.js';
import { clampRaymarchStepCount, refreshRaymarchedIsosurfaceUniforms } from './raymarchedisosurface.js';
import { createCellLineObject } from './viewergeometry.js';

export class StructureViewer extends StructureViewerBase {

    constructor() {
        super();
        this.cell = true;
        this.nx = 1;
        this.ny = 1;
        this.nz = 1;
        this.structureData = null;
        this.baseAtoms = [];
        this.baseLattice = null;
        this.gridCell = null;
        this.isosurfaceValue = null;
        this.currentLattice = null;
        this.shouldAutoFitCamera = true;
    }

    setData(data) {
        this.structureData = data;
        this.baseAtoms = Array.isArray(data.atoms) ? data.atoms.map((atom) => atom.slice()) : [];
        this.baseLattice = Array.isArray(data.lat) ? data.lat.map((vector) => vector.slice()) : null;
        this.gridCell = Array.isArray(data.gridCell) ? data.gridCell.map((vector) => vector.slice()) : this.baseLattice;
        this.atom_numbers = Array.isArray(data.atom_numbers) ? data.atom_numbers.slice() : [];
        this.values = Array.isArray(data.values) ? data.values.slice() : null;
        this.sizex = data.sizex || 1;
        this.sizey = data.sizey || 1;
        this.sizez = data.sizez || 1;
        this.resetIsolevelState(Number.isFinite(data.isolevel) ? data.isolevel : this.getIsolevel());
        this.clearIsosurfacePreviewCache();
        this.initializeBondRulesFromAtoms();
        this.refreshAppearanceControls();
        this.shouldAutoFitCamera = true;
        this.updateStructure();
    }

    setRepetitions(nx, ny, nz) {
        this.nx = Math.max(1, parseInt(nx, 10) || 1);
        this.ny = Math.max(1, parseInt(ny, 10) || 1);
        this.nz = Math.max(1, parseInt(nz, 10) || 1);
        this.updateStructure();
    }

    hasChargeDensity() {
        return Array.isArray(this.values) && this.values.length > 0 && this.gridCell;
    }

    setCellCheckbox(domCheckbox) {
        if (!domCheckbox || !domCheckbox.length) {
            return;
        }
        this.cell = !!domCheckbox.prop('checked');
        domCheckbox.on('change', () => {
            this.cell = !!domCheckbox.prop('checked');
            this.updateStructure();
        });
    }

    setShadingCheckbox(domCheckbox) {
        if (!domCheckbox || !domCheckbox.length) {
            return;
        }
        this.shading = !!domCheckbox.prop('checked');
        domCheckbox.on('change', () => {
            this.shading = !!domCheckbox.prop('checked');
            this.updateStructure();
        });
    }

    getDefaultBondCutoff(atomNumberA, atomNumberB) {
        return getChemicalBondLimit(atomNumberA, atomNumberB, covalent_radii);
    }

    getReplicatedAtoms(nx, ny, nz) {
        if (!this.baseAtoms || !this.baseLattice) {
            return [];
        }

        const atoms = [];
        for (let ix = 0; ix < nx; ix++) {
            for (let iy = 0; iy < ny; iy++) {
                for (let iz = 0; iz < nz; iz++) {
                    for (let i = 0; i < this.baseAtoms.length; i++) {
                        const atom = this.baseAtoms[i];
                        atoms.push([
                            atom[0],
                            atom[1] + ix * this.baseLattice[0][0] + iy * this.baseLattice[1][0] + iz * this.baseLattice[2][0],
                            atom[2] + ix * this.baseLattice[0][1] + iy * this.baseLattice[1][1] + iz * this.baseLattice[2][1],
                            atom[3] + ix * this.baseLattice[0][2] + iy * this.baseLattice[1][2] + iz * this.baseLattice[2][2],
                        ]);
                    }
                }
            }
        }

        return atoms;
    }

    initializeBondRulesFromAtoms() {
        if (!this.baseAtoms || !this.baseAtoms.length || !this.atom_numbers || !this.atom_numbers.length) {
            this.bondRules = {};
            return;
        }

        this.bondRules = buildCrystalBondRules(
            this.getReplicatedAtoms(2, 2, 2),
            this.atom_numbers,
            covalent_radii,
            this.getBondRuleKey.bind(this),
        );
    }

    getSupercellLattice() {
        if (!this.baseLattice) {
            return null;
        }
        return [
            this.baseLattice[0].map((value) => value * this.nx),
            this.baseLattice[1].map((value) => value * this.ny),
            this.baseLattice[2].map((value) => value * this.nz),
        ];
    }

    updateViewParameters() {
        if (!this.atoms || !this.atoms.length) {
            return;
        }

        let maxRadius = 0;
        for (let i = 0; i < this.atoms.length; i++) {
            const atom = this.atoms[i];
            const position = new THREE.Vector3(atom[1], atom[2], atom[3]).sub(this.geometricCenter);
            maxRadius = Math.max(maxRadius, position.length());
        }

        const halfFov = (this.cameraViewAngle * Math.PI / 180) / 2;
        this.cameraDistance = Math.max(45, (maxRadius / Math.tan(halfFov)) * 1.25);
        this.setCameraDirection('z');
    }

    addCell(lat) {
        if (!this.cell || !lat) {
            return;
        }
        this.scene.add(createCellLineObject(lat, this.geometricCenter));
    }

    getMarchingCubesOptions() {
        return {
            periodic: true,
            insideIsAbove: true,
        };
    }

    addMarchingCubes() {
        if (!this.hasChargeDensity()) {
            return;
        }

        this.requestMarchingCubesUpdate();
    }

    addMarchingCubesGeometry(geometry) {
        const material = this.createIsosurfaceMaterial();

        for (let ix = 0; ix < this.nx; ix++) {
            for (let iy = 0; iy < this.ny; iy++) {
                for (let iz = 0; iz < this.nz; iz++) {
                    const mesh = new THREE.Mesh(geometry, material);
                    mesh.name = 'isosurface';
                    mesh.position.set(
                        ix * this.baseLattice[0][0] + iy * this.baseLattice[1][0] + iz * this.baseLattice[2][0] - this.geometricCenter.x,
                        ix * this.baseLattice[0][1] + iy * this.baseLattice[1][1] + iz * this.baseLattice[2][1] - this.geometricCenter.y,
                        ix * this.baseLattice[0][2] + iy * this.baseLattice[1][2] + iz * this.baseLattice[2][2] - this.geometricCenter.z,
                    );
                    this.scene.add(mesh);
                }
            }
        }
    }

    addIsosurfaceObject(object) {
        if (!object) {
            return;
        }

        if (object.userData && object.userData.isRaymarchedIsosurface) {
            const baseMatrix = new THREE.Matrix4().set(
                this.baseLattice[0][0], this.baseLattice[1][0], this.baseLattice[2][0], 0,
                this.baseLattice[0][1], this.baseLattice[1][1], this.baseLattice[2][1], 0,
                this.baseLattice[0][2], this.baseLattice[1][2], this.baseLattice[2][2], 0,
                0, 0, 0, 1,
            );
            const superLattice = this.getSupercellLattice();
            const superMatrix = new THREE.Matrix4().set(
                superLattice[0][0], superLattice[1][0], superLattice[2][0], 0,
                superLattice[0][1], superLattice[1][1], superLattice[2][1], 0,
                superLattice[0][2], superLattice[1][2], superLattice[2][2], 0,
                0, 0, 0, 1,
            );
            const transform = new THREE.Matrix4().copy(superMatrix).multiply(new THREE.Matrix4().copy(baseMatrix).invert());
            const mesh = object.clone();
            mesh.geometry = object.geometry.clone();
            mesh.material = object.material.clone();
            mesh.onBeforeRender = object.onBeforeRender;
            mesh.geometry.applyMatrix4(transform);
            mesh.name = 'isosurface';
            mesh.position.set(-this.geometricCenter.x, -this.geometricCenter.y, -this.geometricCenter.z);
            mesh.userData = {
                ...object.userData,
                gridCell: superLattice.map((vector) => vector.slice()),
                textureRepeat: [this.nx, this.ny, this.nz],
            };
            if (mesh.material && mesh.material.uniforms && mesh.material.uniforms.uTextureRepeat) {
                mesh.material.uniforms.uTextureRepeat.value.set(this.nx, this.ny, this.nz);
            }
            if (mesh.material && mesh.material.uniforms && mesh.material.uniforms.uStepCount) {
                const repeatScale = Math.max(this.nx || 1, this.ny || 1, this.nz || 1);
                mesh.material.uniforms.uStepCount.value = clampRaymarchStepCount(
                    mesh.material.uniforms.uStepCount.value * repeatScale,
                    [this.sizex || 1, this.sizey || 1, this.sizez || 1],
                );
            }
            refreshRaymarchedIsosurfaceUniforms(mesh, mesh.userData.gridCell);
            this.scene.add(mesh);
            return;
        }

        for (let ix = 0; ix < this.nx; ix++) {
            for (let iy = 0; iy < this.ny; iy++) {
                for (let iz = 0; iz < this.nz; iz++) {
                    const mesh = object.clone();
                    mesh.name = 'isosurface';
                    mesh.position.set(
                        ix * this.baseLattice[0][0] + iy * this.baseLattice[1][0] + iz * this.baseLattice[2][0] - this.geometricCenter.x,
                        ix * this.baseLattice[0][1] + iy * this.baseLattice[1][1] + iz * this.baseLattice[2][1] - this.geometricCenter.y,
                        ix * this.baseLattice[0][2] + iy * this.baseLattice[1][2] + iz * this.baseLattice[2][2] - this.geometricCenter.z,
                    );
                    this.scene.add(mesh);
                }
            }
        }
    }

    updateStructure() {
        if (!this.scene || !this.structureData || !this.baseLattice || !this.baseAtoms.length) {
            return;
        }

        this.atoms = this.getReplicatedAtoms(this.nx, this.ny, this.nz);
        this.currentLattice = this.getSupercellLattice();

        this.geometricCenter = new THREE.Vector3(0, 0, 0);
        for (let i = 0; i < this.atoms.length; i++) {
            this.geometricCenter.add(new THREE.Vector3(this.atoms[i][1], this.atoms[i][2], this.atoms[i][3]));
        }
        this.geometricCenter.multiplyScalar(1 / this.atoms.length);

        if (this.shouldAutoFitCamera) {
            this.updateViewParameters();
        }
        this.removeStructure();
        this.addLights();
        this.addMarchingCubes();
        this.getAtomMaterials();
        this.addStructure();
        this.addCell(this.currentLattice);
        this.render();
        this.shouldAutoFitCamera = false;
    }
}
