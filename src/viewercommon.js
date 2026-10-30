import * as THREE from 'three';

import { atomic_symbol, jmol_colors, vesta_colors } from './atomic_data.js';
import { createAtomBadgeHtml } from './atomcolors.js';

export function getSharedLightConfig() {
    return {
        color: 0xdddddd,
        intensity: Math.PI * 1.08,
        position: [1, 1, 2],
        ambient: 0x262626,
        ambientIntensity: 1.08,
        decay: 0,
    };
}

export const sharedViewerMethods = {
    colorToInputHex(colorHex) {
        return `#${Number(colorHex).toString(16).padStart(6, '0')}`;
    },

    normalizeColorHex(value, fallback) {
        if (typeof value === 'number' && Number.isFinite(value)) {
            return value;
        }
        if (typeof value === 'string') {
            let normalized = value.trim();
            if (!normalized) {
                return fallback;
            }
            if (normalized.startsWith('#')) {
                normalized = normalized.slice(1);
            }
            if (/^[0-9a-fA-F]{6}$/.test(normalized)) {
                return parseInt(normalized, 16);
            }
        }
        return fallback;
    },

    getDefaultAtomColor(atomNumber) {
        const palette = this.display === 'vesta' ? vesta_colors : jmol_colors;
        const rgb = palette[atomNumber] || [0.5, 0.5, 0.5];
        return new THREE.Color(rgb[0], rgb[1], rgb[2]).getHex();
    },

    getAtomColorHex(atomNumber) {
        if (Object.prototype.hasOwnProperty.call(this.atomColorOverrides, atomNumber)) {
            return this.atomColorOverrides[atomNumber];
        }
        return this.getDefaultAtomColor(atomNumber);
    },

    getAtomColor(atomNumber) {
        return new THREE.Color(this.getAtomColorHex(atomNumber));
    },

    clearAtomColorOverride(atomNumber) {
        delete this.atomColorOverrides[atomNumber];
    },

    getAtomRadiusScale(atomNumber) {
        if (Object.prototype.hasOwnProperty.call(this.atomRadiusScaleOverrides, atomNumber)) {
            return this.atomRadiusScaleOverrides[atomNumber];
        }
        return this.defaultAtomRadiusScale;
    },

    setAtomRadiusScaleOverride(atomNumber, scale) {
        if (!Number.isFinite(scale)) {
            return;
        }
        this.atomRadiusScaleOverrides[atomNumber] = Math.max(0.1, scale);
    },

    getSelectedAppearanceAtomNumber() {
        if (Number.isFinite(this.appearanceSelectedAtomNumber)) {
            return this.appearanceSelectedAtomNumber;
        }
        if (this.atom_numbers && this.atom_numbers.length) {
            return this.atom_numbers[0];
        }
        return null;
    },

    setSelectedAppearanceAtomNumber(atomNumber) {
        this.appearanceSelectedAtomNumber = atomNumber;
        const atomList = this.domAppearanceAtomList || this.dom_appearance_atom_list;
        if (atomList && atomList.length) {
            atomList.find('button').removeClass('active');
            atomList.find(`button[data-atom-number="${atomNumber}"]`).addClass('active');
        }
    },

    getBondRuleKey(atomNumberA, atomNumberB) {
        const a = Math.min(atomNumberA, atomNumberB);
        const b = Math.max(atomNumberA, atomNumberB);
        return `${a}-${b}`;
    },

    setBondRule(atomNumberA, atomNumberB, cutoff) {
        const key = this.getBondRuleKey(atomNumberA, atomNumberB);
        const a = Math.min(atomNumberA, atomNumberB);
        const b = Math.max(atomNumberA, atomNumberB);
        this.bondRules[key] = {
            a,
            b,
            cutoff: Number.isFinite(cutoff) ? cutoff : this.getDefaultBondCutoff(a, b),
        };
    },

    removeBondRule(atomNumberA, atomNumberB) {
        const key = this.getBondRuleKey(atomNumberA, atomNumberB);
        delete this.bondRules[key];
    },

    hasBondRule(atomNumberA, atomNumberB) {
        const key = this.getBondRuleKey(atomNumberA, atomNumberB);
        return Object.prototype.hasOwnProperty.call(this.bondRules, key);
    },

    refreshBondRulesUI() {
        const domBondRulesList = this.domBondRulesList || this.dom_bond_rules_list;
        if (!domBondRulesList || !domBondRulesList.length) {
            return;
        }

        domBondRulesList.empty();
        const keys = Object.keys(this.bondRules).sort();
        if (!keys.length) {
            domBondRulesList.append('<div>none</div>');
            return;
        }

        for (let i = 0; i < keys.length; i++) {
            const rule = this.bondRules[keys[i]];
            const label =
                '<span class="atom-badge-pair">' +
                createAtomBadgeHtml(atomic_symbol[rule.a], rule.a, this.getAtomColorHex.bind(this)) +
                '<span class="atom-badge-separator">-</span>' +
                createAtomBadgeHtml(atomic_symbol[rule.b], rule.b, this.getAtomColorHex.bind(this)) +
                '</span>';
            const cutoff = Number(rule.cutoff).toFixed(2);
            domBondRulesList.append(
                '<div class="appearance-controls">' +
                '<span>' + label + ' ' + cutoff + '</span>' +
                '<button type="button" data-remove-key="' + keys[i] + '">remove</button>' +
                '</div>'
            );
        }
    },

    addLights() {
        this.scene.add(this.camera);
        if (this.pointLight) {
            this.pointLight.visible = true;
        }
        const lightConfig = getSharedLightConfig();
        this.scene.add(new THREE.AmbientLight(lightConfig.ambient, lightConfig.ambientIntensity));
    },

    updateLightStyle() {
        if (!this.pointLight) {
            return;
        }
        const lightConfig = getSharedLightConfig();
        this.pointLight.color.setHex(lightConfig.color);
        this.pointLight.intensity = lightConfig.intensity;
        this.pointLight.position.set(...lightConfig.position);
    },

    createShadedMaterial(config = {}) {
        if (!this.shading) {
            return new THREE.MeshBasicMaterial(config);
        }
        return new THREE.MeshLambertMaterial({
            blending: THREE.NormalBlending,
            ...config,
        });
    },

    disposeSceneObject(object) {
        if (!object) {
            return;
        }
        if (object.geometry) {
            object.geometry.dispose();
        }
        if (object.material) {
            if (Array.isArray(object.material)) {
                for (let i = 0; i < object.material.length; i++) {
                    if (object.material[i] && object.material[i].dispose) {
                        object.material[i].dispose();
                    }
                }
            } else if (object.material.dispose) {
                object.material.dispose();
            }
        }
    },

    removeNamedSceneObjects(name) {
        if (!this.scene) {
            return;
        }

        for (let i = this.scene.children.length - 1; i >= 0; i--) {
            const child = this.scene.children[i];
            if (child && child.name === name) {
                this.disposeSceneObject(child);
                this.scene.remove(child);
            }
        }
    },

    forEachNamedSceneObject(name, callback) {
        if (!this.scene || typeof callback !== 'function') {
            return;
        }

        for (let i = 0; i < this.scene.children.length; i++) {
            const child = this.scene.children[i];
            if (child && child.name === name) {
                callback(child);
            }
        }
    },
};
