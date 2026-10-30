import $ from 'jquery';
import * as THREE from 'three';
import Detector from '../libs/Detector.js';

import { StructureViewer, StructureWebpage } from './structurewebsite.js';

if (THREE.ColorManagement && typeof THREE.ColorManagement.enabled === 'boolean') {
    THREE.ColorManagement.enabled = false;
}

globalThis.THREE = THREE;
globalThis.$ = $;
globalThis.jQuery = $;

const viewer = new StructureViewer();
viewer.init($('#structureviewer'));

const page = new StructureWebpage(viewer);
page.setTitle($('#name'));
page.setLattice($('#lattice'));
page.setAtomPositions($('#atompos'));
page.setFileInput($('#file-input'));
page.setIsosurfaceModeInput($('#isosurface_mode'));
page.setIsolevelInput($('#isolevel_range'), $('#isolevel_value'));
page.setIsosurfaceOpacityInput($('#isosurface_opacity_range'), $('#isosurface_opacity_value'));
page.setIsosurfaceSection($('#isosurface_section'));
page.setRepetitionControls($('#nx'), $('#ny'), $('#nz'), $('#update_replications'));
page.setCameraDirectionButton($('#camerax'), 'x');
page.setCameraDirectionButton($('#cameray'), 'y');
page.setCameraDirectionButton($('#cameraz'), 'z');
viewer.setCellCheckbox($('#drawcell'));
viewer.setShadingCheckbox($('#drawshading'));
viewer.setDisplayCombo($('#displaystyle'));
viewer.setAppearanceControls(
    $('#appearance_atom_list'),
    $('#atom_color_input'),
    $('#bond_color_input'),
    $('#bond_color_by_atom_checkbox'),
    $('#atom_radius_input'),
    $('#bond_radius_input'),
    $('#bond_rules_list'),
    $('#bond_add_atom_a'),
    $('#bond_add_atom_b'),
    $('#bond_add_cutoff_input'),
    $('#appearance_reset_atom_button'),
    $('#appearance_reset_bonds_button'),
);

page.init();

if (!Detector.webgl) {
    Detector.addGetWebGLMessage();
}
