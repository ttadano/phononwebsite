import * as THREE from 'three';

export function createAtomSphereGeometry(radius, sphereLat, sphereLon) {
    return new THREE.SphereGeometry(radius, sphereLat, sphereLon);
}

export function createBondCylinderGeometry(radius, length, bondSegments, bondVertical) {
    return new THREE.CylinderGeometry(
        radius,
        radius,
        length,
        bondSegments,
        bondVertical,
        true,
    );
}

export function createCellLineObject(lat, shift, color = 0x000000) {
    const material = new THREE.LineBasicMaterial({ color });
    const points = [];
    const zero = new THREE.Vector3(0, 0, 0);
    const cursor = new THREE.Vector3(0, 0, 0);
    const x = new THREE.Vector3(lat[0][0], lat[0][1], lat[0][2]);
    const y = new THREE.Vector3(lat[1][0], lat[1][1], lat[1][2]);
    const z = new THREE.Vector3(lat[2][0], lat[2][1], lat[2][2]);

    cursor.copy(zero);
    cursor.sub(shift); points.push(cursor.clone());
    cursor.add(x); points.push(cursor.clone());
    cursor.add(y); points.push(cursor.clone());
    cursor.sub(x); points.push(cursor.clone());
    cursor.sub(y); points.push(cursor.clone());

    cursor.copy(zero).add(z);
    cursor.sub(shift); points.push(cursor.clone());
    cursor.add(x); points.push(cursor.clone());
    cursor.add(y); points.push(cursor.clone());
    cursor.sub(x); points.push(cursor.clone());
    cursor.sub(y); points.push(cursor.clone());

    cursor.copy(zero);
    cursor.sub(shift); points.push(cursor.clone());
    cursor.add(z); points.push(cursor.clone());

    cursor.add(x); points.push(cursor.clone());
    cursor.sub(z); points.push(cursor.clone());

    cursor.add(y); points.push(cursor.clone());
    cursor.add(z); points.push(cursor.clone());

    cursor.sub(x); points.push(cursor.clone());
    cursor.sub(z); points.push(cursor.clone());

    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    return new THREE.Line(geometry, material);
}
