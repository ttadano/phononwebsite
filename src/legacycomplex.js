function normalizeOperand(value) {
    if (value && value.__rawComplex) {
        return value.__rawComplex;
    }

    if (value && typeof value.real === 'function' && typeof value.imag === 'function') {
        return { re: Number(value.real()), im: Number(value.imag()) };
    }

    if (value && typeof value.real === 'number') {
        let imag = 0;
        if (typeof value.im === 'number') {
            imag = value.im;
        } else if (typeof value.imag === 'number') {
            imag = value.imag;
        }
        return { re: Number(value.real), im: Number(imag) };
    }

    let numeric = Number(value);
    if (Number.isFinite(numeric)) {
        return { re: numeric, im: 0 };
    }

    return { re: 0, im: 0 };
}

function createComplex(real, imag) {
    let raw = {
        re: Number(real) || 0,
        im: Number(imag) || 0,
    };

    return {
        __rawComplex: raw,
        mult(other) {
            let rhs = normalizeOperand(other);
            return createComplex(
                raw.re * rhs.re - raw.im * rhs.im,
                raw.re * rhs.im + raw.im * rhs.re
            );
        },
        multiply(other) {
            return this.mult(other);
        },
        real() {
            return raw.re;
        },
        imag() {
            return raw.im;
        },
    };
}

createComplex.Polar = function(r, phi) {
    let radius = Number(r) || 0;
    let angle = Number(phi) || 0;
    return createComplex(radius * Math.cos(angle), radius * Math.sin(angle));
};

createComplex.fromPolar = createComplex.Polar;

export { createComplex as Complex };
