export const RAYMARCH_PERFORMANCE_OPTION_DESCRIPTIONS = {
    opacityCompensation: 'Multiplies the UI opacity target before converting it into per-hit raymarch alpha.',
    interactive: 'Settings used while the user is actively moving the camera or otherwise interacting.',
    settled: 'Settings used after interaction stops and the viewer returns to its higher-quality resting mode.',
    resolutionDivisor: 'Integer divisor applied to the real canvas pixel ratio; 2 means half-resolution rendering in each dimension.',
    interpolation: 'Sampling mode for the scalar field; use "trilinear", "tricubic", or "selected" to inherit the UI choice.',
    stepScale: 'Multiplier that converts grid size into raymarch step count; lower values are faster but less accurate.',
    minSteps: 'Lower bound for the ray step count after scaling.',
    maxSteps: 'Upper bound for the ray step count after scaling.',
    maxRayHits: 'Maximum number of isosurface crossings accumulated along a ray; fewer hits are faster and lighter-looking.',
};

export const RAYMARCH_PERFORMANCE_SETTINGS = {
    opacityCompensation: 1.7,
    interactive: {
        resolutionDivisor: 1,
        interpolation: 'selected',
        stepScale: 1.0,
        minSteps: 48,
        maxSteps: 320,
        maxRayHits: 1,
    },
    settled: {
        resolutionDivisor: 1,
        interpolation: 'selected',
        stepScale: 1.0,
        minSteps: 48,
        maxSteps: 320,
        maxRayHits: 1,
    },
};

export function getRaymarchPerformanceProfile({ interacting, selectedInterpolation, gridSizeMax }) {
    const profile = interacting
        ? RAYMARCH_PERFORMANCE_SETTINGS.interactive
        : RAYMARCH_PERFORMANCE_SETTINGS.settled;

    const interpolation = profile.interpolation === 'selected'
        ? selectedInterpolation
        : profile.interpolation;

    const rawSteps = Math.round(gridSizeMax * 2.0 * profile.stepScale);
    const stepCount = Math.max(profile.minSteps, Math.min(profile.maxSteps, rawSteps));

    return {
        interpolation,
        stepCount,
        resolutionDivisor: Math.max(1, Math.round(profile.resolutionDivisor || 1)),
        activeRayHits: Math.max(1, Math.min(2, Math.round(profile.maxRayHits || 1))),
    };
}
