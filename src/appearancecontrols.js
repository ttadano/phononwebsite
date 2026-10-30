export function createBondColorInputStateUpdater(viewer, bondColorInput) {
    return () => {
        if (bondColorInput && bondColorInput.length) {
            bondColorInput.prop('disabled', viewer.bondColorByAtom);
        }
    };
}

export function bindAppearanceAtomSelection(viewer, atomList, atomColorInput, atomRadiusInput) {
    if (!atomList || !atomList.length) {
        return;
    }

    atomList.on('click', 'button[data-atom-number]', (event) => {
        const atomNumber = Number(event.currentTarget.getAttribute('data-atom-number'));
        if (!Number.isFinite(atomNumber)) {
            return;
        }
        viewer.setSelectedAppearanceAtomNumber(atomNumber);
        if (atomColorInput && atomColorInput.length) {
            atomColorInput.val(viewer.colorToInputHex(viewer.getAtomColorHex(atomNumber)));
        }
        if (atomRadiusInput && atomRadiusInput.length) {
            atomRadiusInput.val(viewer.getAtomRadiusScale(atomNumber));
        }
    });
}

export function bindEnterToApply(inputs, apply) {
    for (let i = 0; i < inputs.length; i++) {
        const domInput = inputs[i];
        if (domInput && domInput.length) {
            domInput.on('keydown', (event) => {
                if (event.key === 'Enter') {
                    event.preventDefault();
                    apply();
                }
            });
        }
    }
}

export function bindBondRuleControls(
    viewer,
    domBondRulesList,
    domBondAddAtomA,
    domBondAddAtomB,
    domBondAddCutoffInput,
    domBondAddButton,
    onRulesChanged,
) {
    if (typeof domBondAddButton === 'function' && typeof onRulesChanged === 'undefined') {
        onRulesChanged = domBondAddButton;
        domBondAddButton = null;
    }

    if (domBondRulesList && domBondRulesList.length) {
        domBondRulesList.on('click', 'button[data-remove-key]', (event) => {
            const key = event.currentTarget.getAttribute('data-remove-key');
            if (key && viewer.bondRules[key]) {
                delete viewer.bondRules[key];
                onRulesChanged();
            }
        });
    }

    const updateBondCutoffInput = () => {
        if (!domBondAddCutoffInput || !domBondAddCutoffInput.length) {
            return;
        }
        const a = Number(domBondAddAtomA.val());
        const b = Number(domBondAddAtomB.val());
        if (!Number.isFinite(a) || !Number.isFinite(b)) {
            return;
        }
        const key = viewer.getBondRuleKey(a, b);
        const value = viewer.bondRules[key] ? viewer.bondRules[key].cutoff : viewer.getDefaultBondCutoff(a, b);
        domBondAddCutoffInput.val(Number(value).toFixed(2));
    };

    const addBondRuleFromControls = () => {
        const a = Number(domBondAddAtomA.val());
        const b = Number(domBondAddAtomB.val());
        if (!Number.isFinite(a) || !Number.isFinite(b)) {
            return;
        }
        let cutoff = viewer.getDefaultBondCutoff(a, b);
        if (domBondAddCutoffInput && domBondAddCutoffInput.length) {
            const parsed = parseFloat(domBondAddCutoffInput.val());
            if (Number.isFinite(parsed) && parsed > 0) {
                cutoff = parsed;
            }
            domBondAddCutoffInput.val(cutoff.toFixed(2));
        }
        viewer.setBondRule(a, b, cutoff);
        onRulesChanged();
    };

    if (domBondAddAtomA && domBondAddAtomA.length) {
        domBondAddAtomA.on('change', updateBondCutoffInput);
    }
    if (domBondAddAtomB && domBondAddAtomB.length) {
        domBondAddAtomB.on('change', updateBondCutoffInput);
    }
    if (domBondAddCutoffInput && domBondAddCutoffInput.length) {
        domBondAddCutoffInput.on('keydown', (event) => {
            if (event.key === 'Enter') {
                event.preventDefault();
                addBondRuleFromControls();
            }
        });
    }
    if (domBondAddButton && domBondAddButton.length) {
        domBondAddButton.on('click', (event) => {
            event.preventDefault();
            addBondRuleFromControls();
        });
    }
}
