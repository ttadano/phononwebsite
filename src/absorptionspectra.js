import $ from 'jquery';

function getTag(tags, object) {
    for (let i = 0; i < tags.length; i++) {
        const tag = tags[i];
        if (tag in object) {
            return object[tag];
        }
    }
    throw new Error(`${tags.join(', ')} not found in the file.`);
}

export class AbsorptionSpectra {

    constructor(container) {
        this.container = container;
        this.chart = null;
        this.excitonIndex = 0;
        this.onExcitonSelected = null;
        this.excitons = [];
        this.energies = [];
        this.eps = [];
    }

    setSelectionHandler(handler) {
        this.onExcitonSelected = handler;
    }

    getDataObject(data) {
        this.rawData = data;
        this.excitons = data.excitons || [];
        this.energies = getTag(['E/ev[1]'], data);
        this.eps = getTag(['EPS-Im[2]'], data);
        this.sizex = data.nx;
        this.sizey = data.ny;
        this.sizez = data.nz;
        this.nndist = data.nndist;
        this.cell = data.lattice;
        this.gridCell = data.supercell_lattice || data.lattice;
        this.atoms = data.atoms;
        this.natoms = this.atoms.length;
        this.atom_numbers = data.atypes;
        this.excitonIndex = 0;
    }

    getDataFilename(filename) {
        return $.getJSON(filename).then((data) => {
            this.getDataObject(data);
            this.render();
            return data;
        });
    }

    setExcitonIndex(index, notify = false) {
        if (!this.excitons.length) {
            this.excitonIndex = 0;
            return;
        }

        const nextIndex = Math.max(0, Math.min(index, this.excitons.length - 1));
        this.excitonIndex = nextIndex;
        this.highlightSelectedExciton();

        if (notify && this.onExcitonSelected) {
            this.onExcitonSelected(nextIndex);
        }
    }

    render() {
        const options = {
            chart: { type: 'line', zoomType: 'xy' },
            accessibility: { enabled: false },
            title: { text: 'Absorption Spectra' },
            xAxis: {
                title: { text: 'Energy (eV)' },
                plotLines: this.getExcitonPlotLines(),
            },
            yAxis: {
                min: 0,
                title: { text: 'Intensity (arb. units)' },
                plotLines: [{ value: 0, color: '#808080', width: 1 }],
            },
            tooltip: {
                formatter() {
                    return `${this.x.toFixed(3)} eV, ${this.y.toFixed(3)}`;
                }
            },
            legend: { enabled: false },
            plotOptions: {
                line: { animation: false },
                series: {
                    cursor: 'pointer',
                    point: {
                        events: {
                            click: (event) => {
                                const index = this.getClosestExcitonIndex(event.point.x);
                                this.setExcitonIndex(index, true);
                            }
                        }
                    }
                }
            },
            series: [
                {
                    name: 'eps',
                    color: '#0066ff',
                    marker: { radius: 1.5, symbol: 'circle' },
                    data: this.energies.map((energy, index) => [energy, this.eps[index]]),
                },
                {
                    type: 'scatter',
                    name: 'excitons',
                    color: '#ef4444',
                    marker: { radius: 4, symbol: 'circle' },
                    data: this.excitons.map((exciton, index) => ({
                        x: exciton.energy,
                        y: this.getIntensityAtEnergy(exciton.energy),
                        excitonIndex: index,
                    })),
                    point: {
                        events: {
                            click: function () {
                                const index = this.options.excitonIndex;
                                this.series.chart.userOptions.__spectra.setExcitonIndex(index, true);
                            }
                        }
                    }
                }
            ],
            __spectra: this,
        };

        this.chart = globalThis.Highcharts.chart(this.container[0], options);
        if (this.chart && typeof this.chart.reflow === 'function') {
            this.chart.reflow();
        }
        this.highlightSelectedExciton();
    }

    getExcitonPlotLines() {
        return this.excitons.map((exciton) => ({
            value: exciton.energy,
            color: '#9ca3af',
            width: 1,
            zIndex: 4,
        }));
    }

    highlightSelectedExciton() {
        if (!this.chart || !this.chart.xAxis || !this.excitons.length) {
            return;
        }

        const axis = this.chart.xAxis[0];
        axis.removePlotLine('selected-exciton');
        axis.addPlotLine({
            id: 'selected-exciton',
            value: this.excitons[this.excitonIndex].energy,
            color: '#ef4444',
            width: 2,
            zIndex: 5,
        });
    }

    getClosestExcitonIndex(energy) {
        let minDistance = Infinity;
        let closest = 0;

        for (let i = 0; i < this.excitons.length; i++) {
            const distance = Math.abs(this.excitons[i].energy - energy);
            if (distance < minDistance) {
                minDistance = distance;
                closest = i;
            }
        }

        return closest;
    }

    getIntensityAtEnergy(energy) {
        if (!this.energies.length || !this.eps.length) {
            return 0;
        }

        let closestIndex = 0;
        let minDistance = Infinity;

        for (let i = 0; i < this.energies.length; i++) {
            const distance = Math.abs(this.energies[i] - energy);
            if (distance < minDistance) {
                minDistance = distance;
                closestIndex = i;
            }
        }

        return this.eps[closestIndex];
    }
}
