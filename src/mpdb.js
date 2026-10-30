export class MaterialsProjectDB {
    /*
    Interact with the local database of phonons
    Hosted on Github
    */

    constructor() {
        this.name = "mpdb";
        this.year = 2018;
        this.author = "G. Petretto et al.";
        this.url = "https://www.nature.com/articles/sdata201865";
        this.probeUrl = "https://materialsproject-parsed.s3.amazonaws.com/ph-bandstructures/dfpt/mp-1000.json.gz";
    }

    isAvailable() {
        return false;
    }

    checkAvailability(callback) {
        if (MaterialsProjectDB.availabilityState !== undefined) {
            callback(MaterialsProjectDB.availabilityState);
            return;
        }

        if (typeof fetch !== 'function') {
            MaterialsProjectDB.availabilityState = false;
            callback(false);
            return;
        }

        fetch(this.probeUrl, { method: 'HEAD' })
            .then(function(response) {
                let available = response.ok;
                MaterialsProjectDB.availabilityState = available;
                callback(available);
            })
            .catch(function(error) {
                console.log("Materials Project OpenData unavailable from browser:", error);
                MaterialsProjectDB.availabilityState = false;
                callback(false);
            });
    }

    get_materials(callback) {
        /*
        this function load the materials from a certain source and returns then to the callback
        Some pre-processing of the data might be required and can be implemented here
        */
        let reference = this.author+", "+"<a href="+this.url+">"+this.name+"</a> ("+this.year+")";
        let name = this.name;

        function dothings(materials) {

            for (let i=0; i<materials.length; i++) {
                let m = materials[i];
                m.source = name;
                m.type = "json";
                m.reference = reference;
                m.url = "https://materialsproject-parsed.s3.amazonaws.com/ph-bandstructures/dfpt/mp-"+m.id+".json.gz";
                m.name = m.name;
                m.link = "https://materialsproject.org/materials/mp-"+m.id;
            }
            callback(materials);
        }

        $.get('data/mpdb/models.json', dothings);
    }

}

MaterialsProjectDB.availabilityState = undefined;
