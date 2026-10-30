import json

from phononweb.phonondb import PhononDB


def test_phonondb_save_and_load_roundtrip(tmp_path):
    db = PhononDB(url="2015")
    db.materials = [["149", "mp-149", "graphene"]]

    savefile = tmp_path / "phonondb.json"
    db.save_materials(str(savefile))
    assert savefile.exists()

    loaded = PhononDB(url="2015")
    loaded.load_materials(str(savefile))
    assert loaded.materials == db.materials


def test_phonondb_string_representation():
    db = PhononDB(url="2018")
    db.materials = [["149", "mp-149", "graphene"]]

    text = str(db)
    assert "total materials: 1" in text
    assert "mp-149" in text


def test_phonondb_file_format_is_json(tmp_path):
    db = PhononDB(url="2015")
    db.materials = [["149", "mp-149", "graphene"]]
    savefile = tmp_path / "phonondb.json"

    db.save_materials(str(savefile))
    content = json.loads(savefile.read_text())
    assert isinstance(content, list)
