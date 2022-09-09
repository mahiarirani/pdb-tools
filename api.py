from fastapi import FastAPI, File, UploadFile
import os
import json

from modules import PDBTools
from models import Model, Chain, Residue

app = FastAPI()


@app.get("/")
def get_index():
    return {"PDBTools.io Public API"}


@app.post("/pdb")
async def post_pdb(file: UploadFile):
    try:
        contents = await file.read()
        name = file.filename
        with open('./PDB/' + name, 'wb') as f:
            f.write(contents)
    except Exception:
        return {"message": "There was an error uploading the file"}
    finally:
        file.file.close()
    return name


@app.get("/sasa")
def sasa(pdb: str, probe_points: int = 100, probe_radius: float = 1.4):
    PDBTools('./PDB/' + pdb + '.pdb', probe_points, probe_radius).sasa()
    return __load_json(pdb)


@app.get("/residue")
def sasa(pdb: str, chain: str, residue: int, model: int = 0,  probe_points: int = 100, probe_radius: float = 1.4):
    PDBTools('./PDB/' + pdb + '.pdb', probe_points, probe_radius).residue_neighbors(Model(model), Chain(chain), Residue(residue))
    return __load_json(pdb)


@app.get("/chain")
def sasa(pdb: str, chain: str, model: int = 0,  probe_points: int = 100, probe_radius: float = 1.4):
    PDBTools('./PDB/' + pdb + '.pdb', probe_points, probe_radius).chain_neighbors(Model(model), Chain(chain))
    return __load_json(pdb)


@app.get("/critical")
def sasa(pdb: str, chain: str, model: int = 0,  probe_points: int = 100, probe_radius: float = 1.4):
    name, extension = os.path.splitext(pdb)
    PDBTools('./PDB/' + pdb + '.pdb', probe_points, probe_radius).critical_residues(Model(model), Chain(chain))
    return __load_json(pdb)


def __load_json(name):
    with open('./PDBObject/' + name + '.json', 'r') as f:
        return json.load(f)
