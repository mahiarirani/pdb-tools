import typer
import os
from PDB import Residue
from PDBTools import PDBTools

app = typer.Typer()


@app.command()
def sasa(pdb_file: str, probe_points: int = 100, probe_radius: float = 1.4, batch: bool = typer.Option(False)):
    if batch:
        directory = pdb_file
        for file in os.listdir(pdb_file):
            pdb_file = directory + '/' + file
            PDBTools(pdb_file, probe_points, probe_radius).sasa()
    else:
        PDBTools(pdb_file, probe_points, probe_radius).sasa()


@app.command()
def residue_neighbors(pdb_file: str, residue, probe_points: int = 100, probe_radius: float = 1.4, batch: bool = typer.Option(False)):
    if batch:
        directory = pdb_file
        for file in os.listdir(pdb_file):
            pdb_file = directory + '/' + file
            PDBTools(pdb_file, probe_points, probe_radius).residue_neighbors(residue)
    else:
        PDBTools(pdb_file, probe_points, probe_radius).residue_neighbors(residue)


@app.command()
def chain_neighbors(pdb_file: str, chain, probe_points: int = 100, probe_radius: float = 1.4, batch: bool = typer.Option(False)):
    if batch:
        directory = pdb_file
        for file in os.listdir(pdb_file):
            pdb_file = directory + '/' + file
            PDBTools(pdb_file, probe_points, probe_radius).chain_neighbors(chain)
    else:
        PDBTools(pdb_file, probe_points, probe_radius).chain_neighbors(chain)


if __name__ == "__main__":
    app()
