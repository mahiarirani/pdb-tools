import typer
import os
from PDB import Model, Chain, Residue
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
def residue_neighbors(
        pdb_file: str,
        chain: str,
        residue: int,
        model: int = '0',
        probe_points: int = 100,
        probe_radius: float = 1.4,
        batch: bool = typer.Option(False),
        extended: bool = typer.Option(False),
        minimal: bool = typer.Option(True)
):
    if batch:
        directory = pdb_file
        for file in os.listdir(pdb_file):
            pdb_file = directory + '/' + file
            residue_neighbors(pdb_file, chain, residue, model, probe_points, probe_radius, False, extended, minimal)
    else:
        PDBTools(pdb_file, probe_points, probe_radius, extended, minimal)\
            .residue_neighbors(Model(model), Chain(chain), Residue(residue))


@app.command()
def chain_neighbors(
        pdb_file: str,
        chain: str,
        model: int = '0',
        probe_points: int = 100,
        probe_radius: float = 1.4,
        batch: bool = typer.Option(False)
):
    if batch:
        directory = pdb_file
        for file in os.listdir(pdb_file):
            pdb_file = directory + '/' + file
            PDBTools(pdb_file, probe_points, probe_radius)\
                .chain_neighbors(Model(model), Chain(chain))
    else:
        PDBTools(pdb_file, probe_points, probe_radius)\
            .chain_neighbors(Model(model), Chain(chain))


if __name__ == "__main__":
    app()
