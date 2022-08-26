import typer
import os
import csv
from PDB import Model, Chain, Residue
from PDBTools import PDBTools

app = typer.Typer()


@app.command()
def sasa(
        pdb_file: str,
        probe_points: int = 100,
        probe_radius: float = 1.4,
        extended: bool = typer.Option(False),
        minimal: bool = typer.Option(True)
):
    PDBTools(pdb_file, probe_points, probe_radius, extended, minimal).sasa()


@app.command()
def residue_neighbors(
        pdb_file: str,
        chain: str,
        residue: int,
        model: int = 0,
        probe_points: int = 100,
        probe_radius: float = 1.4,
        extended: bool = typer.Option(False),
        minimal: bool = typer.Option(True)
):
    PDBTools(pdb_file, probe_points, probe_radius, extended, minimal)\
        .residue_neighbors(Model(model), Chain(chain), Residue(residue))


@app.command()
def chain_neighbors(
        pdb_file: str,
        chain: str,
        model: int = 0,
        probe_points: int = 100,
        probe_radius: float = 1.4,
        extended: bool = typer.Option(False),
        minimal: bool = typer.Option(True)
):
    PDBTools(pdb_file, probe_points, probe_radius, extended, minimal)\
        .chain_neighbors(Model(model), Chain(chain))


@app.command()
def batch(
        directory: str,
):
    csv_file = directory + '/' + 'batch.csv'
    if not os.path.isfile(csv_file):
        print('BATCH.CSV doesn\'t exists!')
    else:
        with open(csv_file, 'r') as data:
            reader = csv.reader(data)
            for line in reader:
                pdb, command = line[0], line[1]
                pdb_file = directory + '/' + pdb + '.pdb'
                if command == 'sasa':
                    sasa(pdb_file)
                elif command == 'residue-neighbors':
                    chain, residue = line[2], int(line[3])
                    residue_neighbors(pdb_file, chain, residue)
                elif command == 'chain-neighbors':
                    chain = line[2]
                    chain_neighbors(pdb_file, chain)
                else:
                    print('Command not found!')


if __name__ == "__main__":
    app()
