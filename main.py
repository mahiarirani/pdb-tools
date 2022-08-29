import typer
import os
import csv
from PDB import Model, Chain, Residue
from PDBTools import PDBTools

app = typer.Typer(help="PDBTools CLI")


@app.command(help="Calculate SASA")
def sasa(
        pdb_file: str = typer.Argument(..., help="PDB File"),
        probe_points: int = typer.Argument(100, help="Number of points in probe structure"),
        probe_radius:float = typer.Argument(1.4, help="Probe radius size"),
        extended: bool = typer.Option(False, help="Extend to atom detail"),
        minimal: bool = typer.Option(True, help="Keep results minimal")
):
    PDBTools(pdb_file, probe_points, probe_radius, extended, minimal).sasa()


@app.command(help="Get residue neighbors")
def residue_neighbors(
        pdb_file: str = typer.Argument(..., help="PDB File"),
        chain: str = typer.Argument('A', help="Chain name on which the command will execute"),
        residue: int = typer.Argument(0, help="Model number on which the command will execute"),
        model: int = typer.Argument(0, help="Model number on which the command will execute"),
        probe_points: int = typer.Argument(100, help="Number of points in probe structure"),
        probe_radius: float = typer.Argument(1.4, help="Probe radius size"),
        extended: bool = typer.Option(False, help="Extend to atom detail"),
        minimal: bool = typer.Option(True, help="Keep results minimal")
):
    PDBTools(pdb_file, probe_points, probe_radius, extended, minimal)\
        .residue_neighbors(Model(model), Chain(chain), Residue(residue))


@app.command(help="Get chain neighbors")
def chain_neighbors(
        pdb_file: str = typer.Argument(..., help="PDB File"),
        chain: str = typer.Argument('A', help="Chain name on which the command will execute"),
        model: int = typer.Argument(0, help="Model number on which the command will execute"),
        probe_points: int = typer.Argument(100, help="Number of points in probe structure"),
        probe_radius: float = typer.Argument(1.4, help="Probe radius size"),
        extended: bool = typer.Option(False, help="Extend to atom detail"),
        minimal: bool = typer.Option(True, help="Keep results minimal")
):
    PDBTools(pdb_file, probe_points, probe_radius, extended, minimal)\
        .chain_neighbors(Model(model), Chain(chain))


@app.command(help="Batch run commands")
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
