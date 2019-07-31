# FEM based solver for calculating effective permittivity of a heterogeneous
# structure made of inner inclusions with inner_permittivity (mesh.subdomain = 1)
# and outside material matrix with outer_permittivity (mesh.subdomain = 2)

    # Domain defining parameters (permittivity) are hardcoded in main part as:
    #   patch_permittivity = 1
    #   matrix_permittivity = 11.7

# Function call: python3 Yproblem.py mesh_folder mesh_name output_folder
# ie. python3 Yproblem.py mesh hexagonal results

# input = unit_cell mesh with subdomain markers in .h5 format
# output = txt file with 2x2 matrix of effective permittivity

# Using FEniCS 2017.2.0
import petsc4py
import sys
petsc4py.init(sys.argv)
from petsc4py import PETSc
import subprocess

import dolfin as df
import argparse

import yproblem

#TODO: maybe read file?
if df.MPI.rank(df.mpi_comm_world()) == 0:
    # test version like this because of portability chaos...
    dolfin_version = subprocess.run(['dolfin-version'],
                                    stdout=subprocess.PIPE).stdout.decode().strip('\n')
    if dolfin_version != "2017.2.0":
        raise AssertionError("You are using {} ".format(dolfin_version) +
                "FEniCS and code is tested using 2017.2.0 FEniCS version.")

def main(mesh_filename, subdomains):
    y = yproblem.Yproblem(mesh_filename, subdomains)
    V = df.FunctionSpace(y.mesh, "P", 1,
            constrained_domain = yproblem.PeriodicBoundary())
    #TODO tic
    effective = y.get_effective(V)
    #TODO toc

    if df.MPI.rank(df.mpi_comm_world()) == 0:
        print (effective)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "FEM based solver for " +
    "calculating effective permittivity of a heterogeneous\n" +
    "structure made of inner inclusions with inner permittivity and\n" +
    "outside material matrix with outer permittivity")

    parser.add_argument("permittivity", nargs='+',
                        help="Permittivity accordingly to subdomain_id")
    parser.add_argument("-o", "--output", help="Folder for outputing results")
    # don't use -m because that is reserved and there is no bash autocomplete
    parser.add_argument("--mesh", help="Unit cell mesh filename")

    args, petsc_args = parser.parse_known_args()

    print ("Using PETSc args {}".format(petsc_args))

    subdomains = {i: args.permittivity[i-1] for i in
                        range(1, len(args.permittivity)+1)}

    main(args.mesh, subdomains)
