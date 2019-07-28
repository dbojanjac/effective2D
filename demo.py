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
import dolfin as df
import subprocess
import yproblem

# test version like this because of portability chaos...
dolfin_version = subprocess.run(['dolfin-version'],
                                stdout=subprocess.PIPE).stdout.decode().strip('\n')
if dolfin_version != "2017.2.0":
    raise AssertionError("You are using {} FEniCS and code is tested using 2017.2.0 FEniCS version.".format(dolfin_version))

def main(mesh_filename, subdomains):
    y = yproblem.Yproblem(mesh_filename, subdomains)
    V = df.FunctionSpace(y.mesh, "P", 1,
            constrained_domain = yproblem.PeriodicBoundary())
    #TODO tic
    print (y.get_effective(V))
    #TODO toc

if __name__ == "__main__":
    subdomains = {1: 1, 2: 11.8}
    main("mesh/hexagonal.h5", subdomains)
