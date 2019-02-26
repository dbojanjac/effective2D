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
import dolfin as df
import sys


def Y_solver_2D(mesh_folder, mesh_name, inner_permittivity, outer_permittivity):
    """Unit cell solver function"""

    # Read mesh and subdomain markers from mesh_folder/mesh_name
    #---------------------------------------------------------------------------
    mesh_folder = mesh_folder + '/'

    mesh = df.Mesh()
    hdf = df.HDF5File(mesh.mpi_comm(), mesh_folder + mesh_name + '.h5', 'r')

    hdf.read(mesh, mesh_folder + "mesh", False)
    markers = df.MeshFunction('int', mesh)
    hdf.read(markers, mesh_folder + "subdomains")

    #---------------------------------------------------------------------------
    # Periodic boundary condition map
    #---------------------------------------------------------------------------
    class PeriodicBoundary(df.SubDomain):

        def inside(self, x, on_boundary):
            # return True if on left or bottom boundary AND NOT on one of the
            # two corners (0, 1) and (1, 0)
            return bool((df.near(x[0], 0) or df.near(x[1], 0)) and
                    (not ((df.near(x[0], 0) and df.near(x[1], 1)) or
                            (df.near(x[0], 1) and df.near(x[1], 0)))) and on_boundary)

        def map(self, x, y):

            # if on right upper corner copy it to left down corner
            if df.near(x[0], 1) and df.near(x[1], 1):
                y[0] = x[0] - 1.
                y[1] = x[1] - 1.

            # if on right boundary copy it to the left boundary
            elif df.near(x[0], 1):
                y[0] = x[0] - 1.
                y[1] = x[1]

            # if on upper boundary copy it to the lower boundary
            else:
                y[0] = x[0]
                y[1] = x[1] - 1.

    # Function space P1 with periodic boundary conditions
    V = df.FunctionSpace(mesh, "P", 1, constrained_domain = PeriodicBoundary())

    #---------------------------------------------------------------------------
    # Permittivity coefficient for previously defined subdomains
    #---------------------------------------------------------------------------
    class Coeff(df.Expression):

        def __init__(self, mesh, **kwargs):
            self.markers = markers

        def eval_cell(self, values, x, cell):
            if markers[cell.index] == 1:
                values[0] = inner_permittivity
            else:
                values[0] = outer_permittivity

    # Interpolation to zeroth order polynomial
    permittivity = Coeff(mesh, degree = 0)

    # Weak formulation
    #---------------------------------------------------------------------------
    f1 = df.TrialFunction(V); f2 = df.TrialFunction(V)
    v1 = df.TestFunction(V);  v2 = df.TestFunction(V)

    # Variational form for the first corrector (f1)
    a1 = permittivity * df.dot(df.grad(f1), df.grad(v1)) * df.dx
    L1 = - df.Dx(v1, 0) * permittivity * df.dx

    # Variational form for the second corrector (f2)
    a2 = permittivity * df.dot(df.grad(f2), df.grad(v2)) * df.dx
    L2 = - df.Dx(v2, 1) * permittivity * df.dx

    # System assembly
    #---------------------------------------------------------------------------
    # Solution Functions (Correctors)
    f1 = df.Function(V); F1 = f1.vector()
    f2 = df.Function(V); F2 = f2.vector()

    # Assemble LHS, RHS and solve the system A*F=b
    A1 = df.assemble(a1);  b1 = df.assemble(L1);  df.solve(A1, F1, b1)
    A2 = df.assemble(a2);  b2 = df.assemble(L2);  df.solve(A2, F2, b2)

    # Effective permittivity calculation
    #---------------------------------------------------------------------------
    effective_11 = df.assemble(permittivity * (df.Dx(f1, 0) + 1) * df.dx)
    effective_12 = 0
    effective_21 = 0
    effective_22 = df.assemble(permittivity  * (df.Dx(f2, 1) + 1) * df.dx)

    # Write calculated effective parameters to the file effective (2x2 matrix)
    #---------------------------------------------------------------------------
    ofile = open('effective', 'w')

    ofile.write('%12.6e %12.6e \n' %(effective_11, effective_12))
    ofile.write('%12.6e %12.6e \n' %(effective_21, effective_22))

    ofile.close()

    return mesh, f1, f2
#-------------------------------------------------------------------------------

def save_PVD(output_folder, output_name, u):
    """Save function u and coresponding mesh to .pvd file format"""

    # Input Variables:
        # output_folder: folder where .h5 file will be store
        # mesh_name: name of mesh containig .h5 file
        # u: function that will be saved in 'output_folder/output_name.pvd'

    vtkfile = df.File(output_folder + output_name + '.pvd')
    vtkfile << u

    return 0
#-------------------------------------------------------------------------------


def save_HDF5(output_folder, mesh, mesh_name, Field, u):
    """Save function u and coresponding mesh to .h5 file format"""

    # Input Variables:
        # output_folder: folder where .h5 file will be store, format folder/
        # mesh: mesh keeping variable
        # mesh_name: name of .h5 file in which mesh is stored
        # Field:
        # u: function that will be saved in 'output_folder/Field_mesh_name.h5' file

    hdf = df.HDF5File(mesh.mpi_comm(), output_folder + Field + mesh_name + '.h5', 'w')
    hdf.write(mesh, output_folder + 'mesh')
    hdf.write(u, output_folder + 'solution');
    hdf.close()

    return 0
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Main part
#-------------------------------------------------------------------------------


if __name__ == '__main__':

    # Function call: python3 Yproblem.py mesh_folder mesh_name output_folder
    # ie. python3 Yproblem.py mesh hexagonal results

    # Input parameters
    mesh_folder = sys.argv[1]
    mesh_name = sys.argv[2]
    output_folder = sys.argv[3]

    # Domain defining permittivity coefficients
    inner_permittivity = 1
    outer_permittivity = 11.7

    # Call Y_solver_2D
    mesh, F1, F2 = Y_solver_2D(mesh_folder, mesh_name, inner_permittivity, outer_permittivity)

    # Save Correctors to XDMF File
    xdmffile_F1 = df.XDMFFile('results/XDMF/F1_' + mesh_name + '.xdmf');
    xdmffile_F1.write(F1)

    xdmffile_F2 = df.XDMFFile('results/XDMF/F2_' + mesh_name + '.xdmf');
    xdmffile_F2.write(F2)

    # Output files in PVD (for ParaView) and HDF5 (for later processing) format
    save_PVD(output_folder + '/PVD/', 'F1_' + mesh_name, F1);
    save_PVD(output_folder + '/PVD/', 'F2_' + mesh_name, F2)

    save_HDF5(output_folder +'/XDMF/', mesh, mesh_name, 'F1_', F1)
    save_HDF5(output_folder +'/XDMF/', mesh, mesh_name, 'F2_', F2)
