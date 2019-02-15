# FEM based solver for calculating effective permittivity of a heterogeneous
# structure made of inclusions with inner permittivity (markers = 1) and
# outside matrix with outer_permittivity (markers = 2)

# input = unit_cell mesh with subdomain markers
# output = txt file with 2x2 matrix of effective permittivity
# output = txt file with 2x2 matrix of effective permittivity

from dolfin import *

def Y_solver(Mesh_Folder, Mesh_Name):

    # Read Mesh from Mesh_Folder/Mesh
    #---------------------------------------------------------------------------
    Mesh_Folder = Mesh_Folder + '/'

    mesh = Mesh()
    hdf = HDF5File(mesh.mpi_comm(), Mesh_Folder + Mesh_Name, 'r')

    hdf.read(mesh, Mesh_Folder + "mesh", False)
    markers = MeshFunction('int', mesh)
    hdf.read(markers, Mesh_Folder + "subdomains")

    #---------------------------------------------------------------------------
    # Periodic boundary condition map
    #---------------------------------------------------------------------------
    class PeriodicBoundary(SubDomain):

        # Left boundary is "target domain" G
        def inside(self, x, on_boundary):
            # return True if on left or bottom boundary AND NOT on one of the
            # two corners (0, 1) and (1, 0)
            return bool((near(x[0], 0) or near(x[1], 0)) and
                    (not ((near(x[0], 0) and near(x[1], 1)) or
                            (near(x[0], 1) and near(x[1], 0)))) and on_boundary)

        def map(self, x, y):
            # if on right upper corner copy it to left down corner
            if near(x[0], 1) and near(x[1], 1):
                y[0] = x[0] - 1.
                y[1] = x[1] - 1.
            elif near(x[0], 1): # if on right boundary
                y[0] = x[0] - 1.
                y[1] = x[1]
            else:   # if on upper boundary
                y[0] = x[0]
                y[1] = x[1] - 1.

    # Function space with periodic boundary condition
    V = FunctionSpace(mesh, "P", 1, constrained_domain=PeriodicBoundary())

    #---------------------------------------------------------------------------
    # Permittivity coefficient previously defined subdomains
    #---------------------------------------------------------------------------
    class Coeff(Expression):
        def __init__(self, mesh, **kwargs):
            self.markers = markers
        def eval_cell(self, values, x, cell):
            if markers[cell.index] == 1:
                values[0] = inner_permittivity
            else:
                values[0] = outer_permittivity

    permittivity = Coeff(mesh, degree=0)

    # Weak formulation
    #---------------------------------------------------------------------------
    f1 = TrialFunction(V); f2 = TrialFunction(V)
    v1 = TestFunction(V);  v2 = TestFunction(V)

    # Variational form for the first corrector (f1)
    a1 = permittivity * dot(grad(f1), grad(v1)) * dx
    L1 = -Dx(v1, 0) * permittivity * dx

    # Variational form for the second corrector (f2)
    a2 = permittivity * dot(grad(f2), grad(v2)) * dx
    L2 = -Dx(v2, 1) * permittivity * dx

    # System assembly
    #---------------------------------------------------------------------------
    # Solution Functions
    f1 = Function(V); F1 = f1.vector()
    f2 = Function(V); F2 = f2.vector()

    # Assemble RHS, LHS and solve the system
    A1 = assemble(a1);  b1 = assemble(L1);  solve(A1, F1, b1)

    A2 = assemble(a2);  b2 = assemble(L2);  solve(A2, F2, b2)

    # Effective permittivity calculation
    #---------------------------------------------------------------------------

    A11 = assemble(permittivity * (Dx(f1, 0) + 1) * dx)
    A12 = 0
    A21 = 0
    A22 = assemble(permittivity  * (Dx(f2, 1) + 1) * dx)

    # Write calculated effective parameters to the file permittivity (2x2 matrix)
    #---------------------------------------------------------------------------
    ofile = open('effective', 'w')

    ofile.write('%12.6e %12.6e \n' %(A11, A12))
    ofile.write('%12.6e %12.6e \n' %(A21, A22))

    ofile.close()

    return f1, f2


#-------------------------------------------------------------------------------
# Main part
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # Domain defining Coefficients
    inner_permittivity = 1
    outer_permittivity = 11.7

    # Mesh Defining Parameters
    Mesh_Folder = 'mesh'; Mesh_Name = 'Ymesh' + '.h5'

    # Call Y_solver
    F1, F2 = Y_solver(Mesh_Folder, Mesh_Name)

    # Save Correctors to XDMF File
    xdmffile_F1 = XDMFFile('rezultati/XDMF/F1.xdmf');   xdmffile_F1.write(F1)
    xdmffile_F2 = XDMFFile('rezultati/XDMF/F2.xdmf');   xdmffile_F2.write(F2)
