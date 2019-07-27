# Using FEniCS 2017.2.0
import petsc4py
import sys
import os
petsc4py.init(['-log_view'])
from petsc4py import PETSc
import dolfin as df


#---------------------------------------------------------------------------
#                   Periodic boundary condition map
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


class Coeff(df.Expression):
    def __init__(self, markers, subdomains, **kvargs):
        self.markers = markers
        self.subdomains = subdomains

    def eval_cell(self, values, x, cell):
        values[0] = self.subdomains[self.markers[cell.index]]

class Yproblem:
    #TODO: pass dict to constructor: subdomain_id : permittivity
    def __init__(self, mesh_filename, subdomains):
        self.mesh_filename = mesh_filename
        if mesh_filename.endswith(".h5"):
            self._parse_hdf5()
        elif mesh_filename.endswith(".msh"):
            self._parse_gmsh()
        else:
            raise TypeError("Mesh {} has unsupported \
                    type".format(mesh_filename))
        self.subdomains = subdomains
        self.effective = []

    def _parse_hdf5(self):
        # generate relative mesh_filename, mesh_folder and mesh_name
        self._mesh_filename = os.path.relpath(self.mesh_filename)
        self._mesh_folder = "/".join(self.mesh_filename.split("/")[:-1])
        self._mesh_name = self.mesh_filename.split("/")[-1]

        self.mesh = df.Mesh()
        hdf = df.HDF5File(self.mesh.mpi_comm(), self._mesh_filename, 'r')

        hdf.read(self.mesh, self._mesh_folder + "/mesh", False)
        self.mesh_markers = df.MeshFunction('int', self.mesh)
        hdf.read(self.mesh_markers, self._mesh_folder + "/subdomains")

        return self.mesh

    def _parse_gmsh(self):
        pass

    def get_effective(self, V):
        # Interpolation to zeroth order polynomial
        permittivity = Coeff(self.mesh_markers,
                             self.subdomains, degree=0)

        #--------------------------------------------------------------------
        #                       Weak formulation
        #--------------------------------------------------------------------
        f1 = df.TrialFunction(V); f2 = df.TrialFunction(V)
        v1 = df.TestFunction(V);  v2 = df.TestFunction(V)

        # Variational form for the first corrector (f1)
        a1 = permittivity * df.dot(df.grad(f1), df.grad(v1)) * df.dx
        L1 = - df.Dx(v1, 0) * permittivity * df.dx

        # Variational form for the second corrector (f2)
        a2 = permittivity * df.dot(df.grad(f2), df.grad(v2)) * df.dx
        L2 = - df.Dx(v2, 1) * permittivity * df.dx

        #--------------------------------------------------------------------
        #                       System assembly
        #--------------------------------------------------------------------
        # Solution Functions (Correctors)
        self.f1 = df.Function(V); F1 = self.f1.vector()
        self.f2 = df.Function(V); F2 = self.f2.vector()

        # Assemble LHS, RHS and solve the system A*F=b
        A1 = df.assemble(a1);  b1 = df.assemble(L1);  df.solve(A1, F1, b1)
        A2 = df.assemble(a2);  b2 = df.assemble(L2);  df.solve(A2, F2, b2)

        #--------------------------------------------------------------------
        #               Effective permittivity calculation
        #--------------------------------------------------------------------
        permittivity = df.interpolate(permittivity, V)

        self.effective.append(df.assemble(
                                permittivity * (df.Dx(self.f1, 0) + 1)
                                * df.dx))
        self.effective.append(0)
        self.effective.append(0)
        self.effective.append(df.assemble(
                                permittivity * (df.Dx(self.f2, 1) + 1)
                                * df.dx))

        return self.effective


def Y_solver_2D(mesh_folder, mesh_name, inner_permittivity, outer_permittivity):


    # Function space P1 with periodic boundary conditions
    V = df.FunctionSpace(mesh, "P", 1, constrained_domain = PeriodicBoundary())


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
    permittivity = df.interpolate(permittivity, V)

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

    return mesh, markers, f1, f2
#-------------------------------------------------------------------------------


def save_PVD(output_folder, output_name, output_file):
    """Save function output_file to .pvd file format"""
    # Input Variables:
        # output_folder: folder where .h5 file will be store, format folder/
        # mesh_name: name of .h5 file in which mesh is stored
        # output_file: function that will be saved in
        #   'output_folder/Field_mesh_name.pvd' file

    # Output:
        # function output_file will be saved to output_folder/output_file.pvd

    vtkfile = df.File(output_folder + output_name + '.pvd')
    vtkfile << output_file

    return 0
#-------------------------------------------------------------------------------


def save_HDF5(output_folder, output_name, mesh, markers, output_file):
    """Save function output_file and coresponding mesh to .h5 file format"""
    # Input Variables:
        # output_folder: folder where .h5 file will be store
        # output_name: name of the output file
        # mesh: FEniCS mesh function
        # markers: FEniCS subdomains function
        # output_file: function that will be saved in output

    # Output:
        # function output_file will be saved to output_folder/output_file.h5

    # Open file output_folder/output_name.h5 for writing
    hdf = df.HDF5File(mesh.mpi_comm(), output_folder + output_name + '.h5', 'w')

    hdf.write(mesh, output_folder + 'mesh')
    hdf.write(markers, output_folder + 'subdomains')
    hdf.write(output_file, output_folder + 'solution');

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
    mesh, markers, F1, F2 = Y_solver_2D(mesh_folder, mesh_name, \
        inner_permittivity, outer_permittivity)


    # Output files in PVD (for ParaView) and HDF5 (for later processing) format
    save_PVD(output_folder + '/PVD/', 'F1_' + mesh_name, F1)
    save_PVD(output_folder + '/PVD/', 'F2_' + mesh_name, F2)

    save_HDF5(output_folder +'/XDMF/', 'F1_' + mesh_name, mesh, markers, F1)
    save_HDF5(output_folder +'/XDMF/', 'F2_' + mesh_name, mesh, markers, F2)
