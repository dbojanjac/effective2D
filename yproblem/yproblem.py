# Using FEniCS 2017.2.0
import os
import uuid

import numpy as np
import dolfin as df
import dolfin_utils
import dolfin_utils.meshconvert
from dolfin_utils.meshconvert import meshconvert

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
    def __init__(self, mesh_filename, subdomains):
        self.mesh_filename = mesh_filename
        if mesh_filename.endswith(".h5"):
            self._parse_hdf5()
        else:
            if df.MPI.rank(df.mpi_comm_world()) == 0:
                print ("[Y] Calling FEniCS meshconvert util")
                self._convert_mesh()
        self.permittivity = Coeff(self.mesh_markers, subdomains, degree=0)

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

    def _convert_mesh(self):
        mesh_id = str(uuid.uuid4())
        mesh_xml = "/tmp/" + mesh_id + ".xml"
        meshconvert.convert2xml(self.mesh_filename, mesh_xml)

        self.mesh = df.Mesh(mesh_xml)
        without_xml = os.path.splitext(mesh_xml)[0]
        self.mesh_markers = df.MeshFunction("size_t", self.mesh, without_xml + "_physical_region.xml");

    def solve(self, degree=1):
        # Interpolation to zeroth order polynomial
        V = df.FunctionSpace(self.mesh, "P", degree,
                             constrained_domain = PeriodicBoundary())

        #--------------------------------------------------------------------
        #                       Weak formulation
        #--------------------------------------------------------------------
        u1 = df.TrialFunction(V); u2 = df.TrialFunction(V)
        v1 = df.TestFunction(V);  v2 = df.TestFunction(V)

        # Variational form for the first corrector (u1)
        a1 = self.permittivity * df.dot(df.grad(u1), df.grad(v1)) * df.dx
        L1 = - df.Dx(v1, 0) * self.permittivity * df.dx

        # Variational form for the second corrector (u2)
        a2 = self.permittivity * df.dot(df.grad(u2), df.grad(v2)) * df.dx
        L2 = - df.Dx(v2, 1) * self.permittivity * df.dx

        #--------------------------------------------------------------------
        #                       System assembly
        #--------------------------------------------------------------------
        # Solution Functions (Correctors)
        f1 = df.Function(V); F1 = f1.vector()
        f2 = df.Function(V); F2 = f2.vector()

        # Assemble LHS, RHS and solve the system A*F=b
        A1 = df.assemble(a1);  b1 = df.assemble(L1);  df.solve(A1, F1, b1, 'gmres', 'hypre_amg')
        A2 = df.assemble(a2);  b2 = df.assemble(L2);  df.solve(A2, F2, b2, 'gmres', 'hypre_amg')

        #--------------------------------------------------------------------
        #               Effective permittivity calculation
        #--------------------------------------------------------------------
        permittivity = df.interpolate(self.permittivity, V)
        d1 = df.assemble(permittivity * (df.Dx(f1, 0) + 1) * df.dx)
        d2 = df.assemble(permittivity * (df.Dx(f2, 1) + 1) * df.dx)

        return np.array([[d1, 0], [0, d2]]), (f1, f2)
