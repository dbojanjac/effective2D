# Effective2D

FEM based solver for calculating effective permittivity of a heterogeneous structure made of inner inclusions with inner_permittivity and outside material matrix with outer_permittivity.

Using FEniCS 2017.2.0

## Using in docker:

```
$ docker run -ti --rm -v $(pwd):/home/fenics/shared quay.io/fenicsproject/stable:2017.2.0
$ cd shared
```

To generate effective params for hexagonal mesh you should firstly generate
the mesh:
```
$ cd mesh/
$ make hexagonal.msh
```
After generating the mesh, running demo script is as simple as
```
$ python3 demo.py
[Y] Calling FEniCS meshconvert util
Converting from Gmsh format (.msh, .gmsh) to DOLFIN XML format
Expecting 148940 vertices
Found all vertices
Expecting 296538 cells
Found all cells
Conversion done
Geting effective parameters in 38.69988560676575 sec
[[ 5.54646606  0.        ]
 [ 0.          5.84589027]]
```
