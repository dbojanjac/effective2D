meshes = circle.msh cross.msh ellipse.msh hexagonal.msh rectangle.msh square.msh

all: $(meshes)

%.msh: %.geo
	gmsh -2 $<

clean:
	$(RM) -f $(meshes)
