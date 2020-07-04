import os
import time
import numpy as np
import yproblem

def save_effective(output_folder, effective):
    with open(output_folder + "/effective.npy", "wb") as f:
        np.save(f, effective)

mesh = "mesh/hexagonal.msh"
subdomains = {1: 1, 2: 11.8}
output_folder = "results"
y = yproblem.Yproblem(mesh, subdomains)
start = time.time()
effective, (f1, f2) = y.solve()
elapsed_time = time.time() - start

print ("Geting effective parameters in {} sec".format(elapsed_time))
print (effective)

os.makedirs(output_folder + "/plots", exist_ok=True)
os.makedirs(output_folder + "/PVD", exist_ok=True)
yproblem.save_field_plots(output_folder, f1, f2)
save_effective(output_folder, effective)
yproblem.save_pvd(output_folder, f1, f2)
