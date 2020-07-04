import dolfin as df
import matplotlib.pyplot as plt

def save_field_plots(output_folder, f1, f2):
    df.plot(f1)
    plt.savefig(output_folder + "/plots/f1.pdf")
    df.plot(f2)
    plt.savefig(output_folder + "/plots/f2.pdf")


def save_pvd(output_folder, f1, f2):
    f = df.File(output_folder + "/PVD/f1.pvd")
    f << f1

    f = df.File(output_folder + "/PVD/f2.pvd")
    f << f2
