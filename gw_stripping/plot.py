import matplotlib.pyplot as plt
import pandas as pd

def main():
    fd_dist_mas = pd.read_table("stripping_dist_mass.dat", delimiter=",", header=None)
    fd_rad = pd.read_table("stripping_rad.dat", delimiter=",", header=None)

    plt.plot(fd_rad[0], fd_rad[1], '.', label="Neutrino")
    plt.plot(fd_dist_mas[0], fd_dist_mas[4], '.', label="GW")
    plt.ylabel("log(L)")
    plt.xlabel("t (time), s")
    plt.yscale("log")
    plt.grid(True)
    plt.legend()


    _, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.plot(fd_dist_mas[0], fd_dist_mas[1], '.')
    ax1.set(ylabel="a (distanse between NSs), km")
    ax1.grid(True)

    ax2.plot(fd_dist_mas[0], fd_dist_mas[2], '.')
    ax2.set(ylabel="q")
    ax2.grid(True)

    ax3.plot(fd_dist_mas[0], fd_dist_mas[3], '.')
    ax3.grid(True)
    ax3.set(xlabel="t (time), s", ylabel="m2")

    plt.show()

if __name__ == "__main__":
    main()
