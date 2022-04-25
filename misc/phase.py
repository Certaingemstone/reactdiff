# phase.py : For plotting points in the context of phase portraits
import csv
import matplotlib.pyplot as plt
import os

rk4out = r"C:\Users\jade2\source\repos\rk4testclient\rk4testclient"
fnames = [os.path.join(rk4out, f"{i}.csv") for i in range(1,5)]

plt.figure(figsize=(8, 6), dpi=120)
for name in fnames:
    xi = []
    chi = []
    with open(name, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(reader):
            print(row)
            if i % 5 == 0:
                xi.append(float(row[0]))
                chi.append(float(row[1]))
    plt.scatter(xi, chi)
    #plt.xlim(-0.1, 0.1)
    #plt.ylim(-0.1, 0.1)

plt.title("Time Evolution\n k1=k2=k3=1.0, S=1.0, initial data offsets (0.1, 0.1), dt=0.005")
plt.xlabel(r"$\xi$")
plt.ylabel(r"$\chi$")
plt.show()