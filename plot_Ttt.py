#sample python3 file to plot Ttt at a given timestep
import numpy as np
import matplotlib.pyplot as plt

#load data
x   = np.loadtxt('datafiles/x.txt')
t   = np.loadtxt('datafiles/t.txt')
Ttt = np.loadtxt('datafiles/Ttt.txt')

#arbitrarily pick a timestep
timestep = 55

#plot data; notice that because the data is formatted
#such that a row is at a given timestep, plotting at
#a fixed time is the same as plotting a given row
plt.plot(x, Ttt[timestep, :])
plt.title(f'Ttt at t = {t[timestep]}')
plt.xlabel('x')
plt.ylabel('Ttt')
plt.show()
