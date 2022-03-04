# Main Program
import numpy as np
import math
import random
import matplotlib.pyplot as plt
from interpf import interpfx
from interpf import test1
from interpf import test2

#declare variables
xmin = -7
xmax = 7
nx = 100
x = np.linspace(xmin, xmax, nx)
xmin2 = -5
xmax2 = 5
nx2 = 100
x = np.linspace(xmin, xmax, nx)
x2 = np.linspace(xmin2, xmax2, nx2)
xdata = np.array([-2*math.pi,-3/2* math.pi, -1*math.pi, -1/2* math.pi, 0,1/2*math.pi, math.pi, 3/2*math.pi,2* math.pi])
xdata2 = np.array([-4,-2,0,2,4])

#call functions
fdata = test1(xdata)
fdata2 = test1(xdata2)
exact = test1(x)
exact2 = test1(x2)
p = interpfx(x, xdata, fdata)
p2 = interpfx(x2, xdata2, fdata2)
error = test2(p, exact)
error2 = test2(p2, exact2)
print("Percent error b/w interpfx and exact (not really): {:%}".format(error))
print("There is an option commented out in interpolation.py to print plots seperate.")
#Option to save as larger single images to seperate files
"""
plt.figure(1)
plt.plot(xdata, fdata, "s", color = "black")
plt.plot(x, p, color = "red")
plt.plot(x, exact, "r--", color = "black")
plt.title("8th Order Interpolation")
plt.xlabel('X')
plt.ylabel("p/ fdata/ exact data")
plt.savefig('interp1.pdf')

plt.figure(2)
plt.plot(xdata2, fdata2, "s", color = "black")
plt.plot(x2, p2, color = "red")
plt.plot(x2, exact2, "r--", color = "black")
plt.title("4th Order Interpolation")
plt.xlabel('X')
plt.ylabel("p/ fdata/ exact data")
plt.savefig('interp2.pdf')
"""
fig, (ax1, ax2) = plt.subplots(1,2)
ax1.plot(x, p, color = "red")
ax1.plot(x, exact, "r--", color = "black")
ax1.plot(xdata, fdata, "s", color = "black")
ax1.set_title("8th Order Interp")
ax1.set(xlabel = "x/ xdata", ylabel = "p/ exact/ fdata")

ax2.plot(x2, p2, color = "red")
ax2.plot(x2, exact2, "r--", color = "black")
ax2.plot(xdata2, fdata2, "s", color = "black")
ax2.set_title("4th Order Interp")
plt.savefig("interp.pdf")
