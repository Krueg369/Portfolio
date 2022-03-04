from numba import jit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy import ndimage
import time
from statistics import mean

#Create initial state of spin up(+1) or down(-1) in NxN matrix
def initState(N):
    state = 2*np.random.randint(2, size = (N,N)) -1
    pad = np.pad(state, ((1,1),(1,1)), 'constant')

    return state, pad

#Calculate Energy of the Configuration using convolution, nearest neighbors algorithm
def Hamiltonian(lattice, H, J):
    energy = 0
    kernel = ndimage.generate_binary_structure(2,1)
    kernel[1][1] = 0
    ij = ndimage.convolve(lattice, kernel, mode = 'constant', cval=0)
    internal = J*(lattice * ij)
    external = sum(map(sum, lattice))*H
    energy = -1 * sum(map(sum, internal)) - external
    return energy

@jit(nopython=True)
def Metropolis(lattice, beta, energy, J, N, H):
    pad = lattice
    spinNet = pad.sum()
    realNet = spinNet

    #print(pad)

    net_spin = [spinNet/(N**2)]
    net_energy = [energy]
    i = 0

    while i < 500000:
        #Choose random xy coordinate and do a spin flip
        x = np.random.randint(1,N+1)
        y = np.random.randint(1,N+1)
        init_spin = pad[x,y]
        flip_spin = init_spin * -1
        #print(flip_spin)
        #Change in enery Calculation
        Ei = 0
        Ef = 0

        summ = pad[(x+1), y] + pad[x,(y+1)] + pad[(x-1), y] + pad[x,(y-1)]
        Ei = -J*summ*init_spin - init_spin*H
        Ef = -J*summ*flip_spin - flip_spin*H

        deltaE = Ef - Ei

        prob = np.random.random()
        #State change in accordance to probability
        if(deltaE > 0)*(prob < np.exp((-beta*deltaE))):
            pad[x,y] = flip_spin
            energy += deltaE
        elif deltaE <= 0:
            pad[x,y] = flip_spin
            energy += deltaE

        #print(pad)
        spinNet = abs(pad.sum())
        realNet = pad.sum()
        net_spin.append((realNet/N**2))
        net_energy.append((energy))
        i = i + 1

    return i, net_spin, net_energy

def Partitionfx(E,T):
    Z = np.exp(-E/T)
    return Z

def Magnetization(lattice, BJ, energy, J, N, H):
    ms = np.zeros(len(BJ))
    mcsteps = np.zeros(len(BJ))
    for i,bj in enumerate(BJ):
        steps, spins, energies = Metropolis(lattice, bj, energy, J, N, H)
        ms[i] = np.mean(spins)
        mcsteps[i] = steps
        if i % 10 == 0:
            print(i)
    return ms, mcsteps

#______________________________________________________________________________#
" MAIN "

print("Iteration Counter: ")
start = time.time()

#Create grid for lattice
N = 10 # Size of lattice
T = 600 # Temperature
H1 = 0
H2 = 1  # strength of external magnetic field
J = 1  # spin-spin interaction
betaJ = 1/T
temps = np.linspace(50,5000,40)

#BJs = 1/temps
BJs = np.linspace(0.1,2,150)
switch1 = 0
switch2 = 1

# Call functions to generate initial state and cal. energy
state, spad = initState(N)
energy1 = Hamiltonian(state, H1, J)
energy2 = Hamiltonian(state, H2, J)



if switch1 == 1:
    steps, spins, energies = Metropolis(spad, betaJ, energy, J, N, H)
    plt.figure(1)
    plt.plot(spins)
    plt.ylim(-1.25,1.25)
    plt.xlabel("mcSteps")
    plt.ylabel("Average Spin")

    plt.figure(2)
    plt.xlabel("mcSteps")
    plt.ylabel("Energy")
    plt.plot(energies)

if switch2 == 1:
    ms1, mcsteps1 = Magnetization(spad, BJs, energy1, J, N, H1)
    ms2, mcsteps2 = Magnetization(spad, BJs, energy2, J, N, H2)

    plt.figure(3)
    plt.plot(1/BJs, abs(ms1), "+", color = 'black')
    plt.plot(1/BJs, abs(ms2), "o", color = 'red')
    plt.ylim(-1.25,1.25)
    plt.xlabel("T")
    plt.ylabel("|M|")
    plt.title("Absolute Value of Magnetization Vs. Temperature")
    median = mlines.Line2D([], [], color='black', marker='+', markersize=10, label='H = 0')
    function = mlines.Line2D([], [], color='red', marker = 'o', label='H = 1')
    plt.legend(handles=[median, function])

end = time.time()
total = end - start
print("Time to run:", total)
plt.show()
