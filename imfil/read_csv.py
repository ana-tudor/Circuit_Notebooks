import pandas as pd
import numpy

df = pd.read_csv("imfil_1.57rad_out2.csv")
print(df.columns)
initamp1 = df.initamp1
initamp2 = df.initamp2
finalamp1 = df.finalamp1
finalamp2 = df.initamp2

y = b['1.57080rad Optimal Amplitudes']
halfpi_opt_amplitudes_neg = finalamp1
halfpi_opt_amplitudes_pos = finalamp2
halfpi_opt_energies = b['1.57080rad Optimal Energies']

y = b['1.57080rad Initial Amplitudes']
halfpi_init_amplitudes_neg = initamp1
halfpi_init_amplitudes_pos = initamp2
halfpi_init_energies = b['1.57080rad Initial Energies']

fig = pylab.figure(figsize=(20,15))
pylab.subplot(3,2,1)
pylab.plot(halfpi_init_amplitudes_neg, halfpi_init_energies, 'ro')
pylab.xlabel("Initial amplitude at angle pi/2")
pylab.ylabel("Initial energy")
pylab.title("First amplitude element")

pylab.subplot(3,2,2)
pylab.plot(halfpi_init_amplitudes_pos, halfpi_init_energies, 'ro')
pylab.xlabel("Initial amplitude at angle pi/2")
pylab.ylabel("Initial energy")
pylab.title("Second amplitude element")
