import sys
import numpy
from cantera import *
from matplotlib.pylab import *
import csv

# Mechanism used for process
gas = Solution('gri30.cti')

# Number of time steps and the time step size
nt = 100000
dt = 10 ** (-6)

# Amount of iterations
n_points = 11

# Simulation parameters
T_min = 1100
T_max = 2400
P_min = 101325
P_max = 3 * 101325
fi_min = 0.4
fi_max = 3.4

# Storage space
Ti = numpy.zeros(n_points, 'd')
Pi = numpy.zeros(n_points, 'd')
fi = numpy.zeros(n_points, 'd')
tim = numpy.zeros(nt, 'd')
temp = numpy.zeros(nt, 'd')
d_temp = numpy.zeros(nt-1, 'd')
Autoignition_csv = numpy.zeros(n_points**2 * n_points, 'd')
FinalTemp_csv = numpy.zeros(n_points**2 * n_points, 'd')

# Storage space for plots
Autoignition_Temp = numpy.zeros(n_points, 'd')
FinalTemp_Temp = numpy.zeros(n_points, 'd')
Autoignition_Pressure = numpy.zeros(n_points, 'd')
FinalTemp_Pressure = numpy.zeros(n_points, 'd')
Autoignition_Fi = numpy.zeros(n_points, 'd')
FinalTemp_Fi = numpy.zeros(n_points, 'd')

s = 0


for k in range(n_points):
    Ti[k] = T_min + (T_max - T_min) * k / (n_points - 1)

    for l in range(n_points):
        Pi[l] = P_min + (P_max - P_min) * l / (n_points - 1)

        for m in range(n_points):
            fi[m] = fi_min + (fi_max - fi_min) * m / (n_points - 1)
            # initial temperature, pressure and stoichiometry
            n2 = float(0.79 * 0.2857 / (fi[m] * 0.060024))  # amount of N2 moles
            o2 = float(0.21 * 0.2857 / (fi[m] * 0.060024))   # amount of O2 moles
            X = 'C2H6:0.2857 N2:{0} O2:{1}'.format(str(n2), str(o2))
            gas.TPX = Ti[k], Pi[l], X
            r = IdealGasReactor(gas)  # batch reactor
            sim = ReactorNet([r])  # reactor network consisting of single batch reactor

            # initial simulation time
            time = 0.0

            # Simulation
            for n in range(nt):  # loop for nt times steps of dt seconds
                time += dt
                sim.advance(time)
                tim[n] = time
                temp[n] = r.T

            # catching the autoignition timing
            Dtmax = [0, 0.0]
            for n in range(nt - 1):
                d_temp[n] = (temp[n + 1] - temp[n]) / dt
                if d_temp[n] > Dtmax[1]:
                    Dtmax[0] = n
                    Dtmax[1] = d_temp[n]
            Autoignition = (tim[Dtmax[0]] + tim[Dtmax[0] + 1]) / 2.
            # print 'For T=' + str(Ti[k]) + 'K P=' + str(Pi[l]) + 'Pa and fi=' + str(fi[m]) + ', Autoignition time=' + str(Autoignition) + '(s)'
            Autoignition_csv[s] = Autoignition * 1000  # [s] to [ms]
            FinalTemp_csv[s] = temp[nt - 1]
            # print 'FinalTemp =' + str(FinalTemp_cas[s])    # Prints final temperature
            s += 1


            # Saving data for plots
            if Pi[l] == 101325 and fi[m] == 1:
                print 'For T=' + str(Ti[k]) + 'K P=' + str(Pi[l]) + 'Pa and fi=' + str(
                    fi[m]) + ', Autoignition time=' + str(Autoignition) + '(s)'
                print 'FinalTemp =' + str(FinalTemp_csv[s-1])
                FinalTemp_Temp[k] = temp[nt-1]
                Autoignition_Temp[k] = Autoignition*1000

            if Ti[k] == 1100 and fi[m] == 1:
                print 'For T=' + str(Ti[k]) + 'K P=' + str(Pi[l]) + 'Pa and fi=' + str(
                    fi[m]) + ', Autoignition time=' + str(Autoignition) + '(s)'
                print 'FinalTemp =' + str(FinalTemp_csv[s-1])
                FinalTemp_Pressure[l] = temp[nt-1]
                Autoignition_Pressure[l] = Autoignition*1000

            if Pi[l] == 101325.0 and Ti[k] == 1100.0:
                print 'For T=' + str(Ti[k]) + 'K P=' + str(Pi[l]) + 'Pa and fi=' + str(
                    fi[m]) + ', Autoignition time=' + str(Autoignition) + '(s)'
                print 'FinalTemp =' + str(FinalTemp_csv[s-1])
                FinalTemp_Fi[m] = temp[nt-1]
                Autoignition_Fi[m] = Autoignition*1000

# Saving results to .csv file
b = 0
csv_file = 'Autoignition of Ethane and Air mixture.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Initial temperature', 'Pressure', 'Fi', 'Autoignition time', 'Final Temperature'])
    for k in range(n_points):
        writer.writerow([Ti[k]])
        for l in range(n_points):
            writer.writerow(['', Pi[l]])
            for m in range(n_points):
                writer.writerow(['', '', fi[m], Autoignition_csv[b], FinalTemp_csv[b]])
                b += 1




# PLOTS

# Autoignition time(Fi)

plot(fi, Autoignition_Fi, '-', color='red')
xlabel(r'Fi', fontsize=20)
ylabel("Autoignition [ms]")
title(r'Autoignition of $C_{2}H_{6}$ and Air at P=1atm and T=1100K', fontsize=20, horizontalalignment='center')
axis([0.2, 3.5, 35.0, 55.0])
grid()
savefig('AutoignitionTime(Fi).png', bbox_inches='tight')

# Final Temperature(Fi)

plot(fi, FinalTemp_Fi, '-', color='red')
xlabel(r'Fi', fontsize=20)
ylabel("Final Temp[K]")
title(r'Autoignition of $C_{2}H_{6}$ and Air at P=1atm and T=1100K', fontsize=20, horizontalalignment='center')
axis([0.2, 3.5, 1500.0, 3000.0])
grid()
savefig('FinalTemp(Fi).png', bbox_inches='tight')

# Autoignition time(InitialTemp)

plot(Ti, Autoignition_Temp, '-', color='red')
xlabel(r'Temp [K]', fontsize=20)
ylabel("Autoignition time [ms]")
title(r'Autoignition of $C_{2}H_{6}$ and Air at P=1atm and $\Phi$=1', fontsize=20, horizontalalignment='center')
axis([1050, 2450, 0.0, 50.0])
grid()
savefig('Autoignition(InitialTemperature).png', bbox_inches='tight')

# Final Temperature(InitialTemp)

plot(Ti, FinalTemp_Temp, '-', color='red')
xlabel(r'Temp [K]', fontsize=20)
ylabel("Final temperature [K]")
title(r'Autoignition of $C_{2}H_{6}$ and Air at P=1atm and $\Phi$=1', fontsize=20, horizontalalignment='center')
axis([1050.0, 2450.0, 2750.0, 3200.0])
grid()
savefig('FinalTemp(Temp).png', bbox_inches='tight')

# Autoignition time(pressure)

plot(Pi/1000, Autoignition_Pressure, '-', color='red')
xlabel(r'Pressure [kPa]', fontsize=20)
ylabel("Autoignition time [ms]")
title(r'Autoignition of $C_{2}H_{6}$ and Air at T=1100K and $\Phi$=1', fontsize=20, horizontalalignment='center')
axis([100.0, 310.0, 15.0, 50.0])
grid()
savefig('AutoignitionTime(Pressure).png', bbox_inches='tight')

# Final Temp(pressure)

plot(Pi/1000, FinalTemp_Pressure, '-', color='red')
xlabel(r'Pressure [kPa]', fontsize=20)
ylabel("Final Temperature[K]")
title(r'Autoignition of $C_{2}H_{6}$ and T=1100K and Air at $\Phi$=1', fontsize=20, horizontalalignment='center')
axis([100.0, 310.0, 2700.0, 3000.0])
grid()
savefig('FinalTemp(Pressure).png', bbox_inches='tight')

