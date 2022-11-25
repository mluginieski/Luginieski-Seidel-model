## IdcLSm: Id current of Luginieski-Seidel model
#
#  code created by: Msc. Marcos Luginieski
#          contact: mluginieski@ifsc.usp.com
#             date: 05/10/2022

### Modules ###
import os.path                          # To place the output on the right folder
import numpy as np                      # For easy vectorization
import matplotlib.pyplot as plt         # To plot graphs
import constants as cte                 # Physical constants module
from unit_converter import converter    # To convert the physical units to SI

## Input reading sequence and pre-calculation

values = []                                                  # Empty list
with open('IdcLSm.input', 'r') as fi:                        # Read the input file
    parameters = [list(map(str, i.split('\t'))) for i in fi] # Create a list of lists, where "i" is
                                                             #   the string readed from input file
for i in range(0, len(parameters)):
    parameters[i][1] = parameters[i][1]                      # Convert the [i][1] value of the list to a float number
    values.append(converter(parameters[i])[1])               # Unit converter module. The hole [i] list is given, but
                                                             #   just the [1] value is appended to the values list

kappa, T, Ci, D, W, L, Vg, Vt, mi_sat = values               # Definition of every variable with the list values

def criticalThickness(y):

    epson = kappa * cte.epsilon0           # Dielectric constant
    Tc    = T*y                            # Characteristc temperature
    xi    = (4*epson*cte.kb*Tc)/(cte.e*Ci) # Xi constant
    
    return xi/(Vg - Vt)

def Id_current(Dc, a, y, maxExpVd, lenExpData):

    ## Definition of the output folder and the output files
    folder         = 'Output/'                                                   # Output folder
    main_file_name = 'CalculatedData.csv'                                        # Main file name
    os.makedirs(os.path.join(os.path.dirname(__file__),folder), exist_ok = True) # Creates the folder Output/
    output_file    = os.path.join(folder,main_file_name)                         # Output file (Vds, Ids, Idsat, ...)
    data_file      = os.path.join(folder,'DATA_%s' % main_file_name)             # Data file: parameters in use, values calculated, ...

    fo = open(data_file, 'w') # Creates the data file

    fo.write('Variables in use: \n') # Print on data file

    for i in range(0, len(parameters)):
        fo.write(parameters[i][0])                    # Variable name
        fo.write(' \t= ')                             # tab
        fo.write(str(parameters[i][1]))               # Variable value (SI)
        fo.write('   {} \n'.format(parameters[i][2])) # Variable's SI unit

    fo.write('Dc = {} m\n'.format(Dc)) # Fitted Dc
    fo.write('a = {}  .\n'.format(a))  # Fitted alfa
    fo.write('y = {}  .\n'.format(y))  # Fitted gamma
    
    ## Constants calculation
    epson = kappa * cte.epsilon0           # Dielectric constant
    Tc    = T*y                            # Characteristc temperature
    xi    = (4*epson*cte.kb*Tc)/(cte.e*Ci) # Xi constant
    Vtr   = xi/D                           # Transition voltage
    Vl    = Vg - Vt - Vtr                  # V'' voltage
    lc    = xi/(Vg - Vt)                   # Critical thickness
    
    VDS   = np.linspace(0,maxExpVd,lenExpData) # Vds values range
    IDS   = np.empty(lenExpData)          # Empty Ids vector
    dVds  = VDS[1] - VDS[0]        # Integration step
    
    # Data file printings
    fo.write('epson = %e\n' % epson)
    fo.write('Tc    = %e\n' %    Tc)
    fo.write('dVds  = %e\n' %  dVds)
    fo.write('\nConstants calculated:\n')
    fo.write('xi    = %e\n' %    xi)
    fo.write('Vtr   = %e\n' %   Vtr)
    fo.write("V''   = %e\n" %    Vl)
    fo.write('lc    = %e\n' %    lc)

    fo.close() # Close the data file
    ## End of reading sequence and precalculation
    
    ## Beginning of simulation

    i = 0 # For vectors use

    # OFET
    if Vl < 0:
        with open(output_file, 'w') as fo:
            print('#Vds(V);Isat(A);Ids(A);l(m);Vgs(V);mi', file=fo)
    
            # Simulation loop
            I_sat  = 0.0
    
            for Vds in VDS:
                l = xi/(Vg - Vt - Vds)

                # Accumulation regime
                if Vds < (Vg - Vt):
                    mi  = (((D-Dc)/D)**a)
                    Id = ((mi_sat*W*Ci)/L) * (((D-Dc)/D)**a) * ((Vg - Vt) - Vds*0.5)*Vds
                    IDS[i] = Id
                # Saturation regime
                if Vds >= (Vg - Vt):
                    mi  = (((D-Dc)/D)**a)
                    I_sat = ((mi*mi_sat*W*Ci)/(2*L)) * (Vg-Vt)**2
                    IDS[i] = I_sat
                
                i += 1
                print("%.3f;%e;%e;%e;%f;%e" % (Vds, I_sat, Id, l, Vg, mi), file=fo)

    # EGOFET
    if Vl >= 0.0:
        with open(output_file, 'w') as fo:
            print('#Vds(V);Isat(A);Ids(A);l(m);Vgs(V);mi', file=fo)
    
            # Simulation loop
            Id_i   = 0.0
            I_sat  = 0.0
    
            for Vds in VDS:
                # Accumulation regime
                if Vds < (Vg - Vt):

                    if Vds < Vl: #V(x) < V''
                        l = xi/(Vg - Vt - Vds)
                        
                    if Vds >= Vl: #V(x) >= V''
                        l = D

                    mi     = (((l-Dc)/D)**a) * ((l/D)**(1-y))
                    Id     = ((mi*mi_sat*W*Ci)/L) * (Vg - Vt - Vds) * dVds
                    Id_i  += Id
                    IDS[i] = Id_i

                # Saturation regime
                if Vds >= (Vg - Vt):
                    mi     = (((D-Dc)/D)**a)
                    I_sat  = Id_i + ((mi*mi_sat*W*Ci)/(2*L)) * Vtr**2
                    IDS[i] = I_sat

                i       += 1
                print("%.3f;%e;%e;%e;%f;%.8e" % (Vds, I_sat, Id_i, l, Vg, mi), file=fo)
        
    ## End of simulation
    if lc < Dc:
        return np.zeros(lenExpData) # In this case the drain current is set as zero
                                    # so the fitting code must seef for another set of Dc, a and y
    else:
        return IDS
