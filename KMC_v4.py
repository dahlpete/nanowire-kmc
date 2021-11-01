"""
---------------------------------------------------------------------------------------------------------
KINETIC MONTE CARLO SIMULATION OF ELECTRON TRANSFER IN OMCS NANOWIRES

Scheme is adopted from the scheme described in Byun et al. (2014). ChemElectroChem.

This script takes as an input the heme-to-heme electron transfer rates for each of the six heme
pairs in OmcS. The probability of electron transfer from a heme to an adjacent neighbor is 
proportional to the cooresponding rate and is also dependent on the reduced population of the
hemes. Essentially, an electron cannot be transferred to a heme if that heme is already in a 
reduced state. The total electron flux (current) can be determined by monitoring the total number
of electrons that leave the system.

Code written by:
	Peter Dahl (1)(2)
	Malvankar Lab (1)(2)
	Batista Lab (3)

	(1) Yale University, Department of Molecular Biophysics and Biochemistry, New Haven, CT 06511
	(2) Yale University, Microbial Sciences Institute, West Haven, CT 06516
	(3) Yale University, Department of Chemistry, New Haven, CT 06511
  
                                                                                              
    IN   |//////|      |//////|      |//////|      |//////|      |//////|      |//////| OUT      
    ---\ |  H1  | ---> |  H2  | ---> |  H3  | ---> |  H4  | ---> |  H5  | ---> |  H6  | ---\      
    ---/ |RED/OX| <--- |RED/OX| <--- |RED/OX| <--- |RED/OX| <--- |RED/OX| <--- |RED/OX| ---/     
         |//////|      |//////|      |//////|      |//////|      |//////|      |//////|          
                                                                                              
---------------------------------------------------------------------------------------------------------
"""
import sys
import numpy as np
from numpy import random as rand
import marcus_func2 as kET
from decimal import Decimal

T = int(sys.argv[1]) #(K)
bias = float(sys.argv[2]) #(V)

tot_steps = 100000000
number_subunits = 3

der1 = 1.0; #0.12804130749 / 0.128030150754
der2 = 1.0; #0.113391360379 / 0.11032160804
der3 = 1.0; #0.136938372155 / 0.116417085427
der4 = 1.0; #0.0912013822974 / 0.0895979899497
der5 = 1.0; #0.114665406177 / 0.114216080402
der6 = 1.0; #0.0938840547411 / 0.090472361809

habr1 = 1.0; #0.00876963209851 / 0.00877249841624
habr2 = 1.0; #0.0497088764928 / 0.0498598267415
habr3 = 1.0; #0.012046494579 / 0.0116074937589
habr4 = 1.0; #0.044232302219 / 0.0441327613777
habr5 = 1.0; #0.0122095070541 / 0.0122016765867
habr6 = 1.0; #0.0440450731169 / 0.0433170656103

###     Hole transfer rate parameters      ###
#  deltaG (eV)       reorg (eV)      coupling (eV)  
h12dg = -0.098;     h12r = 0.736;    Hab12 = 0.014851183 * habr1
h23dg = -0.418;     h23r = 0.904;    Hab23 = 0.034763377 * habr2
h34dg =  0.175;     h34r = 0.795;    Hab34 = 0.011409582 * habr3
h45dg = -0.040;     h45r = 0.672;    Hab45 = 0.033269693 * habr4
h56dg =  0.197;     h56r = 0.748;    Hab56 = 0.012135786 * habr5
h61dg =  0.184;     h61r = 0.702;    Hab61 = 0.026080843 * habr6

h12de = (h12dg + h12r)**2 / (4*h12r) * der1;     h21de = (-h12dg + h12r)**2 / (4*h12r) * der1
h23de = (h23dg + h23r)**2 / (4*h23r) * der2;     h32de = (-h23dg + h23r)**2 / (4*h23r) * der2
h34de = (h34dg + h34r)**2 / (4*h34r) * der3;     h43de = (-h34dg + h34r)**2 / (4*h34r) * der3
h45de = (h45dg + h45r)**2 / (4*h45r) * der4;     h54de = (-h45dg + h45r)**2 / (4*h45r) * der4
h56de = (h56dg + h56r)**2 / (4*h56r) * der5;     h65de = (-h56dg + h56r)**2 / (4*h56r) * der5
h61de = (h61dg + h61r)**2 / (4*h61r) * der6;     h16de = (-h61dg + h61r)**2 / (4*h61r) * der6


###    HEME-TO-HEME ELECTRON TRANSFER RATES (k_ij)    ###
k12 = kET.marcus_func2(h12de,h12r,Hab12,T); k21 = kET.marcus_func2(h21de,h12r,Hab12,T)
k23 = kET.marcus_func2(h23de,h23r,Hab23,T); k32 = kET.marcus_func2(h32de,h23r,Hab23,T)
k34 = kET.marcus_func2(h34de,h34r,Hab34,T); k43 = kET.marcus_func2(h43de,h34r,Hab34,T)
k45 = kET.marcus_func2(h45de,h45r,Hab45,T); k54 = kET.marcus_func2(h54de,h45r,Hab45,T)
k56 = kET.marcus_func2(h56de,h56r,Hab56,T); k65 = kET.marcus_func2(h65de,h56r,Hab56,T)
k61 = kET.marcus_func2(h61de,h61r,Hab61,T); k16 = kET.marcus_func2(h16de,h61r,Hab61,T)



print "%.3E %.3E %.3E %.3E %.3E %.3E" % (k12,k23,k34,k45,k56,k61)
if bias >= 0:
	k_inv_monomer = np.array([1/k12,1/k23,1/k34,1/k45,1/k56])
	k_inv_poly = np.array([1/k61,1/k12,1/k23,1/k34,1/k45,1/k56])
elif bias < 0:
	k_inv_monomer = np.array([1/k65,1/k54,1/k43,1/k32,1/k21])
        k_inv_poly = np.array([1/k16,1/k65,1/k54,1/k43,1/k32,1/k21]) 
k_inv_tot = np.sum(k_inv_monomer) + number_subunits*np.sum(k_inv_poly)

###    These are used if bias >= 0    ####    These are used if bias < 0    ###     
dV12 = bias * (1/k12) * k_inv_tot**(-1); dV21 = bias * (1/k21) * k_inv_tot**(-1)
dV23 = bias * (1/k23) * k_inv_tot**(-1); dV32 = bias * (1/k32) * k_inv_tot**(-1)
dV34 = bias * (1/k34) * k_inv_tot**(-1); dV43 = bias * (1/k43) * k_inv_tot**(-1)
dV45 = bias * (1/k45) * k_inv_tot**(-1); dV54 = bias * (1/k54) * k_inv_tot**(-1)
dV56 = bias * (1/k56) * k_inv_tot**(-1); dV65 = bias * (1/k65) * k_inv_tot**(-1)
dV61 = bias * (1/k61) * k_inv_tot**(-1); dV16 = bias * (1/k16) * k_inv_tot**(-1)

if bias >= 0:
	k12 = kET.marcus_func2(h12dg - dV12,h12r,Hab12,T); k21 = kET.marcus_func2(-1*h12dg + dV12,h12r,Hab12,T)
	k23 = kET.marcus_func2(h23dg - dV23,h23r,Hab23,T); k32 = kET.marcus_func2(-1*h23dg + dV23,h23r,Hab23,T)
	k34 = kET.marcus_func2(h34dg - dV34,h34r,Hab34,T); k43 = kET.marcus_func2(-1*h34dg + dV34,h34r,Hab34,T)
	k45 = kET.marcus_func2(h45dg - dV45,h45r,Hab45,T); k54 = kET.marcus_func2(-1*h45dg + dV45,h45r,Hab45,T)
	k56 = kET.marcus_func2(h56dg - dV56,h56r,Hab56,T); k65 = kET.marcus_func2(-1*h56dg + dV56,h56r,Hab56,T)
	k61 = kET.marcus_func2(h61dg - dV61,h61r,Hab61,T); k16 = kET.marcus_func2(-1*h61dg + dV61,h61r,Hab61,T)

	k_monomer_f = np.array([k12,k23,k34,k45,k56]); k_poly_f = np.array([k61,k12,k23,k34,k45,k56])
	k_monomer_r = np.array([k21,k32,k43,k54,k65]); k_poly_r = np.array([k16,k21,k32,k43,k54,k65])

elif bias < 0:
	k12 = kET.marcus_func2(h12dg - dV21,h12r,Hab12,T); k21 = kET.marcus_func2(-1*h12dg + dV21,h12r,Hab12,T)
        k23 = kET.marcus_func2(h23dg - dV32,h23r,Hab23,T); k32 = kET.marcus_func2(-1*h23dg + dV32,h23r,Hab23,T)
        k34 = kET.marcus_func2(h34dg - dV43,h34r,Hab34,T); k43 = kET.marcus_func2(-1*h34dg + dV43,h34r,Hab34,T)
        k45 = kET.marcus_func2(h45dg - dV54,h45r,Hab45,T); k54 = kET.marcus_func2(-1*h45dg + dV54,h45r,Hab45,T)
        k56 = kET.marcus_func2(h56dg - dV65,h56r,Hab56,T); k65 = kET.marcus_func2(-1*h56dg + dV65,h56r,Hab56,T)
        k61 = kET.marcus_func2(h61dg - dV16,h61r,Hab61,T); k16 = kET.marcus_func2(-1*h61dg + dV16,h61r,Hab61,T)

	k_monomer_f = np.array([k65,k54,k43,k32,k21]); k_poly_f = np.array([k16,k65,k54,k43,k32,k21])
	k_monomer_r = np.array([k56,k45,k34,k23,k12]); k_poly_r = np.array([k61,k56,k45,k34,k23,k12])

kin = 1e15; kout = 1e15;  # In and Out rates are chosen such that they are greater than any 
                          # of the heme-to-heme rates. Therefore will not be rate limiting.

###    Initialize variables    ###
heme_occupancy = list(np.zeros(number_subunits*6)) # Initialize all hemes in an oxidized state
cumulative_occupation = np.zeros(number_subunits*6)
time = 0
electron_out = 0
timestep = 1

for_rate_array = np.array([])
rev_rate_array = np.array([])
subunit = 1
while subunit <= number_subunits:
	if subunit == 1:
		for_rate_array = np.append(for_rate_array,k_monomer_f)
		rev_rate_array = np.append(rev_rate_array,k_monomer_r)
	else:
		for_rate_array = np.append(for_rate_array,k_poly_f)
		rev_rate_array = np.append(rev_rate_array,k_poly_r)
	subunit += 1

'''
MAIN SIMULATION LOOP
'''
convergenceOUT = open('flux_convergence.txt','w')
print >> convergenceOUT, "step  time    flux (s-1)      current (A)"

while timestep < tot_steps:
	ejection_possibility = False
	ktot = 0; rates = np.array([])
	next_step_occupancy = []

	if heme_occupancy[0] == 0:
		# injection of electron is an option
		rates = np.append(rates, kin)       
		electron_injection = np.zeros(number_subunits*6); electron_injection[0] = 1
		inject_to_system = list(np.array(heme_occupancy)+electron_injection)
		next_step_occupancy.append(inject_to_system)

	if heme_occupancy[-1] == 1:
		# ejection of electron is an option
		ejection_possibility = True
		rates = np.append(rates, kout)
		ejection_index = len(rates) - 1
		electron_ejection = np.zeros(number_subunits*6); electron_ejection[-1] = -1
		eject_from_system = list(np.array(heme_occupancy)+electron_ejection)
		next_step_occupancy.append(eject_from_system)

		if heme_occupancy[-2] == 0:
			# Electron could move back to position -2 (cannot implement in heme loop)
			rates = np.append(rates,rev_rate_array[-1])
			electron_transfer = np.zeros(number_subunits*6)
			electron_transfer[-1] = -1; electron_transfer[-2] = 1
			move_the_electron = list(np.array(heme_occupancy)+electron_transfer)
			next_step_occupancy.append(move_the_electron)

	if heme_occupancy[0] == 1 and heme_occupancy[1] == 0:
		# Electron could move forward from heme 1 subunit 1 (cannot implement in heme loop)
		rates = np.append(rates, for_rate_array[0])
		electron_transfer = np.zeros(number_subunits*6)
		electron_transfer[0] = -1; electron_transfer[1] = 1
		move_the_electron = list(np.array(heme_occupancy)+electron_transfer)
		next_step_occupancy.append(move_the_electron)


	for heme in range(1,number_subunits*6-1):
		if heme_occupancy[heme] == 1:
			# electron transfer to adjacent neighbor is an option
			if heme_occupancy[heme+1] == 0:
				rates = np.append(rates, for_rate_array[heme])
				electron_transfer = np.zeros(number_subunits*6) 
				electron_transfer[heme] = -1; electron_transfer[heme+1] = 1
				move_the_electron = list(np.array(heme_occupancy)+electron_transfer)
				next_step_occupancy.append(move_the_electron)

			if heme_occupancy[heme-1] == 0:
				rates = np.append(rates, rev_rate_array[heme])
				electron_transfer = np.zeros(number_subunits*6)
				electron_transfer[heme] = -1; electron_transfer[heme-1] = 1
				move_the_electron = list(np.array(heme_occupancy)+electron_transfer)
				next_step_occupancy.append(move_the_electron)

	### determine criterium for event selection ###
	ktotal = np.sum(rates)
	random_number1 = rand.random()
	criterium = random_number1 * ktotal

	if criterium < rates[0]:
		event_index = 0
	else:
		i = 1
		while ((np.sum(rates[0:i]) >= criterium) and (criterium >= np.sum(rates[0:i+1]))):
			i += 1
		event_index = i

	heme_occupancy = next_step_occupancy[event_index]
	
	###    Count electrons exiting system    ###
	if ejection_possibility == True and event_index == ejection_index:
		electron_out += 1

	random_number2 = rand.random()
	tau = -np.log(random_number2) * ktotal**(-1)
	time += tau
	
	###    WRITE TEMPORARY RESULT FILE    ###
	if timestep % 1000000 == 0:
		temp_filename = "omcs_kmc_results.temp"
		fileTEMP = open(temp_filename,'w')
		print >> fileTEMP, "TEMPORARY RESULTS FILE"		
		print >> fileTEMP, "This file is written every 1e6 steps"
		print >> fileTEMP, "Current timestep: %s" % (timestep)
		print >> fileTEMP, "Elapsed time: %.3E seconds" % (time)

		temp_occupation = cumulative_occupation * time**(-1)
		print >> fileTEMP, "\nAVERAGE OCCUPANCY"
		for tempvalue in temp_occupation:
			print >> fileTEMP, tempvalue

		tempflux = electron_out * time**(-1)
		if bias >= 0:
			tempI = tempflux * 1.602e-19
		elif bias < 0:
			tempI = -1 * tempflux * 1.602e-19
		
		print >> fileTEMP, "\nTEMPORARY FLUX MEASUREMENT"
		print >> fileTEMP, "electron flux (J) = %.3E s-1" % Decimal(tempflux)
		print >> fileTEMP, "Current (I) = %.3E A" % Decimal(tempI)
		
		print >> convergenceOUT, "%s    %.3E    %.3E    %.3E" % (timestep,time,Decimal(tempflux),Decimal(tempI))		

		fileTEMP.close()

	# move to the next timestep
	timestep += 1
	
	### Determine the time-averated electron occupation density ###
	cumulative_occupation += np.array(heme_occupancy) * tau
average_occupation = cumulative_occupation * time**(-1)

###    Determine the total electron flux and current    ###
electron_flux = electron_out * time**(-1)
if bias >= 0:
	current = electron_flux * 1.602e-19
elif bias < 0:
	current = -1 * electron_flux * 1.602e-19

filename = 'OmcS_KMC_results.out'
fileOUT = open(filename,'w')
print >> fileOUT, "Kinetic Monte Carlo Simulation of ET in OmcS\n"
print >> fileOUT, "results for ET through %s subunits" % (number_subunits)
print >> fileOUT, "total number of timesteps: %s" % (tot_steps)
print >> fileOUT, "total elapsed time: %.3E seconds" % (time)
print >> fileOUT, "\nAVERAGE HEME OCCUPANCY"
for value in average_occupation:
	print >> fileOUT, value

print >> fileOUT, "\nELECTRON FLUX AND CURRENT"
print >> fileOUT, "electron flux (J) = %.3E s-1" % Decimal(electron_flux)
print >> fileOUT, "Current (I) = %.3E A" % Decimal(current)
fileOUT.close()

print "%s electrons left the system in an elapsed time of %.3E seconds" % (electron_out,time)
print "electron flux (J) = %.3E s-1" % Decimal(electron_flux)
print "Current (I) = %.3E A" % Decimal(current)

convergenceOUT.close()
