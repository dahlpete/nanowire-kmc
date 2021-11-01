import numpy as np
def marcus_func2(dE,reorg,Hab,T):

	kb=8.617e-5
	hbar=6.582e-16

	MarcusET = (2*np.pi / hbar) * Hab**2 * (1/np.sqrt(4*np.pi*reorg*kb*T)) * np.exp(-dE/(kb*T))
	
	return MarcusET
