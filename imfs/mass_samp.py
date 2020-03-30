## for sampling from an IMF
## by LTL

import numpy as np
import scipy.integrate as spint
from popstar.imf import imf, multiplicity


if __name__ == '__main__':
	# stellar population
	# defined in terms of the isochrone that is used
	(lAge,feh) = (9.929,4.916e-05)

	## sampling an imf
	imf_multi = multiplicity.MultiplicityUnresolved()

	#massLimits = np.array([0.08,5, 120])
	#powers = np.array([-1.3, -2.3]) 
	massLimits = np.array([0.1,1.8244])
	powers = np.array([-1.3]) 
	my_imf = imf.IMF_broken_powerlaw(massLimits, powers)

	(masses,is_multi,comps,bleh) = my_imf.generate_cluster(1e6)

	np.savetxt("set1.txt",masses)


