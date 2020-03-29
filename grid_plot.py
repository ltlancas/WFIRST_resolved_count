# to plot the results of the star detection calculations
# over a grid in distance vs. surface brightness
# by LTL

import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

	J_outs = np.loadtxt("outputs/J_out.txt")
	J_detect_5sig = 25.939

	arr = np.loadtxt("outputs/feh1_age9005.txt")
	d = np.linspace(1.,25.,51)
	mucorr = J_detect_5sig - J_outs
	dcorr = 1e-5*(10**(mucorr/5))

	darr = np.array([d for i in range(arr.shape[1])]).T

	ax = plt.subplot(211)
	plt.imshow(np.log10(arr),extent=(28.,33.,25.,1.))
	plt.colorbar()
	ax.set_aspect(0.2)
	ax2 = plt.subplot(212)
	plt.hist(dcorr,bins="auto",histtype="step")
	plt.xlim(1.,25.)
	plt.yscale("log")
	plt.show()