# to plot the results of the star detection calculations
# over a grid in distance vs. surface brightness
# by LTL

import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
	arr = np.loadtxt("output.txt")

	plt.imshow(np.log10(arr),extent=(28.,33.,5.,1.))
	plt.colorbar()
	plt.show()