import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from joblib import Parallel, delayed

def main():
	Parallel(n_jobs=4)(delayed(img_process)(i) for i in trange(1000))

def img_process(i):
	arr = np.loadtxt("output/pcount_{0}.dat".format(i))

	plt.figure(frameon=False)

	ax = plt.axes()
	plt.imshow(arr, interpolation='nearest')

	ax.set_xticks([]); ax.set_yticks([]);
	plt.axis('off')

	plt.savefig("frames/check_rand/pcount_{0}.png".format(i), bbox_inches='tight', pad_inches=0)
	plt.close()

if __name__ == '__main__':
	main()
