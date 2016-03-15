import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from joblib import Parallel, delayed

def main():
	Parallel(n_jobs=4)(delayed(img_process)(i) for i in trange(5000))
		
def img_process(i):
	arr = np.loadtxt("output/pcount_{0}.dat".format(i))
	plt.matshow(arr)
	plt.savefig("frames/pcount_{0}.png".format(i))
	plt.close()

if __name__ == '__main__':
	main()