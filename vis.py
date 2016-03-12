import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

for i in tqdm(range(500)):
	arr = np.loadtxt("output/pcount_{0}.dat".format(i))
	plt.matshow(arr)
	plt.savefig("frames/pcount_{0}.png".format(i))
	plt.close()