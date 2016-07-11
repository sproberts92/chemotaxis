import numpy as np
import matplotlib.pyplot as plt

def main():
	arr = np.loadtxt("output/total_visits.dat")

	plt.figure(frameon=False)
	ax = plt.axes()

	plt.imshow(arr)

	ax.set_xticks([]); ax.set_yticks([]);
	plt.axis('off')

	plt.savefig("frames/total_visits.png", bbox_inches='tight', pad_inches=0)
	# plt.show()
	plt.close()

if __name__ == '__main__':
	main()
