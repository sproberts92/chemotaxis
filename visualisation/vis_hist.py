import numpy as np
import matplotlib.pyplot as plt

def main():
	arr = np.loadtxt("output/visited_hist.dat")
	plt.plot(arr)
	plt.show()
	# plt.savefig("frames/total_visits.png")
	plt.close()

if __name__ == '__main__':
	main()
