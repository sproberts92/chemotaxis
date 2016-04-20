import numpy as np
import matplotlib.pyplot as plt

def main():
	arr = np.loadtxt("output/total_visits.dat")
	plt.matshow(arr)
	plt.savefig("frames/total_visits.png")
	plt.close()

if __name__ == '__main__':
	main()
