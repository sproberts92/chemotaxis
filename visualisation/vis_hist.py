import numpy as np
import matplotlib.pyplot as plt
import sys

def main(argv):
	if len(argv) > 1:
		arr = np.loadtxt(argv[1])
		plt.plot(arr)
		plt.show()
		# plt.savefig("frames/total_visits.png")
		plt.close()
	else:
		print("Specify input file.")

if __name__ == '__main__':
	main(sys.argv)
