import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys

fnames = ["modeijk", "modeikj", "modejik", "modejki", "modekij", "modekji", "intrinsic"]
for fname in fnames:
	fl = open("file_to_fit.txt", "w+")
	fl.write(fname)
	fl.close()
	subprocess.call(["gnuplot",'timefit.plt'])
