import numpy as np

ar = []
for line in open('./divergences_low_abund.txt'):
    v = line.split(':')[2]
    ar.append(float(v))

print(np.mean(ar), np.median(ar))

