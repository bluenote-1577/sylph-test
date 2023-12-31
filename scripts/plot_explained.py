import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from collections import defaultdict
import seaborn as sns
from natsort import natsorted
from scipy import stats

t_expl = []
expl = []
for line in open('./results_oct15/holdout_explained.txt','r'):
    expl.append(float(line.strip()))
for line in open('./results_oct15/holdout_explained_true.txt','r'):
    t_expl.append(float(line.strip()) * 100)

plt.scatter(expl,t_expl);
plt.plot([0,100],[0,100])
plt.show()



