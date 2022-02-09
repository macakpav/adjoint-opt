#!/usr/bin/env python3

from math import *

u2orig=1.4408807664587323e+02
v2orig=5.6636378742932301e+01
alpha2orig=21.45819002769148

alphas = [ alpha2orig + a for a in range(-6,9,2) ]
print(alphas)
v2targets = [ tan(alp/180*pi)*u2orig for alp in alphas ]
print(v2targets)

