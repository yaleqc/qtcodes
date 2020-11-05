import numpy as np
import os
import random


decoder_keys = [(d, 1) for d in range(3, 14, 2)]

with open("dsq_script", "w") as f:
    f.write("")

with open("dsq_script", "a") as f:
    for d, T in decoder_keys:
        line = "python3 simulation.py {} {} ".format(
            d, T
        ) + "| tee_simulation_d_{}_T_{}.out\n".format(d, T)
        f.write(line)

