# stark
Stark broadened line profiles for Hydrogen, using the tables of Lemke 1997

## installation
```
git clone https://github.com/StuartLittlefair/stark
cd stark
pip install .
```

## usage

```python
from stark import make_hydrogen_line_profile
import numpy as np
from matplotlib import pyplot as plt

# lower level of line of interest
nlo = 2
# upper levelt
nhi = 3

# temp in K
temp = 2500

# electron density in cm^-3
nelec = 1e10

# macro turbulence and instrumental line broadening in km/s
v_macro = 5
v_inst = 5

# make the line profile
velocity, line_profile = make_hydrogen_line_profile(nlo, nup, temp, nelec, v_inst, v_macro)

# save it
np.savetxt('lprof.txt', np.column_stack((velocity, line_profile)))

# plot it
plt.plot(velocity, line_profile)
plt.show()
```


