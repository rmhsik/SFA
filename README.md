# pySFA

Python library for computing High-Order Harmonic Generation (HHG) spectra relying on the Strong Field Approximation (SFA) with saddle point.


## Building

Create a new `build` folder, execute `cmake` and `make`
```bash
git clone git@github.com:rmhsik/pySFA.git
cd pySFA
mkdir build && cd build
cmake ..
make
```
## API Reference

#### pySFA(Ip, efield, t, nthreads=1)

Compute the dipole response in time domain using the Strong-Field Approximation. 

| Parameter | Type         | Description                       |
| :-------- | :------------| :---------------------------------|
| `Ip`      | `double` | Ionization potential of the atomic target in atomic units |
|`efield`   | `complex` | Complex driving electric field array. It must be the same legnth as `t`. It must be in atomic units |
|`t`        | `double` | Temporal array in atomic units    |
|`nt`       | `int`     | Number of threads to perform the integration. *Optional*|

| Return | Type         | Description   |
| :-------- | :------------| :--------------------------------|
| `dipole_array`      | `complex`| Complex dipole response array in the temporal domain |

## Usage/Examples
Check examples directory.
```python
import numpy as np
import matplotlib.pyplot as plt
from pySFA import pySFA

def sin2(t,tmax):
  '''
  Sin-squared temporal envelope
  '''
  env = 0
  if -tmax/2<t and t < tmax/2:
      env = np.sin(np.pi*(t-tmax/2)/tmax)**2
  return env

w0 = 0.057 # 800 nm wavelength in atomic untis
e0 = 0.067 # ~1.6x10^14 W/cm2 
period = 2*np.pi/w0

t,dt = np.linspace(-10*period,10*period, 8192,retstep=True)
n = t.shape[0]
tmax = 16*period
mask = np.array([sin2(ti,tmax) for ti in t])

efield = e0*mask*np.exp(1j*w0*t)

Ip = 0.578 # Ionization potential for Argon
dipole  = pySFA(Ip, efield, t, nthreads=16)

## Spectrum calculation
dw = 2*np.pi/(n*dt)
wmax = n*dw/2
w = np.arange(0,wmax,dw)
hhg_spectrum = np.abs(np.fft.rfft(dipole.real))**2
```


## License

[MIT](https://choosealicense.com/licenses/mit/)


