**This script has been written by Rémi Khatib.
Its goal is to treat the outputs of the vvaf program.**

# Behaviour
  1. First it does an apodisation of the signal (filtering)
  2. Then it does a Fourier transform

# Compilation
```gfortran -O2 filter.f95 -o filter```

# Using
## Without input file (Just reply to the questions)
```./filter```

## Launch it with an input file (check the example below)
```./filter < input_filter```

# Remark about the outputs
1. The number of signals in the vvaf file (number of columns - 1) is associated with the number of block in the output file

2. Each block contains 6 columns:
  - Frequency
  - Real part of the FT
  - Im part of the FT
  - Im part of a Laplace transform (from 0 to $+\infty$ and from 0 to $-\infty$). This is really important for odd functions. This is what I use for the simulation of the SFG spectra.
  - Norm of the FT
  - Phase of the FT

