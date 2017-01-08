#The goal fo this script is to treat the outputs
#of the vvaf program.
#1) First it does an apodisation of the signal (filtering)
#2) Then it does a Fourier transform

#Compilation
> gfortran -O2 filter.f95 -o filter

#Running without input file (Just reply to the questions)
> ./filter

#Running with input file (check the example)
> ./filter < input_filter

#========================
#Remark about the outputs
#========================
A) The number of signals the vvaf file (number of columns - 1) is
associated with the number of block in the output file

B) Each block contains 6 columns:
  1) Frequency
  2) Real part of the FT
  3) Im part of the FT
  4) Im part of a Laplace transform (from 0 to +inf and from 0 to -inf)
     This is really important for odd functions. This is what I use
     for the simulation of the SFG spectra
  5) norm of the FT
  6) Phase of the FT

