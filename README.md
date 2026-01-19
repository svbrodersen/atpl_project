# ATPL Project

To compile and execute the different tests and part of the program there is a Makefile provided.
* For a direct futhark benchark you can run `make futhark`. Further, if you
wish to configure the run, you can change the values within the Makefile in the
futhark directory. Note, the default target for the futhark is CUDA, so if run
on a non-compatible platform, you must change the target in this file.

* For the different interface tests you can run `make interface`, `make
hybrid`, `make unsupoported`, and `make GPUHQP`

* For the HQP vs Qiskit benchmak you can run `make python`

The futhark backend is inside of the futhark directory, while the interface
(and one benchmark file) is inside of the interface directory, finally the
benchmark files are on root and in the qiskit directory.
