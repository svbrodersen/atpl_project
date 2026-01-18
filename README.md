# ATPL Project
To compile and execute the different tests and part of the program you can go into the futhark directory and run the make file.
* For a direct futhark benchark you can run `make benchmark`
* For the different interface tests you can run `make interface`, `make hybrid`, `make unsupoported`, and `make GPUHQP`
* For the HQP vs Qiskit benchmak you can run `make python`

The futhark backend is inside of the futhark directory, while the interface (and one benchmark file) is inside of the interface directory,
finally the benchmark files are on root and in the qiskit directory.