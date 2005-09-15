Tcov directory:  Torture test for CHOLMOD, with statement coverage.
--------------------------------------------------------------------------------

This test suite is not required to compile and use CHOLMOD.  It is thus
not ported to all architectures.  Linux is assumed; see the Makefile for
running on Solaris.  Use tcov instead of gcov in the "covall" script.  Edit
the Makefile and change the definition of CC.  You may need to change PRETTY
as well.  You will need to edit LIB to reflect the proper LAPACK BLAS
libraries.

Requires all CHOLMOD modules except the Partition Module, which it can
optionally use (and test).  Also acts as a statement coverage test for
AMD, COLAMD, and CCOLAMD.

Type "make" in this directory to compile and run CHOLMOMD with statement
coverage testing.  Every line of CHOLMOD will be exercised, and its results
checked.  The line "All tests passed" should be printed for each test on
stderr.  Some matrices will report NaN as their maximum error; these are the
three singular test matrices.  This test result is expected.

The source code files are first preprocessed with cc -E, and the resulting
file (z_*.c or zz_*.c) is then compiled.  This is to ensure that all lines
within macros and included *.c files are tested (take a look at z_updown.c
and z_solve.c if you'd like to see how loop-unrolling and real/complex
templates are done in CHOLMOD, and compare those files with their source files
in ../Modify and ../Cholesky).

Note that many, many error messages will appear in the test output itself
(tmp/*.output), because all of CHOLMOD's error handling is checked as well.
These errors are expected.  Any unexpected error will cause the test to fail.
The last line of each output file should be "All tests successful".

To remove all but the original source files and output files from
this directory, type "make clean".  To remove all but the
files in the original distribution, type "make distclean".
