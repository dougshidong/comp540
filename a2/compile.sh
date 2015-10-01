gfortran -c prec.f
gfortran -Wall -Wextra -Wconversion -pedantic -fbounds-check -ffpe-trap=zero,overflow,underflow -g -o chol.exe chol.f
