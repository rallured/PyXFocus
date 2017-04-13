f2py -c -m zernsurf zernsurf.f95 --f90flags=-fopenmp -lgomp
f2py -c -m transformationsf transformationsf.f95 --f90flags=-fopenmp -lgomp
f2py -c -m surfacesf surfacesf.f95 --f90flags=-fopenmp -lgomp
f2py -c -m woltsurf woltsurf.f95 --f90flags=-fopenmp -lgomp
f2py -c -m reconstruct reconstruct.f95 --f90flags=-fopenmp -lgomp
