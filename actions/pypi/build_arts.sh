#!/bin/sh

# Get arts source
mkdir arts && cd arts && git init .
git remote add origin https://github.com/${GITHUB_REPOSITORY}
git fetch --prune --depth=1 origin

# Build pyarts
git checkout --force ${GITHUB_REF#refs/heads/}; mkdir build; cd build;
cmake3 -DENABLE_FORTRAN=1 -DBLAS_blas_LIBRARY=/usr/lib64/atlas/libtatlas.so -DLAPACK_lapack_LIBRARY=/usr/lib64/atlas/libtatlas.so ..
make pyarts

# Packaging
cd python
python3 setup.py sdist bdist_wheel
auditwheel repair dist/pyarts*.whl
python3 -m twine upload wheelhouse/pyarts*.whl -u __token__ -p ${INPUT_pypi_access}
