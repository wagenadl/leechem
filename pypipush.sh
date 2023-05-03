#!/bin/sh

python3 setup.py sdist
twine upload dist/leechem-1.2.tar.gz
