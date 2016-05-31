#!/bin/bash


## Installation comamnds
python setup.py install --user



## Clean up
rm -rf PyDAIR.egg-info
rm -rf build
rm -rf dist
rm ./PyDAIR/test/data/test_output_*
rm PyDAIR/*.pyc
rm PyDAIR/app/*.pyc
rm PyDAIR/io/*.pyc
rm PyDAIR/seq/*.pyc
rm PyDAIR/test/*.pyc
rm PyDAIR/bin/*.pyc
rm PyDAIR/plot/*.pyc
rm PyDAIR/stats/*.pyc
rm PyDAIR/utils/*.pyc



## PyPI registration commands
# python setup.py sdist
# python setup.py bdist_wheel --universal
# twine upload dist/*


