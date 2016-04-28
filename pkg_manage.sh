#!/bin/bash


## Installation comamnds
python setup.py install --user



## Clean up
rm -rf PyDAIR.egg-info
rm -rf build
rm -rf dist
rm ./PyDAIR/test/data/test_output_*




## PyPI registration commands
# python setup.py egg_info
# python setup.py sdist
# twine register dist/*

