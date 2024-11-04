#!/bin/bash 

# clean up caches files
# rm -rf dist build

find . \
  -type d \
  \( \
    -name "*cache*" \
    -o -name "*.dist-info" \
    -o -name "*.egg-info" \
  \) \
  -not -path "./.venv/*" \
  -exec rm -r {} +

# python3 -m build --sdist --wheel ./
# cd dist 
# unzip *.whl
# cd ..
