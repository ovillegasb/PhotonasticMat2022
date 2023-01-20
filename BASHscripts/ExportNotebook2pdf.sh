#!/bin/bash

# To PDF
jupyter nbconvert --to pdf --template hidecode Example.ipynb

# TO slides
jupyter nbconvert --to slides --no-input --post serve Example.ipynb