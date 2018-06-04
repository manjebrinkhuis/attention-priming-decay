#!/bin/bash
jupyter nbconvert \
    --to python $1 \
    --output out.py
~/.conda/neuro/bin/python out.py
