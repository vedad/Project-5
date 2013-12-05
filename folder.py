#!/usr/bin/env python
"""
Created on Thu 5 Dec 2013

Creates data-structure directories.
"""
import os

if __name__ == '__main__':
    if not os.path.exists('data'):
        os.mkdir('data')

    if not os.path.exists('data/energy'):
        os.mkdir('data/energy')

    if not os.path.exists('data/energy/cluster'):
        os.mkdir('data/energy/cluster')

    if not os.path.exists('data/objects'):
        os.mkdir('data/objects')

    if not os.path.exists('data/parameters'):
        os.mkdir('data/parameters')
