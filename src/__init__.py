#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import os
import sys

# Makes relative imports work consistently across Python versions.
here = os.path.dirname(os.path.realpath(__file__))
sys.path.append(here)

# Single source of truth for version.
try:
    version = open(os.path.join(here, '..', 'VERSION'), 'r').readlines()[0].strip()
except IOError:
    version = 'unknown'
