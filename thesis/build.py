#!/usr/bin/env python3
"""Build thesis and clean auxiliary files afterward."""

import os
import subprocess
import sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))

subprocess.run([sys.executable, "-c", ""], check=True, capture_output=True)  # sanity check
result = subprocess.run(["latexmk", "-pdf", "main.tex"])
if result.returncode != 0:
    sys.exit(result.returncode)

subprocess.run(["latexmk", "-c"])
