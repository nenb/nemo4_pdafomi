#!/usr/bin/env python2

import sys, importlib

for argv in sys.argv[1:]:
	try:
		importlib.import_module(argv)
	except ImportError:
		print("Package %s is missing in Python" % argv)

