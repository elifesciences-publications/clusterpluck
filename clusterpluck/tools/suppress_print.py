#!/usr/bin/env Python

import sys
import os
from contextlib import contextmanager


@contextmanager
def suppress_stdout():
	with open(os.devnull, 'w') as devnull:
		old_stdout = sys.stdout
		sys.stdout = devnull
		try:
			yield
		finally:
			sys.stdout = old_stdout


def main():
	suppress_stdout()


if __name__ == '__main__':
	main()
