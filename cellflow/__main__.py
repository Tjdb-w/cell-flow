"""Module entry point: python -m cellflow."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
