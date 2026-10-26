"""Compare protein sequence alignments to identify contrasting sites."""

import logging

__version__ = "0.3.0.dev0"


def configure_logging(quiet=False):
    """Log to standard error for the command-line programs.

    By default, status messages are shown with the module and line that
    logged them; with quiet, only warnings and errors, with the module.
    """
    if quiet:
        logging.basicConfig(level=logging.WARNING, format="%(module)s: %(message)s")
    else:
        logging.basicConfig(
            level=logging.INFO, format="%(module)s [@%(lineno)s]: %(message)s"
        )
