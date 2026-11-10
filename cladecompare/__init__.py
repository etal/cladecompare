"""Compare protein sequence alignments to identify contrasting sites."""

import contextlib
import logging

__version__ = "0.3.0.dev0"


class InputError(ValueError):
    """An input file or option that CladeCompare cannot use.

    The command-line programs print the message without a traceback, and the
    web form returns it as the response.
    """


class ToolError(RuntimeError):
    """An external program that CladeCompare cannot run.

    As for InputError, the cladecompare command prints the message without a
    traceback, and cladeweb returns it as the response.
    """


@contextlib.contextmanager
def reading(fname, form=None):
    """Report a failure to read an input file as InputError naming the file.

    An OSError, or bytes that are not UTF-8 text, always become InputError;
    given the form the file should have, such as "a fasta alignment", so
    does a ValueError raised while reading, by the parser of that form or by
    the caller. Of its message, only the first paragraph is kept: Biopython
    follows it with advice on calling Biopython.
    """
    try:
        yield
    except InputError:
        raise
    except OSError as exc:
        raise InputError(f"Cannot read {fname}: {exc.strerror or exc}") from exc
    except UnicodeDecodeError:
        raise InputError(f"Cannot read {fname}: it is not a text file") from None
    except ValueError as exc:
        if form is None:
            raise
        detail = " ".join(str(exc).split("\n\n")[0].split())
        raise InputError(f"Cannot read {fname} as {form}: {detail}") from exc


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
