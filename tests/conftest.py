"""Fixtures shared by the test modules."""

import os
import shutil

import pytest


@pytest.fixture
def require():
    """Return a function that skips the test unless programs are on PATH.

    The function takes program names. If CLADECOMPARE_REQUIRE_TOOLS is
    non-empty, as ``pixi run test-integration`` sets it, a missing program
    fails the test instead of skipping it.
    """

    def require(*programs):
        missing = [name for name in programs if shutil.which(name) is None]
        if missing:
            reason = "not on PATH: " + ", ".join(missing)
            if os.environ.get("CLADECOMPARE_REQUIRE_TOOLS"):
                pytest.fail(reason)
            pytest.skip(reason)

    return require
