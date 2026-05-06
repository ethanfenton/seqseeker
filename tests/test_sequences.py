"""Tests for built-in sequence registry."""

import pytest
from seqseeker._sequences import get_sequence, BUILTIN_SEQUENCES


def test_all_builtins_loadable():
    for name in BUILTIN_SEQUENCES:
        seq = get_sequence(name)
        assert isinstance(seq, str)
        assert len(seq) > 100
        assert all(c in "ACGTN" for c in seq)


def test_egfp_length():
    assert len(get_sequence("egfp")) == 720


def test_unknown_raises():
    with pytest.raises(ValueError, match="Unknown sequence"):
        get_sequence("notasequence")


def test_case_insensitive_not_applicable():
    # get_sequence expects exact case; the C binary handles -n case-insensitively
    with pytest.raises(ValueError):
        get_sequence("EGFP")
