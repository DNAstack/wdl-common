#!/usr/bin/env python3
"""Tests for detect_duplicate_rgs."""

import tempfile
from pathlib import Path

import pytest

from detect_duplicate_rgs import detect, is_valid_pacbio_rgid, write_fixed_header


def build_sam(rg_list):
    """Build SAM content from a list of RG dicts."""
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:248956422",
    ]
    for rg in rg_list:
        fields = "\t".join(f"{k}:{v}" for k, v in rg.items())
        lines.append(f"@RG\t{fields}")
    return "\n".join(lines) + "\n"


def run_detect(rg_list):
    """Run detect() on a temporary SAM file built from RG dicts."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sam", delete=False) as f:
        f.write(build_sam(rg_list))
        sam_path = f.name
    try:
        return detect(sam_path)
    finally:
        Path(sam_path).unlink()


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

RGS_NO_ISSUE = [
    {"ID": "94c47df0/36--36", "PL": "PACBIO", "LB": "lib1", "SM": "sample1"},
]

RGS_SINGLE_ISSUE = [
    {"ID": "94c47df0/36--36", "PL": "PACBIO", "LB": "lib1", "SM": "sample1"},
    {"ID": "94c47df0/36--36-35ED621", "PL": "PACBIO", "LB": "lib1", "SM": "sample1"},
    {"ID": "94c47df0/36--36-21D71C81", "PL": "PACBIO", "LB": "lib1", "SM": "sample1"},
    {"ID": "94c47df0/36--36-C911267", "PL": "PACBIO", "LB": "lib1", "SM": "sample1"},
]

RGS_ISSUE_PLUS_UNRELATED = [
    {"ID": "94c47df0/36--36", "PL": "PACBIO", "LB": "lib1", "SM": "sample1"},
    {"ID": "94c47df0/36--36-35ED621", "PL": "PACBIO", "LB": "lib1", "SM": "sample1"},
    {"ID": "a1b2c3d4/10--10", "PL": "PACBIO", "LB": "lib2", "SM": "sample2"},
]

RGS_NON_PACBIO = [
    {"ID": "sample1", "PL": "ILLUMINA", "LB": "lib1", "SM": "sample1"},
    {"ID": "sample1-ABC123", "PL": "ILLUMINA", "LB": "lib1", "SM": "sample1"},
]

RGS_NO_BARCODES = [
    {"ID": "94c47df0", "PL": "PACBIO", "LB": "lib1", "SM": "sample1"},
    {"ID": "94c47df0-ABC1234", "PL": "PACBIO", "LB": "lib1", "SM": "sample1"},
]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestDetect:
    def test_no_issue_single_rg(self):
        code, result = run_detect(RGS_NO_ISSUE)
        assert code == 0
        assert result is None

    def test_single_rg_issue_fixable(self):
        code, result = run_detect(RGS_SINGLE_ISSUE)
        assert code == 1
        assert result == "94c47df0/36--36"

    def test_issue_plus_unrelated_rgs(self):
        code, result = run_detect(RGS_ISSUE_PLUS_UNRELATED)
        assert code == 2
        assert "Other unrelated RGs present" in result

    def test_non_pacbio_rg_format(self):
        code, result = run_detect(RGS_NON_PACBIO)
        assert code == 2
        assert "does not match PacBio spec" in result

    def test_pacbio_no_barcodes(self):
        code, result = run_detect(RGS_NO_BARCODES)
        assert code == 1
        assert result == "94c47df0"


class TestPacBioRgIdValidation:
    def test_valid_8hex(self):
        assert is_valid_pacbio_rgid("abcd1234")

    def test_valid_8hex_with_barcodes(self):
        assert is_valid_pacbio_rgid("abcd1234/0--0")
        assert is_valid_pacbio_rgid("abcd1234/36--36")

    def test_invalid_uppercase_hex(self):
        assert not is_valid_pacbio_rgid("ABCD1234")

    def test_invalid_short_hex(self):
        assert not is_valid_pacbio_rgid("abc123")

    def test_invalid_non_hex(self):
        assert not is_valid_pacbio_rgid("sample1")


class TestWriteFixedHeader:
    def test_removes_suffixed_rgs(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sam", delete=False) as f:
            f.write(build_sam(RGS_SINGLE_ISSUE))
            sam_path = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sam", delete=False) as f:
            header_path = f.name

        try:
            write_fixed_header(sam_path, header_path, "94c47df0/36--36")
            header = Path(header_path).read_text()
            assert "ID:94c47df0/36--36\t" in header
            assert "ID:94c47df0/36--36-" not in header
        finally:
            Path(sam_path).unlink()
            Path(header_path).unlink()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
