#!/usr/bin/env python3
"""
Detects PacBio BAMs with duplicate RGs from chunk merging.

Exit codes:
  0 - no issue detected
  1 - issue detected and safe to fix (base ID printed to stdout, fixed header written)
  2 - issue detected but NOT safe to fix (other unrelated RGs present)
"""

import argparse
import pysam
import re
import sys
from collections import defaultdict

# Valid PacBio RG ID: 8 lowercase hex chars, optionally followed by /bcFwd--bcRev
# See: https://pacbiofileformats.readthedocs.io/en/latest/BAM.html
PACBIO_RGID_PATTERN = re.compile(r"^[0-9a-f]{8}(/\d+--\d+)?$")


def is_valid_pacbio_rgid(rgid):
    """Check if RG ID matches PacBio spec: 8 hex chars, optionally with barcode pair."""
    return bool(PACBIO_RGID_PATTERN.match(rgid))


def find_base_id(ids):
    """
    Find common base ID if all IDs match pattern: base or base-HEXSUFFIX.

    Returns (base_id, error_msg) where:
    - (base_id, None) = valid PacBio base found
    - (base_id, error_msg) = base found but doesn't match PacBio spec
    - (None, None) = no suffix pattern detected
    """
    if len(ids) < 2:
        return None, None
    base = min(ids, key=len)
    pattern = re.compile(r"^" + re.escape(base) + r"(-[0-9A-Fa-f]+)?$")
    if not all(pattern.match(i) for i in ids):
        return None, None
    # Suffix pattern detected - validate base against PacBio spec
    if not is_valid_pacbio_rgid(base):
        return base, f"Base RG ID '{base}' does not match PacBio spec (expected: 8 lowercase hex chars, optionally with /bcFwd--bcRev)"
    return base, None


def detect(bam_path):
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        rgs = bam.header.get("RG", [])

    if len(rgs) <= 1:
        return 0, None

    # Group RGs by non-ID fields
    groups = defaultdict(list)
    for rg in rgs:
        key = tuple(sorted((k, v) for k, v in rg.items() if k != "ID"))
        groups[key].append(rg["ID"])

    # Find groups with duplicate RGs (same base ID with different suffixes)
    affected = []
    affected_count = 0
    errors = []
    for ids in groups.values():
        base, error = find_base_id(ids)
        if base:
            affected.append(base)
            affected_count += len(ids)
            if error:
                errors.append(error)

    if not affected:
        return 0, None
    if errors:
        return 2, "; ".join(errors)
    if len(affected) > 1:
        return 2, f"Multiple affected RG groups: {affected}"
    if affected_count < len(rgs):
        return 2, f"Other unrelated RGs present (base: {affected[0]})"
    return 1, affected[0]


def write_fixed_header(bam_path, output_path, base_id):
    """Write header with duplicate suffixed RGs removed."""
    suffix_pattern = re.compile(r"@RG\t.*ID:" + re.escape(base_id) + r"-[0-9A-Fa-f]+")
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        header_text = str(bam.header)
    lines = [ln for ln in header_text.splitlines() if not suffix_pattern.match(ln)]
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Detect PacBio BAMs with duplicate RGs from chunk merging"
    )
    parser.add_argument(
        "-i", "--input-bam", type=str, required=True, help="Input BAM/SAM file"
    )
    parser.add_argument(
        "-o",
        "--fixed-header",
        type=str,
        required=True,
        help="Output path for fixed header (only written if issue detected and safe to fix)",
    )
    args = parser.parse_args()

    code, result = detect(args.input_bam)
    if code == 1:
        print(result)
        write_fixed_header(args.input_bam, args.fixed_header, result)
    elif code == 2:
        print(f"ERROR: {result}", file=sys.stderr)
    sys.exit(code)


if __name__ == "__main__":
    main()
