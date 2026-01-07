"""Detect PacBio BAMs with duplicate RGs from chunk merging."""

from .detect_duplicate_rgs import (
    detect,
    find_base_id,
    is_valid_pacbio_rgid,
    main,
    write_fixed_header,
)

__all__ = [
    "detect",
    "find_base_id",
    "is_valid_pacbio_rgid",
    "main",
    "write_fixed_header",
]
