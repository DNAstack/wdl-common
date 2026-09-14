"""Guards the samtools sort behaviour that pbmm2_align_wgs depends on.

pbmm2 writes reads that share an alignment start coordinate in an order that
varies between runs, because alignment is multithreaded and a coordinate sort
keys only on reference, position and strand. That order reaches DeepVariant,
which fills a bounded pileup image by walking the BAM, and any duplicate marker
that breaks a base-quality tie on encounter order. Both change variant calls.

The task therefore re-sorts the aligned BAM by read name and then by
coordinate. That works only because a samtools coordinate sort is stable, so
the name order survives inside a coordinate tie. Stability is real behaviour
rather than a documented guarantee, so these tests fail the image build if a
samtools upgrade ever changes it.

They also pin the choice of -N over -n. The natural ordering used by -n ignores
leading zeros within a digit run, so it reports two distinct read names as
equal and leaves their order unfixed.
"""

import itertools
import subprocess

SEQUENCE = "ACGT" * 25
QUALITY = "I" * len(SEQUENCE)
TIED_NAMES = ["read_c", "read_a", "read_d", "read_b"]
TIED_POSITION = 5000
DISTANT = ("read_elsewhere", 9000)


def write_bam(tmp_path, name, records):
    """Build a BAM whose records are stored in exactly the given order."""
    sam_path = tmp_path / (name + ".sam")
    bam_path = tmp_path / (name + ".bam")
    lines = ["@HD\tVN:1.6\tSO:unknown", "@SQ\tSN:chr1\tLN:200000"]
    for read_name, position in records:
        lines.append(
            "\t".join(
                [
                    read_name,
                    "0",
                    "chr1",
                    str(position),
                    "60",
                    "{}M".format(len(SEQUENCE)),
                    "*",
                    "0",
                    "0",
                    SEQUENCE,
                    QUALITY,
                ]
            )
        )
    sam_path.write_text("\n".join(lines) + "\n")
    subprocess.run(
        ["samtools", "view", "-b", "-o", str(bam_path), str(sam_path)], check=True
    )
    return bam_path


def coordinate_sort(source, destination):
    subprocess.run(
        ["samtools", "sort", "-o", str(destination), str(source)], check=True
    )
    return destination


def canonical_sort(source, destination, tmp_path, name_flag="-N"):
    """The two-stage sort pbmm2_align_wgs runs on the aligned BAM."""
    name_sorted = tmp_path / (destination.name + ".name_sorted.bam")
    subprocess.run(
        ["samtools", "sort", name_flag, "-o", str(name_sorted), str(source)], check=True
    )
    return coordinate_sort(name_sorted, destination)


def read_names(path):
    result = subprocess.run(
        ["samtools", "view", str(path)], check=True, stdout=subprocess.PIPE
    )
    return [
        line.split("\t")[0]
        for line in result.stdout.decode().splitlines()
        if line.strip()
    ]


def tied(order):
    return [(name, TIED_POSITION) for name in order] + [DISTANT]


class TestCoordinateSortIsStable:
    """The property the whole approach rests on."""

    def test_tied_records_keep_their_input_order(self, tmp_path):
        for index, order in enumerate([TIED_NAMES, list(reversed(TIED_NAMES))]):
            source = write_bam(tmp_path, "stable{}".format(index), tied(order))
            result = coordinate_sort(
                source, tmp_path / "stable{}.out.bam".format(index)
            )
            assert read_names(result) == order + [DISTANT[0]], (
                "samtools coordinate sort reordered records tied at one position. "
                "pbmm2_align_wgs relies on that sort being stable."
            )


class TestCanonicalSortIsOrderInvariant:
    def test_every_permutation_converges_on_one_order(self, tmp_path):
        results = set()
        for index, order in enumerate(itertools.permutations(TIED_NAMES)):
            source = write_bam(tmp_path, "perm{}".format(index), tied(list(order)))
            result = canonical_sort(
                source, tmp_path / "perm{}.out.bam".format(index), tmp_path
            )
            results.add(tuple(read_names(result)))
        assert len(results) == 1, "canonical sort produced {} distinct orders".format(
            len(results)
        )

    def test_canonical_order_is_the_name_order(self, tmp_path):
        source = write_bam(tmp_path, "named", tied(TIED_NAMES))
        result = canonical_sort(source, tmp_path / "named.out.bam", tmp_path)
        assert read_names(result) == sorted(TIED_NAMES) + [DISTANT[0]]


class TestAsciiNameSortIsRequired:
    """Why the task uses -N and not -n."""

    COLLIDING = ["read007", "read7"]

    def test_ascii_name_sort_fixes_the_order(self, tmp_path):
        results = set()
        for index, order in enumerate([self.COLLIDING, self.COLLIDING[::-1]]):
            source = write_bam(tmp_path, "ascii{}".format(index), tied(order))
            result = canonical_sort(
                source, tmp_path / "ascii{}.out.bam".format(index), tmp_path, "-N"
            )
            results.add(tuple(read_names(result)))
        assert len(results) == 1, (
            "ASCII name sort left the order of two distinct read names unfixed; "
            "the canonical sort cannot be relied on"
        )

    def test_natural_name_sort_conflates_them(self, tmp_path):
        """Documents the defect that -N avoids, so the choice is not undone."""
        results = set()
        for index, order in enumerate([self.COLLIDING, self.COLLIDING[::-1]]):
            source = write_bam(tmp_path, "natural{}".format(index), tied(order))
            result = canonical_sort(
                source, tmp_path / "natural{}.out.bam".format(index), tmp_path, "-n"
            )
            results.add(tuple(read_names(result)))
        assert len(results) == 2, (
            "samtools -n no longer conflates names differing only by leading "
            "zeros; the reason for preferring -N may no longer hold"
        )
