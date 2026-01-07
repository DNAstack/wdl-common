from setuptools import setup, find_packages

setup(
    name="detect_duplicate_rgs",
    version="0.1.0",
    description="Detect PacBio BAMs with duplicate RGs from chunk merging",
    author="DNAstack",
    author_email="bioinformatics@dnastack.com",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "detect_duplicate_rgs = detect_duplicate_rgs:main",
        ],
    },
    python_requires=">=3.8",
    install_requires=["pysam"],
)
