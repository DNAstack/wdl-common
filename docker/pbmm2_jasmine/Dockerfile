FROM quay.io/pacbio/pbmm2@sha256:b58eef0645dd9adca850c9b4811f289ddcd57f630596d9f524701bc787467f30

LABEL MAINTAINER="Bioinformatics <bioinformatics@dnastack.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME="${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG="${IMAGE_TAG}"

RUN apt-get -qq update && \
	apt-get -qq install \
        wget

ARG JASMINE_VERSION
RUN wget https://github.com/PacificBiosciences/jasmine/releases/download/v${JASMINE_VERSION}/jasmine.tar.xz && \
    tar xf jasmine.tar.xz -C /usr/local && \
    chown root:root /usr/local/bin/jasmine && \
    chown root:root /usr/local/lib/libonnxruntime.so.* && \
    chmod +x /usr/local/bin/jasmine

# Install detect_duplicate_rgs package and run tests
COPY detect_duplicate_rgs/ /opt/detect_duplicate_rgs/
RUN python3 -m pip install --quiet pytest /opt/detect_duplicate_rgs/ && \
    pytest /opt/detect_duplicate_rgs/tests/ -v && \
    rm -rf \
        /opt/detect_duplicate_rgs/tests/__pycache__/ \
        /opt/detect_duplicate_rgs/.pytest_cache/ \
        /opt/detect_duplicate_rgs/*.egg-info/

# pbmm2_align_wgs fixes the order of coordinate-tied reads by sorting on name
# and then on coordinate, which is only correct while this image's samtools
# coordinate sort remains stable. Fail the build rather than let a samtools
# upgrade silently reintroduce non-deterministic variant calls.
COPY tests/ /opt/pbmm2_jasmine_tests/
RUN pytest /opt/pbmm2_jasmine_tests/ -v && \
    rm -rf \
        /opt/pbmm2_jasmine_tests/__pycache__/ \
        /opt/pbmm2_jasmine_tests/.pytest_cache/
