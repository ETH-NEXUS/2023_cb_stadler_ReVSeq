FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7d48f7cc7ab6d592a6c8034e4c75c1362d0b57a0c34c9e1bb6f610421e014611"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/bwa.yaml
#   prefix: /conda-envs/ce9ce66279cf2d937f4d9f921103fbc2
#   name: bwa
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - bwa
RUN mkdir -p /conda-envs/ce9ce66279cf2d937f4d9f921103fbc2
COPY workflow/envs/bwa.yaml /conda-envs/ce9ce66279cf2d937f4d9f921103fbc2/environment.yaml

# Conda environment:
#   source: workflow/envs/consensus.yaml
#   prefix: /conda-envs/23d56b9ee5e851cd8b51167844647120
#   name: svvc
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - python==3.6*
#   - pip
#   - pip:
#     - numpy
#     - matplotlib
#     - seaborn
#     - pandas
#     - biopython
RUN mkdir -p /conda-envs/23d56b9ee5e851cd8b51167844647120
COPY workflow/envs/consensus.yaml /conda-envs/23d56b9ee5e851cd8b51167844647120/environment.yaml

# Conda environment:
#   source: workflow/envs/dehuman.yaml
#   prefix: /conda-envs/8cf35ad5a110f7ed1b9a4254ff7aa732
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bwa=0.7.17
#     - samtools>1.11
#     - xxhash
#     - coreutils # [not linux]
RUN mkdir -p /conda-envs/8cf35ad5a110f7ed1b9a4254ff7aa732
COPY workflow/envs/dehuman.yaml /conda-envs/8cf35ad5a110f7ed1b9a4254ff7aa732/environment.yaml

# Conda environment:
#   source: workflow/envs/fetch_primers.yaml
#   prefix: /conda-envs/ad2d76c5aa9984f517f4fd1a4708e08a
#   name: fetch_primers
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - python==3.6*
#   - pip
#   - pip:
#     - pandas
#     - biopython
RUN mkdir -p /conda-envs/ad2d76c5aa9984f517f4fd1a4708e08a
COPY workflow/envs/fetch_primers.yaml /conda-envs/ad2d76c5aa9984f517f4fd1a4708e08a/environment.yaml

# Conda environment:
#   source: workflow/envs/pair_frequencies.yaml
#   prefix: /conda-envs/d12165c4cb6fe427ebb059d043880f4a
#   name: svvc
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - python==3.6*
#   - pip
#   - pip:
#     - numpy
#     - pysam
#     - pandas
RUN mkdir -p /conda-envs/d12165c4cb6fe427ebb059d043880f4a
COPY workflow/envs/pair_frequencies.yaml /conda-envs/d12165c4cb6fe427ebb059d043880f4a/environment.yaml

# Conda environment:
#   source: workflow/envs/samtools.yaml
#   prefix: /conda-envs/a35995b0a5d83c423649af8d9c31d29d
#   name: svvc
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - samtools
RUN mkdir -p /conda-envs/a35995b0a5d83c423649af8d9c31d29d
COPY workflow/envs/samtools.yaml /conda-envs/a35995b0a5d83c423649af8d9c31d29d/environment.yaml

# Conda environment:
#   source: workflow/envs/trim_galore.yaml
#   prefix: /conda-envs/1e28605d97ec109bf58f11e2f91e5485
#   name: trim_galore
#   channels:
#   - conda-forge
#   - bioconda
#   - defaults
#   dependencies:
#   - trim-galore
RUN mkdir -p /conda-envs/1e28605d97ec109bf58f11e2f91e5485
COPY workflow/envs/trim_galore.yaml /conda-envs/1e28605d97ec109bf58f11e2f91e5485/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/ce9ce66279cf2d937f4d9f921103fbc2 --file /conda-envs/ce9ce66279cf2d937f4d9f921103fbc2/environment.yaml && \
    mamba env create --prefix /conda-envs/23d56b9ee5e851cd8b51167844647120 --file /conda-envs/23d56b9ee5e851cd8b51167844647120/environment.yaml && \
    mamba env create --prefix /conda-envs/8cf35ad5a110f7ed1b9a4254ff7aa732 --file /conda-envs/8cf35ad5a110f7ed1b9a4254ff7aa732/environment.yaml && \
    mamba env create --prefix /conda-envs/ad2d76c5aa9984f517f4fd1a4708e08a --file /conda-envs/ad2d76c5aa9984f517f4fd1a4708e08a/environment.yaml && \
    mamba env create --prefix /conda-envs/d12165c4cb6fe427ebb059d043880f4a --file /conda-envs/d12165c4cb6fe427ebb059d043880f4a/environment.yaml && \
    mamba env create --prefix /conda-envs/a35995b0a5d83c423649af8d9c31d29d --file /conda-envs/a35995b0a5d83c423649af8d9c31d29d/environment.yaml && \
    mamba env create --prefix /conda-envs/1e28605d97ec109bf58f11e2f91e5485 --file /conda-envs/1e28605d97ec109bf58f11e2f91e5485/environment.yaml && \
    mamba clean --all -y
