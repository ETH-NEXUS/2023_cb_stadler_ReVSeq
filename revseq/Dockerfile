# Docker image for deployment in production

FROM alpine

RUN groupadd -g 999 appuser && \
    useradd -r -u 999 -g appuser appuser
USER appuser

WORKDIR /app

COPY revseq sync_sftp secrets viollier.conf run_revseq snakemake.yaml /app/
RUN mkdir /home/appuser/.ssh
COPY ../../secrets/config /home/appuser/.ssh/
COPY ../../secrets/id_rsa_belfry* /home/appuser/.ssh/
RUN ["chmod", "+x", "/app/run_revseq"]

RUN apk update
RUN apk add --no-cache bash openssh lftp wget

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/appuser/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version
ENV PATH="/home/appuser/miniconda3/bin:${PATH}"
ARG PATH="/home/appuser/miniconda3/bin:${PATH}"
RUN conda env create -f snakemake.yaml

ENTRYPOINT [ "sleep", "10d" ]
#ENTRYPOINT [ "run_revseq" ]