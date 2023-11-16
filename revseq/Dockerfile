# Docker image for deployment in production

#Official miniconda Ubuntu
FROM continuumio/miniconda3:23.5.2-0

RUN addgroup --gid 6666 appgroup && adduser --ingroup appgroup --uid 6666 appuser

WORKDIR /root

RUN mkdir /home/appuser/.ssh

RUN chown -R appuser:appgroup /home/appuser/.ssh
RUN chown -R appuser:appgroup /opt/conda

RUN apt-get update
RUN apt-get install -y ssh lftp wget

WORKDIR /data/revseq

USER appuser

RUN export PATH=/data/revseq/:$PATH

RUN conda install -y -n base -c conda-forge mamba libarchive
RUN mamba install -y -c conda-forge -c bioconda -n base snakemake=7.32.3

ENTRYPOINT [ "/data/revseq/entrypoint.sh" ]
