# Docker image for deployment in production

#Official miniconda Ubuntu
FROM continuumio/miniconda3:23.5.2-0

RUN addgroup --gid 6666 appgroup && adduser --ingroup appgroup --uid 6666 appuser

WORKDIR /root

RUN mkdir /home/appuser/.ssh

#RUN --mount=type=secret,required,id=config && cat /run/secrets/config > /home/appuser/.ssh/config
#RUN --mount=type=secret,required,id=prkey && cat /run/secrets/prkey > /home/appuser/.ssh/id_rsa_revseq
#RUN --mount=type=secret,required,id=pukey && cat /run/secrets/pukey > /home/appuser/.ssh/id_rsa_revseq.pub
#RUN --mount=type=secret,required,id=khost && cat /run/secrets/khost > /home/appuser/.ssh/known_hosts

RUN chown -R appuser:appgroup /home/appuser/.ssh
RUN chown -R appuser:appgroup /opt/conda

RUN apt-get update
RUN apt-get install -y ssh lftp wget

WORKDIR /data/revseq

USER appuser

RUN export PATH=/data/revseq/:$PATH

RUN conda install -y -n base -c conda-forge mamba libarchive
RUN mamba install -y -c conda-forge -c bioconda -n base snakemake=7.32.3

#ENTRYPOINT [ "sleep", "10d" ]
ENTRYPOINT [ "conda", "run", "-n", "base", "/data/revseq/run_revseq" ]
