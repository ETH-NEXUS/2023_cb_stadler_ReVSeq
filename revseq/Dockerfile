# Docker image for deployment in production

FROM alpine

WORKDIR /app

COPY revseq sync_sftp secrets viollier.conf run_revseq /app/
RUN ["chmod", "+x", "/app/run_revseq"]
#COPY SVVC SVVC/

### HERE RUN CONDA FOR THE BASE
#ENV PATH="/root/miniconda3/bin:${PATH}"
#ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apk update
RUN apk add --no-cache bash
RUN apk add --no-cache openssh
RUN apk add --no-cache lftp
#RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

#RUN wget \
#    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#    && mkdir /root/.conda \
#    && bash Miniconda3-latest-Linux-x86_64.sh -b \
#    && rm -f Miniconda3-latest-Linux-x86_64.sh 
#RUN conda --version

ENTRYPOINT [ "sleep", "10d" ]