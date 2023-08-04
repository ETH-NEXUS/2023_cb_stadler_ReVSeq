FROM python:3.11-bookworm

# Optimise python and debian in docker
ENV PYTHONUNBUFFERED 1
ENV PIP_DISABLE_PIP_VERSION_CHECK 1
ENV DEBIAN_FRONTEND noninteractive


# Used for LDAP authentication
# Update repositories and install required packages
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y libsasl2-dev python3-dev libldap2-dev libssl-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Create the required directories
RUN mkdir -p /app /vol/web/media /vol/web/static /vol/backups

# Global python requirements
RUN pip install --upgrade pip
COPY ./requirements.txt /
RUN pip install -r /requirements.txt

# Set working directory
WORKDIR /app