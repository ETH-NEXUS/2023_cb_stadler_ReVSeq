#!/usr/bin/env bash

######
echo "Setting up the credentials"
cp /run/secrets/config /home/appuser/.ssh/config
cp /run/secrets/prkey /home/appuser/.ssh/id_rsa_revseq
cp /run/secrets/pukey /home/appuser/.ssh/id_rsa_revseq.pub
cp /run/secrets/khost /home/appuser/.ssh/known_hosts

sleep 10d
#conda run -n base /data/revseq/run_revseq