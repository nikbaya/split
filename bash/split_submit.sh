#! /usr/bin/env bash

set -e

#### USE PREEMPTIBLES
# cluster start ukbb-nb -m n1-standard-16 --worker-machine-type n1-highmem-8 --num-workers 2 --num-preemptible-workers 200

cluster submit ukbb-nb split.py
# cluster submit ukbb-nb split_temp.py

# eof
