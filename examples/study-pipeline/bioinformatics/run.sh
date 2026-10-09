#!/bin/sh
set -eu
if [ "$#" -ne 1 ]; then
    echo 'Usage: sh run.sh /absolute/path/to/jlinalg.jar' >&2
    exit 2
fi
case "$1" in /*) jar_path=$1 ;; *) jar_path="$(pwd)/$1" ;; esac
cd "$(dirname "$0")"
java -jar "$jar_path" study auto study.yaml --offline
java -jar "$jar_path" study verify --project .
