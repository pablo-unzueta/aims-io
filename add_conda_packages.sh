#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 package1 package2 package3 ..."
    exit 1
fi

for pkg in "$@"; do
    version=$(conda list | grep "^$pkg " | awk '{print $2}')
    if [ -n "$version" ]; then
        echo "Adding $pkg version $version"
        poetry add "${pkg}@^${version}"
    else
        echo "Package $pkg not found in the current Conda environment"
    fi
done