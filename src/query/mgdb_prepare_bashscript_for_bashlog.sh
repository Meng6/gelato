#!/bin/bash

curl --data-binary @"${snakemake_input[datalog]}" https://www.thomasrebele.org/projects/bashlog/api/datalog\?query= > "${snakemake_output[0]}"
