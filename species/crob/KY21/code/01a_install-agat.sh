#!/usr/bin/env bash

module load miniconda

conda create -n agat
conda activate agat

mamba install -c bioconda agat

conda deactivate
