#!/usr/bin/env bash

# from: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

cd /home/arc78/applications

wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1686152135&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2ODYxNTIxMzV9fX1dfQ__&Signature=XaUrHQ1rb1C0xc9oAm05S9hwwssBlP6lAfoHBsSE8oQ24fD0nzgq9XFGFwpWJ0gnhtgCOIGRrTVryZe8mkm3i-YlwXn1LdPogEakzasjW8Y3SI-OO3WvgtUNiaVpcnONcNjUozhpKgTI72ZPmAQhxhO2UdE9jruNZJTwChDzPBO886MXhZau6tRsiyfrCcOjoGj3h2vTKp1Zi8I0qlqLYwrgK5fXgMg8i1eG-1eDwlZL1bJ0O6QbPiF0~7RvA5M~Lz~ljf~r9OLmK5XDX3TDBOiebdYvqkHJCA~6ofZb1S5IldTdFRcc9irWFl6NFbgZUUPqXYGazytCHs0LshNjLg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar -xzvf cellranger-7.1.0.tar.gz

# check installation
~/applications/cellranger-7.1.0/cellranger sitecheck > sitecheck.txt

~/applications/cellranger-7.1.0/cellranger testrun --id=tiny
