#!/usr/bin/env bash

# from: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

cd /home/arc78/applications

wget -O cellranger-7.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.0.tar.gz?Expires=1656070463&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NTYwNzA0NjN9fX1dfQ__&Signature=FfYBAjw~150aT7jTVfm7Y-Pawx0mdSU9pTqdr7npdSo9-nLsuY~PZyVt6GBEY371D~~DW~Y09wOmS7xTadQWN5Yx5j0F1NvZIzIUZUXxVkMhgLsZ9vNLyH9D7UDQSWD4mNbGBycEQF3BeNDhHo~vlfsJ49vEhZ9YDLqTvuQff8jtI9TGWH2QBrZVIiC~O59j8Jrcxdonp8Lzccov7GfSLUXmXLR2e1oORo-10xkLLT0BoOYx0c9lIalsW6ryaqFrao5s~uZPh9swG5RVs7RsXWDOexkBZbzRdJXCo5iG5KSIUt9T8ruXJMQsZW5X8pbL1HU7JJTp6VX-argPOLEuhw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar -xzvf cellranger-7.0.0.tar.gz

# check installation
~/arc78/applications/cellranger-7.0.0/cellranger sitecheck > sitecheck.txt

~/arc78/applications/cellranger-7.0.0/cellranger testrun --id=tiny
