#! /bin/bash

# The downloaded fa file has an empty first line, which is tripping up cellranger mkref. Remove this first line using sed and write to a new file. 

cd /home/arc78/scratch60/macrophages/bflo/genome/BfsuperV5/files

sed '/^[[:space:]]*$/d' BfsuperV5_fix.fa > BfsuperV5_fix_empty-line-fixed.fa



