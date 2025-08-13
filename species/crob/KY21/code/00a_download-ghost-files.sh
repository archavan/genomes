#!/usr/bin/env bash

# Here we download genome file for Ciona robusta (KY21 version)
# from http://ghost.zool.kyoto-u.ac.jp/download_ht.html


cd ~/Dropbox/research/projects/endoderm/species/crob/genome/KY21/data/ghost-db/

# genome sequence
curl http://ghost.zool.kyoto-u.ac.jp/datas/HT.Ref.fasta.zip -O

# gene models
curl http://ghost.zool.kyoto-u.ac.jp/datas/HT.KY21Gene.2.fasta.zip -O
curl http://ghost.zool.kyoto-u.ac.jp/datas/HT.KY21Gene.2.gff3.zip -O

# proteome
curl http://ghost.zool.kyoto-u.ac.jp/datas/HT.KY21Gene.protein.2.fasta.zip -O

# These files are zipped; unzip them and gzip them
unzip HT.Ref.fasta.zip
unzip HT.KY21Gene.2.fasta.zip
unzip HT.KY21Gene.2.gff3.zip
unzip HT.KY21Gene.protein.2.fasta.zip

gzip HT.Ref.fasta
gzip HT.KY21Gene.2.fasta
gzip HT.KY21Gene.2.gff3
gzip HT.KY21Gene.protein.2.fasta

rm HT.Ref.fasta.zip
rm HT.KY21Gene.2.fasta.zip
rm HT.KY21Gene.2.gff3.zip
rm HT.KY21Gene.protein.2.fasta.zip


