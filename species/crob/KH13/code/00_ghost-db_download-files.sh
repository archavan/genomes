#!/usr/bin/env bash

# Here we download genome file for Ciona robusta (KH 2013 version)
# from http://ghost.zool.kyoto-u.ac.jp/download_kh.html


cd ~/research/endoderm/species/crob/genome/KH13/data/ghost-db

# genome sequence
curl http://ghost.zool.kyoto-u.ac.jp/datas/JoinedScaffold.zip -O

# gene models
curl http://ghost.zool.kyoto-u.ac.jp/datas/KH.KHGene.2013.fasta.zip -O
curl http://ghost.zool.kyoto-u.ac.jp/datas/KH.KHGene.2013.gff3.zip -O

# proteome
curl http://ghost.zool.kyoto-u.ac.jp/datas/KH.KHGene.2013.Longest.protein.zip -O
curl http://ghost.zool.kyoto-u.ac.jp/datas/KH.KHGene.2013.Met.protein.zip -O

