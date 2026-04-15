#!/bin/sh

query='/mnt/4TB/giovanna/foldseek/query/TCR3d/MHCI_groove_aligned/db_samuel'
target='/mnt/4TB/giovanna/foldseek/target/pdb100/pdb'
alignment_file='/mnt/4TB/giovanna/foldseek/version_02/pdb/dbf_pdb_aln'
tmp_file='/mnt/4TB/giovanna/foldseek/version_02/pdb/tmpFolder'
fmt_output='query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,lddtfull,prob'

foldseek easy-search $query $target $alignment_file $tmp_file \
    --cluster-search 1 \
    --exhaustive-search 1 \
    --format-output $fmt_output
