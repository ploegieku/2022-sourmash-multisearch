#!/bin/bash


for SCALED in 1 2 5 ; do
  for KSIZE in $(seq 5 1 30); do
    for ALPHABET in hp dayhoff; do
      SUFFIX=$ALPHABET.k$KSIZE.scaled$SCALED

      for FASTA in *fa ; do
        SIG=$FASTA.$SUFFIX.sig
        CSV=$FASTA.$SUFFIX.describe.csv
        KMERS=$FASTA.$SUFFIX.kmers.csv
        sourmash sketch protein \
          --output $SIG \
          --singleton --$ALPHABET \
          -p k=$KSIZE,scaled=$SCALED,abund \
          $FASTA
        sourmash sig describe \
          --csv $CSV \
          --quiet \
          $SIG
#        python sig2kmer.py \
#          --output-kmers $KMERS \
#          --input-is-protein \
#          --ksize $KSIZE \
#          --$ALPHABET \
#          --no-dna \
#          $SIG $FASTA
      done

      C_ELEGANS=Caenorhabditis_elegans.WBcel235.pep*$SUFFIX.sig
      HUMAN=gencode.v38.basic.annotation.protein*$SUFFIX.sig

      echo "C_ELEGANS" $C_ELEGANS
      echo "HUMAN" $HUMAN

      # --num-results 0  to report all results
      OUTDIR=multisearch-$ALPHABET-k$KSIZE-scaled$SCALED
      mkdir $OUTDIR
      python sourmash-multisearch.py \
        --threshold 0.0001 \
        --$ALPHABET --no-dna -k $KSIZE \
        --num-results 10 \
        --query $C_ELEGANS \
        --databases $HUMAN \
        --add-kmer-stats \
        --output-dir $OUTDIR \
        > $OUTDIR/stdout.log
    done
  done
done
