#!/bin/bash -e

# Figure out how to call bwa based on our name.
case $0 in
    *match-bwa-mem.sh) alg=mem;;
    *match-bwa-aln.sh) alg=aln; flags=;;
    *match-bwa-aln-l.sh) alg=aln; flags='-l 1000';;
    *) echo "I don't recognize my own name: $0" >&2; exit 1;;
esac

# Allow a seqgen config file to be passed, else use config.json.
case $# in
    0) config=config.json;;
    1) config=$1;;
    *) echo "Usage $(basename $0) [config-file.json]" >&2; exit 1;;
esac

rm -f query.fasta subject.fasta
seq-gen.py --spec $config | fasta-split-by-id.py --force

if [ ! -f query.fasta -o ! -f subject.fasta ]
then
    echo "$(basename $0): seqgen followed by the fasta split did not produce query.fasta and subject.fasta" >&2
    exit 1
fi

bwa index subject.fasta 2>/dev/null

if [ $alg = mem ]
then
    bwa mem subject.fasta query.fasta 2>/dev/null > bwa.sam
elif [ $alg = aln ]
then
    bwa aln $flags subject.fasta query.fasta 2>/dev/null > bwa.sai
    bwa samse -f bwa.sam subject.fasta bwa.sai query.fasta 2>/dev/null
    rm bwa.sai
fi

egrep -v '^@' < bwa.sam |
    awk '{if ($4 == 0) {print "NO MATCH"} else {print "POS:", $4, "CIGAR:", $6}}'
