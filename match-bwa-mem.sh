#!/bin/bash -e

# How many times to run.
n=10

# Figure out which matching algorithm to call based on our name.
case $0 in
    *bwa-mem*) alg=bwa; variant=mem;;
    *bwa-aln-l*) alg=bwa; variant=aln; flags='-l 1000';;
    *bwa-aln*) alg=bwa; variant=aln; flags=;;
    *blastn*) alg=blastn;;
    *bowtie2*) alg=bowtie2;;
    *) echo "I don't recognize my own name: $0" >&2; exit 1;;
esac

# Allow a seqgen config file to be passed, else use config.json.
case $# in
    0) config=config.json;;
    1) config=$1;;
    *) echo "Usage $(basename $0) [config-file.json]" >&2; exit 1;;
esac



for i in $(seq $n)
do
    rm -f query.fasta subject.fasta
    seq-gen.py --spec $config | fasta-split-by-id.py --force

    if [ ! -f query.fasta -o ! -f subject.fasta ]
    then
        echo "$(basename $0): seqgen followed by the fasta split did not produce query.fasta and subject.fasta" >&2
        exit 1
    fi

    if [ $alg = bwa ]
    then
        bwa index subject.fasta 2>/dev/null
        if [ $variant = mem ]
        then
            bwa mem subject.fasta query.fasta 2>/dev/null > bwa.sam
        elif [ $variant = aln ]
        then
            bwa aln $flags subject.fasta query.fasta 2>/dev/null > bwa.sai
            bwa samse -f bwa.sam subject.fasta bwa.sai query.fasta 2>/dev/null
            rm bwa.sai
        fi

        egrep -v '^@' < bwa.sam |
            awk '{if ($4 == 0) {print "NO MATCH"} else {print "POS:", $4, "CIGAR:", $6}}'
    elif [ $alg = blastn ]
    then
        makeblastdb -in subject.fasta -dbtype nucl -out subject-blastn >/dev/null
        blastn -db subject-blastn -query query.fasta -outfmt '6 sstart btop' -task blastn > blast.out
        test -s blast.out && (cat blast.out | tr '\t' ' ') || echo 'NO MATCH'
    elif [ $alg = bowtie2 ]
    then
        bowtie2-build -f --quiet subject.fasta subject-bowtie2
        bowtie2 -f -U query.fasta --local -x subject-bowtie2 --no-hd --xeq --quiet -S bowtie.sam
        awk '{if ($4 == 0) {print "NO MATCH"} else {print "POS:", $4, "CIGAR:", $6}}' < bowtie.sam
    fi
done
