# MazF cleavage specificity

## Usage

### Extracting randomized region from Illumina reads

```
gzip -cd example_input/input_reads.fastq.gz > example_input/input_reads.fastq

extract.py \
    --preset CAA,7,TAAAGATC \
    --output-file example_results/read_counts.csv \
    --output-fasta example_results/trimmed_reads.fa \
    example_input/input_reads.fastq
```

### Determining enriched sequences

```
bin/enrich.py \
    --bootstrap 100 \
    --quantile 0.95 \
    --gen-fasta \
    --output-fasta example_results/enriched_sequences.fasta \
    --output-csv example_results/enriched_sequences.csv \
    example_input/R347.csv example_results/read_counts.csv > example_results/enriched_sequences.log
```

### Generating MEME motif (requires installed [MEME](https://anaconda.org/channels/bioconda/packages/meme/overview))

```
meme example_results/enriched_sequences.fasta \
    -oc example_results/meme \
    -rna \
    -nostatus \
    -time 14400 \
    -mod zoops \
    -nmotifs 3 \
    -minw 3 \
    -maxw 8 \
    -objfun classic \
    -markov_order 0
```
