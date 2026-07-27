---
title: CLI Usage
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/terminal
---

```text
usage: pyfgs <seq> [options]

🔗🐍⏭️	PyO3 bindings and Python interface to FragGeneScanRs,
	a gene prediction model for short and error-prone reads.

Input options 💽:

  seq                 Sequence file (or '-' for stdin)
  -m, --model         Sequence error model (default: complete)
                       - short1: Illumina sequencing reads with about 0.1% error rate
                       - short5: Illumina sequencing reads with about 0.5% error rate
                       - short10: Illumina sequencing reads with about 1% error rate
                       - sanger5: Sanger sequencing reads with about 0.5% error rate
                       - sanger10: Sanger sequencing reads with about 1% error rate
                       - pyro5: 454 pyrosequencing reads with about 0.5% error rate
                       - pyro10: 454 pyrosequencing reads with about 1% error rate
                       - pyro30: 454 pyrosequencing reads with about 3% error rate
                       - complete: Complete genomic sequences or short sequence reads without sequencing error
  -r, --reads         Force FASTQ parsing (Overrides auto-detection)
  -w, --whole-genome  Strict contiguous ORFs. Disables error-tolerant frameshift detection.

Output options ⚙️:
  Provide a PATH to save to a file, or use the flag alone to print to stdout.

  --faa [PATH]        Output protein FASTA
  --fna [PATH]        Output nucleotide FASTA
  --bed [PATH]        Output BED6+1 format
  --gff [PATH]        Output GFF3 format
  --vcf [PATH]        Output VCF v4.2 format

Other options 🚧:

  -t, --threads       Number of threads (default: optimal)
  -v, --version       Print version and exit
  -h, --help          Print help and exit

```
