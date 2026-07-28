---
title: Performance & Benchmarks
author: Tom Stanton
comments: true
tags: [markdown, documentation, web]
icon: lucide/timer
---

`pyfgs` is continuously benchmarked against NCBI RefSeq ground-truth datasets on every commit to `main` to ensure we
never introduce performance regressions.

When configured in `whole_genome=True` mode, `pyfgs` uses zero-copy memory operations and Rayon-based multi-threading to out-perform even SIMD-accelerated tools like `pyrodigal`—especially on fragmented datasets (e.g. multi-contig draft genomes).

## Continuous Accuracy & Speed Tracking

<div style="border: 1px solid #e0e0e0; border-radius: 8px; overflow: hidden; margin-top: 1em;">
    <iframe src="benchmarks/index.html" width="100%" height="800px" style="border:none;"></iframe>
</div>