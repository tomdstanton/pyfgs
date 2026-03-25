window.BENCHMARK_DATA = {
  "lastUpdate": 1774435330696,
  "repoUrl": "https://github.com/tomdstanton/pyfgs",
  "entries": {
    "Benchmark": [
      {
        "commit": {
          "author": {
            "email": "tomdstanton@gmail.com",
            "name": "Tom Stanton",
            "username": "tomdstanton"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "37e54db0d235850271957ad2d43e281c757b01d3",
          "message": "Change branch for benchmark data from gh-pages",
          "timestamp": "2026-03-25T21:32:12+11:00",
          "tree_id": "c58f06e981cf77ab2f3ca159b2e4a4270f2411de",
          "url": "https://github.com/tomdstanton/pyfgs/commit/37e54db0d235850271957ad2d43e281c757b01d3"
        },
        "date": 1774435330250,
        "tool": "pytest",
        "benches": [
          {
            "name": "benchmarks/bench_accuracy.py::test_pyrodigal_accuracy[E_coli]",
            "value": 0.5055964071435928,
            "unit": "iter/sec",
            "range": "stddev: 0.006129775083433936",
            "extra": "mean: 1.9778621562000012 sec\nrounds: 5"
          },
          {
            "name": "benchmarks/bench_accuracy.py::test_pyrodigal_accuracy[S_aureus]",
            "value": 1.3243771200455143,
            "unit": "iter/sec",
            "range": "stddev: 0.0030285541614447214",
            "extra": "mean: 755.0719389999983 msec\nrounds: 5"
          },
          {
            "name": "benchmarks/bench_accuracy.py::test_pyrodigal_accuracy[P_aeruginosa]",
            "value": 0.36534813186911647,
            "unit": "iter/sec",
            "range": "stddev: 0.014294190298677312",
            "extra": "mean: 2.737115405200001 sec\nrounds: 5"
          },
          {
            "name": "benchmarks/bench_accuracy.py::test_pyfgs_accuracy[E_coli]",
            "value": 0.5655532329249502,
            "unit": "iter/sec",
            "range": "stddev: 0.003937772728353946",
            "extra": "mean: 1.7681801495999963 sec\nrounds: 5"
          },
          {
            "name": "benchmarks/bench_accuracy.py::test_pyfgs_accuracy[S_aureus]",
            "value": 0.906743582412125,
            "unit": "iter/sec",
            "range": "stddev: 0.0010379956301260539",
            "extra": "mean: 1.1028476180000013 sec\nrounds: 5"
          },
          {
            "name": "benchmarks/bench_accuracy.py::test_pyfgs_accuracy[P_aeruginosa]",
            "value": 0.4951477308866112,
            "unit": "iter/sec",
            "range": "stddev: 0.06614835655890672",
            "extra": "mean: 2.0195992784000056 sec\nrounds: 5"
          }
        ]
      }
    ]
  }
}