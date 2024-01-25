# Crystal Metric

This is set of 21 benchmarks for Crystal language, in one file. Program were taken from my benchmarks: [jit-benchmarks](https://github.com/kostya/jit-benchmarks), [crystal-benchmarks-game](https://github.com/kostya/crystal-benchmarks-game), [benchmarks](https://github.com/kostya/benchmarks). Benchmarks implemented as separated classes, this is good for test multi module compilation.

## Run

```
crystal build run.cr --release
./run
```

Output

```
Binarytrees: ok in 0.815s
Brainfuck: ok in 4.504s
Brainfuck2: ok in 1.547s
Fannkuchredux: ok in 2.313s
Fasta: ok in 0.961s
Knuckeotide: ok in 0.817s
Mandelbrot: ok in 0.881s
Matmul: ok in 0.844s
Nbody: ok in 0.886s
Pidigits: ok in 0.788s
RegexDna: ok in 0.701s
Revcomp: ok in 0.981s
Spectralnorm: ok in 0.812s
Threadring: ok in 0.867s
Base64Encode: ok in 0.890s
Base64Decode: ok in 0.904s
Primes: ok in 1.076s
JsonGenerate: ok in 0.876s
JsonParsePure: ok in 0.942s
JsonParseSerializable: ok in 0.897s
JsonParsePull: ok in 0.889s
----
OK 21
24.1921s
```

## Crystal releases:

![Crystal releases](https://github.com/kostya/crystal-metric/blob/master/releases.png)
