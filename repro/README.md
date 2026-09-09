# Reproducing the experiments

The final artifact will contain two explicitly separated campaigns.

## Suite A: unrestricted MWM

Suite A compares Edmonds and Gabow. It contains selected real matrices from
SuiteSparse together with synthetic Erdos--Renyi and Barabasi--Albert graphs.
Both algorithms receive the same normalized graph and stored weights. Scaling
is excluded when the graph has no verified perfect optimum because its current
entry point solves MWPM.

## Suite B: common MWM/MWPM objective

Suite B contains controlled synthetic families. Candidate graphs are generated
from frozen topology and weight seeds, then independently screened with LEMON
for both unrestricted MWM and MWPM. A graph enters the three-way timing set only
when the two exact objective values are equal. Edmonds, Gabow, and Scaling then
receive identical graph bytes and weights.

The final scripts will record:

- source commit and source-tree digest;
- compiler, flags, linked libraries, CPU, memory, OS, and Slurm version;
- graph source, checksum, normalized digest, dimensions, family, and seeds;
- algorithm order, warmup policy, trial, timeout, wall time, peak memory, and
  correctness certificate; and
- Scaling configuration and operation counters.

Raw datasets and multi-gigabyte result archives will be deposited separately;
this repository will retain manifests, compact summaries, and retrieval
instructions.

See [`HPC_JOB_LOG.md`](HPC_JOB_LOG.md) for the current campaign bookkeeping.
