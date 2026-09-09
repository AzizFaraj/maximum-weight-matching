# Public-release checklist

## Required before making the repository public

- [ ] Freeze the final Edmonds, Gabow, and Scaling source revision.
- [ ] Confirm that the release source matches the final HPC binaries by source
      manifest/tree digest.
- [ ] Complete the direct-MTX LEMON validation tracked as Slurm job `2323302`.
- [ ] Import the final Suite A and Suite B manifests, scripts, compact summaries,
      and reproduction commands.
- [ ] Verify the exact upstream Matchbox revision, retrieval date, attribution,
      and license provenance in `NOTICE.md`.
- [ ] Add verified human author names, affiliations, and optional ORCIDs to
      `AUTHORS.md` and `CITATION.cff`.
- [ ] Add modification notices to materially changed inherited source files.
- [ ] Confirm that no unlicensed dataset, downloaded paper, credential,
      machine-specific secret, or personal path is present.
- [ ] Run Release, Debug, ASan, and UBSan correctness gates.
- [ ] Build and test an exported `git archive` in a fresh directory.
- [ ] Confirm all documented commands on GNU/Linux.
- [ ] Create the approved GitHub remote, enable branch protection, and tag the
      archived paper release.

## Intentionally deferred until final experiments

- final result tables and operation-count summaries;
- archival dataset/result DOI;
- paper DOI and preferred citation;
- final release version and date.
