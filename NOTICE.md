# Provenance and third-party notice

This is a **pre-release provenance record**. The repository must remain private
until the unresolved fields below are verified.

## Matchbox base

This software is derived from the open-source Matchbox codebase distributed
under the GNU General Public License version 3. The local research tree was
named `matchbox-10.12` and was imported into the current Git history in 2026.

Candidate upstream repository:
<https://github.com/CSCsw/matchbox>

Before public release, record all of the following from primary evidence:

- exact upstream revision, tag, or source archive;
- retrieval date and archive checksum;
- upstream attribution exactly as supplied by upstream; and
- verified names of the human authors responsible for the 2026 modifications.

Do not infer or invent a copyright holder from repository ownership.

## Local modifications

The research version adds and modifies weighted-matching implementations,
including the Gabow and Scaling paths, correctness gates, graph generators,
instrumented experiment drivers, and HPC reproduction scripts. The final
release will identify the exact frozen commit and source-tree digest used for
the published experiments.

## External validation software

LEMON 1.3.1 is used only as an optional independently built validation oracle.
LEMON is not included in this repository and retains its own copyright and
license terms.

## Graph datasets

Large real graph inputs are not distributed here. The final dataset manifest
must give each source URL, citation, license or redistribution status,
checksum, and preprocessing/normalization rule. Small files under
`tests/fixtures/` are project-generated regression fixtures.
