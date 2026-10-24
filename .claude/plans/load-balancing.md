# Load-balancing PR stack: working plan

Living plan for the parallel load-balancing work (see `.claude/issues/101.md`).
Keep this file updated as work progresses so it can be resumed after a break.

_Last updated: 2026-10-24_

## Goal

Get the load-balancing work (and the bug fixes found along the way) into a set
of pull requests against upstream `philchalmers/SimDesign` that:

1. are compatible with upstream `main` (now **2.28, released to CRAN 2026-10-06**);
2. are in line with the #101 discussion with the owner;
3. merge in the documented order **with zero conflicts**, with the final merged
   tree passing the full test suite.

## Branch stack

The branches are **stacked**: each one is built on the one above it, so the
PRs must be merged in this order. (Stacking is needed because every branch adds
entries to the same new `## Changes in SimDesign 2.29` NEWS section, which
conflicts if the branches are independent.)

| # | Branch | Content |
|---|--------|---------|
| 1 | `claude/remove-leftover-browser` | Removes leftover `browser()` in `Analysis()` (still present on upstream main). Opens the 2.29 NEWS section; dev version `2.28.1`. |
| 2 | `claude/fix-array-progress-default` | **New.** `runArraySimulation()` signature still says `progress = TRUE` although the 2.28 NEWS/docs (and #101) say `FALSE`; the mismatch shipped in CRAN 2.28. Fixes signature + Rd, adds test. |
| 3 | `claude/fix-parallel-export-lookup` | Fixes "object not found" when exporting non-global user functions through `runArraySimulation()`/`SimSolve()` in parallel. |
| 4 | `claude/fix-max-time-partial-results` | `max_time` crash fixes (NULL results, all-timed-out condition) + wall-clock deadlines on workers + no unbounded retry after an interrupt (SLURM safety). |
| 5 | `claude/simulation-load-balancing-06jgq4` | The load-balancing feature + hardening from the adversarial review. Contains `.claude/` notes. |

## Status

### Done
- [x] Merged upstream `main` (2.28) into the bottom of the stack; rebuilt the
      stack 1 → 5 with merge commits (no history rewritten).
- [x] Resolved conflicts: test-03 appended blocks (kept both); upstream's
      `missing(verbose)` handling (#105) next to the wall-clock `max_time.start`;
      DESCRIPTION version; NEWS restructured so all entries sit in 2.29.
- [x] LB docs now say "versions prior to 2.29" (was 2.27, which is released).
- [x] Reinstalled deps after container recycle: **mirai 2.7.3, pbapply 1.7.5**.
      Full suite passes on every branch (LB branch: 197 pass, 0 fail, 0 warn).
      Re-verified on mirai 2.7.3: crash watchdog, PSOCK crash teardown,
      export-clobber safety, re-entrant `runSimulation()`.
- [x] #101 experiments (see Findings).
- [x] Committed the bit-identity test (`test-07-load-balancing.R`):
      runArraySimulation ≡ runSimulation with `genSeeds(iseed, arrayID)`,
      serial and parallel. Full suite on branch 5: 199 pass, 0 fail, 0 warn.

### In progress
- (nothing)

### To do
- [x] Pushed all five branches (2026-10-24).
- [ ] Re-simulate the merge order 1 → 5 onto `origin/main` in a scratch clone;
      confirm zero conflicts; install the merged tree and run the full suite.
- [ ] Update the merge-order artifact page (now 5 PRs, stacked, 2.29).
- [ ] Before opening the upstream PR for branch 5: decide whether to drop the
      `.claude/` directory (issue notes + this plan) from the PR.

## Findings for the #101 discussion

- **`pboptions(use_lb = TRUE)` is reproducible, as the owner said** — but only
  because pbapply splits the work into chunks of exactly `ncores` items, so each
  worker still runs one replication per chunk. That is the same structure that
  causes the "waves", so `use_lb` cannot balance across them.
- **`use_lb` has no effect when `progress = FALSE`**: the legacy path then calls
  `parallel::parLapply()` directly. Since `runArraySimulation()` now defaults to
  `progress = FALSE`, the suggested workaround does nothing for array jobs.
- Benchmark (4 cores, 80 reps, 14 irregular 1 s stragglers, startup subtracted):
  upstream 6.8–7.2 s of compute in every configuration (`use_lb` included);
  LB branch 4.8–5.2 s; ideal 4.3 s.
- Bit-identity under the LB branch: runArraySimulation ≡ runSimulation (serial
  and parallel); parallel identical across core counts and mirai/PSOCK.
  **Serial ≠ parallel** (one continuous stream vs one stream per replication).

## Open decisions (for the user / owner)

- Make the serial path use per-replication streams too? This would make results
  bit-identical across all four [runSimulation, runArraySimulation] ×
  [serial, parallel] combinations, but changes serial results for every seeded
  user. Not done; needs a decision.

## Resuming after a break

```sh
# 1. dependencies (the container may have been recycled)
Rscript -e 'options(repos=c(P3M="https://packagemanager.posit.co/cran/__linux__/noble/latest"), Ncpus=8); install.packages(c("parallelly","dplyr","sessioninfo","beepr","pbapply","future","future.apply","progressr","R.utils","clipr","e1071","qs2","mirai","extraDistr","testthat"))'
# 2. install the checked-out branch and run the suite
R CMD INSTALL --no-docs --no-byte-compile .
Rscript -e "library('testthat',quietly=TRUE);library('SimDesign',quietly=TRUE);options(warn=2);test_dir('tests/tests')"
# 3. merge-order check: clone, checkout origin/main, merge branches 1..5 in order
```
