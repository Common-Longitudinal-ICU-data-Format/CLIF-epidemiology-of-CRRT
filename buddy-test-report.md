<!--
Buddy test report template. See guides/buddy-testing-guide.md for what each check means.
To use: copy this file to the project root as BUDDY_TEST_REPORT.md, fill it in, and commit it.
In the Result column, bold the outcome that applies (or delete the other), then add a note.
Pass with notes is only for the Overall verdict at the bottom.
-->

# Buddy Test Report

| | |
|---|---|
| **Buddy site / institution** | *Emory* |
| **Tester** | *Navya Ramesh* |
| **Date** | *2026-07-08* |


## Checks

| # | Check | Result | Notes |
|:-:|-------|--------|-------|
| 1 | Environment reproduces (`uv sync` / `00_renv_restore.R`, nothing by hand) | Pass  | |
| 2 | Configuration works from `config/README.md` alone; no hardcoding | Pass  | |
| 3 | Required tables/fields match what the code reads (mCIDE-valid) | Pass  | |
| 4 | Runs end to end with no manual edits between steps | Pass  | |
| 5 | Outputs in `output/final_no_phi/` with right naming/type, no raw dumps | Pass  | |
| 6 | **Data security**: no PHI, every stat n ≥ 10, no raw data *(blocking)* | Pass l | |
| 7 | Clinical sanity: aggregates plausible for the cohort | Pass  | |
| 8 | Documentation usable: could run from the README alone | Pass  | |

## Overall verdict

**Verdict:** Pass / Pass with notes / Fail  *(keep one)*

- **Pass** — runs cleanly, output sane and secure, docs followable. Ready to distribute.
- **Pass with notes** — works; address the notes below.
- **Fail** — at least one blocking issue (doesn't run, output wrong/insecure, docs unfollowable). A
  data-security failure (check 6) is always a Fail.

### Blocking issues (must fix before distribution)
1. *…*

### Non-blocking notes / suggestions
1. *…*

