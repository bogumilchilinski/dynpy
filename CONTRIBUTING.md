# Contributing to DynPy

## Workflow

1. **Issue first.** Describe the change in an issue, using the usual title form,
   e.g. "Implementation of `X` class that ... (defined in `dynpy.y` module)".
2. **Branch from the current `master`.** Name it `<issue-id>-<short-description>`,
   e.g. `1208-fix-report-cache`. One issue per branch.
3. **Commit messages** say what was done and reference the issue,
   e.g. "Fix of `ReportCache` serving stale figures (within #1208 issue)".
4. **Test locally** before opening a pull request (see below).
5. **Open a pull request to `master`** and fill in the template. At least one other person
   with write access reviews and approves it, then it is merged and the branch is deleted.

   This is an agreed practice. GitHub does not enforce it yet; it will once branch
   protection is enabled for `master`.

## Tests

Tests live in `tests/` and use pytest. From the parent folder of `dynpy`:

    python -m pytest dynpy/tests -v

There is no CI yet, so state the result in the pull request, together with any manual
checks (notebooks, scripts, generated documents).

## Do not commit

- `__pycache__/`, `.pytest_cache/`, build outputs (`output/`, `tikzplots/`, PDFs)
- notebook outputs and local data
- unrelated changes: open a separate pull request for each topic
