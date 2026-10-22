# Contributing

Install with `python -m pip install -e ".[dev]"`, then run pytest, Ruff, the
parity report, the example and a build. Changes to numerical behavior need a
meaningful regression case and an explicit support/provenance update. Preserve
input order and reject unsupported statistical settings; do not silently fall
back to another method. Review fixture changes against the R generator.

Keep runtime Python code independent of R. New upstream adaptations must retain
source revisions and compatible licenses. In particular, do not translate
locfdr's GPL-2-only routines into this GPL-3.0-or-later distribution.
