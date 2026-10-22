# Authors, licenses, and modifications

Saturna modifications and original Python code: Copyright 2026 PyRverse
contributors. Distributed under GNU GPL version 3 or any later version; see
LICENSE. Saturna is independently maintained and is not an official upstream
satuRn or limma release.

## satuRn

The supported model/test workflow is adapted from satuRn 1.20.0, by Jeroen Gilis,
with contributors Kristoffer Vitting-Seerup, Koen Van den Berge, and Lieven Clement.
The standard version is available at
https://github.com/bioc/satuRn/tree/51f8364f4f1a0ca606bd5f7c6ad4caac47d3cfae
and https://bioconductor.org/packages/satuRn . Upstream is Artistic License 2.0;
the complete license is retained in licenses/Artistic-2.0.txt. This modified
version is named Saturna and made available under GPL-3.0-or-later, with source
and modifications available to recipients. Changes include a new Python API,
NumPy/SciPy implementation, explicit failure checks, dense bulk-only inputs,
and a separately authored empirical-null optimizer. Upstream plotting,
containers, and locfdr-derived routines are not distributed here.

## limma

src/saturna/_moderation.py is adapted from fitFDist.R, fitFDistUnequalDF1.R,
and squeezeVar.R in limma, by Gordon Smyth, Belinda Phipson, and Lizhong Chen
(as attributed in these source files). limma is licensed GPL (>= 2), and these
adaptations are distributed under GPL-3.0-or-later. The audited source was
limma 3.68.5, revision 825d1c8330d024adb6c54af66d79b0f4ec258530:
https://github.com/bioc/limma/tree/825d1c8330d024adb6c54af66d79b0f4ec258530/R .
The actual reference runtime was limma 3.68.4, revision
7c11dae9f0568f17d879280450e557f8851743d5. These routines were checked against
that runtime. Changes include NumPy arrays, SciPy special functions and bounded
optimization, explicit supported-input validation, and omission of robust and
covariate-trend estimators. Further upstream authors are listed in limma's
DESCRIPTION at the cited revision.

## Empirical-null implementation boundary

locfdr 1.1-8 (Bradley Efron, Brit Turnbull, Balasubramanian Narasimhan) declares
GPL-2 without a later-version grant in the inspected distribution. Its Newton
iteration and covariance code, including copies named locfdr_locmle and
locfdr_loccov in satuRn, have **not** been copied or translated into this package.
The original code in _empirical.py instead optimizes the standard conditional
truncated-normal log likelihood with SciPy. Upstream source was inspected during
the audit; no clean-room development claim is made. R fixtures execute the
installed reference packages separately; neither R nor locfdr is a runtime
dependency or bundled with Saturna.

The Saturna IRLS implementation uses standard grouped-binomial equations and
matches R GLM initialization and stopping conventions. It is newly written
Python code. Mathematical methods and third-party package names retain their
original attribution; this package makes no upstream endorsement claim.
