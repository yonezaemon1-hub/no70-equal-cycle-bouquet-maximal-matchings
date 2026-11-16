# No.70 - Exact Nonreal-Root Counts for Maximal-Matching Polynomials of Equal-Cycle Bouquets

Ryutaro Yonezu - Independent Researcher. Manuscript version v1.0, 16 November 2026.
Unrefereed preprint. Human peer review and formal verification have not been performed.
NOVELTY INCONCLUSIVE / REFERENCES PARTIAL.
Copyright (c) 2026 Ryutaro Yonezu. All rights reserved. No reuse license is granted.

## Claim and scope

We study the polynomial enumerating inclusion-maximal matchings in a bouquet of equal-length cycles sharing one hub. For every integer q >= 5 and r >= 1, the proposed number of nonreal roots, counted with algebraic multiplicity, is 2(floor((q-1)/2)-floor(q/3)) floor((r-1)/2). The manuscript presents an analytical proof candidate based on a three-residue interlacing induction and a rational-function root count, accompanied by finite exact checks. This preprint has received the AI reviews described in its evidence statement; it has not undergone human peer review or formal verification. Novelty assessment: inconclusive; full-text comparison with Burnham-Erey (2025) incomplete.

The matching polynomial counts inclusion-maximal matchings. All-q/all-r statements rest on the proposed analytical proof; finite tests check finite cases. No general cut-vertex decomposition, uniqueness or priority claim is made.

## Files

- [Preprint PDF](No70_v1.0_PREPRINT.pdf), [LaTeX source](No70_v1.0_SOURCE.tex), and [Markdown source](No70_v1.0_SOURCE.md). Sections 1-4 preserve the adopted v0.14 mathematics byte-for-byte.
- [Reproducibility ZIP](No70_v1.0_REPRODUCIBILITY.zip): self-contained code/data layout, exact rational-interval certificates, provenance and strict checksums.
- [Release asset SHA256](No70_v1.0_SHA256SUMS.txt): checksums for PDF, TeX and reproducibility ZIP.
- Three unchanged finite-check source files are also exposed at the repository root for inspection. Use the ZIP's code/data layout for the complete reproduction driver.
- [Review scope](REVIEW_SCOPE.md), [prior-art limits](PRIOR_ART.md), and [rights](LICENSE.txt).

## Reproduce

Use Python 3 with the pinned SymPy version in requirements.txt. Run from the extracted package root:

```sh
python -m pip install -r requirements.txt
python -B verify_release.py
```

The default run checks every distributed file hash, checks all 196 stored q=5..200 certificates against the exact path polynomials, reruns the v0.9 finite algebra/interlacing test through q=45 (including 90 exact nonreal-root-count cases q=5..19, r=1..6), and reruns the v0.10 guidepost test through q=48. Regenerated results are compared with the frozen summaries; only SymPy version metadata is excluded from that summary comparison. Regenerated output stays in a temporary directory.

To rerun root isolation for the full certificate range, which can take longer:

```sh
python -B verify_release.py --full-intervals
```

The full option reproduces the three original interval batches q=5..140, q=141..180 and q=181..200 in a temporary directory and compares the complete mathematical certificate records. The default certificate check instead uses exact endpoint signs or exact point roots in disjoint rational intervals; the number of intervals exhausts each polynomial's degree, proving that these real roots are simple. It also checks the prescribed root order. Neither route is a formal proof of the all-q induction.


## Review and operating-system evidence

The author supplied separate AI reports for v0.10 proof review, v0.11-v0.13 editorial deltas, and v0.14 mathematical-delta PASS. The v0.14 report confirms the r=1 pole-evaluation correction, verifies the exact patch, and reports no new defect within its stated review scope. The old v0.13 sentence P(c_i)=a(c_i)^r failed at r=1; this error and its repair are disclosed in Section 5. Prior AI reviews missed that sentence. AI review is not a guarantee of correctness.

Native Windows v0.12 PASS (exit code 0) is attributed to the author's handoff. No fresh native Windows v0.14/v1.0 run is claimed. Preparation independently reran finite checks on Linux. Raw review and Windows execution logs are excluded from this distribution. No human peer review or formal proof certificate is claimed.

## PDF build and rights

Compile No70_v1.0_SOURCE.tex twice with XeLaTeX. This PDF was newly generated from v0.14 mathematical text; the superseded v1.0 PDF is not used. PDF metadata carries title, author, subject and keywords.

Public availability does not grant redistribution, modification or software reuse rights. LICENSE.txt retains the original All rights reserved statement; no Creative Commons or open-source license is granted.
