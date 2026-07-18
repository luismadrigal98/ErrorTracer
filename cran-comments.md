This is a resubmission of an archived package (1.2.1). ErrorTracer was
archived on 2026-06-08 because a test failure on macOS arm64 was not
corrected. This version fixes that failure. No reverse dependencies were
detected on CRAN.

## The archival issue, and how it is fixed

**Cause.** In the uncertainty decomposition, the "no environmental noise"
case was detected by testing a computed variance difference for exact
equality with zero (`v_perturbed - v_param_sub == 0`), and the two column
variances were computed with different `stats::var()` `na.rm` settings.
Those settings select different C code paths in `cov.c`
("na.or.complete" vs "everything"), which round identically on x86_64 but
not on aarch64. A last-bit difference therefore flipped the reported Monte
Carlo standard error of `env_var` from 0 to 0.20, failing
`test-v-env-stability.R:53` on the two `r-*-macos-arm64` flavours.

**Fix.** Both variances are now computed the same way, and the zero test
uses a relative tolerance (`sqrt(.Machine$double.eps)`, the `all.equal()`
convention) rather than exact floating-point equality, so the result no
longer depends on bit-level floating-point reproducibility across
platforms. Two regression tests were added, one of which reproduces the
aarch64 behaviour on any platform by introducing a ULP-scale difference
between the two inputs.

The third ERROR on the archived check page
(r-devel-linux-x86_64-debian-gcc) was `Package required but not
available: 'brms'`, i.e. a missing dependency in that check environment
rather than a defect in this package. `brms` is declared in Imports and
the package checks cleanly where it is installed.

## Other changes since the last CRAN release (1.1.0)

See NEWS.md for the full 1.2.0 entry. In summary: the variance budget is
now internally consistent (`total_var` is defined as the sum of its
components, so shares sum to 100%), covariate error can be propagated
pathwise via a driver ensemble, the autocorrelation/temporal variance
component was corrected and made tail-robust, and `shelf_life()` gained a
null-model skill gate.

## Test environments

* local Linux install (Ubuntu), R 4.5.3
* win-builder (devel and release)

The archival failure was specific to aarch64, which I cannot test
directly. It was addressed by making the affected code path independent
of bit-exact floating-point equality, and by a regression test that
reproduces the divergence on x86_64.

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the expected resubmission note:

```
New submission
Package was archived on CRAN
```

`R CMD check --as-cran --run-donttest` executes all \donttest{}
examples successfully.

## Reverse dependencies

None on CRAN.
