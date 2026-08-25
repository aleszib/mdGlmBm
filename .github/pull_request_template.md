## Summary

<!-- What changed? -->

## Task / scope

<!-- Link or name the OAP task file, e.g. instructions/tasks/001-package-foundation.md -->

## Statistical impact

- [ ] No statistical behavior changed
- [ ] Statistical behavior changed and is documented
- [ ] Objective scale / likelihood / deviance behavior changed and is tested
- [ ] PPML pseudo-likelihood interpretation is preserved where relevant

## Tests and checks

Commands run:

```text

```

Results:

```text

```

## Evidence

- [ ] Tests added or updated
- [ ] R CMD check run or reason documented
- [ ] CI expected to run
- [ ] No skipped test is reported as passed

## Safety / scope control

- [ ] Read `AGENTS.md`
- [ ] PR is small and task-scoped
- [ ] No unrelated broad refactoring
- [ ] No unrelated files changed
- [ ] No production secrets or credentials touched
- [ ] No direct merge to main
- [ ] Documentation updated if behavior changed

## Completion checklist

- [ ] Tests added or updated for every behavior change
- [ ] Documentation updated
- [ ] `Rscript -e "devtools::document()"` passed
- [ ] `Rscript -e "devtools::test()"` passed
- [ ] `_R_CHECK_FORCE_SUGGESTS_=false Rscript -e 'rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")'` passed
- [ ] Limitations documented
- [ ] Next task proposed

## Known limitations

<!-- Be explicit. -->

## Follow-up

<!-- Recommended next PR-sized task. -->
