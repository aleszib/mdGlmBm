# Task 000: Install OAP instructions and CI

## Goal

Install the OAP project instructions and GitHub CI files into the repository without changing package source code.

## Scope

- Extract the provided instruction ZIP into `instructions/`.
- Copy or synchronize the root `AGENTS.md` from `instructions/AGENTS.md` to repository root.
- Install GitHub Actions workflow files from `instructions/.github/workflows/` into `.github/workflows/`.
- Install the pull request template from `instructions/.github/pull_request_template.md` into `.github/pull_request_template.md`.
- Commit these instruction and CI files on a feature branch.
- Open a pull request if GitHub authentication permits.

## Non-goals

- Do not modify R package source code.
- Do not import `mdsbm` or `GLM_blockmodeling` code yet.
- Do not implement model functions.
- Do not run long package checks if no package exists yet.

## Required checks

Run:

```sh
git status
find instructions -maxdepth 3 -type f | sort
find .github -maxdepth 3 -type f | sort
```

If the repository already has `.github/workflows` or `AGENTS.md`, do not overwrite silently. Preserve existing files or merge carefully, and report what happened.

## Final report

Report branch, commit, PR URL if created, installed files, and any authentication issues.
