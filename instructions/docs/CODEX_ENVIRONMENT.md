# Codex Environment and Workflow

## Runtime context

This repository is worked on with Codex CLI using `--yolo` and the `gpt-5.4-Mini` model inside a hardened WSL Ubuntu guest.

Passwordless `sudo` may be used only to install local dependencies inside the hardened guest. Do not access host secrets, production credentials, cloud credentials, or unrelated files. Never invent credentials.

## Required workflow

- Read the repository-root `AGENTS.md` before doing work.
- Start from current `main`, check status and remotes, and check `gh auth status`.
- Use small, task-scoped branches and commits.
- Keep PRs small and do not merge PRs.
- Do not broaden the task or rewrite existing behavior without explicit instruction.
- If authentication fails, stop after safe authentication checks and ask the human to authenticate.
- Record any locally installed dependencies in the final report.

## Required final report

Every execution task reports the branch, commit SHA(s), PR URL if opened, files changed, exact tests and checks with results, local dependencies installed, limitations or skipped/blocked checks, known risks, and the recommended next task.
