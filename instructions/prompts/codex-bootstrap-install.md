You are working in my GitHub-connected repository for the dynGLMbm R package project.

This is a bootstrap/instructions installation task, not the package implementation task yet.

Context:
- I added an OAP instruction ZIP file to this repository.
- Find the instruction ZIP. It is likely named something like `dynGLMbm_oap_repo_instructions_v2.zip`.
- Extract it into an `instructions/` folder.
- Read `instructions/AGENTS.md`, `instructions/docs/PROJECT_MEMORY.md`, `instructions/docs/MODEL_STRATEGY.md`, `instructions/docs/DATA_REPRESENTATION.md`, `instructions/docs/IMPLEMENTATION_PLAN.md`, and `instructions/tasks/000-bootstrap-instructions-ci.md`.
- Then execute Task 000 only.

Environment:
- You are running in a hardened WSL Ubuntu guest.
- You may use passwordless sudo for local system dependencies inside this guest.
- Do not access host credentials, production secrets, cloud credentials, or files outside this repository except for normal dependency installation inside the guest.

GitHub authentication:
- First inspect the remote with `git remote -v` and current auth state.
- If GitHub authentication fails, first try to refresh authentication using the available GitHub tooling, for example `gh auth status` and `gh auth refresh` if `gh` is installed.
- If authentication still fails, stop and ask me to authenticate. Do not invent credentials or switch to an unsafe workflow.

Required actions:
1. Confirm the working tree state with `git status`.
2. Locate the instruction ZIP.
3. Create or cleanly update `instructions/`.
4. Extract the ZIP into `instructions/`.
5. Copy `instructions/AGENTS.md` to repository-root `AGENTS.md`, unless an existing root `AGENTS.md` must be merged. Do not overwrite silently.
6. Copy `instructions/.github/workflows/R-CMD-check.yaml` to `.github/workflows/R-CMD-check.yaml`.
7. Copy `instructions/.github/pull_request_template.md` to `.github/pull_request_template.md`.
8. Create a branch named `chore/oap-instructions-ci` or similar.
9. Commit only the instruction and CI files.
10. Open a pull request if authentication allows it. Do not merge.

Required checks:
- `git status`
- `find instructions -maxdepth 3 -type f | sort`
- `find .github -maxdepth 3 -type f | sort`

Final report:
- Branch name
- Commit SHA
- PR URL if opened
- Installed files
- Whether GitHub authentication worked
- Any conflicts or files preserved
- Next recommended action, which should be to execute `instructions/tasks/001-package-foundation.md`
