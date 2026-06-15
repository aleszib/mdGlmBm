You are working in my GitHub-connected dynGLMbm R package repository.

Read these files first:
- `AGENTS.md`
- `instructions/docs/PROJECT_MEMORY.md`
- `instructions/docs/MODEL_STRATEGY.md`
- `instructions/docs/DATA_REPRESENTATION.md`
- `instructions/docs/IMPLEMENTATION_PLAN.md`
- `instructions/docs/TESTING_STRATEGY.md`
- `instructions/tasks/001-package-foundation.md`

Then execute Task 001 exactly as a PR-sized package-foundation task.

Important constraints:
- Use the existing `mdsbm` package as the package skeleton if the source material is available.
- Import selected static GLM blockmodeling functions, but do not implement the full dynamic GLM-Markov optimizer yet.
- Add tests.
- Preserve the long-term doctrine: mdsbm-style dynamic logic, independent time-specific GLMs, deviance-scale local scoring with `-2 log` Markov penalties, binomial and PPML first, PPML uses common result names and `pseudo = TRUE`.
- You are running in a hardened WSL Ubuntu guest and may use passwordless sudo for local system dependencies inside this guest.
- If GitHub authentication fails, first try to refresh authentication. If that fails, ask me to authenticate.

Workflow:
1. Start from current `main` unless the repo policy indicates otherwise.
2. Create a feature branch such as `feature/package-foundation`.
3. Make only Task 001 changes.
4. Run roxygen and package checks where possible.
5. Commit related files.
6. Open a PR if authentication allows it. Do not merge.

Final report must include:
- Branch
- Commit SHA(s)
- PR URL if opened
- Source material detected
- Files imported or created
- Functions exported
- Tests added
- Exact commands run and results
- Dependencies installed via sudo or R package installation
- Skipped or blocked checks
- Known risks
- Recommended next task
