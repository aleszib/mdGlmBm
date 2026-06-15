# dynGLMbm OAP repository instruction pack

This ZIP is meant to be added to your GitHub-connected repository and installed by Codex CLI.

Recommended bootstrap workflow:

1. Add this ZIP to the repository root.
2. Start Codex CLI in your hardened WSL Ubuntu working folder.
3. Paste the prompt from `prompts/codex-bootstrap-install.md`.
4. Review the PR that installs instructions and CI.
5. After that, run the first implementation task using `prompts/codex-task-001-package-foundation.md`.

The project doctrine encoded here is:

- mdsbm-style dynamic membership logic;
- GLM-based static/within-time observation models;
- independent GLMs per time point in the first dynamic model;
- actor entry/exit support through actor-time lineage;
- binomial and PPML first;
- same result names across model families;
- `pseudo = TRUE` for PPML;
- R reference backend first;
- C++ scoring/update backend later;
- tests required for implementation PRs.
