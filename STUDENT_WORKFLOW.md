# Student Workflow for Editing Code on GitHub

**BENG207 • Winter/Spring 2026 • Version 2026-02-17**

This guide describes the required workflow for contributing code changes to the course repository using GitHub. The goal is to keep changes reviewable, reproducible, and safe for the project history.

## Core rules (read first)

- Do not commit directly to `main`. Always work in a branch and open a Pull Request (PR).
- Keep PRs small and focused: one logical change per PR.
- Before requesting review, run the project’s checks (tests / lint / a basic run) and confirm outputs.
- Never commit secrets (tokens, passwords), large raw datasets, or patient/private data.

## Workflow (local development)

1. Clone the repository (once).
2. Create a branch for each task.
3. Edit code locally; keep changes scoped.
4. Run checks (tests / lint / minimal example run).
5. Commit with a clear message.
6. Push your branch and open a PR.

## Commands (copy/paste)

```bash
git clone <REPO_URL>
cd <REPO_NAME>
git checkout main
git pull
git checkout -b <initials>/<short-description>
# ...edit files...
git status
git diff
# run project checks (examples)
# python -m pytest
# python -m ruff check .
# python -m black .
git add <files>
git commit -m "Concise description of the change"
git push -u origin <initials>/<short-description>
```
