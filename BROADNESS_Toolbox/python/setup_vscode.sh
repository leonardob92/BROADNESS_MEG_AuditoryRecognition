#!/usr/bin/env bash
# Run from any directory: bash setup_vscode.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

python_command=""
for candidate in "${BROADNESS_PYTHON:-}" .venv/bin/python python3.11 python3.12 python3.13 python3.14 python3 python; do
    if [ -n "$candidate" ] && command -v "$candidate" >/dev/null 2>&1 && \
       "$candidate" -c 'import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))' >/dev/null 2>&1; then
        python_command="$candidate"
        break
    fi
done

if [ -z "$python_command" ]; then
    printf '%s\n' 'No compatible Python was found. BROADNESS currently supports Python 3.11 through 3.14.'
    printf '%s\n' 'Install Python 3.14 from https://www.python.org/downloads/release/python-31410/'
    printf '%s\n' 'Python 3.15 is currently too new for an easy h5py installation.'
    if command -v open >/dev/null 2>&1; then
        open 'https://www.python.org/downloads/release/python-31410/'
    elif command -v xdg-open >/dev/null 2>&1; then
        xdg-open 'https://www.python.org/downloads/release/python-31410/' >/dev/null 2>&1 || true
    fi
    printf '%s\n' 'After installation, reopen the terminal and run this script again.'
    printf '%s\n' 'If you already have a compatible Python, set BROADNESS_PYTHON to its full path.'
    exit 1
fi

printf 'Using %s: ' "$python_command"
"$python_command" --version

if [ -x .venv/bin/python ] && \
   ! .venv/bin/python -c 'import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))' >/dev/null 2>&1; then
    backup=".venv-previous-$(date +%Y%m%d%H%M%S)"
    mv .venv "$backup"
    printf 'Moved the incompatible environment to %s\n' "$backup"
fi

if [ ! -x .venv/bin/python ]; then
    "$python_command" -m venv .venv
fi

.venv/bin/python -m pip install --only-binary=:all: -e '.[vscode]'
.venv/bin/python -c 'import broadness, ipykernel; print("BROADNESS and the VS Code kernel are ready")'
printf '%s\n' 'In VS Code, select the Python interpreter in this folder: .venv/bin/python'
