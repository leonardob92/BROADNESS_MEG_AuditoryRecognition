#!/usr/bin/env bash
# Run from any directory: bash setup_vscode.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

python_command=""
for candidate in python3.11 python3 python; do
    if command -v "$candidate" >/dev/null 2>&1 && \
       "$candidate" -c 'import sys; sys.exit(sys.version_info < (3, 11))' >/dev/null 2>&1; then
        python_command="$candidate"
        break
    fi
done

if [ -z "$python_command" ]; then
    printf '%s\n' 'Python 3.11 or newer was not found.'
    printf '%s\n' 'Download and install Python from https://www.python.org/downloads/'
    if command -v open >/dev/null 2>&1; then
        open 'https://www.python.org/downloads/'
    elif command -v xdg-open >/dev/null 2>&1; then
        xdg-open 'https://www.python.org/downloads/' >/dev/null 2>&1 || true
    fi
    printf '%s\n' 'After installation, reopen the terminal and run this script again.'
    exit 1
fi

printf 'Using %s: ' "$python_command"
"$python_command" --version

if [ ! -x .venv/bin/python ]; then
    "$python_command" -m venv .venv
fi

.venv/bin/python -m pip install -e '.[vscode]'
.venv/bin/python -c 'import broadness, ipykernel; print("BROADNESS and the VS Code kernel are ready")'
printf '%s\n' 'In VS Code, select the Python interpreter in this folder: .venv/bin/python'
