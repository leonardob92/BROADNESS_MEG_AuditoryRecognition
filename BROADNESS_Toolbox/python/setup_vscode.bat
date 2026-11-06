@echo off
setlocal
cd /d "%~dp0"

set "PYTHON_COMMAND="
py -3.11 -c "import sys; sys.exit(sys.version_info < (3, 11))" >nul 2>&1
if not errorlevel 1 set "PYTHON_COMMAND=py -3.11"
if not defined PYTHON_COMMAND (
    py -3 -c "import sys; sys.exit(sys.version_info < (3, 11))" >nul 2>&1
    if not errorlevel 1 set "PYTHON_COMMAND=py -3"
)
if not defined PYTHON_COMMAND (
    python -c "import sys; sys.exit(sys.version_info < (3, 11))" >nul 2>&1
    if not errorlevel 1 set "PYTHON_COMMAND=python"
)

if not defined PYTHON_COMMAND (
    echo Python 3.11 or newer was not found.
    echo Download and install Python from https://www.python.org/downloads/
    start "" "https://www.python.org/downloads/"
    echo After installation, reopen the terminal and run this script again.
    exit /b 1
)

echo Using:
%PYTHON_COMMAND% --version
if not exist ".venv\Scripts\python.exe" (
    %PYTHON_COMMAND% -m venv .venv
    if errorlevel 1 exit /b 1
)

".venv\Scripts\python.exe" -m pip install -e ".[vscode]"
if errorlevel 1 exit /b 1
".venv\Scripts\python.exe" -c "import broadness, ipykernel; print('BROADNESS and the VS Code kernel are ready')"
if errorlevel 1 exit /b 1
echo In VS Code, select the Python interpreter in this folder: .venv\Scripts\python.exe
