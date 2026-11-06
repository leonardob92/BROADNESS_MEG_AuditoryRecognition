@echo off
setlocal EnableDelayedExpansion
cd /d "%~dp0"

set "PYTHON_COMMAND="
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" -c "import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))" >nul 2>&1
    if not errorlevel 1 set "PYTHON_COMMAND=.venv\Scripts\python.exe"
)
if not defined PYTHON_COMMAND (
    py -3.11 -c "import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))" >nul 2>&1
    if not errorlevel 1 set "PYTHON_COMMAND=py -3.11"
)
if not defined PYTHON_COMMAND (
    py -3.12 -c "import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))" >nul 2>&1
    if not errorlevel 1 set "PYTHON_COMMAND=py -3.12"
)
if not defined PYTHON_COMMAND (
    py -3.13 -c "import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))" >nul 2>&1
    if not errorlevel 1 set "PYTHON_COMMAND=py -3.13"
)
if not defined PYTHON_COMMAND (
    py -3.14 -c "import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))" >nul 2>&1
    if not errorlevel 1 set "PYTHON_COMMAND=py -3.14"
)
if not defined PYTHON_COMMAND (
    py -3 -c "import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))" >nul 2>&1
    if not errorlevel 1 set "PYTHON_COMMAND=py -3"
)
if not defined PYTHON_COMMAND (
    python -c "import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))" >nul 2>&1
    if not errorlevel 1 set "PYTHON_COMMAND=python"
)

if not defined PYTHON_COMMAND (
    echo No compatible Python was found. BROADNESS currently supports Python 3.11 through 3.14.
    echo Install Python 3.14 from https://www.python.org/downloads/release/python-31410/
    start "" "https://www.python.org/downloads/release/python-31410/"
    echo After installation, reopen the terminal and run this script again.
    exit /b 1
)

echo Using:
%PYTHON_COMMAND% --version
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" -c "import sys; sys.exit(not ((3, 11) <= sys.version_info < (3, 15)))" >nul 2>&1
    if errorlevel 1 (
        set "BACKUP=.venv-previous-%RANDOM%"
        ren ".venv" "!BACKUP!"
        if errorlevel 1 exit /b 1
        echo Moved the incompatible environment to !BACKUP!
    )
)
if not exist ".venv\Scripts\python.exe" (
    %PYTHON_COMMAND% -m venv .venv
    if errorlevel 1 exit /b 1
)

".venv\Scripts\python.exe" -m pip install --only-binary=:all: -e ".[vscode]"
if errorlevel 1 exit /b 1
".venv\Scripts\python.exe" -c "import broadness, ipykernel; print('BROADNESS and the VS Code kernel are ready')"
if errorlevel 1 exit /b 1
echo In VS Code, select the Python interpreter in this folder: .venv\Scripts\python.exe
