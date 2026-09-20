@echo off
setlocal
cd /d "%~dp0"

where python >nul 2>nul
if errorlevel 1 (
  echo Python was not found.
  echo Install Python 3.13.9 from https://www.python.org/downloads/ and try again.
  pause
  exit /b 1
)

if not exist ".venv\Scripts\python.exe" (
  echo Creating the TEMPO virtual environment...
  python -m venv .venv
  if errorlevel 1 goto :failed
)

".venv\Scripts\python.exe" -c "import numpy, scipy, pandas, matplotlib" >nul 2>nul
if errorlevel 1 (
  echo Installing the tested TEMPO dependencies...
  ".venv\Scripts\python.exe" -m pip install -r requirements.txt
  if errorlevel 1 goto :failed
)

echo Starting TEMPO Design Explorer...
".venv\Scripts\python.exe" app.py
if errorlevel 1 goto :failed
exit /b 0

:failed
echo.
echo TEMPO could not start. Review the message above or run test_installation.py.
pause
exit /b 1
