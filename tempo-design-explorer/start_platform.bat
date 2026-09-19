@echo off
cd /d "%~dp0"
python server.py --open
if errorlevel 1 (
  echo.
  echo TEMPO could not start. Check that Python, NumPy, SciPy, pandas and matplotlib are installed.
  pause
)
