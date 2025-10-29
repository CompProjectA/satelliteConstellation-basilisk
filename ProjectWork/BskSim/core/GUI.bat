@echo off
REM === Auto-launch Basilisk Fault Simulator GUI ===
cd /d "%~dp0"
call "..\..\..\..\basilisk-2.2.1\basilisk-env\Scripts\activate.bat"
python run_fault_simulator.py
pause
