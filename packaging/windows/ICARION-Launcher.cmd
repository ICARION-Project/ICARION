@echo off
setlocal
set "SCRIPT_DIR=%~dp0"
start "" powershell.exe -NoProfile -ExecutionPolicy Bypass -STA -File "%SCRIPT_DIR%ICARION-Launcher.ps1"
