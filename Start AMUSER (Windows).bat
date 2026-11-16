@echo off
rem ===================================================================
rem  Start the AMUSER GUI on Windows.
rem  Double-click this file. A browser tab opens automatically.
rem  Keep this window open while you work; close it to stop the server.
rem ===================================================================
cd /d "%~dp0"
where perl >nul 2>nul
if errorlevel 1 (
  echo Perl was not found on this PC. Install Strawberry Perl from
  echo https://strawberryperl.com/ and then double-click this file again.
  echo.
  pause
  exit /b 1
)
perl amuser_gui.pl
echo.
echo The AMUSER server has stopped.
pause
