@echo off
echo =============================================
echo   qPCR App - Run Locally (for testing)
echo =============================================
echo.
echo Starting Shiny app...
echo The app will open in your browser.
echo Close this window or press Ctrl+C to stop.
echo.
cd /d "%~dp0"
call :find_rscript
if "%RSCRIPT%"=="" (
    echo ERROR: Rscript.exe not found. Install R or add it to PATH.
    pause
    exit /b 1
)
"%RSCRIPT%" -e "shiny::runApp('.', launch.browser=TRUE)"
pause
exit /b 0

:find_rscript
set "RSCRIPT="
for %%R in (Rscript.exe) do set "RSCRIPT=%%~$PATH:R"
if not "%RSCRIPT%"=="" exit /b 0
for /d %%D in ("C:\Program Files\R\R-*") do set "RSCRIPT=%%~D\bin\Rscript.exe"
if exist "%RSCRIPT%" exit /b 0
set "RSCRIPT="
exit /b 0
