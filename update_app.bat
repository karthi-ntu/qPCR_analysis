@echo off
echo =============================================
echo   qPCR App - Update and Deploy to GitHub
echo =============================================
echo.
cd /d "%~dp0"
call :find_rscript
if "%RSCRIPT%"=="" (
    echo ERROR: Rscript.exe not found. Install R or add it to PATH.
    pause
    exit /b 1
)

echo [1/4] Running tests...
"%RSCRIPT%" -e "testthat::test_dir('tests/testthat', reporter='summary', stop_on_failure=TRUE)"
if errorlevel 1 (
    echo ERROR: Tests failed. Fix them before deploying.
    pause
    exit /b 1
)
echo.

echo [2/4] Exporting app to static site (this takes ~30 seconds)...
"%RSCRIPT%" -e "shinylive::export(appdir = '.', destdir = 'docs')"
if errorlevel 1 (
    echo ERROR: Export failed. Is the shinylive package installed?
    pause
    exit /b 1
)
echo Export complete.
echo.

echo [3/4] Staging changes in git...
git add docs/ app.R R/ tests/ README.md
echo.

echo [4/4] Committing and pushing to GitHub...
set /p MSG="Commit message (Enter for default): "
if "%MSG%"=="" set "MSG=update: refresh app and Shinylive export"
git commit -m "%MSG%"
git push
echo.

echo =============================================
echo   Done! Your app is live at:
echo   https://karthi-ntu.github.io/qPCR_analysis/
echo   (wait ~60 seconds for GitHub to update)
echo =============================================
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
