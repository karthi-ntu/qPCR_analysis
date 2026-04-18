@echo off
echo =============================================
echo   qPCR App - Update and Deploy to GitHub
echo =============================================
echo.

echo [1/3] Exporting app to static site (this takes ~30 seconds)...
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" -e "shinylive::export(appdir = 'E:/qPCR_analysis', destdir = 'E:/qPCR_analysis/docs')"
if errorlevel 1 (
    echo ERROR: Export failed. Is R installed?
    pause
    exit /b 1
)
echo Export complete.
echo.

echo [2/3] Staging changes in git...
cd /d "E:\qPCR_analysis"
git add docs/ app.R R/ tests/
echo.

echo [3/3] Committing and pushing to GitHub...
git commit -m "update: refresh app and Shinylive export"
git push
echo.

echo =============================================
echo   Done! Your app is live at:
echo   https://karthi-ntu.github.io/qpcr-analysis/
echo   (wait ~60 seconds for GitHub to update)
echo =============================================
pause
