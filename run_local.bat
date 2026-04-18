@echo off
echo =============================================
echo   qPCR App - Run Locally (for testing)
echo =============================================
echo.
echo Starting Shiny app...
echo The app will open in your browser.
echo Close this window or press Ctrl+C to stop.
echo.
"C:\Program Files\R\R-4.5.3\bin\Rscript.exe" -e "shiny::runApp('E:/qPCR_analysis/app.R', launch.browser=TRUE)"
pause
