RD /S /Q "%HOME%\HybridCheck"
RD /S /Q "%APPDATA%\Microsoft\Windows\Start Menu\Programs\HybridCheck
set SCRIPT="%TEMP%\%RANDOM%-%RANDOM%-%RANDOM%-%RANDOM%.vbs"
echo MsgBox "Uninstalled Windows HybridCheck Launcher." >> %SCRIPT%
C:\Windows\System32\cscript /nologo %SCRIPT%
del %SCRIPT%
