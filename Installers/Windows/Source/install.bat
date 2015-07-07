mkdir "%HOME%\HybridCheck"
mkdir "%APPDATA%\Microsoft\Windows\Start Menu\Programs\HybridCheck"

C:\Windows\System32\robocopy "%cd%" "%HOME%\HybridCheck"

set SCRIPT="%TEMP%\%RANDOM%-%RANDOM%-%RANDOM%-%RANDOM%.vbs"
echo Set oWS = WScript.CreateObject("WScript.Shell") >> %SCRIPT%
echo sLinkFile = "%APPDATA%\Microsoft\Windows\Start Menu\Programs\HybridCheck\HybridCheck.lnk" >> %SCRIPT%
echo Set oLink = oWS.CreateShortcut(sLinkFile) >> %SCRIPT%
echo oLink.TargetPath = "%HOME%\HybridCheck\HybridCheck.bat" >> %SCRIPT%
echo oLink.WorkingDirectory = "%HOME\HybridCheck" >> %SCRIPT%
echo oLink.IconLocation = "%HOME%\HybridCheck\HybridCheck.icl, 0" >> %SCRIPT%
echo oLink.Save >> %SCRIPT%
C:\Windows\System32\cscript /nologo %SCRIPT%
del %SCRIPT%

set SCRIPT="%TEMP%\%RANDOM%-%RANDOM%-%RANDOM%-%RANDOM%.vbs"
echo Set oWS = WScript.CreateObject("WScript.Shell") >> %SCRIPT%
echo sLinkFile = "%APPDATA%\Microsoft\Windows\Start Menu\Programs\HybridCheck\Update.lnk" >> %SCRIPT%
echo Set oLink = oWS.CreateShortcut(sLinkFile) >> %SCRIPT%
echo oLink.TargetPath = "%HOME%\HybridCheck\update.bat" >> %SCRIPT%
echo oLink.WorkingDirectory = "%HOME%\HybridCheck" >> %SCRIPT%
echo oLink.IconLocation = "%HOME%\HybridCheck\HybridCheck.icl, 0" >> %SCRIPT%
echo oLink.Save >> %SCRIPT%
C:\Windows\System32\cscript /nologo %SCRIPT%
del %SCRIPT%

set SCRIPT="%TEMP%\%RANDOM%-%RANDOM%-%RANDOM%-%RANDOM%.vbs"
echo Set oWS = WScript.CreateObject("WScript.Shell") >> %SCRIPT%
echo sLinkFile = "%APPDATA%\Microsoft\Windows\Start Menu\Programs\HybridCheck\Uninstall.lnk" >> %SCRIPT%
echo Set oLink = oWS.CreateShortcut(sLinkFile) >> %SCRIPT%
echo oLink.TargetPath = "%HOME%\HybridCheck\remove.bat" >> %SCRIPT%
echo oLink.WorkingDirectory = "%HOME%\HybridCheck" >> %SCRIPT%
echo oLink.IconLocation = "%HOME%\HybridCheck\HybridCheck.icl, 0" >> %SCRIPT%
echo oLink.Save >> %SCRIPT%
C:\Windows\System32\cscript /nologo %SCRIPT%
del %SCRIPT%
