@echo off
echo Building 'ExternalViewContext' sample...
devenv ExternalViewContext.sln /useenv /build "Debug|Win32" > BuildDebug.log
if ERRORLEVEL 1 goto error
devenv ExternalViewContext.sln /useenv /build "Release|Win32" > BuildRelease.log
if ERRORLEVEL 1 goto error
goto end
:error
ECHO *** Error building 'ExternalViewContext' sample - Exiting ***
:end
