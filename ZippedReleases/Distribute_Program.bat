@echo off
echo Copying MSGFPlus.jar to C:\DMS_Programs\ and \\Proto-3\DMS_Programs_Dist\AnalysisToolManagerDistribution\
pause

@echo on
xcopy ..\target\msgfplus.jar C:\DMS_Programs\MSGFPlus\ /Y
xcopy ..\target\msgfplus.jar \\Proto-3\DMS_Programs_Dist\AnalysisToolManagerDistribution\MSGFPlus\ /Y

@echo off
pause
