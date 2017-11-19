@echo off
REM automatically generated
ECHO generated on host: megatron
ECHO generated on date: 2017-11-19
ECHO didehpc version: 0.1.1
ECHO context version: 0.1.2
ECHO running on: %COMPUTERNAME%
set CONTEXT_WORKDRIVE=Q:
set CONTEXT_WORKDIR=zika
set CONTEXT_ROOT=Q:\zika\contexts_test
set CONTEXT_ID=d16f65c692c98b612ca816016ba1f0eb
set CONTEXT_PROPAGATE_ERROR=TRUE
set CONTEXT_BOOTSTRAP=TRUE
call setr64_3_3_2.bat
ECHO mapping Q: -^> \\fi--san03\homes\jah113
net use Q: \\fi--san03\homes\jah113 /y
set REDIS_HOST=11.0.0.1
set REDIS_URL=redis://11.0.0.1:6379
%CONTEXT_WORKDRIVE%
cd \%CONTEXT_WORKDIR%
ECHO working directory: %CD%
ECHO this is a single task
set CONTEXT_TASK_ID=8a2655d3f1cd8103c514922764b5a602
set CONTEXT_LOGFILE=Q:\zika\contexts_test\logs\%CONTEXT_TASK_ID%
ECHO logfile: %CONTEXT_LOGFILE%
@REM The quoting here is necessary for paths with spaces.
ECHO on
Rscript "Q:\zika\contexts_test\bin\task_run" "%CONTEXT_ROOT%" %CONTEXT_TASK_ID% > "%CONTEXT_LOGFILE%" 2>&1
@ECHO off
if %ERRORLEVEL% neq 0 (
  ECHO Error running task
  EXIT /b %ERRORLEVEL%
)
@ECHO Quitting
