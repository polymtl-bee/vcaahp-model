rem   Run this script to copy the files from TRNSYS directories to this repo,
rem   after making modifications but before staging and committing.

rem   Fortran 90 source code
xcopy "C:\TRNSYS18\SourceCode\UserTypes\Type3223.f90" "%CD%" /D
xcopy "C:\TRNSYS18\SourceCode\UserTypes\Type3254.f90" "%CD%" /D

rem   tmf proformas
xcopy "C:\TRNSYS18\Studio\Proformas\UserTypes\Type3223.tmf" "%CD%" /D
xcopy "C:\TRNSYS18\Studio\Proformas\UserTypes\Type3254.tmf" "%CD%" /D
