cd corvo-main
call gradlew.bat shadowJar
copy build\libs\*-all.jar ..\src\corvolauncher\resources
cd ..
python -m build
if errorlevel 1 exit 1
