cd corvo-main
call gradlew.bat shadowJar
move build\libs\*-all.jar ..\src\corvolauncher\resources
cd ..
python -m pip install --no-deps --ignore-installed .
if errorlevel 1 exit 1
