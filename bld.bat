cd corvo-main
call gradlew.bat shadowJar
copy build\libs\*-all.jar ..\corvolauncher\resources
cd ..
python setup.py install
if errorlevel 1 exit 1
