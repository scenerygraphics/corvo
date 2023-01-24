cd corvo-main
gradlew.bat shadowJar
copy build\libs\*-all.jar corvolauncher\resources
cd ..
"%PYTHON%" setup.py install
if errorlevel 1 exit 1