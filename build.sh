#!/usr/bin/env bash
cd corvo-main
./gradlew shadowJar
cd ..
cp corvo-main/build/libs/corvo-0.1.0-SNAPSHOT-all.jar src/corvolauncher/resources/corvo-0.1.0-SNAPSHOT-all.jar
python -m build
