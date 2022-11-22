import sys
from setuptools import setup


if sys.version_info < (3, 7):
    sys.stderr.write(
        f'You are using Python '
        + "{'.'.join(str(v) for v in sys.version_info[:3])}.\n\n"
        + 'copylot only supports Python 3.9 and above.\n\n'
        + 'Please install Python 3.9 using:\n'
        + '  $ pip install python==3.9\n\n'
    )
    sys.exit(1)

setup(
    use_scm_version={"write_to": "corvolauncher/_version.py"},
    setup_requires=['setuptools_scm'],
    entry_points={'console_scripts': ['corvo=corvolauncher.gui.qt_main_window:main']},
    data_files=[("resources", ["corvolauncher/resources/corvo-0.1.0-SNAPSHOT-all.jar"])],
    packages=['corvolauncher',
              'corvolauncher.gui',
              'corvolauncher.gui.job_runners',
              'corvolauncher.gui.test',
              'corvolauncher.utilities',
              'corvolauncher.utilities.test'],
)
