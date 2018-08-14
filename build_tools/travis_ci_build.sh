#!/bin/sh
python setup.py install
coverage run test/data_handling_test.py
pylint --rcfile build_tools/pylint.rc -f parseable reflectivity_ui/interfaces/data_handling reflectivity_ui/interfaces/event_handlers
