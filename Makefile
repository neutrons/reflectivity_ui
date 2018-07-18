all:
	@echo "Run make install to install the reflectivity reduction application"

check:
	# Check dependencies

local: compile_ui
	python setup.py install --user

install: compile_ui
	python setup.py install

compile_ui:
	python setup.py pyuic
	python setup.py pyrcc

rpm:
	python setup.py bdist_rpm

test:
	cd $(prefix)/app; python manage.py test

.PHONY: check
.PHONY: install
