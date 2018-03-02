all:
	@echo "Run make install to install the reflectivity reduction application"

check:
	# Check dependencies

install: compile_ui
	python setup.py install --user

compile_ui:
	python setup.py pyuic
	python setup.py pyrcc

test:
	cd $(prefix)/app; python manage.py test

.PHONY: check
.PHONY: install
