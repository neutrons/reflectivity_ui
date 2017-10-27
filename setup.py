#pylint: disable=invalid-name
"""
    Setup script for magnetic reflectivity reduction application
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
from setuptools import setup, find_packages
import os

if 'pyuic' in sys.argv[:]:
    indir = 'designer'
    outdir = 'reflectivity_ui/interfaces/generated'
    files = os.listdir(indir)
    files = [os.path.join('designer', item) for item in files]
    files = [item for item in files if item.endswith('.ui')]

    done = 0
    for inname in files:
        outname = inname.replace('.ui', '.py')
        outname = outname.replace(indir, outdir)
        print("Converting '%s' to '%s'" % (inname, outname))
        command = "pyuic5 %s -o %s"  % (inname, outname)
        os.system(command)
        done += 1
    if not done:
        print("Did not convert any '.ui' files")
    sys.exit(0)

if 'pyrcc' in sys.argv[:]:
    infile = './icons/icons.qrc'
    assert os.path.isfile(infile)
    outfile = './reflectivity_ui/interfaces/generated/icons_rc.py'
    print("Converting icons_rc file:")
    command = "pyrcc5  %s -o %s" % (infile, outfile)
    print("> %s" %command)
    os.system(command)
    sys.exit(0)

setup(name="reflectivity_ui",
      version='0.0',
      description = "Magnetic Reflectivity Reduction",
      url = "https://github.com/mdoucet/reflectivity_ui",
      long_description = """Desktop application for magnetic reflectivity reduction""",
      license = "Apache License 2.0",
      scripts=["bin/RefRedM"],
      zip_safe=False,
      packages=find_packages(),
      package_dir={},
      install_requires=['numpy','matplotlib'],
      setup_requires=[],
)
