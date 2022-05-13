#pylint: disable=invalid-name
"""
    Setup script for magnetic reflectivity reduction application
"""

import sys
from setuptools import setup, find_packages
import os
import reflectivity_ui

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

package_data = {"reflectivity_ui.interfaces.data_handling": ["genx_templates/*.gx",], "reflectivity_ui.config": ["./settings.json"]}

setup(name="reflectivity_ui",
      version=reflectivity_ui.__version__,
      description = "Magnetic Reflectivity Reduction",
      url = "https://github.com/mdoucet/reflectivity_ui",
      long_description = """Desktop application for magnetic reflectivity reduction""",
      license = "Apache License 2.0",
      scripts=["bin/RefRedM", "bin/quicknxs2"],
      zip_safe=False,
      packages=find_packages(),
      package_dir={},
      package_data=package_data,
      install_requires=['numpy','matplotlib'],
      setup_requires=[],
)
