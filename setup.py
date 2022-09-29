# pylint: disable=invalid-name
"""
    Setup script for magnetic reflectivity reduction application
"""

import sys
from setuptools import setup, find_packages
import os
from versioningit import get_cmdclasses

package_data = {
    "reflectivity_ui.interfaces.data_handling": [
        "genx_templates/*.gx",
    ],
    "reflectivity_ui.config": ["./settings.json"],
}

setup(
    name="reflectivity_ui",
    description="Magnetic Reflectivity Reduction",
    cmdclass=get_cmdclasses(),
    url="https://github.com/neutrons/reflectivity_ui",
    long_description="""Desktop application for magnetic reflectivity reduction""",
    license="Apache License 2.0",
    scripts=["bin/RefRedM", "bin/quicknxs2"],
    zip_safe=False,
    packages=find_packages(),
    package_dir={},
    package_data=package_data,
    install_requires=[
        "numpy",
        "matplotlib",
        "mantidworkbench",
    ],
    setup_requires=[],
)
