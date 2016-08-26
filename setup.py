from __future__ import (print_function, absolute_import)
from setuptools import setup

setup(
    name="stark",
    version="0.1",
    description="Stark broadened line profiles from Lemke 1997",
    author_email="s.littlefair@shef.ac.uk",
    packages=['stark'],
    package_data={'': ['data/*']},
    include_package_data=True,
    zip_safe=False
)