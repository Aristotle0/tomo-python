from setuptools import setup, find_packages
import os
# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
setup(
    name = "tomopy",
    version = "0.1",
    packages = find_packages(),
    scripts = ['say_hello.py'],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires = ['docutils>=0.3'],

    # metadata for upload to PyPI
    author = "Lei Pan",
    author_email = "panlei7@gmail.com",
    description = "An additional package for tomo7",
    long_description = read('README'),
    keywords = "tomo7, FWI, model",

    # could also include long_description, download_url, classifiers, etc.
)