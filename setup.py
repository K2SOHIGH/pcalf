import os
import glob
from setuptools import setup, find_packages


setup(
    name='pycalf',
    version='0.1',
    description='search calcyanin within a set of amino acid sequences',
    url='',
    author='Maxime Millet',
    author_email='maxime.luc.millet@gmail.com',
    license='MIT',
    packages=find_packages(),
    install_requires=[
         "pyhmmer",
         "biopython",
         "numpy==1.22.4",
         "pyyaml==6.0",
         "pandas",
    ],    
    
    include_package_data=True,
    scripts = [script for script in glob.glob("bin/*") if not os.path.basename(script).startswith("_") ],
    py_modules = ["bin"],
    zip_safe=False
)
