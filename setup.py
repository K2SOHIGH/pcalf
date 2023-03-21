import os
import glob
from setuptools import setup, find_packages


setup(
    name='pcalf',
    packages=['pcalf'],
    package_dir={'':'src'},
    version='2.0',
    description='Search calcyanin in a set of amino acid sequences',
    url='',
    author='Maxime Millet',
    author_email='maxime.luc.millet@gmail.com',
    license='MIT',
    # packages=find_packages(),
    install_requires=[
         "pyhmmer==0.7.1",
         "biopython==1.81",
         "numpy==1.22.4",
         "pyyaml==6.0",
         "snakemake==7.22",
         "pandas==1.5.3",
         "tqdm==4.64.1",
    ],    
    include_package_data=True,
    scripts = [script for script in glob.glob("scripts/*") if not os.path.basename(script).startswith("_") ],
    py_modules = ["bin"],
    zip_safe=False
)
