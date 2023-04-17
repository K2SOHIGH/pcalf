import os
import glob
from setuptools import setup, find_packages

with open("README.md", "r", encoding = "utf-8") as fh:
    long_description = fh.read()

setup(
    name='pcalf',
    version='1.0.0',
    description='Search calcyanin in a set of amino acid sequences',
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url='https://github.com/K2SOHIGH/pcalf',
    author='Maxime Millet',
    author_email='maxime.luc.millet@gmail.com',
    license='MIT',
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
         "pyhmmer==0.7.1",
         "biopython==1.81",
         "numpy>=1.21",
         "pyyaml>=6",
         "snakemake==7.22",
         "pandas==1.5.3",
         "tqdm==4.64.1",
    ],
    python_requires = ">=3.9",
    packages = find_packages(where="src"),
    package_dir = {"": "src"},
    include_package_data=True,
    scripts = [script for script in glob.glob("scripts/*") if not os.path.basename(script).startswith("_") ],
    py_modules = ["bin"],
    zip_safe=False
)