

python3 -m build
python3 -m pip install --upgrade twine
python setup.py sdist bdist_wheel
python3 -m twine upload --skip-existing --repository testpypi dist/*