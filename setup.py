from setuptools import setup

setup(
    name='Graph of Frags',
    version='0.1.0',
    packages=['chemicalgof'
              ],
    url='https://github.com/f48r1/chemicalgof',
    license='MIT',
    author='Fabrizio Mastrolorito',
    author_email='fabrizio.mastrolorito@uniba.it',
    description='Python package for molecular graph reduction to molecular fragment-based graph',
    install_requires=[
                'networkx',
                'numpy',
                'pandas',
                'rdkit',
        ],
    python_requires='>=3.10',
)