from setuptools import setup, find_packages
from pathlib import Path

requirements = Path("requirements.txt").read_text().splitlines()

readme = Path("README.md").read_text() if Path("README.md").exists() else ""

setup(
    name='chemicalgof',
    version='0.2.0',
    packages=find_packages(),
    url='https://github.com/f48r1/chemicalgof',
    license='MIT',
    author='Fabrizio Mastrolorito',
    author_email='fabrizio.mastrolorito@uniba.it',
    description='fragSMILES, fragment-based chemical notation from graph reduction process.',
    long_description=readme,
    long_description_content_type='text/markdown',
    install_requires=requirements,
    python_requires='>=3.10',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: AI :: Chemistry',
    ],
)
