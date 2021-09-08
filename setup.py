"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
from os import path
from io import open
from GENCOV import __version__

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='GENCOV',
      version=__version__,
      packages=find_packages(),
      package_data={"":["GENCOV/*", ]},
      description='SARS-CoV-2 snp-calling and consensus generation - amplicon - Illumina paired-ends reads',
      url='https://github.com/metagenlab/GENCOV',
      author='Trestan Pillonel',
      author_email='trestan.pillonel@gmail.com',
      entry_points="""
      [console_scripts]
      gencov = GENCOV.gencov:cli
      """,
      include_package_data=True,
      keywords=[],
      zip_safe=False,
      python_requires='>=3.6',
      install_requires=['snakemake>=5.2', 
                        'biopython>=1.77', 
                        'strictyaml',
                        'pyvcf', 
                        'biopython'],
    )


