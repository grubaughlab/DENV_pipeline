from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from denv_pipeline import __version__, _program

setup(name='denv_pipeline',
      version=__version__,
      packages=find_packages(),
      scripts=[
            "denv_pipeline/scripts/denv_pipeline.smk",
            "denv_pipeline/scripts/DENV_MAPPER.sh",
            "denv_pipeline/scripts/mapper_done_slurm.sh"],
      description='Bioinformatic pipeline to generate reads and consensus sequences for DENV',
      package_data={"denv_pipeline":["primers/*"]},
      install_requires=["biopython>=1.70"],
      url='https://github.com/ViralVerity/DENV_pipeline',
      author='Verity Hill',
      author_email='verity.hill@yale.edu',
      entry_points="""
      [console_scripts]
      {program} = denv_pipeline.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
