"""Setup script for the obsdb package."""
import glob

from setuptools import setup

PACKAGES = ['obsdb']

REQUIRES = ['sqlalchemy>=1.2',
            'pymysql',
            'astropy',
            'configobj',
            ]

# Get the version string
__version__ = None
with open('obsdb/version.py') as f:
    exec(f.read())  # Should set __version__

setup(name='obsdb',
      version=__version__,
      description='Observation Database API for GOTO',
      url='http://github.com/GOTO/goto-obsdb',
      author='Martin Dyer',
      author_email='martin.dyer@sheffield.ac.uk',
      install_requires=REQUIRES,
      packages=PACKAGES,
      package_data={'': ['data/*']},
      include_package_data=True,
      scripts=glob.glob('scripts/*'),
      zip_safe=False,
      )
