from setuptools import setup
import glob
from obsdb import __version__

setup(name='obsdb',
      version=__version__,
      description='Observation Database API for GOTO',
      url='http://github.com/GOTO/goto-obsdb',
      author='Martin Dyer',
      author_email='martin.dyer@sheffield.ac.uk',
      packages=['obsdb'],
      package_data={'': ['data/*']},
      install_requires=['sqlalchemy>=1.2', 'astropy'],
      scripts=glob.glob('scripts/*'),
      include_package_data=True,
      zip_safe=False)
