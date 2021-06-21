"""Setup script for the gtecs-obs package."""
import glob

from setuptools import setup, find_namespace_packages

REQUIRES = ['sqlalchemy>=1.2',
            'pymysql',
            'astropy',
            ]

setup(name='gtecs-obs',
      version='0',
      description='G-TeCS functions for observation scheduling',
      url='http://github.com/GOTO/goto-obsdb',
      author='Martin Dyer',
      author_email='martin.dyer@sheffield.ac.uk',
      install_requires=REQUIRES,
      packages=find_namespace_packages(include=['gtecs*']),
      package_data={'gtecs': ['obs/data/*']},
      scripts=glob.glob('scripts/*'),
      zip_safe=False,
      )
