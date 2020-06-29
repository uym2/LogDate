from setuptools import setup, find_packages
import logdate
from os import walk, listdir
from os.path import join,normpath,isfile

param = {
    'name': logdate.PROGRAM_NAME,
    'version': logdate.PROGRAM_VERSION,
    'description': logdate.PROGRAM_DESCRIPTION,
    'author': logdate.PROGRAM_AUTHOR,
    'author_email': ['umai@eng.ucsd.edu'],
    'license': logdate.PROGRAM_LICENSE,
    'packages': find_packages(),
    'include_package_data': True,
    'scripts' : ['launch_LogDate.py'],
    'zip_safe': True,
    'install_requires': ['dendropy>=4.3.0','scipy>=1.3.1','bitsets>=0.7.15'],
    'keywords': 'Phylogenetics Evolution Biology',
    'long_description': """A Python implementation of the wLogDate algorithm""",
    'classifiers': ["Environment :: Console",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "License :: OSI Approved :: GNU General Public License (GPL)",
                    "Natural Language :: English",
                    "Operating System :: OS Independent",
                    "Programming Language :: Python",
                    "Topic :: Scientific/Engineering :: Bio-Informatics",
                    ],
    }
    
setup(**param)
