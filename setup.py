from setuptools import find_packages, setup

NAME = 'primertool'
DESCRIPTION = 'Generating primers for sanger sequencing'
# URL = 'https://github.com/me/myproject'
EMAIL = 'ddey@ukaachen.de'
AUTHOR = 'Daniela Dey'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = None

with open('requirements.txt') as f:
    REQUIRED = f.read().splitlines()

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=open('README.md').read(),
      # url=URL,
      author=AUTHOR,
      author_email=EMAIL,
      license='',
      packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
      install_requires=REQUIRED,
      python_required=REQUIRES_PYTHON,
      )