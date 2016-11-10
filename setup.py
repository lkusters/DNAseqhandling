from distutils.core import setup

setup(name='DNAseqhandling',
      version='0.1dev',
      description='handling of DNA sequences',
      author='Lieneke Kusters',
      long_description=open('README.md').read(),
      packages=['DNAseqhandling'],
      install_requires=['BioPython'],
      )
