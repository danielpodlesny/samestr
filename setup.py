from setuptools import setup, find_packages
setup(name='samestr',
      version='2020.10',
      description='SameStr identifies shared strains between pairs of '
                  'metagenomic samples based on the similarity of SNV profiles.',
      author='Daniel Podlesny',
      license='GPLv3',
      packages=find_packages(),
      package_data={'': ['*.r', '*.R']},
      include_package_data=True)
