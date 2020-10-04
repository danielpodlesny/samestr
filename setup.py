from setuptools import setup, find_packages

setup(name='samestr',
      version='2020.10',
      description='SameStr identifies shared strains between pairs of '
                  'metagenomic samples based on the similarity of SNV profiles.',
      author='Daniel Podlesny',
      author_email='daniel.podlesny@uni-hohenheim.de',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      url='https://github.com/danielpodlesny/samestr/',
      license=open('LICENSE').read(),
      packages=find_packages(),
      package_data={'': ['samestr','kpileup.pl', '*.R', 'LICENSE']},
      include_package_data=True)
