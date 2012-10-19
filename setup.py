#from distutils.core import setup
from setuptools import setup

setup(
    name='calc_LZC',
    version='0.1.0',
    author='N. E. Sanders',
    author_email='nsanders@cfa.harvard.edu',
    packages=['LZC'],
    package_dir={'LZC':'LZC'},
    #package_data={'LZC':['LZC/LZCpar.p']},
    data_files=[('LZC',['LZC/LZCpar.p'])],
    include_package_data=True,
    scripts=['LZC/LZC.py'],
    url='https://www.cfa.harvard.edu/~nsanders/papers/LZC/summary.htm',
    license='LICENSE.txt',
    description='Calculate photometric galaxy metallicities.',
    long_description=open('README.md').read(),
    install_requires=[
    ],
)