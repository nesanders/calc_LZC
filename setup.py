from distutils.core import setup

setup(
    name='calc_LZC',
    version='0.1.0',
    author='N. E. Sanders',
    author_email='nsanders@cfa.harvard.edu',
    packages=['calc_LZC'],
    scripts=['calc_LZC/calc_LZC.py'],
    url='https://www.cfa.harvard.edu/~nsanders/papers/LZC/summary.htm',
    license='LICENSE.txt',
    description='Calculate photometric galaxy metallicities.',
    long_description=open('README.txt').read(),
    install_requires=[
    ],
)