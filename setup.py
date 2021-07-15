from setuptools import setup, find_packages
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

setup(
    name='spacings',
    version='0.0.3',
    packages=find_packages(),
    license='Apache 2.0',
    author='Philipp Elleri, Lolian Shtembari',
    author_email='peller.phys@gmail.com',
    url='https://github.com/philippeller/spacings',
    description='Spacings',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={'rps_tables':['rps_tables/*.csv',]},
    python_requires='>=3.6',
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'scipy>=0.17',
        'numpy>=1.11',
    ],
)
