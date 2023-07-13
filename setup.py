"""Minimal setup file for tasks project."""

from setuptools import setup, find_packages

setup(
    name='tasks',
    version='0.1.0',
    license='proprietary',
    description='Module Experiment',

    author='greenteabiscuit',
    author_email='hogehoge@gmail.com',
    url='None.com',

    packages=find_packages(where='src'),
    package_dir={'': 'src'},
)
