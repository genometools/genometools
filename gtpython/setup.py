#!/usr/bin/python
# -*- coding: utf-8 -*-

from distutils.core import setup

setup(name='GenomeTools', version='0.1', description=
      'Python bindings for GenomeTools', author='Sascha Steinbiss',
      author_email='steinbiss@zbh.uni-hamburg.de', url=
      'http://www.genometools.org ', packages=['gt', 'gt.core',
      'gt.annotationsketch', 'gt.extended'])
