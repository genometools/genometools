#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Parse getopt-style help texts for options
# and generate zsh(1) completion function.
# http://github.com/RobSis/zsh-completion-generator

# Usage: program --help | ./help2comp.py program_name

# Changelog:
# 9th September, 2014:
#   - Copied from http://github.com/RobSis/zsh-completion-generator.
#   - Print argument list instead of completion function.

import os
import sys
import re
import argparse
from string import Template


URL = 'http://github.com/RobSis/zsh-completion-generator'
STRIP_CHARS = " \t\n,="

COMPLETE_FUNCTION_TEMPLATE = """
#compdef $program_name

# zsh completions for '$program_name'
# automatically generated with $url
local arguments

arguments=(
$argument_list
    '*:filename:_files'
)

_arguments -s $arguments
"""

ARGUMENT_TEMPLATE = """    {$opts}'[$description]$style'"""
SINGLE_ARGUMENT_TEMPLATE = """    '$opt[$description]$style'"""


def cut_option(line):
    """
    Cuts out the first option (short or long) and its argument.
    """
    # TODO: dare to make it regex-free?
    newline = line.strip(STRIP_CHARS)
    opt = re.findall(r'^(-[a-zA-Z0-9\-]+(?:[\[\ =][^\-\ ][a-zA-Z\<\>\[\|\:\]\-\_\?#]*\]?)?)', line)
    if len(opt) > 0:
        newline = line.replace(opt[0], "", 1).strip(STRIP_CHARS)
        # return without parameter
        return newline, re.split('[\ \[=]', opt[0])[0]
    else:
        return newline, None


def parse_options(help_text):
    """
    Parses the options line by line.
    When description is missing and options are missing on
    consecutive line, link them together.
    """
    all_options = []
    previous_description_missing = False
    for line in help_text:
        line = line.strip(STRIP_CHARS)
        if re.match(r'^--?[a-zA-Z0-9]+', line) != None:  # starts with option
            previous_description_missing = False
            options = []
            while True:
                line, opt = cut_option(line)
                if opt == None:
                    break

                options.append(opt)

            if (len(line) == 0):
                previous_description_missing = True

            options.append(line)
            all_options.append(options)
        elif previous_description_missing:
            all_options[-1][-1] = line
            previous_description_missing = False

    return all_options


def _escape(line):
    """
    Escape the syntax-breaking characters.
    """
    line = line.replace('[','\[').replace(']','\]')
    line = re.sub('\'', '', line)  # ' is unescapable afaik
    return line


def generate_argument_list(options):
    """
    Generate list of arguments from the template.
    """
    argument_list = []
    for opts in options:
        model = {}
        # remove unescapable chars.

        model['description'] = _escape(opts[-1])
        model['style'] = ""
        if (len(opts) > 2):
            model['opts'] = ",".join(opts[:-1])
            argument_list.append(Template(ARGUMENT_TEMPLATE).safe_substitute(model))
        elif (len(opts) == 2):
            model['opt'] = opts[0]
            argument_list.append(Template(SINGLE_ARGUMENT_TEMPLATE).safe_substitute(model))
        else:
            pass

    return "\n".join(argument_list)


def generate_completion_function(options, program_name):
    """
    Generate completion function from the template.
    """
    model = {}
    model['program_name'] = program_name
    model['argument_list'] = generate_argument_list(options)
    model['url'] = URL
    return Template(COMPLETE_FUNCTION_TEMPLATE).safe_substitute(model).strip()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        options = parse_options(sys.stdin.readlines())
        if (len(options) == 0):
            sys.exit(2)

        #print generate_completion_function(options, sys.argv[1])
        print generate_argument_list(options)
    else:
        print "Please specify program name."
