#!/usr/bin/env python
# -*- coding: utf-8 -*-

def gtbytes(data):
    data = str(data)
    try:
        data = bytes(data, 'utf8')
    except:
        data = bytes(data)
    return data
