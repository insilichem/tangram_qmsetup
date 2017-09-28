#!/usr/bin/env python
# encoding: utf-8


from __future__ import print_function, division 
# Python stdlib
# Chimera stuff
# Additional 3rd parties
# Own


class Controller(object):

    def __init__(self, gui=None, model=None, *args, **kwargs):
        self.gui = gui
        self.model = model

    def set_mvc(self):
        pass

class Model(object):

    def __init__(self, *args, **kwargs):
        return