#!/usr/bin/env python
# encoding: utf-8


from __future__ import print_function, division
import chimera.extension


class CauchianExtension(chimera.extension.EMO):

    def name(self):
        return 'Tangram Cauchian'

    def description(self):
        return "Prepare input files for Gaussian jobs"

    def categories(self):
        return ['InsiliChem']

    def icon(self):
        return

    def activate(self):
        self.module('gui').showUI()


chimera.extension.manager.registerExtension(CauchianExtension(__file__))
