# -*- coding: utf-8 -*-
"""
scratch work
"""





##############################################################################
"""
Inheritance / overloading
"""
class Parent:
    def overridable(self):
        print('Parent')
        

class Child(Parent):
    def overridable(self):
        super().overridable()
        print('Child')
        

class GrandChild(Child):
    def overridable(self):
        super().overridable()
        print('Grand Child')
##############################################################################