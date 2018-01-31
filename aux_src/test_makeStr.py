"""Tests of the makeStr.py module."""
import unittest as ut
import os
import pytest

class TestMakeStructures(ut.TestCase):
    """Tests of the _make_structures subroutine."""

    def _compare_files(self,file1,file2):
        out1 = []
        out2 = []
        with open(file1,"r") as o1:
            for line in o1:
                out1.append(line.strip().split())
        with open(file2,"r") as o2:
            for line in o2:
                out2.append(line.strip().split())
        self.assertEqual(out1,out2)
                
    def test_str1(self):
        from makeStr import _make_structures
        from os import system
        args = {"structures":[1],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"test_files/struct_enum.out_quaternary",
                "mink":True,
                "species":None,
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0,
                "remove_zeros":False,
                "config":"f",
                "mapping":None
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),
                            "test_files/quaternary_vasp.1")
        system("rm vasp*")

    def test_str2(self):
        from makeStr import _make_structures
        from os import system
        args = {"structures":[11],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"test_files/struct_enum.out_quaternary",
                "mink":True,
                "species":["Al","Ni","Ti","Cu"],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0,
                "remove_zeros":False,
                "config":"f",
                "mapping":None
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),
                            "test_files/quaternary_vasp.2")
        system("rm vasp*")        

    def test_str3(self):
        from makeStr import _make_structures
        from os import system
        args = {"structures":[11],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"test_files/struct_enum.out_quaternary",
                "mink":True,
                "species":["Al","Ni","Ti","Cu"],
                "verbose":None,
                "outfile":"vasp.{}",
                "rattle":0.0,
                "remove_zeros":True,
                "config":"f",
                "mapping":None
                }
        _make_structures(args)
        self._compare_files("vasp.{}".format(args["structures"][0]),
                            "test_files/quaternary_vasp.3")
        system("rm vasp*")        

    def test_str4(self):
        from makeStr import _make_structures
        from os import system
        args = {"structures":[1,2,3,4,5,6,7,8,9,10,11],
                "debug":False,
                "examples":False,
                "displace":0.0,
                "input":"test_files/struct_enum.out_quaternary",
                "mink":True,
                "species":["Al","Ni","Ti","Cu"],
                "verbose":None,
                "outfile":"train.cfg",
                "rattle":0.0,
                "remove_zeros":False,
                "config":"t",
                "mapping":None
                }
        _make_structures(args)
        self._compare_files("train.cfg","test_files/quaternary.cfg.1")
        system("rm train.cfg")        
