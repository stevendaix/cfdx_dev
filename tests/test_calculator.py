#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from calculator import add, multiply


def test_add():
    assert add(2, 3) == 5
    assert add(-1, 1) == 0
    assert add(0, 0) == 0
    assert add(10, -5) == 5
    assert add(1.5, 2.5) == 4.0


def test_multiply():
    assert multiply(2, 3) == 6
    assert multiply(-1, 5) == -5
    assert multiply(0, 10) == 0
    assert multiply(3, 4) == 12
    assert multiply(1.5, 2) == 3.0


if __name__ == "__main__":
    test_add()
    test_multiply()
    print("All tests passed!")