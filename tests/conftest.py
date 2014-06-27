import os
import os.path
import sys

sys.path.append(os.path.join(os.getcwd(), '.'))
sys.path.append(os.path.join(os.getcwd(), '..'))
sys.path.append(os.path.join(os.getcwd(), '../wrf_microphys/kessler/'))


def pytest_addoption(parser):
    parser.addoption("--libname", action="append", default=[],
        help="name of the tested library")

def pytest_generate_tests(metafunc):
    if 'libname' in metafunc.fixturenames:
        metafunc.parametrize("libname", metafunc.config.option.libname)
