from subprocess import call

call(["py.test", "test_common_eqs.py"])
call(["py.test", "--libname=libcloudphxx_blk_1m_pytest", "blk_1mom_test.py"])

