from oop_adv_hor import Micro

Test_micro = Micro(nx=5, dx=1, time_adv_tot=1, dt=1, C=1, aerosol={}, RHenv=.5, sl_sg=slice(1,3), apr="trad", setup="rhoconst",  n_intrp=1, sl_act_time=1, dir_name="none", test=True, it_output_l=[], scheme="sd")

def test_supersat():
    rv0 = Test_micro.state["rv"].copy()
    Test_micro.rv2absS()
    Test_micro.absS2rv()
    assert (rv0 == Test_micro.state["rv"]).all()
