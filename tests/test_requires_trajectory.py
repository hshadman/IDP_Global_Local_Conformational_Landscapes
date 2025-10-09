import pytest

def test_methods_raise_before_set_trajectory():
    from pyheteromap import PyHeteroMap
    h = PyHeteroMap("TEST")

    with pytest.raises(RuntimeError):
        h.calculate_nu_KLL_from_seq_name(1, 5)

    with pytest.raises(RuntimeError):
        h.plot_subchain_RSA(4, 3)

    with pytest.raises(RuntimeError):
        h.mod_RSA_Rs_compute_3dplot_from_seq_name()
