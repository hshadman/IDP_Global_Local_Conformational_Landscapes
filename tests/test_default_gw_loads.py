import pandas as pd

def test_default_gw_reference_loads_nonempty_df():
    from pyheteromap import PyHeteroMap
    h = PyHeteroMap("TEST")
    h.reinitialize_gw_reference()
    assert h.gw_df is not None
    assert isinstance(h.gw_df, pd.DataFrame)
    assert not h.gw_df.empty
