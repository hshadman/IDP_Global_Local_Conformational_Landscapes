import pandas as pd

def test_custom_gw_reference_path(tmp_path):
    # create a tiny CSV
    custom_csv = tmp_path / "mini.csv"
    pd.DataFrame({"dummy": [1]}).to_csv(custom_csv, index=False)

    from pyheteromap import PyHeteroMap
    h = PyHeteroMap("TEST", gw_reference_csv=str(custom_csv))
    h.reinitialize_gw_reference()

    assert h.gw_df is not None
    assert not h.gw_df.empty
