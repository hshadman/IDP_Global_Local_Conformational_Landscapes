def test_import_and_load():
    from pyheteromap import PyHeteroMap
    h = PyHeteroMap("TEST")
    assert h.gw_reference_csv.endswith("reference_GW_chainlength_100.csv")
    assert h.gw_df is None
