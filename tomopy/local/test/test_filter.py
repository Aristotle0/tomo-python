from tomopy.local.filterseism import seism_filter

def test_filterseism():
    seism_filter(200, 'data')
    