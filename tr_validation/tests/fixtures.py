import pytest

@pytest.fixture
def cigar_tuple():
    return ([(0, 40), (1, 3), (0, 55), (2, 2)])