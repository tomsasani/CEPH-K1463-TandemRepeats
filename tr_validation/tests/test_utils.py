import pytest

from utils import count_indel_in_read, bp_overlap


@pytest.mark.parametrize(
    "ct,rs,vs,ve,slop,exp",
    [
        pytest.param(
            [(0, 10), (1, 2), (0, 5), (2, 3)],
            0,
            9,
            20,
            0,
            -1,
        ),
        pytest.param(
            [(0, 10), (1, 2), (0, 5), (2, 3)],
            0,
            2,
            14,
            0,
            2,
        ),
    ],
)
def test_count_indel_in_read(ct, rs, vs, ve, slop, exp):
    assert count_indel_in_read(ct, rs, vs, ve, slop=slop) == exp


@pytest.mark.parametrize(
    "s1,e1,s2,e2,exp",
    [
        pytest.param(
            0, 10, 5, 11, 5
        ),
        pytest.param(
            6, 7, 8, 10, 0
        ),
        pytest.param(
            1, 11, 5, 6, 1
        )
    ],
)
def test_bp_overlap(s1, e1, s2, e2, exp):
    assert bp_overlap(s1, e1, s2, e2) == exp

