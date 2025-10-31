import pytest
import LebwohlLasher_numpy as LL
import numpy as np

def test_initdat():
    '''
    Checks that random values in array range from 0 to 2*pi
    '''
    import numpy as np
    nmax=10
    array = LL.initdat(nmax)
    assert (array>=0).all
    assert (array<=2*np.pi).all
    assert array.size==nmax*nmax

def test_order():
    '''
    Checks that the order calculation is correct for a fully ordered system
    '''
    import numpy as np
    nmax=5
    array = np.zeros((nmax, nmax))
    assert LL.get_order(array,nmax)==pytest.approx(1.0, abs=1e-6)
    
@pytest.mark.parametrize("array, nmax, expected", [
    (np.zeros((3,3)),   3,    -36)
])

def test_totenergy(array, nmax, expected):
    '''
    Checks that the energy sum is functionng correctly
    '''
    energy = LL.all_energy(array,nmax)
    assert energy == pytest.approx(expected, abs=1e-6)

@pytest.mark.parametrize("ix, iy, nmax, expected", [
    (0,  2,  3,    -2.5),
    (1,  1,  3,    -2.5),
    (2,  1,  3,    -2.5),
    (1,  2,  3,   -4.0)
])

def test_oneenergy(ix, iy, nmax, expected):
    '''
    Checks that individual energy sums are functioning as intended
    taking into account periodic boundaries
    '''
    array = np.zeros((nmax,nmax))
    array[0,1] = np.pi/2
    energy = LL.one_energy(array, ix, iy, nmax)
    assert energy == pytest.approx(expected, abs=1e-6)
  