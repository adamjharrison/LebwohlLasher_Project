import pytest
import LebwohlLasher as LL
import numpy as np
def test_initdat():
    import numpy as np
    nmax=10
    array = LL.initdat(nmax)
    assert (array>=0).all
    assert (array<=2*np.pi).all

def test_order():
    import numpy as np
    nmax=5
    array = np.zeros((nmax, nmax))
    assert LL.get_order(array,nmax)==pytest.approx(1.0, abs=1e-6)
    
@pytest.mark.parametrize("array, nmax, energy", [
    (np.zeros((3,3)),   3,    -36)
])

def test_oneenergy(array, nmax, energy):
    assert LL.all_energy(array,nmax) == pytest.approx(energy, abs=1e-6)

    