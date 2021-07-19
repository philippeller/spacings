# spacings

This python package implements several test statistics based on spacings, i.e. the spaces between ordered samples.
The API is similar to scipy's tests, such as the [kstest](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html).

## Installation

```
pip install spacings
```

## Usage Example

```python
import numpy as np
from spacings import rps, moran
from scipy.stats import kstest, cramervonmises
```
Generate random samples
```python
rand = np.random.RandomState(0)
uniform = rand.rand(100)
```

* RPS test:
```python
>>> rps(uniform, "uniform")
RPStestResult(statistic=0.9861622857523498, pvalue=0.6747286655166371)

```

* Moran test:
```python
>>> moran(np.random.rand(100), "uniform")
MorantestResult(statistic=525.7712608675467, pvalue=0.3920410695917047)
```

* KS test (scipy):
```python
>>> kstest(uniform, "uniform")
KstestResult(statistic=0.0736727958814345, pvalue=0.6226614157861845)
```

* Cramer von Mises test (scipy):
```python
>>> cramervonmises(uniform, "uniform")
CramerVonMisesResult(statistic=0.1203440927515015, pvalue=0.4947678804693505)
```
