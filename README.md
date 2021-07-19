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
```

```python
rps(np.random.rand(100), "uniform")
```
```python
moran(np.random.rand(100), "uniform")
```
