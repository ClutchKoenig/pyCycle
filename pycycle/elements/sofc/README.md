# pyCycle SOFC

Solid-oxide fuel-cell (SOFC) element for [pyCycle](https://github.com/ClutchKoenig/pyCycle).

This repository is mounted as a git submodule at `pycycle/elements/sofc` inside pyCycle
and depends on pyCycle's thermo and element infrastructure.

```python
from pycycle.elements.sofc import SOFC, SegmentSOFC
```

- `api.py` – public API of all SOFC sub-components
- `sofc.py` – full `SOFC` element; `sofc_segment.py` – single discretisation segment
- `tests/` – test scripts and notebooks
- `docs/` – model notes and change logs
