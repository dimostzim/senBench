"""Frozen quadratic classifier artifact and scoring interface."""
from dataclasses import asdict, dataclass
from pathlib import Path

import joblib
import numpy as np
from scipy.special import expit

from representation import OrdinalState, QuadraticState, LABEL_MAPPING, SCHEMA, validate_quadratic_state

MODEL_NAME = 'model.joblib'
ARTIFACT_VERSION = 17


@dataclass(frozen=True)
class LinearClassifier:
    """Linear in the joint ordinal and quadratic design."""
    w: np.ndarray
    b: float

    def decision_function(self, z):
        return np.asarray(z, dtype=np.float64) @ self.w + self.b

    def predict_proba(self, z):
        s = self.decision_function(z)
        return np.column_stack((expit(-s), expit(s)))


def save_pipeline(path: str | Path, model: LinearClassifier, state: QuadraticState):
    joblib.dump({'artifact_version': ARTIFACT_VERSION, 'representation_schema': SCHEMA,
                 'label_mapping': LABEL_MAPPING, 'quadratic_state': asdict(state),
                 'w': model.w, 'b': model.b}, path, compress=3)


def load_pipeline(path: str | Path) -> tuple[LinearClassifier, QuadraticState]:
    # Load only trusted artifacts.
    payload = joblib.load(path)
    if (payload.get('artifact_version') != ARTIFACT_VERSION
            or payload.get('representation_schema') != SCHEMA
            or payload.get('label_mapping') != LABEL_MAPPING):
        raise ValueError('Unsupported model artifact or label mapping')
    fields = dict(payload['quadratic_state'])
    fields['ordinal'] = OrdinalState(**fields['ordinal'])
    state = QuadraticState(**fields)
    validate_quadratic_state(state)
    w = np.asarray(payload['w'], dtype=np.float64)
    b = float(payload['b'])
    if w.shape != (state.n_features,) or not np.isfinite(w).all() or not np.isfinite(b):
        raise ValueError('Invalid binary classifier parameters')
    return LinearClassifier(w, b), state
