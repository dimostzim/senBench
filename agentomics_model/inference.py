"""Score a label-free expression table with the fitted ordinal–quadratic classifier."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import tempfile

import numpy as np
from scipy.special import expit
from threadpoolctl import threadpool_limits

from pipeline import MODEL_NAME, load_pipeline
from representation import DEFAULT_BATCH_ROWS, iter_quadratic_features

DEFAULT_ARTIFACTS = Path(__file__).resolve().parent / 'artifacts'
COLUMNS = ('id', 'prediction', 'probability_0', 'probability_1')


def _scorer(model, use_cuda):
    if use_cuda:
        try:
            import torch
        except ImportError:
            pass
        else:
            if torch.cuda.is_available():
                try:
                    with torch.inference_mode():
                        weights = torch.as_tensor(model.w, dtype=torch.float64, device='cuda')
                except torch.cuda.OutOfMemoryError:
                    return model.decision_function

                def cuda_score(x):
                    nonlocal weights
                    if weights is not None:
                        try:
                            with torch.inference_mode():
                                rows = torch.as_tensor(x, dtype=torch.float64, device=weights.device)
                                return (rows @ weights + model.b).cpu().numpy()
                        except torch.cuda.OutOfMemoryError:
                            weights = None
                    return model.decision_function(x)

                return cuda_score
    return model.decision_function


def predict(input_dir, output, artifacts_dir=DEFAULT_ARTIFACTS, *, batch_rows=DEFAULT_BATCH_ROWS,
            use_cuda=True):
    """Write every input row in order, committing the CSV only after a successful scan."""
    model, state = load_pipeline(Path(artifacts_dir) / MODEL_NAME)
    score = _scorer(model, use_cuda)
    output = Path(output)
    if output.resolve() == (Path(input_dir) / 'data.csv').resolve():
        raise ValueError('Output must not overwrite the input table')
    output.parent.mkdir(parents=True, exist_ok=True)
    temp_path = None
    n = 0
    try:
        with (threadpool_limits(limits=8, user_api='blas'),
              tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='', dir=output.parent,
                                          prefix='.predictions-', suffix='.tmp', delete=False) as f):
            temp_path = Path(f.name)
            writer = csv.writer(f)
            writer.writerow(COLUMNS)
            for batch, z, _ in iter_quadratic_features(input_dir, state, batch_rows):
                logits = score(z)
                decision = (logits >= 0).astype(np.int64)
                writer.writerows(zip(batch.ids, decision, expit(-logits), expit(logits)))
                n += len(batch.ids)
        temp_path.replace(output)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--artifacts-dir', type=Path, default=DEFAULT_ARTIFACTS)
    args = parser.parse_args()
    predict(args.input, args.output, args.artifacts_dir)


if __name__ == '__main__':
    main()
