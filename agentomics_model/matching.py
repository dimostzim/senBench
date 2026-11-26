"""Training-only reciprocal opposite-class matching in partially whitened context."""
from dataclasses import dataclass

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.utils.extmath import randomized_svd


@dataclass(frozen=True)
class PairGraph:
    positive: np.ndarray
    control: np.ndarray
    weight: np.ndarray
    distance: np.ndarray

    def __len__(self):
        return len(self.weight)


def reciprocal_graph(coords, y, k=10):
    """Exact Euclidean neighbors; ties follow the original training row order."""
    pos, ctrl = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    if not len(pos) or not len(ctrl) or coords.shape[1] == 0:
        empty = np.empty(0, dtype=np.intp)
        return PairGraph(empty, empty, np.empty(0), np.empty(0))

    def neighbors(src, dst):
        take = min(k, len(dst))
        ix = np.empty((len(src), take), dtype=np.intp)
        dist = np.empty((len(src), take), dtype=np.float64)
        dst_coords = coords[dst]
        for start in range(0, len(src), 256):
            sl = slice(start, start + 256)
            distances = cdist(coords[src[sl]], dst_coords, metric='euclidean')
            order = np.argsort(distances, axis=1, kind='stable')[:, :take]
            ix[sl] = order
            dist[sl] = np.take_along_axis(distances, order, axis=1)
        return ix, dist

    forward, distance = neighbors(pos, ctrl)
    reverse, _ = neighbors(ctrl, pos)
    mutual = (reverse[forward] == np.arange(len(pos))[:, None, None]).any(axis=2)
    i, slot = np.nonzero(mutual)
    if not len(i):
        return PairGraph(pos[i], ctrl[i], np.empty(0), np.empty(0))
    j = forward[i, slot]
    dp, dc = np.bincount(i, minlength=len(pos)), np.bincount(j, minlength=len(ctrl))
    w = 0.5 / (np.count_nonzero(dp) * dp[i]) + 0.5 / (np.count_nonzero(dc) * dc[j])
    graph = PairGraph(pos[i], ctrl[j], w, distance[i, slot])
    assert np.all(y[graph.positive] == 1) and np.all(y[graph.control] == 0)
    assert np.all(w > 0)
    np.testing.assert_allclose(w.sum(), 1., rtol=0., atol=1e-14)
    return graph


def build_graph(b, y):
    rank = min(24, len(y) - 1, b.shape[1])
    coords = np.empty((len(y), 0), dtype=np.float64)
    if rank > 0:
        b = np.asarray(b, dtype=np.float64)
        centered = b - b.mean(axis=0)
        # Eliminate arithmetic residuals from exactly constant columns.
        centered[:, np.all(b == b[:1], axis=0)] = 0.
        if np.any(centered):
            _, s, vt = randomized_svd(centered, n_components=rank, n_iter=7, random_state=2026)
            projected = centered @ vt[s > 1e-8 * s[0]].T
            coords = projected / np.sqrt(projected.std(axis=0, ddof=0))
    return reciprocal_graph(coords, y), coords.shape[1]


def audit_graph(graph, rank, y, cells):
    active = np.zeros(len(y), dtype=bool)
    active[graph.positive] = active[graph.control] = True
    def counts(mask):
        return {'active': int((active & mask).sum()), 'total': int(mask.sum())}
    return {'rank': rank, 'edges': len(graph), 'weight_sum': float(graph.weight.sum()),
            'distance_quantiles': dict(zip(('min', 'p25', 'p50', 'p75', 'max'),
                                          np.quantile(graph.distance, [0, .25, .5, .75, 1]).tolist()))
            if len(graph) else {},
            'endpoint_coverage': {c: {'overall': counts(y == c),
                                      'by_cell': {g: counts((y == c) & (cells == g)) for g in range(-1, 8)}}
                                  for c in (0, 1)}}


def test_graph():
    # The eleventh tied positive is rejected by each control's reverse top ten.
    y = np.array([0, 1, 0] + [1] * 10)
    graph = reciprocal_graph(np.zeros((len(y), 1)), y)
    p = np.flatnonzero(y)[:10]
    np.testing.assert_array_equal(graph.positive, np.repeat(p, 2))
    np.testing.assert_array_equal(graph.control, np.tile([0, 2], 10))
    np.testing.assert_allclose(graph.weight, np.full(20, 1 / 20))
    # Unequal degrees: control 3 prefers positives 0 and 2; control 4
    # prefers positives 2 and 5. Positive 0 also takes control 1 by row tie.
    coords = np.array([0., -2., 4., 2., 6., 8.])[:, None]
    y = np.array([1, 0, 1, 0, 0, 1])
    graph = reciprocal_graph(coords, y, k=2)
    expected = {(0, 1): 1/6, (0, 3): 1/6, (2, 3): 1/6,
                (2, 4): 1/6, (5, 4): 1/4}
    # Control 1 chooses 0 and 2, but 2 does not reciprocate; its
    # single endpoint share is 1/6, in addition to 0's 1/12.
    expected[0, 1] = 1/4
    assert set(zip(graph.positive, graph.control)) == set(expected)
    np.testing.assert_allclose(graph.weight, [expected[p, c] for p, c in zip(graph.positive, graph.control)])
    g, rank = build_graph(np.full((len(y), 2), 0.1), y)
    assert rank == len(g) == 0
    g, rank = build_graph(coords, y)
    assert rank == 1 and len(g) == 9
    np.testing.assert_allclose(g.distance, np.abs(coords[g.positive, 0] - coords[g.control, 0])
                               / np.sqrt(coords[:, 0].std()))
