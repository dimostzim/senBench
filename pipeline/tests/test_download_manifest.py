import pytest

import download_manifest


def fake_attempts(monkeypatch, attempts):
    """Let each call of curl give the next (exit code, HTTP status, bytes written); return the waits."""
    waits = []
    monkeypatch.setattr(download_manifest.time, "sleep", waits.append)

    def fetch(url, target, resume):
        code, status, data = attempts.pop(0)
        with target.open("ab") as handle:
            handle.write(data)
        return code, status

    monkeypatch.setattr(download_manifest, "_fetch", fetch)
    return waits


def test_a_refusal_for_too_many_requests_waits_the_rate_limit_delay(monkeypatch, tmp_path):
    waits = fake_attempts(monkeypatch, [(22, "429", b""), (18, "200", b"ab"), (0, "200", b"cd")])
    download_manifest._fetch_until_complete("https://x", tmp_path / "f", 4, retries=3, retry_delay=10,
                                            rate_limit_delay=300)
    assert waits == [300, 10]
    assert (tmp_path / "f").read_bytes() == b"abcd"


def test_attempts_without_new_bytes_give_up_after_the_retries(monkeypatch, tmp_path):
    fake_attempts(monkeypatch, [(22, "429", b"")] * 3)
    with pytest.raises(download_manifest._FileProblem, match="HTTP 429"):
        download_manifest._fetch_until_complete("https://x", tmp_path / "f", 4, retries=1, retry_delay=10,
                                                rate_limit_delay=300)
