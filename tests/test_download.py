"""Size limit of download_file, with a fake server response (no network). Run: python -m unittest discover tests"""
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from belgium_lidar_metrics.download import SizeLimitReached, download_file


def fake_response(chunks):
    response = mock.MagicMock()
    response.__enter__.return_value = response
    response.iter_content.return_value = chunks
    return response


class TestDownloadLimit(unittest.TestCase):

    def test_limit(self):
        chunks = [b"x" * 10] * 3   # a 30-byte "file"
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "a.zip"
            with mock.patch("requests.get", return_value=fake_response(chunks)):
                self.assertEqual(download_file("url", out, max_bytes=30), 30 / 1e6)   # exactly the limit: kept
            self.assertTrue(out.exists())

            out = Path(tmp) / "b.zip"
            with mock.patch("requests.get", return_value=fake_response(chunks)):
                with self.assertRaises(SizeLimitReached):
                    download_file("url", out, max_bytes=29)
            self.assertFalse(out.exists())                                           # never looks finished
            self.assertEqual((Path(tmp) / "partial_b.zip").stat().st_size, 20)       # stopped before passing


if __name__ == "__main__":
    unittest.main()
