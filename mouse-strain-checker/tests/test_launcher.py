"""Windows launcher checks use stubs and never run Docker or open a browser."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
POWERSHELL = shutil.which("powershell.exe") or shutil.which("pwsh")


class LauncherFileTests(unittest.TestCase):
    def test_batch_preserves_script_path_and_failure_pause(self):
        batch = (ROOT / "start.bat").read_text(encoding="ascii")
        self.assertIn('-File "%~dp0start.ps1"', batch)
        self.assertRegex(batch, r"(?s)if errorlevel 1\s*\(.*pause.*exit /b 1")

    def test_powershell_source_is_ascii_for_windows_powershell(self):
        (ROOT / "start.ps1").read_text(encoding="ascii")


@unittest.skipUnless(POWERSHELL, "PowerShell is required for launcher execution checks")
class LauncherPowerShellTests(unittest.TestCase):
    def run_powershell(self, body):
        environment = os.environ.copy()
        environment["MOUSE_STRAIN_LAUNCHER"] = str(ROOT / "start.ps1")
        result = subprocess.run(
            [POWERSHELL, "-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass", "-Command", body],
            capture_output=True, text=True, timeout=15, env=environment,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout

    def test_windows_powershell_parser_accepts_launcher(self):
        self.run_powershell(r"""
            $tokens = $null
            $errors = $null
            [void][System.Management.Automation.Language.Parser]::ParseFile(
                $env:MOUSE_STRAIN_LAUNCHER, [ref]$tokens, [ref]$errors)
            if ($errors.Count -gt 0) { throw ($errors | Out-String) }
        """)

    def run_stubbed_launcher(self, running=True, app="mouse-strain-checker", compose_failure=False):
        body = r"""
            . $env:MOUSE_STRAIN_LAUNCHER
            $script:commands = @()
            $script:opened = @()
            function Test-DockerReady { return $true }
            function docker {
                $script:commands += ($args -join ' ')
                $global:LASTEXITCODE = 0
                if (($args -join ' ') -eq 'compose ps --status running --quiet strain-checker') {
                    if ($env:LAUNCHER_TEST_RUNNING -eq 'true') { 'existing-container' }
                }
                if (($args -join ' ') -eq 'compose up -d --build strain-checker' -and $env:LAUNCHER_TEST_FAIL -eq 'true') {
                    $global:LASTEXITCODE = 1
                }
            }
            function Invoke-WebRequest {
                param($Uri, [switch]$UseBasicParsing, $TimeoutSec)
                if ($Uri -cne 'http://127.0.0.1:8502/api/health') { throw 'Incorrect health URL' }
                return [pscustomobject]@{ StatusCode = 200; Content = ('{"app":"' + $env:LAUNCHER_TEST_APP + '"}') }
            }
            function Start-Process { param($FilePath) $script:opened += $FilePath }
            $failure = $null
            try { Invoke-MouseStrainLauncher | Out-Null } catch { $failure = $_.Exception.Message }
            [pscustomobject]@{ commands = @($script:commands); opened = @($script:opened); failure = $failure } |
                ConvertTo-Json -Compress
        """
        # Test settings use process-local environment variables, never a native command.
        with patch.dict(os.environ, {
            "LAUNCHER_TEST_RUNNING": str(running).lower(),
            "LAUNCHER_TEST_APP": app,
            "LAUNCHER_TEST_FAIL": str(compose_failure).lower(),
        }):
            output = self.run_powershell(body)
        return json.loads(output.strip().splitlines()[-1])

    def test_existing_container_is_reused_without_build_or_recreate(self):
        result = self.run_stubbed_launcher()
        self.assertIsNone(result["failure"])
        self.assertFalse(any(" up " in command for command in result["commands"]))
        self.assertEqual(result["opened"], ["http://localhost:8502"])

    def test_stopped_service_is_built_and_started(self):
        result = self.run_stubbed_launcher(running=False)
        self.assertIsNone(result["failure"])
        self.assertIn("compose up -d --build strain-checker", result["commands"])
        self.assertEqual(result["opened"], ["http://localhost:8502"])

    def test_unrelated_health_server_does_not_open_browser(self):
        result = self.run_stubbed_launcher(app="some-other-app")
        self.assertIn("unexpected app", result["failure"])
        self.assertEqual(result["opened"], [])
        self.assertFalse(any(" up " in command for command in result["commands"]))

    def test_compose_failure_does_not_open_browser(self):
        result = self.run_stubbed_launcher(running=False, compose_failure=True)
        self.assertIn("Docker Compose failed", result["failure"])
        self.assertEqual(result["opened"], [])


if __name__ == "__main__":
    unittest.main()
