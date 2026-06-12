#!/usr/bin/env python3
import json
import os
import re
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

try:
    import paramiko
except Exception:
    paramiko = None


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PROFILE = SCRIPT_DIR / "cluster_profile.local.json"
DEFAULT_EXAMPLE_PROFILE = SCRIPT_DIR / "cluster_profile.example.json"


@dataclass
class CommandResult:
    returncode: int
    stdout: str
    stderr: str
    command: list[str]


def load_profile(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def save_profile(path: Path, profile: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(profile, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def default_profile() -> dict:
    return load_profile(DEFAULT_EXAMPLE_PROFILE)


def _ssh_target(profile: dict) -> str:
    user = profile.get("user", "").strip()
    host = profile.get("host", "").strip()
    if not user or not host:
        raise ValueError("Cluster profile requires both user and host")
    return f"{user}@{host}"


def _base_ssh_args(profile: dict) -> list[str]:
    args = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10"]
    port = profile.get("port")
    if port:
        args.extend(["-p", str(port)])
    ssh_key = str(profile.get("ssh_key", "")).strip()
    if ssh_key:
        args.extend(["-i", os.path.expanduser(ssh_key)])
    args.extend(str(x) for x in profile.get("ssh_extra_args", []))
    args.append(_ssh_target(profile))
    return args


def _base_scp_args(profile: dict) -> list[str]:
    args = ["scp", "-r", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10"]
    port = profile.get("port")
    if port:
        args.extend(["-P", str(port)])
    ssh_key = str(profile.get("ssh_key", "")).strip()
    if ssh_key:
        args.extend(["-i", os.path.expanduser(ssh_key)])
    args.extend(str(x) for x in profile.get("ssh_extra_args", []))
    return args


def _run(command: list[str], timeout: int = 60) -> CommandResult:
    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    return CommandResult(completed.returncode, completed.stdout, completed.stderr, command)


def _uses_password_auth(profile: dict) -> bool:
    return str(profile.get("auth_method", "")).lower() == "password" or bool(profile.get("password"))


def _require_paramiko():
    if paramiko is None:
        raise RuntimeError(
            "Password login requires the Python package 'paramiko'. "
            "Install it in the Python environment used by the GUI, then restart the GUI."
        )


def _connect_paramiko(profile: dict):
    _require_paramiko()
    if _uses_password_auth(profile) and not profile.get("password"):
        raise ValueError("Password authentication selected, but no password was provided in the GUI.")
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    connect_kwargs = {
        "hostname": profile.get("host", "").strip(),
        "port": int(profile.get("port") or 22),
        "username": profile.get("user", "").strip(),
        "timeout": 15,
        "look_for_keys": False,
        "allow_agent": False,
    }
    password = profile.get("password")
    if password:
        connect_kwargs["password"] = password
    client.connect(**connect_kwargs)
    return client


def _run_paramiko(profile: dict, command: str, timeout: int = 60) -> CommandResult:
    client = _connect_paramiko(profile)
    wrapped = f"bash -lc {shlex.quote(command)}"
    command_label = ["paramiko", _ssh_target(profile), wrapped]
    try:
        stdin, stdout, stderr = client.exec_command(wrapped, timeout=timeout)
        channel = stdout.channel
        stdout_chunks = []
        stderr_chunks = []
        deadline = time.monotonic() + timeout
        while True:
            while channel.recv_ready():
                stdout_chunks.append(channel.recv(65536))
            while channel.recv_stderr_ready():
                stderr_chunks.append(channel.recv_stderr(65536))
            if channel.exit_status_ready():
                break
            if time.monotonic() > deadline:
                channel.close()
                stdout_text = b"".join(stdout_chunks).decode("utf-8", errors="replace")
                stderr_text = b"".join(stderr_chunks).decode("utf-8", errors="replace")
                stderr_text += f"\nTimed out after {timeout} seconds.\n"
                return CommandResult(124, stdout_text, stderr_text, command_label)
            time.sleep(0.1)
        while channel.recv_ready():
            stdout_chunks.append(channel.recv(65536))
        while channel.recv_stderr_ready():
            stderr_chunks.append(channel.recv_stderr(65536))
        exit_code = channel.recv_exit_status()
        stdout_text = b"".join(stdout_chunks).decode("utf-8", errors="replace")
        stderr_text = b"".join(stderr_chunks).decode("utf-8", errors="replace")
        return CommandResult(exit_code, stdout_text, stderr_text, command_label)
    finally:
        client.close()


def _remote_path_expr(path: str) -> str:
    text = str(path).strip()
    if not text:
        raise ValueError("Remote path is empty")
    if text == "$HOME":
        return '"$HOME"'
    if text.startswith("$HOME/"):
        return '"$HOME"/' + shlex.quote(text[len("$HOME/"):])
    return shlex.quote(text)


def run_remote(profile: dict, command: str, timeout: int = 60) -> CommandResult:
    if _uses_password_auth(profile):
        return _run_paramiko(profile, command, timeout=timeout)
    ssh_command = _base_ssh_args(profile)
    ssh_command.extend(["bash", "-lc", command])
    return _run(ssh_command, timeout=timeout)


def test_connection(profile: dict) -> CommandResult:
    return run_remote(profile, "echo connected:$(hostname); command -v bjobs || true; command -v bsub || true; command -v bkill || true")


def ensure_remote_dir(profile: dict, remote_dir: str) -> CommandResult:
    return run_remote(profile, f"mkdir -p {_remote_path_expr(remote_dir)}")


def expand_remote_path(profile: dict, remote_path: str) -> str:
    result = run_remote(profile, f"printf '%s' {_remote_path_expr(remote_path)}", timeout=30)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or f"Failed to expand remote path: {remote_path}")
    return result.stdout.strip()


def upload_workflow(profile: dict, local_workflow_dir: Path) -> CommandResult:
    local_dir = Path(local_workflow_dir).resolve()
    if not local_dir.is_dir():
        raise FileNotFoundError(f"Local workflow directory not found: {local_dir}")

    parent = str(profile.get("remote_upload_parent") or "").strip()
    if not parent:
        remote_workflow_dir = str(profile.get("remote_workflow_dir") or "").rstrip("/")
        if not remote_workflow_dir:
            raise ValueError("Profile requires remote_upload_parent or remote_workflow_dir")
        parent = remote_workflow_dir.rsplit("/", 1)[0]

    ensure = ensure_remote_dir(profile, parent)
    if ensure.returncode != 0:
        return ensure

    if _uses_password_auth(profile):
        return _upload_workflow_paramiko(profile, local_dir, expand_remote_path(profile, parent))

    destination = f"{_ssh_target(profile)}:{parent.rstrip('/')}/"
    command = _base_scp_args(profile) + [str(local_dir), destination]
    return _run(command, timeout=600)


def _upload_workflow_paramiko(profile: dict, local_dir: Path, remote_parent: str) -> CommandResult:
    client = _connect_paramiko(profile)
    command_label = ["paramiko-sftp-upload", str(local_dir), remote_parent]
    try:
        sftp = client.open_sftp()
        remote_root = f"{remote_parent.rstrip('/')}/{local_dir.name}"
        _sftp_mkdirs(sftp, remote_root)
        uploaded = []
        for path in local_dir.rglob("*"):
            relative = path.relative_to(local_dir).as_posix()
            remote_path = f"{remote_root}/{relative}"
            if path.is_dir():
                _sftp_mkdirs(sftp, remote_path)
            else:
                _sftp_mkdirs(sftp, str(Path(remote_path).parent).replace("\\", "/"))
                sftp.put(str(path), remote_path)
                uploaded.append(remote_path)
        return CommandResult(0, f"Uploaded {len(uploaded)} files to {remote_root}\n", "", command_label)
    except Exception as exc:
        return CommandResult(1, "", f"{type(exc).__name__}: {exc}\n", command_label)
    finally:
        try:
            sftp.close()
        except Exception:
            pass
        client.close()


def _sftp_mkdirs(sftp, remote_dir: str) -> None:
    parts = [part for part in remote_dir.replace("\\", "/").split("/") if part]
    current = "/" if remote_dir.startswith("/") else ""
    for part in parts:
        current = f"{current.rstrip('/')}/{part}" if current else part
        try:
            sftp.stat(current)
        except IOError:
            sftp.mkdir(current)


def submit_stage(profile: dict, stage_name: str) -> CommandResult:
    scripts = profile.get("stage_scripts", {})
    script = scripts.get(stage_name)
    if not script:
        raise ValueError(f"No script configured for stage: {stage_name}")
    remote_dir = profile.get("remote_workflow_dir")
    if not remote_dir:
        raise ValueError("Profile requires remote_workflow_dir")
    command = f"cd {_remote_path_expr(remote_dir)} && bash {shlex.quote(script)}"
    return run_remote(profile, command, timeout=120)


def list_jobs(profile: dict) -> tuple[CommandResult, list[dict]]:
    user = str(profile.get("bjobs_user") or "$USER").strip()
    if user == "$USER":
        command = 'bjobs -u "$USER" -w || true'
    elif user:
        command = f"bjobs -u {shlex.quote(user)} -w || true"
    else:
        command = "bjobs -w || true"
    result = run_remote(profile, command, timeout=30)
    return result, parse_bjobs(result.stdout)


def describe_jobs(profile: dict, job_ids: list[str]) -> CommandResult:
    clean_ids = [
        str(job_id).strip()
        for job_id in job_ids
        if re.fullmatch(r"\d+(?:\[[^\]\s]+\])?", str(job_id).strip())
    ]
    if not clean_ids:
        raise ValueError("No valid job ids selected")
    command = "bjobs -l " + " ".join(shlex.quote(job_id) for job_id in clean_ids) + " || true"
    return run_remote(profile, command, timeout=60)


def _lsf_user_option(profile: dict) -> str:
    user = str(profile.get("bjobs_user") or "$USER").strip()
    if user == "$USER":
        return '-u "$USER"'
    if user:
        return "-u " + shlex.quote(user)
    return ""


def parse_bjobs(output: str) -> list[dict]:
    jobs = []
    lines = [line.rstrip() for line in output.splitlines() if line.strip()]
    for line in lines:
        if line.upper().startswith("JOBID "):
            continue
        parts = line.split()
        if len(parts) < 6 or not re.fullmatch(r"\d+(?:\[[^\]]+\])?", parts[0]):
            continue
        jobs.append({
            "jobid": parts[0],
            "user": parts[1],
            "stat": parts[2],
            "queue": parts[3],
            "from_host": parts[4],
            "exec_host": parts[5],
            "job_name": parts[6] if len(parts) > 6 else "",
            "submit_time": " ".join(parts[7:]) if len(parts) > 7 else "",
            "raw": line,
        })
    return jobs


def kill_jobs(profile: dict, job_ids: list[str]) -> CommandResult:
    clean_ids = [
        str(job_id).strip()
        for job_id in job_ids
        if re.fullmatch(r"\d+(?:\[[^\]\s]+\])?", str(job_id).strip())
    ]
    if not clean_ids:
        raise ValueError("No valid job ids selected")
    command = "bkill " + " ".join(shlex.quote(job_id) for job_id in clean_ids)
    return run_remote(profile, command, timeout=30)


def _scan_remote_files(profile: dict, roots: list[str], names: list[str], max_depth: int) -> tuple[CommandResult, list[dict]]:
    if not roots:
        raise ValueError("No remote roots configured for scanning")
    quoted_roots = " ".join(_remote_path_expr(root) for root in roots)
    name_expr = " -o ".join(f"-name {shlex.quote(pattern)}" for pattern in names)
    command = (
        f"for root in {quoted_roots}; do "
        f"if [ -e \"$root\" ]; then "
        f"find \"$root\" -maxdepth {max_depth} -type f \\( {name_expr} \\) "
        f"-printf '%p\\t%s\\t%TY-%Tm-%Td %TH:%TM\\n'; "
        f"fi; "
        f"done"
    )
    result = run_remote(profile, command, timeout=120)
    return result, parse_result_scan(result.stdout)


def download_results(profile: dict, local_dir: Path | None = None) -> list[CommandResult]:
    remote_workflow_dir = str(profile.get("remote_workflow_dir") or "").rstrip("/")
    if not remote_workflow_dir:
        raise ValueError("Profile requires remote_workflow_dir")
    result_paths = profile.get("result_paths", [])
    if not result_paths:
        raise ValueError("Profile result_paths is empty")

    destination = Path(local_dir or profile.get("local_download_dir") or "downloads").expanduser()
    if not destination.is_absolute():
        destination = SCRIPT_DIR / destination
    destination.mkdir(parents=True, exist_ok=True)

    results = []
    for relative_path in result_paths:
        remote_path = f"{remote_workflow_dir}/{relative_path}".replace("//", "/")
        if _uses_password_auth(profile):
            results.append(_download_file_paramiko(profile, expand_remote_path(profile, remote_path), destination))
            continue
        source = f"{_ssh_target(profile)}:{remote_path}"
        command = _base_scp_args(profile) + [source, str(destination)]
        results.append(_run(command, timeout=600))
    return results


def scan_results(profile: dict) -> tuple[CommandResult, list[dict]]:
    roots = profile.get("scratch_scan_roots") or []
    if not roots:
        remote_workflow_dir = str(profile.get("remote_workflow_dir") or "").strip()
        roots = [remote_workflow_dir] if remote_workflow_dir else []
    if not roots:
        raise ValueError("Profile requires scratch_scan_roots or remote_workflow_dir")

    max_depth = int(profile.get("scratch_scan_max_depth", 8))
    names = profile.get("result_file_patterns") or [
        "*.csv",
        "*.pdb",
        "*.out",
        "*.err",
        "*.json",
        "*.txt",
    ]
    return _scan_remote_files(profile, roots, names, max_depth)


def scan_merged_csvs(profile: dict) -> tuple[CommandResult, list[dict]]:
    roots = profile.get("scratch_scan_roots") or []
    if not roots:
        remote_workflow_dir = str(profile.get("remote_workflow_dir") or "").strip()
        roots = [remote_workflow_dir] if remote_workflow_dir else []
    if not roots:
        raise ValueError("Profile requires scratch_scan_roots or remote_workflow_dir")

    max_depth = int(profile.get("scratch_scan_max_depth", 8))
    result, entries = _scan_remote_files(profile, roots, ["*.csv"], max_depth)
    merged_entries = [
        entry
        for entry in entries
        if "/merged/" in entry["path"].replace("\\", "/").lower() or "merged" in entry["name"].lower()
    ]
    return result, merged_entries


def scan_job_logs(profile: dict) -> tuple[CommandResult, list[dict]]:
    roots = profile.get("job_log_scan_roots") or profile.get("scratch_scan_roots") or []
    if not roots:
        remote_workflow_dir = str(profile.get("remote_workflow_dir") or "").strip()
        roots = [remote_workflow_dir] if remote_workflow_dir else []
    max_depth = int(profile.get("job_log_scan_max_depth") or profile.get("scratch_scan_max_depth", 8))
    return _scan_remote_files(profile, roots, ["*.out", "*.err"], max_depth)


def parse_result_scan(output: str) -> list[dict]:
    entries = []
    for line in output.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        path, size, mtime = parts[0], parts[1], parts[2]
        info = classify_result_path(path)
        entry = {
            "path": path,
            "name": Path(path).name,
            "size": int(size) if size.isdigit() else 0,
            "mtime": mtime,
            **info,
        }
        entries.append(entry)
    entries.sort(key=lambda item: (item["target"], item["pilot"], item["stage"], item["shard"], item["path"]))
    return entries


def classify_result_path(path: str) -> dict:
    normalized = path.replace("\\", "/")
    parts = [part for part in normalized.split("/") if part]
    joined = "/" + "/".join(parts)
    lower = joined.lower()

    target = ""
    for part in parts:
        if part.endswith("_test"):
            target = part[:-5]
        elif re.match(r"Test\d+_", part):
            target = part.split("_", 1)[1]
    if not target:
        for part in parts:
            if re.search(r"[A-Za-z].*\d*_workflow$", part):
                target = part.replace("_workflow", "")
                break

    pilot_match = re.search(r"/(pilot[^/]+)/", joined, flags=re.IGNORECASE)
    pilot = pilot_match.group(1) if pilot_match else ""

    file_name = Path(path).name.lower()
    if ("example_outputs" in lower or "rfdiffusion" in lower) and file_name.endswith("_px0_traj.pdb"):
        stage = "RFDiffusion pX0 traj"
    elif ("example_outputs" in lower or "rfdiffusion" in lower) and file_name.endswith("_xt-1_traj.pdb"):
        stage = "RFDiffusion Xt-1 traj"
    elif "example_outputs" in lower or "rfdiffusion" in lower:
        stage = "RFDiffusion"
    elif "mpnn" in lower or "proteinmpnn" in lower:
        stage = "ProteinMPNN"
    elif "afcyc" in lower or "rmsd" in lower:
        stage = "AfCycDesign"
    elif "pyrosetta" in lower or "pyro_" in lower:
        stage = "PyRosetta"
    elif "/logs/" in lower:
        stage = "Logs"
    else:
        stage = "Other"

    shard_match = re.search(r"(?:shard[_-]?|runlist_)(\d+)", joined, flags=re.IGNORECASE)
    shard = shard_match.group(1) if shard_match else ""
    file_type = Path(path).suffix.lstrip(".").lower() or "file"

    return {
        "target": target or "unknown",
        "pilot": pilot or "unknown",
        "stage": stage,
        "shard": shard or "all",
        "type": file_type,
    }


def download_destination_for_entry(
    profile: dict,
    entry: dict,
    local_dir: Path | None = None,
    layout: str = "target",
) -> Path:
    destination_root = Path(local_dir or profile.get("local_download_dir") or "downloads").expanduser()
    if not destination_root.is_absolute():
        destination_root = SCRIPT_DIR / destination_root

    target = safe_local_name(entry.get("target", "unknown"))
    pilot = safe_local_name(entry.get("pilot", "unknown"))
    stage = safe_local_name(entry.get("stage", "Other"))
    shard = safe_local_name(entry.get("shard", "all"))

    if layout == "stage":
        return destination_root / stage
    if layout == "stage_target":
        return destination_root / stage / target / pilot / shard
    return destination_root / target / pilot / stage / shard


def local_download_path_for_entry(
    profile: dict,
    entry: dict,
    local_dir: Path | None = None,
    layout: str = "target",
) -> Path:
    return download_destination_for_entry(profile, entry, local_dir=local_dir, layout=layout) / Path(entry["path"]).name


def download_files(
    profile: dict,
    entries: list[dict],
    local_dir: Path | None = None,
    layout: str = "target",
) -> list[CommandResult]:
    if not entries:
        raise ValueError("No result files selected")

    results = []
    for entry in entries:
        destination = download_destination_for_entry(profile, entry, local_dir=local_dir, layout=layout)
        destination.mkdir(parents=True, exist_ok=True)
        if _uses_password_auth(profile):
            results.append(_download_file_paramiko(profile, entry["path"], destination))
            continue
        source = f"{_ssh_target(profile)}:{entry['path']}"
        command = _base_scp_args(profile) + [source, str(destination)]
        results.append(_run(command, timeout=600))
    return results


def _download_file_paramiko(profile: dict, remote_path: str, destination: Path) -> CommandResult:
    client = _connect_paramiko(profile)
    command_label = ["paramiko-sftp-download", remote_path, str(destination)]
    try:
        sftp = client.open_sftp()
        destination.mkdir(parents=True, exist_ok=True)
        local_path = destination / Path(remote_path).name
        sftp.get(remote_path, str(local_path))
        return CommandResult(0, f"Downloaded {remote_path} -> {local_path}\n", "", command_label)
    except Exception as exc:
        return CommandResult(1, "", f"{type(exc).__name__}: {exc}\n", command_label)
    finally:
        try:
            sftp.close()
        except Exception:
            pass
        client.close()


def safe_local_name(value: str) -> str:
    text = str(value or "unknown")
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._") or "unknown"


def scratch_migration_paths(profile: dict) -> dict:
    migration = profile.get("scratch_migration", {})
    scratch_root = str(migration.get("scratch_root") or "/scratch").rstrip("/")
    source_date = str(migration.get("source_date") or "").strip()
    target_date = str(migration.get("target_date") or "auto").strip()
    user_dir = str(migration.get("user_dir") or profile.get("user") or "").strip()
    if not source_date:
        raise ValueError("scratch_migration.source_date is required")
    if not user_dir:
        raise ValueError("scratch_migration.user_dir is required")
    if target_date.lower() == "auto":
        target_date = "$(date +%F)"
    source = f"{scratch_root}/{source_date}/{user_dir}"
    target_parent = f"{scratch_root}/{target_date}"
    target = f"{target_parent}/{user_dir}"
    return {
        "scratch_root": scratch_root,
        "source_date": source_date,
        "target_date": target_date,
        "user_dir": user_dir,
        "source": source,
        "target_parent": target_parent,
        "target": target,
    }


def _scratch_path_expr(path: str) -> str:
    if "$(date +%F)" in path:
        return '"' + path.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return _remote_path_expr(path)


def preview_scratch_migration(profile: dict) -> CommandResult:
    paths = scratch_migration_paths(profile)
    source_expr = _scratch_path_expr(paths["source"])
    target_parent_expr = _scratch_path_expr(paths["target_parent"])
    target_expr = _scratch_path_expr(paths["target"])
    command = (
        "set -e; "
        f"echo SOURCE={source_expr}; "
        f"echo TARGET_PARENT={target_parent_expr}; "
        f"echo TARGET={target_expr}; "
        f"if [ -d {source_expr} ]; then echo SOURCE_EXISTS=yes; du -sh {source_expr} 2>/dev/null || true; else echo SOURCE_EXISTS=no; fi; "
        f"if [ -e {target_expr} ]; then echo TARGET_EXISTS=yes; ls -ld {target_expr}; else echo TARGET_EXISTS=no; fi; "
        f"if [ -d {source_expr} ]; then echo TOP_LEVEL_ITEMS:; find {source_expr} -mindepth 1 -maxdepth 1 -printf '%f\\n' | sort | head -50; fi"
    )
    return run_remote(profile, command, timeout=120)


def move_scratch_user_dir(profile: dict) -> CommandResult:
    paths = scratch_migration_paths(profile)
    source_expr = _scratch_path_expr(paths["source"])
    target_parent_expr = _scratch_path_expr(paths["target_parent"])
    target_expr = _scratch_path_expr(paths["target"])
    command = (
        "set -euo pipefail; "
        f"source_dir={source_expr}; "
        f"target_parent={target_parent_expr}; "
        f"target_dir={target_expr}; "
        "echo moving:$source_dir; "
        "echo to:$target_dir; "
        "if [ ! -d \"$source_dir\" ]; then echo 'ERROR: source directory does not exist' >&2; exit 2; fi; "
        "if [ -e \"$target_dir\" ]; then echo 'ERROR: target directory already exists; refusing to merge or overwrite' >&2; exit 3; fi; "
        "mkdir -p \"$target_parent\"; "
        "mv \"$source_dir\" \"$target_parent\"; "
        "echo 'MOVE_DONE'; "
        "ls -ld \"$target_dir\""
    )
    return run_remote(profile, command, timeout=3600)
