#!/usr/bin/env python3
"""Search for stable Harris Lundquist production envelopes over S/tlim grids."""

import argparse
import csv
import math
import shlex
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass
class EnvelopeCase:
    group_index: int
    s_values: str
    tlim: float
    nlim: int
    output_dt: float

    @property
    def case_id(self) -> str:
        tlim_tag = str(self.tlim).replace(".", "p")
        return f"g{self.group_index:02d}_t{tlim_tag}"

    @property
    def n_s(self) -> int:
        return len(parse_s_values(self.s_values))

    @property
    def max_s(self) -> float:
        values = parse_s_values(self.s_values)
        return max(values) if values else math.nan


def parse_s_values(raw: str) -> list[float]:
    values = []
    for token in raw.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        value = float(stripped)
        if value <= 0.0:
            raise ValueError(f"S values must be > 0. Got {value}")
        values.append(value)
    if not values:
        raise ValueError(f"No valid S values found in '{raw}'")
    return values


def parse_s_groups(raw: str) -> list[str]:
    groups = []
    for token in raw.split(";"):
        stripped = token.strip()
        if not stripped:
            continue
        # Validate now so failures happen early.
        parse_s_values(stripped)
        groups.append(stripped)
    if not groups:
        raise ValueError("At least one S-value group must be provided.")
    return groups


def parse_float_csv(raw: str, label: str) -> list[float]:
    values = []
    for token in raw.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        value = float(stripped)
        values.append(value)
    if not values:
        raise ValueError(f"At least one {label} value is required.")
    return values


def parse_kv_stdout(text: str) -> dict[str, str]:
    out = {}
    for line in text.splitlines():
        if "," not in line:
            continue
        key, value = line.split(",", 1)
        out[key.strip()] = value.strip()
    return out


def parse_failure_reason(stdout_text: str, stderr_text: str) -> str:
    stderr_lines = [line.strip() for line in stderr_text.splitlines() if line.strip()]
    if stderr_lines:
        return stderr_lines[-1]
    stdout_lines = [line.strip() for line in stdout_text.splitlines() if line.strip()]
    for line in reversed(stdout_lines):
        if "failed" in line.lower() or "error" in line.lower():
            return line
    return "unknown_failure"


def rel_nlim(base_nlim: int, base_tlim: float, tlim: float) -> int:
    scaled = int(math.ceil(base_nlim * (tlim / base_tlim)))
    return max(1, scaled)


def rel_output_dt(base_output_dt: float, base_tlim: float, tlim: float) -> float:
    return base_output_dt * (tlim / base_tlim)


def build_cases(
    s_groups: list[str],
    tlims: list[float],
    base_tlim: float,
    base_nlim: int,
    base_output_dt: float,
) -> list[EnvelopeCase]:
    cases = []
    for gi, s_values in enumerate(s_groups):
        for tlim in tlims:
            if tlim <= 0.0:
                raise ValueError(f"tlim must be > 0. Got {tlim}")
            cases.append(
                EnvelopeCase(
                    group_index=gi,
                    s_values=s_values,
                    tlim=tlim,
                    nlim=rel_nlim(base_nlim, base_tlim, tlim),
                    output_dt=rel_output_dt(base_output_dt, base_tlim, tlim),
                )
            )
    return cases


def run_case(
    script_path: Path,
    binary: Path,
    controlled_input: Path,
    full_input: Path,
    campaign_workdir: Path,
    case: EnvelopeCase,
    skip_plots: bool,
    skip_topology: bool,
    extra_args: list[str],
) -> dict[str, str]:
    case_workdir = campaign_workdir / "runs" / case.case_id
    case_outdir = case_workdir / "outputs"
    case_workdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(script_path),
        "--binary",
        str(binary),
        "--workdir",
        str(case_workdir),
        "--controlled-input",
        str(controlled_input),
        "--full-input",
        str(full_input),
        "--output-dir",
        str(case_outdir),
        "--s-values",
        case.s_values,
        "--tlim",
        f"{case.tlim:.6g}",
        "--nlim",
        str(case.nlim),
        "--output-dt",
        f"{case.output_dt:.6g}",
    ]
    if skip_plots:
        cmd.append("--skip-plots")
    if skip_topology:
        cmd.append("--skip-topology")
    cmd.extend(extra_args)

    proc = subprocess.run(cmd, cwd=campaign_workdir, text=True, capture_output=True)
    stdout_log = case_workdir / "production_stdout.log"
    stderr_log = case_workdir / "production_stderr.log"
    stdout_log.write_text(proc.stdout, encoding="utf-8")
    stderr_log.write_text(proc.stderr, encoding="utf-8")

    parsed = parse_kv_stdout(proc.stdout)
    overall = parsed.get("overall_status", "UNKNOWN")
    status = "PASS" if (proc.returncode == 0 and overall == "PASS") else "FAIL"
    if status == "PASS":
        reason = "none"
    elif proc.returncode != 0:
        reason = parse_failure_reason(proc.stdout, proc.stderr)
    else:
        reason = f"overall_status={overall}"

    return {
        "case_id": case.case_id,
        "status": status,
        "returncode": str(proc.returncode),
        "overall_status": overall,
        "failure_reason": reason,
        "group_index": str(case.group_index),
        "s_values": case.s_values,
        "n_s": str(case.n_s),
        "max_s": f"{case.max_s:.6e}",
        "tlim": f"{case.tlim:.6e}",
        "nlim": str(case.nlim),
        "output_dt": f"{case.output_dt:.6e}",
        "summary_csv": parsed.get("summary_csv", ""),
        "production_report": parsed.get("production_report", ""),
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
        "cmd": " ".join(shlex.quote(piece) for piece in cmd),
    }


def sort_pass_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    def key(row: dict[str, str]):
        return (
            int(row["n_s"]),
            float(row["max_s"]),
            float(row["tlim"]),
        )

    return sorted(rows, key=key, reverse=True)


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise RuntimeError("No envelope rows to write.")
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    rows: list[dict[str, str]],
    pass_rows: list[dict[str, str]],
    selected: dict[str, str] | None,
) -> None:
    lines = []
    lines.append("# Harris Lundquist Envelope Campaign")
    lines.append("")
    lines.append(f"- Generated (UTC): {datetime.now(timezone.utc).isoformat()}")
    lines.append(f"- Total cases: **{len(rows)}**")
    lines.append(f"- Passing cases: **{len(pass_rows)}**")
    if selected is None:
        lines.append("- Selected envelope: **none**")
    else:
        lines.append(
            "- Selected envelope:"
            f" **{selected['case_id']}**"
            f" (`S={selected['s_values']}`, `tlim={selected['tlim']}`, `nlim={selected['nlim']}`)"
        )
    lines.append("")
    lines.append("## Case Matrix")
    lines.append("")
    lines.append(
        "| case | status | overall | S values | tlim | nlim | n_s | max_s | failure_reason |"
    )
    lines.append("| :--- | :---: | :---: | :--- | ---: | ---: | ---: | ---: | :--- |")
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{row['case_id']}`",
                    f"`{row['status']}`",
                    f"`{row['overall_status']}`",
                    f"`{row['s_values']}`",
                    f"`{row['tlim']}`",
                    f"`{row['nlim']}`",
                    f"`{row['n_s']}`",
                    f"`{row['max_s']}`",
                    row["failure_reason"],
                ]
            )
            + " |"
        )
    lines.append("")
    lines.append("## Ranked Passing Envelopes")
    lines.append("")
    lines.append("| rank | case | S values | tlim | n_s | max_s |")
    lines.append("| ---: | :--- | :--- | ---: | ---: | ---: |")
    for idx, row in enumerate(pass_rows, start=1):
        lines.append(
            "| "
            + " | ".join(
                [
                    str(idx),
                    f"`{row['case_id']}`",
                    f"`{row['s_values']}`",
                    f"`{row['tlim']}`",
                    f"`{row['n_s']}`",
                    f"`{row['max_s']}`",
                ]
            )
            + " |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--workdir", default=".")
    parser.add_argument("--controlled-input", default="inputs/harris_4d_controlled.in")
    parser.add_argument("--full-input", default="inputs/harris_4d_full.in")
    parser.add_argument("--output-dir", default="harris_lundquist_envelope_outputs")
    parser.add_argument(
        "--s-groups",
        default="250,500,1000,2000;250,500,1000,2000,4000",
        help="Semicolon-separated S value groups. Example: 250,500;250,500,1000",
    )
    parser.add_argument("--tlims", default="0.15,0.20,0.25,0.30")
    parser.add_argument("--base-tlim", type=float, default=0.20)
    parser.add_argument("--base-nlim", type=int, default=2000)
    parser.add_argument("--base-output-dt", type=float, default=0.02)
    parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="Extra argument forwarded to run_harris_lundquist_production.sh",
    )
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--skip-topology", action="store_true", default=True)
    parser.add_argument("--max-cases", type=int, default=0)
    parser.add_argument("--fail-on-no-pass", action="store_true")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    script_path = repo_root / "scripts" / "run_harris_lundquist_production.sh"
    binary = Path(args.binary).resolve()
    workdir = Path(args.workdir).resolve()
    controlled_input = Path(args.controlled_input).resolve()
    full_input = Path(args.full_input).resolve()
    output_dir_arg = Path(args.output_dir)
    output_dir = (
        output_dir_arg.resolve()
        if output_dir_arg.is_absolute()
        else (workdir / output_dir_arg).resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.base_tlim <= 0.0:
        raise ValueError("--base-tlim must be > 0")
    if args.base_nlim <= 0:
        raise ValueError("--base-nlim must be > 0")
    if args.base_output_dt <= 0.0:
        raise ValueError("--base-output-dt must be > 0")

    s_groups = parse_s_groups(args.s_groups)
    tlims = parse_float_csv(args.tlims, "tlim")
    cases = build_cases(
        s_groups=s_groups,
        tlims=tlims,
        base_tlim=args.base_tlim,
        base_nlim=args.base_nlim,
        base_output_dt=args.base_output_dt,
    )
    if args.max_cases > 0:
        cases = cases[: args.max_cases]

    rows = []
    for case in cases:
        row = run_case(
            script_path=script_path,
            binary=binary,
            controlled_input=controlled_input,
            full_input=full_input,
            campaign_workdir=workdir,
            case=case,
            skip_plots=args.skip_plots,
            skip_topology=args.skip_topology,
            extra_args=args.extra_arg,
        )
        rows.append(row)
        print(
            ",".join(
                [
                    "envelope_case",
                    row["case_id"],
                    row["status"],
                    row["overall_status"],
                    row["s_values"].replace(",", ";"),
                    row["tlim"],
                    row["failure_reason"],
                ]
            )
        )

    pass_rows = sort_pass_rows([row for row in rows if row["status"] == "PASS"])
    selected = pass_rows[0] if pass_rows else None

    summary_csv = output_dir / "envelope_summary.csv"
    report_md = output_dir / "envelope_report.md"
    write_csv(summary_csv, rows)
    write_report(report_md, rows, pass_rows, selected)

    print(f"envelope_summary_csv,{summary_csv}")
    print(f"envelope_report,{report_md}")
    if selected is None:
        print("selected_envelope,none")
        if args.fail_on_no_pass:
            raise SystemExit("No passing envelope cases found.")
    else:
        print(
            "selected_envelope,"
            + ",".join(
                [
                    selected["case_id"],
                    selected["s_values"],
                    selected["tlim"],
                    selected["nlim"],
                ]
            )
        )


if __name__ == "__main__":
    main()
