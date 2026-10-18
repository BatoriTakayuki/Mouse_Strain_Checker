"""Local Streamlit UI. Run via python -m mouse_strain web."""
import os
from pathlib import Path
import re
import signal
import subprocess
import sys

import streamlit as st

st.set_page_config(page_title="Mouse Strain Checker", layout="wide")
st.title("Mouse Strain Checker")
st.write("SRXを入力し、HISAT2でマッピングしたRNA-seqのSNPから系統候補を調べます。")
st.warning("探索用です。発現していない領域は判定できません。交雑率やF1、DNA上の系統を確定するものではありません。")
root = Path(os.environ.get("STRAIN_DATA", "./data")).resolve()
resources = st.text_input("参照データのディレクトリ", str(root / "resources"))
output = st.text_input("結果の保存先", str(root / "results"))
threads = st.number_input("CPUスレッド数（マッピングとソートを並行するため追加CPUも使用します）", min_value=1, max_value=128, value=4)
raw = st.text_area("SRX番号（空白・改行・カンマで区切り、各SRXを別々に解析）", placeholder="SRX123456\nSRX234567")
single = st.checkbox("単一系統候補のみ比較する", value=False)
st.caption("初回準備: python -m mouse_strain setup --out data/resources --threads 4（大容量データの取得とインデックス作成を行います）")
job = st.session_state.get("job")
running = job is not None and job["process"].poll() is None
if st.button("解析開始", disabled=running):
    accessions = [s for s in re.split(r"[\s,]+", raw.strip()) if s]
    if not accessions or any(not re.fullmatch(r"[SED]RX\d+", x) for x in accessions):
        st.error("SRX/ERX/DRX番号を入力してください。")
    elif not (Path(resources) / "config.json").exists():
        st.error("参照データを先にsetupで準備してください。")
    else:
        directory = Path(output).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        log_path = directory / "job.log"
        command = [sys.executable, "-m", "mouse_strain", "run", "--resources", resources,
                   "--out", str(directory), "--threads", str(threads), "--srx", *accessions]
        if single:
            command.append("--single-only")
        with open(log_path, "w") as handle:
            process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
        st.session_state.job = dict(process=process, log=log_path, output=directory, accessions=accessions)
        st.rerun()
if job:
    if running:
        st.info("解析中です。「状態を更新」で進捗を表示します。このタブを開いたままにしてください。")
        if st.button("解析を中断"):
            os.killpg(job["process"].pid, signal.SIGTERM)
            job["process"].wait()
            st.rerun()
    elif job["process"].returncode == 0:
        st.success("解析が完了しました。")
    else:
        st.error("解析が中断・失敗しました。ログを確認してください。同じ入力で再実行すると検証済みデータを再利用します。")
    st.button("状態を更新")
    if job["log"].exists():
        with open(job["log"], "rb") as handle:
            handle.seek(max(0, job["log"].stat().st_size - 20000))
            st.code(handle.read().decode(errors="replace"), language="text")
    for acc in job["accessions"]:
        directory = job["output"] / acc
        if not (directory / "complete.json").exists():
            continue
        st.subheader(acc)
        import csv
        with open(directory / "chromosomes.csv") as f:
            st.dataframe(list(csv.DictReader(f)), use_container_width=True)
        for filename, mime in (("report.html", "text/html"), ("chromosomes.csv", "text/csv"), ("windows.csv", "text/csv")):
            st.download_button(filename, (directory / filename).read_bytes(), file_name=acc + "_" + filename,
                               mime=mime, key=acc + filename)
