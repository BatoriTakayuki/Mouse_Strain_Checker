/* Local, dependency-free interface. API strings are always rendered as text. */
(() => {
  "use strict";

  const STATUS = {
    queued: "待機中", running: "実行中", pausing: "中断処理中", paused: "中断中",
    completed: "完了", failed: "失敗",
  };
  const FINDING = {
    candidate_supported: "単一系統候補を支持", mixed_or_heterozygous: "2系統モデルを支持",
    regional_discordance: "区間間で不一致", localized_evidence: "根拠が局所的",
    ambiguous: "判別困難", poor_fit: "モデル適合不良", insufficient_sites: "SNP不足",
    no_observed_snps: "SNP未観測",
  };
  const STAGE = {
    queued: "待機", setup: "参照の準備", reference: "参照の準備", reference_setup: "参照の準備",
    download: "配列の取得", downloading: "配列の取得", alignment: "配列の整列", align: "配列の整列",
    inference: "候補系統の評価", infer: "候補系統の評価", analyze: "RNA-seqの解析",
    reference_download: "参照データの取得", reference_extract: "参照FASTAの展開",
    reference_ready: "参照の準備完了", panel: "系統識別用SNPパネルの作成",
    hisat2_index: "ゲノム索引の作成", metadata: "実験・run情報の確認",
    sra_download: "SRAデータの取得（prefetch）", sra_validate: "SRAの検証", fastq_convert: "FASTQへの変換（fasterq-dump）",
    read_preprocess: "リードの品質処理（fastp）",
    mapping: "配列の整列・BAMの作成", merge: "runのBAM統合", sample_completed: "実験の解析完了",
    sex_inference: "性別マーカーの評価",
    completed: "完了", paused: "中断", pausing: "中断処理", failed: "処理の失敗",
  };
  const PALETTE = ["#437db1", "#35876a", "#b9719e", "#bc8a38", "#8b78b4", "#4097a4", "#c66b5e", "#839c4d", "#767db0", "#ac7758", "#48968a", "#a86e88"];
  const UNOBSERVED = "#dce1e0";
  const UNCALLED = "#f0dfb6";

  function parseAccessions(value, required = true) {
    const tokens = String(value || "").toUpperCase().split(/[\s,;、，]+/).filter(Boolean);
    const accessions = [...new Set(tokens)];
    if (required && !accessions.length) throw new Error("SRX / ERX / DRX番号を1件以上入力してください。");
    if (accessions.length > 30) throw new Error("1回に解析できる実験番号は、重複を除いて30件までです。");
    const invalid = accessions.filter(value => !/^[SED]RX\d+$/.test(value));
    if (invalid.length) throw new Error("実験番号を確認してください: " + invalid.slice(0, 5).join(", ") + "（SRX / ERX / DRX + 数字）");
    return accessions;
  }

  function parseThreads(value, fallback = 4) {
    const number = Number(String(value).trim() || fallback);
    if (!Number.isInteger(number) || number < 1 || number > 128) throw new Error("CPUスレッド数は1〜128の整数で指定してください。");
    return number;
  }

  function numeric(value) {
    return value !== "" && value !== null && value !== undefined && Number.isFinite(Number(value)) ? Number(value) : null;
  }

  function numberText(value, digits = 0) {
    const number = numeric(value);
    return number === null ? "—" : number.toLocaleString("ja-JP", { maximumFractionDigits: digits });
  }

  function formatBytes(value) {
    const size = numeric(value);
    if (size === null || size < 0) return "—";
    if (size < 1024) return numberText(size) + " B";
    const exponent = Math.min(4, Math.floor(Math.log(size) / Math.log(1024)));
    return numberText(size / (1024 ** exponent), 1) + " " + ["B", "KiB", "MiB", "GiB", "TiB"][exponent];
  }

  function modelParts(model) {
    return String(model || "").split(/\s+\+\s+/).map(value => value.trim()).filter(Boolean);
  }

  function buildPalette(sample) {
    const names = new Set();
    for (const row of [...(sample.chromosomes || []), ...(sample.windows || [])]) {
      for (const model of [row.call, row.background_call]) modelParts(model).forEach(name => names.add(name));
    }
    return new Map([...names].sort().map((name, index) => [name, PALETTE[index] || `hsl(${Math.round(index * 137.5 % 360)} 45% 48%)`]));
  }

  function windowAppearance(row, palette) {
    if (!(numeric(row.n_sites) > 0)) return { background: UNOBSERVED, kind: "未観測" };
    if (!row.call && !row.background_call) return { background: UNCALLED, kind: "観測あり・判定保留" };
    const reference = !row.call;
    const parts = modelParts(row.call || row.background_call);
    const first = palette.get(parts[0]) || "#657d73";
    const second = palette.get(parts[1]) || "#a77c9c";
    const background = parts.length < 2 ? first : `repeating-linear-gradient(135deg, ${first} 0, ${first} 5px, ${second} 5px, ${second} 10px)`;
    if (reference) return { background: `repeating-linear-gradient(45deg, #ffffffcc 0, #ffffffcc 4px, #ffffff66 4px, #ffffff66 7px), ${parts.length < 2 ? `linear-gradient(${first}, ${first})` : background}`, kind: "参考背景（亜系統は未確定）", reference: true };
    return { background, kind: parts.length < 2 ? "単一系統モデル" : "2系統モデル" };
  }

  function strictReason(row) {
    if (row.call) return "—";
    if (row.reason) return row.reason;
    const reasons = { ambiguous: "候補間の識別SNP数またはスコア差が基準未満です。", poor_fit: "候補モデルとの適合が基準を満たしません。", insufficient_sites: "評価に使えるSNPが不足しています。", no_observed_snps: "評価に使えるSNPが観測されていません。", regional_discordance: "区間ごとの判定が一致しません。", localized_evidence: "判定を支える区間数が不足しています。" };
    return reasons[row.status] || "保留理由の詳細は保存された結果を確認してください。";
  }

  function resultCall(row) {
    if (row.call) return String(row.call);
    return row.candidate_options ? "候補：" + row.candidate_options + "（未確定）" : "保留";
  }

  function candidateSummary(sample) {
    const groups = new Map();
    for (const row of sample.chromosomes || []) {
      if (row.call || !row.candidate_options) continue;
      const names = groups.get(row.candidate_options) || [];
      names.push(String(row.chromosome));
      groups.set(row.candidate_options, names);
    }
    return [...groups].map(([options, chromosomes]) => "染色体 " + chromosomes.join(", ") + "：" + options + "（未確定）").join("\n");
  }

  function sexDescription(sample) {
    const sex = sample.sex;
    if (!sex || typeof sex !== "object" || Array.isArray(sex)) return { label: "未評価（旧結果）", reason: "この保存結果には性別の推定がありません。", status: "unavailable", evidence: [], markers: [] };
    const labels = { male_likely: "オス候補（推定）", female_likely: "メス候補（推定）", inconclusive: "判定保留", unsupported_reference: "未評価（参照未対応）", unavailable: "未評価" };
    const evidenceLabels = { y_fragments: "Yマーカーの断片数", xist_fragments: "Xist領域の断片数", y_genes_supported: "基準を満たすY遺伝子数", xist_spliced_fragments: "Xist固有の接合をまたぐ断片数", xist_junctions: "観測したXist固有の接合数", control_genes_supported: "基準を満たす対照遺伝子数", control_fragments: "対照マーカー断片数（集計上限あり）", autosomal_fragments: "常染色体の断片数", total_fragments: "全断片数", mapped_fragments: "整列した断片数", y_marker_count: "観測されたYマーカー数", detected_y_markers: "観測されたYマーカー数", y_to_autosome_ratio: "Y / 常染色体の比", xist_to_autosome_ratio: "Xist / 常染色体の比", mapping_quality: "最小マッピング品質", min_mapq: "最小マッピング品質" };
    const evidence = [];
    if (sex.evidence && typeof sex.evidence === "object" && !Array.isArray(sex.evidence)) {
      for (const [key, value] of Object.entries(sex.evidence)) {
        if (value === null || value === undefined) continue;
        const readable = typeof value === "number" ? numberText(value, 6) : typeof value === "object" ? JSON.stringify(value) : String(value);
        evidence.push([evidenceLabels[key] || key, readable]);
      }
    }
    return { label: String(sex.label || labels[sex.status] || "判定保留"), reason: String(sex.reason || "推定理由の詳細は性別推定JSONを確認してください。"), status: sex.status, evidence, markers: Array.isArray(sex.markers) ? sex.markers : [] };
  }

  function sampleExplanation(sample) {
    const chromosomes = sample.chromosomes || [];
    const used = numeric(sample.diagnosis?.summary?.used_sites) ?? chromosomes.reduce((sum, row) => sum + (numeric(row.n_sites) || 0), 0);
    const strict = chromosomes.filter(row => row.call).length;
    const background = chromosomes.filter(row => row.background_call).length;
    let text = "判定に使えるSNP " + numberText(used) + "か所を評価しました。厳密判定 " + strict + "染色体、参考背景の基準を満たすもの " + background + "染色体です。";
    const unresolvedCount = chromosomes.filter(row => !row.call && row.candidate_options).length;
    if (unresolvedCount) text += " 未確定の " + unresolvedCount + "染色体は、区別しきれない候補を or で併記します。";
    if (!strict && used > 0) text += " 厳密判定がすべて保留でも、候補と保留理由から背景の手がかりを確認できます。";
    const unresolved = sample.diagnosis?.summary?.unresolved_substrain_pairs;
    for (const pair of Array.isArray(unresolved) ? unresolved : []) {
      if (Array.isArray(pair?.strains) && pair.strains.includes("C57BL_6J") && pair.strains.includes("C57BL_6NJ")) text += " C57BL/6の6J/6NJは識別不能（識別に使えるSNP " + numberText(pair.diagnostic_sites) + "か所）です。";
    }
    if (Array.isArray(sample.workflow?.stages)) text += "\n実行した処理：" + sample.workflow.stages.join(" → ");
    if (Array.isArray(sample.workflow?.quality_notes)) text += sample.workflow.quality_notes.map(note => "\n" + String(note)).join("");
    if (sample.cleanup?.scope === "cache") text += "\n中間ファイルは削除済みです。結果と性別評価は保存済みの解析結果です。";
    return text;
  }

  function canConfirmCleanup(plan, agreed, busy) {
    return Boolean(plan && typeof plan.token === "string" && plan.token.length === 64 && plan.files > 0 && agreed === true && !busy);
  }

  function selectedStillCurrent(expectedId, expectedGeneration, currentId, currentGeneration) {
    return expectedId === currentId && expectedGeneration === currentGeneration;
  }

  // Export pure validation and display rules for the Node-based regression checks.
  if (typeof document === "undefined") {
    if (typeof module !== "undefined") module.exports = { parseAccessions, parseThreads, windowAppearance, buildPalette, modelParts, formatBytes, selectedStillCurrent, strictReason, sampleExplanation, windowDescription, resultCall, candidateSummary, sexDescription, canConfirmCleanup };
    return;
  }

  const $ = id => document.getElementById(id);
  const state = { id: null, generation: 0, storageVersion: 0, pollBusy: false, pollAgain: false, job: null, history: [], health: null, healthAt: 0, results: null, resultsKey: null, sample: null, actionBusy: false, submitting: false };
  const storageState = { busy: false, preview: null };

  function element(tag, content, className) {
    const node = document.createElement(tag);
    if (content !== undefined && content !== null) node.textContent = String(content);
    if (className) node.className = className;
    return node;
  }

  function setText(id, value) { $(id).textContent = value === undefined || value === null ? "" : String(value); }
  function setError(id, error) { setText(id, error instanceof Error ? error.message : error); }
  function current(id, generation) { return selectedStillCurrent(id, generation, state.id, state.generation); }
  function jobUrl(id = state.id) { return "/api/jobs/" + encodeURIComponent(id); }
  function fileUrl(accession, name) { return jobUrl() + "/files/" + encodeURIComponent(accession) + "/" + encodeURIComponent(name); }

  async function api(path, options = {}, timeoutMs = 20000) {
    const controller = new AbortController();
    const timer = window.setTimeout(() => controller.abort(), timeoutMs);
    try {
      const response = await fetch(path, { cache: "no-store", credentials: "same-origin", ...options, signal: controller.signal });
      const body = await response.text();
      let data;
      try { data = body ? JSON.parse(body) : {}; } catch (_) { data = {}; }
      if (!response.ok) throw new Error(typeof data.error === "string" ? data.error : typeof data.message === "string" ? data.message : "サーバーがエラーを返しました（" + response.status + "）。");
      return data;
    } catch (error) {
      if (error.name === "AbortError") throw new Error("サーバーからの応答がありません。接続を確認しています。");
      if (error instanceof TypeError) throw new Error("サーバーに接続できません。起動状態を確認してください。接続は自動で再試行します。");
      throw error;
    } finally { window.clearTimeout(timer); }
  }

  function post(path, body = {}, timeoutMs = 20000) {
    return api(path, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) }, timeoutMs);
  }

  function dateText(value) {
    const date = new Date(value);
    return !value || Number.isNaN(date.getTime()) ? "—" : date.toLocaleString("ja-JP", { month: "2-digit", day: "2-digit", hour: "2-digit", minute: "2-digit", second: "2-digit", hour12: false });
  }

  function elapsedText(job) {
    const start = new Date(job.attempt_started_at || job.started_at || job.created_at).getTime();
    const terminal = ["completed", "failed", "paused"].includes(job.status);
    const end = terminal ? new Date(job.updated_at).getTime() : Date.now();
    if (!Number.isFinite(start) || !Number.isFinite(end)) return "経過時間 —";
    const seconds = Math.max(0, Math.floor((end - start) / 1000));
    const hours = Math.floor(seconds / 3600);
    const minutes = Math.floor(seconds % 3600 / 60);
    return "今回の経過 " + (hours ? hours + "時間 " : "") + minutes + "分 " + seconds % 60 + "秒";
  }

  function requestAccessions(job) {
    const value = job.accessions || job.request?.srx || [];
    return Array.isArray(value) ? value : String(value).split(/[\s,]+/).filter(Boolean);
  }

  function jobTitle(job) {
    if (job.request?.mode === "setup") return "参照データの準備";
    const accessions = requestAccessions(job);
    return accessions.length ? accessions[0] + (accessions.length > 1 ? " ほか" + (accessions.length - 1) + "件" : "") : "RNA-seqの解析";
  }

  function renderHistory(jobs) {
    state.history = jobs;
    const fragment = document.createDocumentFragment();
    if (!jobs.length) fragment.append(element("p", "まだ解析履歴はありません。", "muted"));
    for (const job of jobs) {
      const button = element("button", null, "history-item" + (job.id === state.id ? " selected" : ""));
      button.type = "button";
      button.setAttribute("aria-pressed", String(job.id === state.id));
      button.append(element("strong", jobTitle(job)));
      const meta = element("span", null, "history-meta");
      meta.append(element("span", STATUS[job.status] || job.status, "badge " + (Object.hasOwn(STATUS, job.status) ? job.status : "")), element("small", dateText(job.created_at)));
      button.append(meta);
      button.addEventListener("click", () => selectJob(job.id));
      fragment.append(button);
    }
    $("history").replaceChildren(fragment);
  }

  function renderHealth(data) {
    state.health = data;
    const missing = data.tools && typeof data.tools === "object" ? Object.entries(data.tools).filter(([, found]) => !found).map(([name]) => name) : [];
    const details = [];
    details.push(data.reference_ready ? "参照データ：準備済み" : "参照データ：初回解析時に自動準備");
    if (numeric(data.logical_cpus) !== null) details.push("論理CPU " + numberText(data.logical_cpus));
    if (numeric(data.default_threads) !== null) {
      details.push("既定 " + numberText(data.default_threads) + "スレッド");
      $("threads").placeholder = "既定 " + data.default_threads + "スレッド";
    }
    if (numeric(data.disk_free_gb) !== null) details.push("ディスク空き " + numberText(data.disk_free_gb, 1) + " GB");
    if (missing.length) details.push("不足ツール：" + missing.join(", ") + "（Docker環境を確認してください）");
    setText("health", details.join("　/　"));
    $("health").classList.toggle("warning", missing.length > 0);
  }

  function updateActionButtons(job) {
    $("pause").hidden = !["queued", "running", "pausing"].includes(job.status);
    $("resume").hidden = job.status !== "paused";
    $("retry").hidden = job.status !== "failed";
    for (const id of ["pause", "resume", "retry"]) $(id).disabled = state.actionBusy || job.status === "pausing";
  }

  function renderJob(job) {
    state.job = job;
    $("empty").hidden = true;
    $("job").hidden = false;
    setText("job-id", "JOB / " + job.id);
    setText("job-title", jobTitle(job));
    const badge = $("status");
    badge.textContent = STATUS[job.status] || job.status;
    badge.className = "badge " + (Object.hasOwn(STATUS, job.status) ? job.status : "");
    setText("job-message", job.message || "処理状態を確認しています…");
    const active = ["queued", "running", "pausing"].includes(job.status);
    const rawProgress = numeric(job.progress);
    const progress = job.status === "completed" ? 100 : active && rawProgress === 0 ? null : rawProgress;
    if (progress === null && active) $("progress").removeAttribute("value");
    else if (progress === null) $("progress").value = 0;
    else $("progress").value = Math.max(0, Math.min(100, progress));
    setText("percent", progress === null ? (job.status === "queued" ? "順番待ち" : active ? "工程を処理中" : STATUS[job.status] || "") : numberText(Math.max(0, Math.min(100, progress))) + "%（目安）");
    $("progress").setAttribute("aria-label", progress === null && active ? "工程を処理中。全体の残り時間は未計測です。" : "工程の進捗目安");
    setText("stage", STAGE[job.stage] || job.stage || STATUS[job.status] || "");
    setText("elapsed", elapsedText(job));
    setText("last-update", "最終更新 " + dateText(job.updated_at));
    const doneAccessions = job.completed_accessions || job.done_accessions;
    const completed = Array.isArray(doneAccessions) ? doneAccessions.length : numeric(doneAccessions) || 0;
    const accessions = requestAccessions(job);
    const metrics = [];
    if (job.request?.mode !== "setup" && accessions.length) metrics.push("完了 " + completed + " / " + accessions.length + "実験");
    if (numeric(job.runs_total) !== null) metrics.push("run " + numberText(job.runs_done || 0) + " / " + numberText(job.runs_total));
    setText("run-count", metrics.join("　/　"));
    const activity = [];
    if (job.current_accession || job.accession) activity.push("処理中の実験：" + (job.current_accession || job.accession));
    if (job.activity) activity.push(job.activity);
    setText("activity", activity.join("\n"));
    const transfer = job.transfer;
    $("transfer").hidden = !transfer;
    if (transfer) {
      const done = numeric(transfer.done_bytes), total = numeric(transfer.total_bytes);
      if (done !== null && total > 0) $("transfer-progress").value = Math.max(0, Math.min(100, done / total * 100));
      else $("transfer-progress").removeAttribute("value");
      const rate = numeric(transfer.rate_bps);
      setText("transfer-text", (transfer.name || "ダウンロード") + "　" + formatBytes(done) + (total > 0 ? " / " + formatBytes(total) : "（総量未取得）") + (rate !== null && rate >= 0 ? "　" + formatBytes(rate) + "/秒" : ""));
    }
    updateActionButtons(job);
    const info = {
      pausing: "現在の処理を停止しています。中断完了後に再開できます。",
      paused: "中断しています。再開時は保存済みデータを利用し、未完了の工程を実行します。",
      failed: "ログを確認してください。再試行では設定と保存済みデータを利用します。",
      completed: job.request?.mode === "setup" ? "参照の準備が完了しました。実験番号を入力して解析を開始できます。" : "解析が完了しました。結果はこの下から確認できます。",
    };
    setText("pause-info", info[job.status] || "ブラウザを閉じても処理は続きます。履歴から状態を確認できます。");
    $("log-download").href = jobUrl() + "/files/pipeline.log";
  }

  function renderTable(id, rows, columns) {
    const table = $(id);
    const caption = table.querySelector("caption");
    const head = element("thead"), header = element("tr"), body = element("tbody");
    for (const column of columns) {
      const cell = element("th", column.label);
      cell.scope = "col";
      header.append(cell);
    }
    head.append(header);
    for (const row of rows.slice(0, 1000)) {
      const tr = element("tr");
      for (const column of columns) {
        const cell = element("td");
        if (column.key === "status") {
          cell.className = "status-cell " + (row.call ? "status-supported" : "status-hold");
          cell.append(element("span", FINDING[row.status] || row.status || "—"), element("small", row.status || ""));
        } else if (column.key === "call") {
          cell.className = "result-call-cell " + (row.call ? "status-supported" : "status-hold");
          cell.append(element("strong", resultCall(row)));
        } else if (column.key === "candidate_options") {
          cell.className = "candidate-options-cell";
          cell.append(element("strong", row.candidate_options || "—"));
          if (row.candidate_options && numeric(row.candidate_count) !== null) cell.append(element("small", numberText(row.candidate_count) + (row.call ? "比較候補" : "候補・未確定")));
          if (row.candidate_options_reason) cell.append(element("small", row.candidate_options_reason));
        } else if (column.key === "background_call") {
          cell.className = "background-cell";
          cell.append(element("strong", row.background_call || (row.background_status ? "保留" : "未評価（旧結果）")));
          if (row.background_candidate) cell.append(element("small", "候補：" + row.background_candidate));
          if (row.background_runner_up) cell.append(element("small", "次の背景候補：" + row.background_runner_up));
          if (row.background_reason) cell.append(element("small", row.background_reason));
          if (row.background_status) cell.append(element("small", "群間の判別SNP " + numberText(row.background_diagnostic_sites) + " / 差 " + numberText(row.background_margin_per_site, 4) + " / 適合損失 " + numberText(row.background_fit_loss_per_site, 4)));
          if (row.background_members) cell.append(element("small", "候補内訳：" + row.background_members));
        } else if (column.format) cell.textContent = column.format(row[column.key], row);
        else cell.textContent = row[column.key] === "" || row[column.key] === null || row[column.key] === undefined ? "—" : String(row[column.key]);
        tr.append(cell);
      }
      body.append(tr);
    }
    if (!rows.length) {
      const cell = element("td", "表示できるデータはありません。", "no-data");
      cell.colSpan = columns.length;
      const tr = element("tr"); tr.append(cell); body.append(tr);
    }
    table.replaceChildren(...(caption ? [caption] : []), head, body);
  }

  const EVIDENCE_COLUMNS = [
    { key: "status", label: "厳密評価" }, { key: "call", label: "結果（未確定なら候補）" },
    { key: "candidate_options", label: "絞り込んだ候補と理由" },
    { key: "background_call", label: "参考背景（6J/6NJをまとめた比較）" },
    { key: "reason", label: "厳密判定の保留理由", format: (_, row) => strictReason(row) },
    { key: "candidate", label: "厳密比較の最上位候補" }, { key: "n_sites", label: "使えるSNP", format: value => numberText(value) },
    { key: "diagnostic_sites", label: "厳密比較の判別SNP", format: value => numberText(value) },
    { key: "margin_per_site", label: "厳密比較の差/SNP", format: value => numberText(value, 4) },
    { key: "fit_loss_per_site", label: "厳密比較の適合損失/SNP", format: value => numberText(value, 4) },
  ];

  function windowDescription(row) {
    let text = row.chromosome + ":" + numberText(row.start) + "–" + numberText(row.end) + "　" + (FINDING[row.status] || row.status) + "\n厳密判定：" + (row.call || "保留") + "　最上位候補：" + (row.candidate || "—") + "　使えるSNP：" + numberText(row.n_sites) + "　判別SNP：" + numberText(row.diagnostic_sites);
    if (!row.call) text += "\n厳密判定の保留理由：" + strictReason(row);
    if (!row.call && row.candidate_options) text += "\n候補：" + row.candidate_options + "（未確定）" + (row.candidate_options_reason ? "\n候補の選定理由：" + row.candidate_options_reason : "");
    text += "\n参考背景：" + (row.background_call || (row.background_status ? "保留" : "未評価（旧結果）"));
    if (row.background_candidate) text += "　候補：" + row.background_candidate;
    if (row.background_runner_up) text += "　次の背景候補：" + row.background_runner_up;
    if (row.background_reason) text += "\n参考背景の評価理由：" + row.background_reason;
    if (row.background_status) text += "　群間の判別SNP：" + numberText(row.background_diagnostic_sites) + "　差/SNP：" + numberText(row.background_margin_per_site, 4);
    return text;
  }

  function renderBars(sample) {
    const palette = buildPalette(sample), legend = document.createDocumentFragment();
    const legendItem = (label, background) => {
      const item = element("span", null, "legend-item"), swatch = element("span", null, "swatch");
      swatch.style.background = background;
      item.append(swatch, document.createTextNode(label)); legend.append(item);
    };
    for (const [name, color] of palette) legendItem(name, color);
    legendItem("SNP未観測", UNOBSERVED);
    legendItem("観測あり・判定保留", UNCALLED);
    legendItem("2系統モデル", "repeating-linear-gradient(135deg,#657d73 0,#657d73 5px,#a9bba3 5px,#a9bba3 10px)");
    legendItem("参考背景（淡色・細い斜線）", windowAppearance({ n_sites: 1, background_call: "参考" }, new Map()).background);
    $("legend").replaceChildren(legend);
    const grouped = new Map();
    for (const chromosome of sample.chromosomes || []) grouped.set(String(chromosome.chromosome), []);
    for (const row of (sample.windows || []).slice(0, 1000)) {
      const key = String(row.chromosome);
      if (!grouped.has(key)) grouped.set(key, []);
      grouped.get(key).push(row);
    }
    const fragment = document.createDocumentFragment();
    for (const [name, rows] of grouped) {
      const container = element("div", null, "chr-row");
      container.append(element("span", name, "chr-label"));
      const track = element("div", null, "chr-track");
      track.setAttribute("aria-label", "染色体 " + name + " の固定窓");
      rows.sort((a, b) => Number(a.start) - Number(b.start));
      for (const row of rows) {
        const appearance = windowAppearance(row, palette), cell = element("button", null, "window-cell");
        if (appearance.reference) cell.classList.add("reference-window");
        cell.type = "button";
        cell.style.flexGrow = String(Math.max(1, (numeric(row.end) || 0) - (numeric(row.start) || 0) + 1));
        cell.style.flexBasis = "0";
        cell.style.background = appearance.background;
        const description = windowDescription(row);
        cell.title = description;
        cell.setAttribute("aria-label", description);
        cell.addEventListener("click", () => setText("window-detail", description));
        track.append(cell);
      }
      if (!rows.length) track.append(element("span", "窓データなし", "muted"));
      container.append(track); fragment.append(container);
    }
    if (!grouped.size) fragment.append(element("p", "表示できる窓データはありません。", "muted"));
    $("chromosome-bars").replaceChildren(fragment);
    setText("window-detail", "区間を選択して根拠を確認してください。");
  }

  function renderSample(sample) {
    state.sample = sample.accession;
    const chromosomes = sample.chromosomes || [], windows = sample.windows || [], scores = sample.scores || [];
    const metrics = [
      [chromosomes.filter(row => row.call).length, "厳密判定のある染色体"],
      [chromosomes.filter(row => row.background_call).length, "参考背景の基準を満たす染色体"],
      [numeric(sample.diagnosis?.summary?.used_sites) ?? chromosomes.reduce((sum, row) => sum + (numeric(row.n_sites) || 0), 0), "判定に使えるSNP"],
      [windows.filter(row => numeric(row.n_sites) > 0).length, "観測のある窓"],
    ];
    const summary = document.createDocumentFragment();
    for (const [value, label] of metrics) { const node = element("div"); node.append(element("strong", numberText(value)), document.createTextNode(label)); summary.append(node); }
    $("summary").replaceChildren(summary);
    setText("result-explanation", sampleExplanation(sample));
    const candidateText = candidateSummary(sample);
    setText("candidate-summary", candidateText);
    $("candidate-summary-panel").hidden = !candidateText;
    const sex = sexDescription(sample);
    setText("sex-label", sex.label);
    setText("sex-reason", sex.reason);
    $("sex-label").className = "sex-label " + (["male_likely", "female_likely"].includes(sex.status) ? "sex-supported" : "sex-hold");
    const sexEvidence = document.createDocumentFragment();
    for (const [label, value] of sex.evidence) sexEvidence.append(element("dt", label), element("dd", value));
    $("sex-evidence").replaceChildren(sexEvidence);
    $("sex-details").hidden = !sex.evidence.length && !sex.markers.length;
    renderTable("sex-markers-table", sex.markers, [{ key: "gene", label: "マーカー遺伝子" }, { key: "chromosome", label: "染色体" }, { key: "start", label: "開始", format: value => numberText(value) }, { key: "end", label: "終了", format: value => numberText(value) }, { key: "fragments", label: "観測断片数", format: (value, row) => (row.count_capped ? "≥ " : "") + numberText(value) }]);
    const assumption = sample.diagnosis?.reference_assumption;
    setText("reference-assumption", assumption ? "参照の前提：" + (sample.diagnosis.added_reference || "追加候補（既定はC57BL_6J）") + "は参照配列のREFアレルと仮定した候補で、実測遺伝型ではありません。 " + assumption : "");
    $("reference-assumption").hidden = !assumption;
    const count = sample.counts || {};
    const limited = (count.chromosomes > chromosomes.length || count.windows > windows.length || count.scores > scores.length || state.results?.truncated);
    setText("result-limit", limited ? "画面の表は各1000行までです。上の件数とバーは取得できた行の範囲を表示しています。省略分を含む全データはCSVをダウンロードしてください。" : "画面の表は各1000行まで表示します。元データの全行はCSVで保存できます。");
    renderBars(sample);
    renderTable("chromosomes-table", chromosomes, [{ key: "chromosome", label: "染色体" }, ...EVIDENCE_COLUMNS, { key: "observed_windows", label: "観測窓", format: value => numberText(value) }, { key: "called_windows", label: "判定窓", format: value => numberText(value) }, { key: "runner_up", label: "次点候補" }]);
    renderTable("windows-table", windows, [{ key: "chromosome", label: "染色体" }, { key: "start", label: "開始", format: value => numberText(value) }, { key: "end", label: "終了", format: value => numberText(value) }, ...EVIDENCE_COLUMNS]);
    renderTable("scores-table", scores, [{ key: "chromosome", label: "染色体" }, { key: "rank", label: "順位", format: value => numberText(value) }, { key: "model", label: "候補モデル" }, { key: "n_sites", label: "観測SNP", format: value => numberText(value) }, { key: "mean_log_score", label: "平均対数スコア", format: value => numberText(value, 4) }, { key: "fit_loss_per_site", label: "適合損失/SNP", format: value => numberText(value, 4) }]);
    const downloads = document.createDocumentFragment();
    const names = { "chromosomes.csv": "染色体別 CSV", "windows.csv": "固定窓 CSV", "candidate_scores.csv": "候補スコア CSV", "report.html": "HTMLレポート", "evidence.json": "解析条件・根拠 JSON", "diagnosis.json": "保留理由・背景評価 JSON", "sex.json": "性別推定・根拠 JSON", "workflow.json": "RNA-seq処理・品質管理 JSON", "cleanup.json": "中間ファイル削除記録 JSON", "sites.tsv.gz": "観測SNP TSV.gz", "ena_runs.json": "run情報 JSON", "mpileup.log": "SNP取得ログ", "complete.json": "完了情報 JSON" };
    for (const file of sample.files || []) {
      const name = typeof file === "string" ? file : file.name;
      if (!name) continue;
      const link = element("a", names[name] || name);
      link.href = fileUrl(sample.accession, name);
      link.download = name;
      downloads.append(link);
    }
    if (!downloads.childNodes.length) downloads.append(element("p", "ダウンロードできる結果はありません。", "muted"));
    $("downloads").replaceChildren(downloads);
  }

  function renderResults(data) {
    state.results = data;
    const samples = Array.isArray(data.samples) ? data.samples : [];
    if (!samples.length) { $("results").hidden = true; return; }
    const selected = samples.find(sample => sample.accession === state.sample) || samples[0];
    const fragment = document.createDocumentFragment();
    for (const sample of samples) {
      const option = element("option", sample.accession); option.value = sample.accession;
      fragment.append(option);
    }
    $("sample-select").replaceChildren(fragment);
    $("sample-select").value = selected.accession;
    $("results").hidden = false;
    renderSample(selected);
  }

  function hashId() {
    try { return decodeURIComponent(location.hash.replace(/^#(?:job=)?/, "")) || null; }
    catch (_) { return null; }
  }

  function selectJob(id, updateHash = true) {
    if (!id) return;
    state.id = String(id);
    state.generation += 1;
    state.job = null;
    state.results = null;
    state.resultsKey = null;
    state.sample = null;
    $("empty").hidden = true;
    $("job").hidden = false;
    $("results").hidden = true;
    setText("job-id", "JOB / " + state.id);
    setText("job-title", "解析状態を取得しています…");
    setText("job-message", "保存済みの解析に接続しています。");
    setText("job-error", "");
    setText("log", "ログを読み込んでいます…");
    $("log-download").href = jobUrl() + "/files/pipeline.log";
    for (const key of ["pause", "resume", "retry"]) $(key).hidden = true;
    for (const key of ["status", "stage", "percent", "elapsed", "run-count", "last-update", "activity", "pause-info"]) setText(key, "");
    $("transfer").hidden = true;
    $("progress").removeAttribute("value");
    if (updateHash) history.replaceState(null, "", "#job=" + encodeURIComponent(state.id));
    renderHistory(state.history);
    poll();
  }

  async function poll() {
    if (state.pollBusy) { state.pollAgain = true; return; }
    state.pollBusy = true;
    state.pollAgain = false;
    const id = state.id, generation = state.generation;
    const storageVersion = state.storageVersion;
    try {
      const historyPromise = api("/api/jobs");
      const healthNeeded = !state.health || Date.now() - state.healthAt > 30000;
      const healthPromise = healthNeeded ? api("/api/health") : Promise.resolve(null);
      const jobPromise = id ? api(jobUrl(id)) : Promise.resolve(null);
      const [jobsResult, healthResult, jobResult] = await Promise.allSettled([historyPromise, healthPromise, jobPromise]);
      const errors = [];
      if (jobsResult.status === "fulfilled") {
        const jobs = Array.isArray(jobsResult.value) ? jobsResult.value : jobsResult.value.jobs || [];
        if (storageVersion === state.storageVersion) renderHistory(jobs);
      } else errors.push(jobsResult.reason);
      if (healthResult.status === "fulfilled" && healthResult.value) { renderHealth(healthResult.value); state.healthAt = Date.now(); }
      else if (healthResult.status === "rejected") errors.push(healthResult.reason);
      if (id && current(id, generation)) {
        if (jobResult.status === "fulfilled" && jobResult.value) {
          const job = jobResult.value;
          renderJob(job);
          const key = [job.id, job.status, ...(job.completed_accessions || [])].join("|");
          const needsResults = job.request?.mode !== "setup" && job.status === "completed";
          const extra = await Promise.allSettled([
            needsResults && key !== state.resultsKey ? api(jobUrl(id) + "/results") : Promise.resolve(null),
            $("log-panel").open ? api(jobUrl(id) + "/log") : Promise.resolve(null),
          ]);
          if (current(id, generation)) {
            if (extra[0].status === "fulfilled" && extra[0].value) { renderResults(extra[0].value); state.resultsKey = key; }
            if (extra[1].status === "fulfilled" && extra[1].value) setText("log", extra[1].value.text || "まだログはありません。");
            for (const result of extra) if (result.status === "rejected") errors.push(result.reason);
          }
        } else if (jobResult.status === "rejected") errors.push(jobResult.reason);
      }
      if (errors.length) { setError("connection", errors[0]); $("connection").hidden = false; }
      else { setText("connection", ""); $("connection").hidden = true; }
    } catch (error) {
      setError("connection", error); $("connection").hidden = false;
    } finally {
      state.pollBusy = false;
      if (state.pollAgain || generation !== state.generation) { state.pollAgain = false; window.setTimeout(poll, 0); }
    }
  }

  async function submit(mode) {
    if (state.submitting) return;
    setText("form-error", "");
    let request;
    try {
      const accessions = mode === "setup" ? [] : parseAccessions($("srx").value);
      request = { mode, srx: accessions.join("\n"), threads: parseThreads($("threads").value, state.health?.default_threads || 4), single_only: $("single-only").checked, pause_after_reference: mode === "analyze" && $("pause-after-reference").checked };
    } catch (error) { setError("form-error", error); return; }
    state.submitting = true;
    $("submit").disabled = $("setup").disabled = true;
    try {
      const job = await post("/api/jobs", request);
      selectJob(job.id);
    } catch (error) { setError("form-error", error); }
    finally { state.submitting = false; $("submit").disabled = $("setup").disabled = false; }
  }

  async function jobAction(action) {
    if (!state.id || state.actionBusy || !state.job) return;
    const id = state.id;
    state.actionBusy = true;
    // Invalidate any GET started before this mutation so it cannot roll back the UI.
    state.generation += 1;
    const generation = state.generation;
    setText("job-error", "");
    updateActionButtons(state.job);
    try {
      const job = await post(jobUrl(id) + "/" + action);
      if (current(id, generation)) {
        if (job.id) renderJob(job);
        state.resultsKey = null;
      }
    } catch (error) { if (current(id, generation)) setError("job-error", error); }
    finally { state.actionBusy = false; if (state.job) updateActionButtons(state.job); poll(); }
  }

  function renderStorage(data) {
    setText("storage-space", "ディスク空き " + formatBytes(data.free_bytes) + " / " + formatBytes(data.total_bytes));
    const head = element("thead"), header = element("tr"), body = element("tbody");
    for (const label of ["解析", "状態", "全体の容量", "削除できる中間ファイル", "操作"]) header.append(element("th", label));
    head.append(header);
    for (const row of data.jobs || []) {
      const tr = element("tr"), title = element("td"), actions = element("td", null, "storage-actions");
      const select = element("button", row.accessions.join(", ") || "参照準備", "text-button");
      select.type = "button"; select.addEventListener("click", () => selectJob(row.id));
      title.append(select, element("small", row.id.slice(0, 12)));
      tr.append(title, element("td", STATUS[row.status] || row.status), element("td", formatBytes(row.total_bytes)), element("td", formatBytes(row.cache_bytes)));
      for (const [scope, label, size] of [["cache", "中間ファイルを削除", row.cache_bytes], ["job", "解析全体を削除", row.total_bytes]]) {
        const button = element("button", label, scope === "job" ? "danger-outline" : "secondary");
        button.type = "button"; button.disabled = !row.deletable || !(size > 0);
        button.addEventListener("click", () => previewCleanup(scope, row.id)); actions.append(button);
      }
      if (row.reason) actions.append(element("small", row.reason));
      tr.append(actions); body.append(tr);
    }
    $("storage-table").replaceChildren(head, body);
    const reference = data.reference || {};
    setText("reference-storage", (reference.deletable ? "削除対象 " + formatBytes(reference.bytes) : "削除対象なし、または現在は削除できません。") + (reference.reason ? " " + reference.reason : ""));
    $("cleanup-reference").disabled = !reference.deletable;
  }

  async function loadStorage() {
    if (storageState.busy) return;
    storageState.busy = true; $("refresh-storage").disabled = true; setText("storage-error", "");
    try { renderStorage(await api("/api/storage", {}, 60000)); }
    catch (error) { setError("storage-error", error); }
    finally { storageState.busy = false; $("refresh-storage").disabled = false; }
  }

  async function previewCleanup(scope, jobId = null) {
    if (storageState.busy) return;
    storageState.busy = true; setText("storage-error", ""); setText("storage-feedback", "");
    try {
      const request = jobId ? { scope, job_id: jobId } : { scope };
      const plan = await post("/api/storage/preview", request, 60000);
      storageState.preview = plan;
      setText("cleanup-target", plan.title + "：" + formatBytes(plan.bytes) + "、" + numberText(plan.files) + "ファイル");
      setText("cleanup-effect", plan.effect); setText("cleanup-error", "");
      $("cleanup-paths").replaceChildren(...(plan.paths || []).map(path => element("li", path)));
      $("cleanup-agree").checked = false; $("confirm-cleanup").disabled = true;
      $("cleanup-dialog").showModal();
    } catch (error) { setError("storage-error", error); }
    finally { storageState.busy = false; }
  }

  async function executeCleanup(event) {
    event.preventDefault();
    const plan = storageState.preview;
    if (!canConfirmCleanup(plan, $("cleanup-agree").checked, storageState.busy)) return;
    storageState.busy = true; $("confirm-cleanup").disabled = true; $("cancel-cleanup").disabled = true;
    setText("confirm-cleanup", "削除しています…");
    state.generation += 1; state.storageVersion += 1;
    try {
      const result = await post("/api/storage/cleanup", { token: plan.token, confirm: true }, 300000);
      setText("storage-feedback", result.message + " 削除したファイル容量：" + formatBytes(result.bytes));
      if (plan.job_id === state.id) {
        state.results = null; state.resultsKey = null; $("results").hidden = true;
        if (plan.scope === "job") {
          state.id = null; state.job = null; state.sample = null;
          $("job").hidden = true; $("empty").hidden = false;
          history.replaceState(null, "", location.pathname + location.search);
        }
      }
      $("cleanup-dialog").close();
    } catch (error) {
      storageState.preview = null; setError("cleanup-error", error);
      state.resultsKey = null;
      if (plan.job_id === state.id) { state.results = null; $("results").hidden = true; }
    } finally {
      storageState.busy = false; $("cancel-cleanup").disabled = false;
      setText("confirm-cleanup", "削除する");
      $("confirm-cleanup").disabled = true; state.healthAt = 0;
      await loadStorage(); poll();
    }
  }

  $("analysis-form").addEventListener("submit", event => { event.preventDefault(); submit("analyze"); });
  $("storage-panel").addEventListener("toggle", () => { if ($("storage-panel").open) loadStorage(); });
  $("refresh-storage").addEventListener("click", loadStorage);
  $("cleanup-reference").addEventListener("click", () => previewCleanup("reference_downloads"));
  $("cleanup-form").addEventListener("submit", executeCleanup);
  $("cleanup-agree").addEventListener("change", () => { $("confirm-cleanup").disabled = !canConfirmCleanup(storageState.preview, $("cleanup-agree").checked, storageState.busy); });
  $("cancel-cleanup").addEventListener("click", () => $("cleanup-dialog").close());
  $("cleanup-dialog").addEventListener("close", () => { storageState.preview = null; });
  $("cleanup-dialog").addEventListener("cancel", event => { if (storageState.busy) event.preventDefault(); });
  $("setup").addEventListener("click", () => submit("setup"));
  $("pause").addEventListener("click", () => jobAction("pause"));
  $("resume").addEventListener("click", () => jobAction("resume"));
  $("retry").addEventListener("click", () => jobAction("resume"));
  $("refresh-history").addEventListener("click", poll);
  $("log-panel").addEventListener("toggle", () => { if ($("log-panel").open) poll(); });
  $("sample-select").addEventListener("change", () => {
    const sample = state.results?.samples?.find(value => value.accession === $("sample-select").value);
    if (sample) renderSample(sample);
  });
  window.addEventListener("hashchange", () => { const id = hashId(); if (id && id !== state.id) selectJob(id, false); });
  window.addEventListener("online", poll);
  window.setInterval(() => { if (state.job) setText("elapsed", elapsedText(state.job)); }, 1000);
  window.setInterval(poll, 2000);
  const initialId = hashId();
  if (initialId) selectJob(initialId, false);
  else poll();
})();
