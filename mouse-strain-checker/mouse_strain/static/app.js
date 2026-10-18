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
    sra_download: "SRAデータの取得", fastq_convert: "FASTQへの変換",
    mapping: "配列の整列・BAMの作成", merge: "runのBAM統合", sample_completed: "実験の解析完了",
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
      for (const model of [row.call, row.candidate]) modelParts(model).forEach(name => names.add(name));
    }
    return new Map([...names].sort().map((name, index) => [name, PALETTE[index] || `hsl(${Math.round(index * 137.5 % 360)} 45% 48%)`]));
  }

  function windowAppearance(row, palette) {
    if (!(numeric(row.n_sites) > 0)) return { background: UNOBSERVED, kind: "未観測" };
    if (!row.call) return { background: UNCALLED, kind: "観測あり・判定保留" };
    const parts = modelParts(row.call);
    const first = palette.get(parts[0]) || "#657d73";
    if (parts.length < 2) return { background: first, kind: "単一系統モデル" };
    const second = palette.get(parts[1]) || "#a77c9c";
    return { background: `repeating-linear-gradient(135deg, ${first} 0, ${first} 5px, ${second} 5px, ${second} 10px)`, kind: "2系統モデル" };
  }

  function selectedStillCurrent(expectedId, expectedGeneration, currentId, currentGeneration) {
    return expectedId === currentId && expectedGeneration === currentGeneration;
  }

  // Export pure validation and display rules for the Node-based regression checks.
  if (typeof document === "undefined") {
    if (typeof module !== "undefined") module.exports = { parseAccessions, parseThreads, windowAppearance, buildPalette, modelParts, formatBytes, selectedStillCurrent };
    return;
  }

  const $ = id => document.getElementById(id);
  const state = { id: null, generation: 0, pollBusy: false, pollAgain: false, job: null, history: [], health: null, healthAt: 0, results: null, resultsKey: null, sample: null, actionBusy: false, submitting: false };

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

  async function api(path, options = {}) {
    const controller = new AbortController();
    const timer = window.setTimeout(() => controller.abort(), 20000);
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

  function post(path, body = {}) {
    return api(path, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
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
    { key: "status", label: "評価状態" }, { key: "call", label: "判定", format: value => value || "保留" },
    { key: "candidate", label: "最上位候補" }, { key: "n_sites", label: "観測SNP", format: value => numberText(value) },
    { key: "diagnostic_sites", label: "判別SNP", format: value => numberText(value) },
    { key: "margin_per_site", label: "スコア差/SNP", format: value => numberText(value, 4) },
    { key: "fit_loss_per_site", label: "適合損失/SNP", format: value => numberText(value, 4) },
  ];

  function windowDescription(row) {
    return row.chromosome + ":" + numberText(row.start) + "–" + numberText(row.end) + "　" + (FINDING[row.status] || row.status) + "\n判定：" + (row.call || "保留") + "　最上位候補：" + (row.candidate || "—") + "　観測SNP：" + numberText(row.n_sites) + "　判別SNP：" + numberText(row.diagnostic_sites);
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
      [chromosomes.filter(row => row.call).length, "判定のある染色体"],
      [chromosomes.length, "評価した染色体"],
      [chromosomes.reduce((sum, row) => sum + (numeric(row.n_sites) || 0), 0), "観測SNP"],
      [windows.filter(row => numeric(row.n_sites) > 0).length, "観測のある窓"],
    ];
    const summary = document.createDocumentFragment();
    for (const [value, label] of metrics) { const node = element("div"); node.append(element("strong", numberText(value)), document.createTextNode(label)); summary.append(node); }
    $("summary").replaceChildren(summary);
    const count = sample.counts || {};
    const limited = (count.chromosomes > chromosomes.length || count.windows > windows.length || count.scores > scores.length || state.results?.truncated);
    setText("result-limit", limited ? "画面の表は各1000行までです。上の件数とバーは取得できた行の範囲を表示しています。省略分を含む全データはCSVをダウンロードしてください。" : "画面の表は各1000行まで表示します。元データの全行はCSVで保存できます。");
    renderBars(sample);
    renderTable("chromosomes-table", chromosomes, [{ key: "chromosome", label: "染色体" }, ...EVIDENCE_COLUMNS, { key: "observed_windows", label: "観測窓", format: value => numberText(value) }, { key: "called_windows", label: "判定窓", format: value => numberText(value) }, { key: "runner_up", label: "次点候補" }]);
    renderTable("windows-table", windows, [{ key: "chromosome", label: "染色体" }, { key: "start", label: "開始", format: value => numberText(value) }, { key: "end", label: "終了", format: value => numberText(value) }, ...EVIDENCE_COLUMNS]);
    renderTable("scores-table", scores, [{ key: "chromosome", label: "染色体" }, { key: "rank", label: "順位", format: value => numberText(value) }, { key: "model", label: "候補モデル" }, { key: "n_sites", label: "観測SNP", format: value => numberText(value) }, { key: "mean_log_score", label: "平均対数スコア", format: value => numberText(value, 4) }, { key: "fit_loss_per_site", label: "適合損失/SNP", format: value => numberText(value, 4) }]);
    const downloads = document.createDocumentFragment();
    const names = { "chromosomes.csv": "染色体別 CSV", "windows.csv": "固定窓 CSV", "candidate_scores.csv": "候補スコア CSV", "report.html": "HTMLレポート", "evidence.json": "解析条件・根拠 JSON", "sites.tsv.gz": "観測SNP TSV.gz", "ena_runs.json": "run情報 JSON", "mpileup.log": "SNP取得ログ", "complete.json": "完了情報 JSON" };
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
    try {
      const historyPromise = api("/api/jobs");
      const healthNeeded = !state.health || Date.now() - state.healthAt > 30000;
      const healthPromise = healthNeeded ? api("/api/health") : Promise.resolve(null);
      const jobPromise = id ? api(jobUrl(id)) : Promise.resolve(null);
      const [jobsResult, healthResult, jobResult] = await Promise.allSettled([historyPromise, healthPromise, jobPromise]);
      const errors = [];
      if (jobsResult.status === "fulfilled") {
        const jobs = Array.isArray(jobsResult.value) ? jobsResult.value : jobsResult.value.jobs || [];
        renderHistory(jobs);
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

  $("analysis-form").addEventListener("submit", event => { event.preventDefault(); submit("analyze"); });
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
