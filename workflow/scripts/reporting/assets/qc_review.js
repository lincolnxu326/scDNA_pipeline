/* QC review report — behaviour. ES5, no dependencies, runs from file://.
   Reads window.QC (see design/handoff/README.md § Data contract) and renders into the template's hooks.

   This is design/handoff/qc_review.js with four agreed deltas — `diff` against design/handoff/ to see them:
     1. buildGrid   — sets the 96 cell size (spec §3.3); 384 keeps the CSS default.
     2. copy-cmd    — writes to QC.decisions_path instead of the terminal's CWD.
     3. profile note— bin size from QC.bins_label, not a 500 kb literal.
     4. applyMode() — the cn viewer: hides the writing controls, shows the heatmap. */
(function () {
  'use strict';

  var QC = window.QC || {};
  var WELLS = QC.wells || [];
  var GATE = QC.gate || { warn: 50000, pass: 100000 };
  var IS_REVIEW = QC.mode !== 'cn';
  /* Not every plate was dispensed on a CellenONE, and a run folder can be missing or
     unmapped. When there is no cell image ANYWHERE on the plate the report drops the
     whole image layer rather than showing 96/384 empty frames: no image panel, no
     droplet call on the hover line. The evidence a reviewer gets is then exactly the
     copy-number profile and the bin-count histogram. Set by generate_qc_review.py,
     which counts the images it actually resolved. */
  var HAS_IMAGES = !!QC.cell_images;
  var KEY = 'qc_review:' + (QC.plate || 'plate') + ':' + (QC.size || 384);
  var DECISIONS = ['PASS', 'EXCLUDE', 'REVIEW', 'REPEAT'];
  var REASON_GROUPS = [
    { title: 'Read counts',   items: ['low_read_count', 'missing_qc_metric'] },
    { title: 'Profile shape', items: ['noisy_profile', 'poor_bin_distribution', 'low_complexity'] },
    { title: 'Well contents', items: ['suspected_doublet_or_mixed_well', 'sample_swap_suspected'] },
    { title: 'Escape hatch',  items: ['manual_exception', 'other'] }
  ];
  var CALL_WORD = {
    SINGLE: 'one object, will be ejected',
    PASS: 'extras up the capillary',
    CONTAMINATION: 'second object may be ejected',
    FAIL: 'objects past the ejection line',
    NO_OBJECT: 'nothing detected',
    NO_IMAGE: 'no image'
  };
  /* merge first: transmission with every fluorescence channel screen-blended over
     it, the same composite the cellenONE PDF report shows. */
  /* Button order only. `S.channel` below still defaults to merge, so the report
     still OPENS on the composite — this just puts TRANS first in the strip. */
  var CHANNEL_ORDER = ['trans', 'merge', 'blue', 'green', 'orange', 'red'];
  var CHANNEL_LABEL = { merge:'MERGE', trans:'TRANS', blue:'BLUE', green:'GREEN',
                        orange:'ORANGE', red:'RED' };
  var CHROM = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'];
  var CHROM_W = [249,243,198,190,182,171,159,145,138,134,135,133,114,107,102,90,83,80,59,64,47,51];
  var CHROM_TOTAL = CHROM_W.reduce(function (a, b) { return a + b; }, 0);

  var BY_ID = {};
  for (var i = 0; i < WELLS.length; i++) BY_ID[WELLS[i].id] = WELLS[i];

  /* ---------- state, persisted on every mutation ---------- */
  var S = { sel: null, sub: null, query: '', channel: 'merge', drawer: false,
            imgMore: (typeof window !== 'undefined' && window.innerWidth >= 1680),
            decisions: {}, seen: {}, decided: {} };

  function defaultDecision(w) { return w.status === 'PASS' ? 'PASS' : w.status === 'WARN' ? 'REVIEW' : 'EXCLUDE'; }

  function load() {
    for (var i = 0; i < WELLS.length; i++) {
      var w = WELLS[i];
      S.decisions[w.id] = { decision: defaultDecision(w), reasons: w.status === 'FAIL' ? ['low_read_count'] : [], notes: '' };
    }
    var raw = null;
    try { raw = window.localStorage.getItem(KEY); } catch (e) { raw = null; }
    if (!raw) return;
    var saved;
    try { saved = JSON.parse(raw); } catch (e) { return; }
    if (!saved) return;
    for (var id in saved.decisions || {}) if (BY_ID[id]) S.decisions[id] = saved.decisions[id];
    S.seen = saved.seen || {};
    S.decided = saved.decided || {};
  }

  function save() {
    try {
      window.localStorage.setItem(KEY, JSON.stringify({ decisions: S.decisions, seen: S.seen, decided: S.decided, at: Date.now() }));
    } catch (e) { /* storage disabled: the report still works, it just will not restore */ }
  }

  /* ---------- helpers ---------- */
  function el(id) { return document.getElementById(id); }
  function fmt(n) { return n == null ? '\u2013' : Number(n).toLocaleString('en-US'); }
  function short(n) {
    if (n == null) return '\u2013';
    if (n >= 1e6) return (n / 1e6).toFixed(2) + 'M';
    if (n >= 1e3) return (n / 1e3).toFixed(n >= 1e5 ? 0 : 1) + 'k';
    return String(n);
  }
  function pct(n, d) { return (100 * n / (d || 1)).toFixed(1) + '%'; }
  function gateClass(w) { return 'f-' + w.status; }

  function queue() {
    var list = WELLS.slice();
    if (!IS_REVIEW) list = list.filter(function (w) { return w.included; });
    if (S.sub) list = list.filter(function (w) { return w.subplate === S.sub; });
    var q = S.query.trim().toLowerCase();
    if (q) return list.filter(function (w) { return (w.id + ' ' + (w.pos || '') + ' ' + w.well).toLowerCase().indexOf(q) !== -1; });
    /* Every well, in plate order. There was once a "needs you" subset that combined the
       read gate with the droplet call, but no metric has earned that authority yet —
       while this is exploratory the report shows the plate as it is and lets the
       reviewer decide what deserves attention. */
    return list;
  }

  /* ---------- mutations ---------- */
  function select(id) { S.sel = id; S.seen[id] = true; save(); renderAll(); }
  function step(d) {
    var list = queue(); if (!list.length) return;
    var idx = -1;
    for (var i = 0; i < list.length; i++) if (list[i].id === S.sel) idx = i;
    idx = Math.max(0, Math.min(list.length - 1, (idx === -1 ? 0 : idx) + d));
    select(list[idx].id);
  }
  function setDecision(v) {
    if (!S.sel) return;
    S.decisions[S.sel].decision = v;
    S.decided[S.sel] = true;
    save(); flash('autosaved ' + S.sel); renderAll();
  }
  function toggleReason(r) {
    if (!S.sel) return;
    var rs = S.decisions[S.sel].reasons, at = rs.indexOf(r);
    if (at === -1) rs.push(r); else rs.splice(at, 1);
    S.decided[S.sel] = true;
    save(); flash('autosaved ' + S.sel); renderAll();
  }
  function confirmWell() {
    if (!S.sel) return;
    S.decided[S.sel] = true; S.seen[S.sel] = true;
    save(); flash('confirmed ' + S.sel); step(1);
  }
  var flashTimer = null;
  function flash(msg) {
    var f = el('flash'); if (!f) return;
    f.textContent = msg;
    if (flashTimer) window.clearTimeout(flashTimer);
    flashTimer = window.setTimeout(function () { f.textContent = ''; }, 2500);
  }

  /* ---------- CSV ---------- */
  function csv() {
    var sub = QC.size === 384;
    var lines = [sub ? 'sample_id,subplate,well,decision,reason,notes' : 'sample_id,well,decision,reason,notes'];
    for (var i = 0; i < WELLS.length; i++) {
      var w = WELLS[i], d = S.decisions[w.id], row = [w.id];
      if (sub) row.push(w.subplate);
      var notes = d.notes || '';
      row.push(w.well, d.decision, d.reasons.join(';'), notes.indexOf(',') !== -1 ? '"' + notes.replace(/"/g, '""') + '"' : notes);
      lines.push(row.join(','));
    }
    return lines.join('\n') + '\n';
  }
  function copyText(text, msg) {
    if (navigator.clipboard) navigator.clipboard.writeText(text).then(function () { flash(msg); }, function () { flash('copy blocked'); });
    else flash('copy blocked \u2014 select and copy manually');
  }
  function download() {
    var blob = new Blob([csv()], { type: 'text/csv' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'qc_decisions.csv';
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    flash('qc_decisions.csv downloaded');
  }

  /* ---------- render: stats ---------- */
  function renderStats() {
    var n = WELLS.length || 1, g = { PASS: 0, WARN: 0, FAIL: 0 };
    for (var i = 0; i < WELLS.length; i++) g[WELLS[i].status]++;
    var defs = [['pass', g.PASS], ['review', g.WARN], ['fail', g.FAIL]];
    for (var j = 0; j < defs.length; j++) {
      var box = el('stat-' + defs[j][0]); if (!box) continue;
      box.querySelector('.big').textContent = pct(defs[j][1], n);
      box.querySelector('.sub').textContent = defs[j][1] + '/' + n + ' wells';
      box.querySelector('.track i').style.width = Math.min(100, 100 * defs[j][1] / n) + '%';
    }
  }

  /* ---------- render: plate map ---------- */
  function renderMap() {
    var wrap = el('wells'); if (!wrap) return;
    var cells = wrap.children;
    for (var c = 0; c < cells.length; c++) {
      var node = cells[c], id = node.getAttribute('data-id');
      if (!id) continue;
      var w = BY_ID[id];
      var d = IS_REVIEW ? S.decisions[id].decision : (w.included ? 'PASS' : 'EXCLUDE');
      /* Colour is the decision, full stop. Opacity now only marks the subplate
         filter — it no longer doubles as a "needs attention" signal. */
      node.className = 'well d-' + d
        + (id === S.sel ? ' sel' : '')
        + (S.sub && w.subplate !== S.sub ? ' dim' : '');
      node.style.opacity = '';
      var mark = node.querySelector('.mark');
      mark.className = 'mark' + (S.decided[id] ? ' decided' : S.seen[id] ? ' seen' : '');
      mark.style.display = (S.decided[id] || S.seen[id]) ? '' : 'none';
    }
  }

  /* ---------- render: inspector ---------- */
  function renderInspector() {
    var w = S.sel ? BY_ID[S.sel] : null;
    if (!w) return;
    var d = S.decisions[w.id];

    el('well-id').textContent = w.id;
    var dc = el('chip-decision');
    dc.className = 'chip c-' + (IS_REVIEW ? d.decision : (w.included ? 'PASS' : 'EXCLUDE'));
    dc.textContent = IS_REVIEW ? d.decision : (w.included ? 'KEPT' : 'DROPPED');
    var sc = el('chip-source');
    sc.className = 'chip ' + (S.decided[w.id] ? 'c-manual' : 'c-auto');
    sc.textContent = S.decided[w.id] ? 'MANUAL' : 'AUTO';
    el('chip-pos').textContent = w.pos || w.well;

    var btns = el('dec-grid').children;
    for (var i = 0; i < btns.length; i++) {
      btns[i].className = 'dec-btn' + (btns[i].getAttribute('data-dec') === d.decision ? ' on' : '');
    }
    el('confirm').textContent = '\u2713';
    el('caret').textContent = S.drawer ? '\u25b2' : '\u25bc';
    el('caret').className = 'caret' + (S.drawer ? ' on' : '');
    el('drawer').hidden = !S.drawer;

    var chips = el('drawer').querySelectorAll('.rchip');
    for (var k = 0; k < chips.length; k++) {
      chips[k].className = 'rchip' + (d.reasons.indexOf(chips[k].getAttribute('data-reason')) !== -1 ? ' on' : '');
    }
    el('notes').value = d.notes || '';

    /* usable reads, fixed 0-500k scale */
    el('usable-val').textContent = fmt(w.reads);
    var fill = el('usable-fill');
    fill.className = 'fill ' + gateClass(w);
    fill.style.width = Math.min(100, 100 * w.reads / 500000).toFixed(1) + '%';

    /* The eight diagnostics, in pipeline order: the read funnel first
       (demux -> unique -> dup -> umi kept -> mapped), then the library-quality
       measures (gc, dimer, avg q). `reads` is deliberately absent — it is the
       headline `usable-fill` figure above and would only repeat it here. */
    var m = w.metrics || {};
    var rows = [
      ['Demux', short(m.demux_reads), ''], ['Unique', short(m.unique_reads), ''],
      ['Dup', num(m.dedup_rate), '%'], ['UMI kept', num(m.umi_retention), '%'],
      ['Mapped', num(m.mapping_rate), '%'], ['GC', num(m.gc_content), '%'],
      ['Dimer', num(m.dimer_rate), '%'], ['Avg Q', num(m.average_quality), '']
    ];
    var grid = el('diag');
    grid.innerHTML = '';
    for (var r = 0; r < rows.length; r++) {
      var cell = document.createElement('div');
      cell.innerHTML = '<span class="k"></span><span class="v"></span><span class="u"></span>';
      cell.querySelector('.k').textContent = rows[r][0];
      cell.querySelector('.v').textContent = rows[r][2];
      cell.querySelector('.v').textContent = rows[r][1];
      cell.querySelector('.u').textContent = rows[r][2];
      grid.appendChild(cell);
    }

    /* image. Skipped entirely on a plate with no CellenONE images — applyImageLayer()
       has already removed the row, so there is nothing here to fill. */
    if (HAS_IMAGES) renderImages(w);

    /* plots */
    setPlot('profile', w.profile_src, 'profile');
    setPlot('hist', w.hist_src, 'histogram');
    /* bin size comes from the pipeline config, not a literal — this report is not always 500 kb */
    el('profile-note').textContent = w.status === 'FAIL' ? 'sparse — few reads per bin'
      : 'segmented' + (QC.bins_label ? ', ' + QC.bins_label + ' bins' : '');
  }

  function renderImages(w) {
    var plate = el('plate');
    /* Channel buttons are built from the channels this well actually has: a run may
       not have used every LED, and Red exists for only a handful of wells. Falls back
       to merge when the selected channel is missing for this well. */
    var have = [];
    for (var ci = 0; ci < CHANNEL_ORDER.length; ci++) {
      if (w.images && w.images[CHANNEL_ORDER[ci]]) have.push(CHANNEL_ORDER[ci]);
    }
    if (have.length && have.indexOf(S.channel) === -1) S.channel = have[0];
    var segHtml = '';
    for (var cj = 0; cj < have.length; cj++) {
      segHtml += '<button type="button" data-ch="' + have[cj] + '"' +
                 (have[cj] === S.channel ? ' class="on"' : '') + '>' +
                 CHANNEL_LABEL[have[cj]] + '</button>';
    }
    el('seg').innerHTML = segHtml;

    plate.className = 'plate' + (S.channel === 'trans' ? ' trans' : '');
    var src = w.images ? w.images[S.channel] : null;
    var img = el('plate-img');
    if (src) {
      img.src = src; img.hidden = false;
      /* Against a fluorescence channel, show that channel's measured intensity for the
         isolated cell — the number the reviewer wants while looking at the stain. On
         merge/trans there is no single channel to report, so the caption stays off. */
      var fi = (w.flu && w.flu[S.channel] != null) ? w.flu[S.channel] : null;
      if (fi != null) {
        el('plate-cap').hidden = false;
        el('plate-cap').textContent = CHANNEL_LABEL[S.channel] + ' intensity ' + num(fi);
      } else {
        el('plate-cap').hidden = true;
      }
    } else { img.hidden = true; el('plate-cap').hidden = false; el('plate-cap').textContent = 'no image for this well'; }

    /* Primary column: what OUR detection found in the whole frame \u2014 these are the
       four facts that explain the call. */
    var im = [['Objects', w.nObj], ['In iso window', w.nIso], ['Rightmost x', num(w.rightmost_x)], ['Rightmost \u00f8', num(w.rightmost_dia) + ' \u00b5m']];
    /* Secondary column: CellenONE's own measurements of the cell it isolated. Useful
       for explaining a call (why is nIso 0 when nObj is 3?), but never predictive of
       sequencing yield \u2014 kept behind the disclosure so it cannot crowd the primaries. */
    var more = [['Diameter', w.cell_dia == null ? null : num(w.cell_dia) + ' \u00b5m'],
                ['Elongation', num2(w.cell_elong)],
                ['Circularity', num2(w.cell_circ)],
                ['Intensity', num(w.cell_int)]];
    fillRail(el('img-metrics-body'), im);
    fillRail(el('img-metrics-more'), more);
  }
  function num(v) { return v == null ? '\u2013' : Number(v).toFixed(1); }
  function num2(v) { return v == null ? null : Number(v).toFixed(2); }
  function fillRail(rail, rows) {
    if (!rail) return;
    rail.innerHTML = '';
    for (var t = 0; t < rows.length; t++) {
      var mm = document.createElement('div');
      mm.className = 'm';
      mm.innerHTML = '<div class="micro"></div><div class="v"></div>';
      mm.querySelector('.micro').textContent = rows[t][0];
      mm.querySelector('.v').textContent = rows[t][1] == null ? '\u2013' : rows[t][1];
      rail.appendChild(mm);
    }
  }
  /* Width drives this column, not the reviewer.
       wide enough  -> the extra space goes to the second column rather than to
                       padding around the image
       too narrow   -> it folds away again, so the image never gets squeezed to
                       make room for numbers
     The caret is an override for the current width only: any resize re-asserts the
     width rule, which is what stops a manual "open" from permanently shrinking the
     image once the window comes back down. */
  var IMG_MORE_AT = 1680;
  var imgMoreTouched = false;
  function syncImgMore() {
    var fits = window.innerWidth >= IMG_MORE_AT;
    var open = imgMoreTouched ? S.imgMore : fits;
    S.imgMore = open;
    var col = el('img-metrics-more'), btn = el('img-more'), box = el('img-metrics');
    if (col) col.hidden = !open;
    if (btn) {
      btn.innerHTML = open ? '&#9664;' : '&#9654;';
      btn.setAttribute('aria-expanded', open ? 'true' : 'false');
    }
    if (box) box.className = 'img-metrics' + (open ? ' wide' : '');
  }
  function setPlot(id, src, what) {
    var img = el(id + '-img'), ph = el(id + '-placeholder');
    if (src) { img.src = src; img.hidden = false; if (ph) ph.hidden = true; }
    else { img.hidden = true; if (ph) { ph.hidden = false; ph.textContent = 'no ' + what + ' for this well'; } }
  }

  /* chromosome axis: suppress a label whose slot is under 18px */
  function renderAxis() {
    var axis = el('profile-axis'); if (!axis) return;
    var width = axis.clientWidth || 600;
    axis.innerHTML = '';
    for (var i = 0; i < CHROM.length; i++) {
      var span = document.createElement('span');
      var slot = CHROM_W[i] / CHROM_TOTAL * width;
      span.style.flex = String(CHROM_W[i]);
      span.textContent = slot >= 18 ? CHROM[i] : '';
      axis.appendChild(span);
    }
  }

  /* ---------- render: header + export ---------- */
  function renderHeader() {
    var done = 0;
    for (var i = 0; i < WELLS.length; i++) if (S.decided[WELLS[i].id]) done++;
    var p = el('progress');
    if (p) p.textContent = done + ' of ' + WELLS.length + ' done \u00b7 ' + (WELLS.length - done) + ' left';
    var counts = { PASS: 0, EXCLUDE: 0, REVIEW: 0, REPEAT: 0 }, decided = 0;
    for (var j = 0; j < WELLS.length; j++) {
      counts[S.decisions[WELLS[j].id].decision]++;
      if (S.decided[WELLS[j].id]) decided++;
    }
    var sum = el('export-summary');
    if (sum) sum.textContent = 'PASS ' + counts.PASS + ' \u00b7 EXCLUDE ' + counts.EXCLUDE + ' \u00b7 REVIEW ' + counts.REVIEW
      + ' \u00b7 REPEAT ' + counts.REPEAT + '  \u2014  ' + decided + ' of ' + WELLS.length + ' decided by you';
    var tabs = document.querySelectorAll('#subtabs .tab');
    for (var t = 0; t < tabs.length; t++) {
      var v = tabs[t].getAttribute('data-sub') || null;
      tabs[t].className = 'tab' + (v === S.sub ? ' on' : '');
    }
  }

  function renderAll() { renderHeader(); renderStats(); renderMap(); renderInspector(); renderAxis(); syncImgMore(); }

  /* ---------- wiring ---------- */
  function buildGrid() {
    var wrap = el('wells'), rowhdr = el('rowhdr'), colhdr = el('colhdr');
    if (!wrap) return;
    var rows = QC.size === 384 ? 16 : 8, cols = QC.size === 384 ? 24 : 12;
    /* spec §3.3: the 96 cell is 48x34 and carries its own label; 384 keeps the 26x24 default */
    if (QC.size !== 384) {
      document.documentElement.style.setProperty('--cell-w', '48px');
      document.documentElement.style.setProperty('--cell-h', '34px');
    }
    var byPos = {};
    for (var i = 0; i < WELLS.length; i++) byPos[WELLS[i].row + ':' + WELLS[i].col] = WELLS[i];
    wrap.style.gridTemplateColumns = 'repeat(' + cols + ',var(--cell-w))';
    colhdr.style.gridTemplateColumns = 'repeat(' + cols + ',var(--cell-w))';
    rowhdr.style.gridTemplateRows = 'repeat(' + rows + ',var(--cell-h))';
    colhdr.innerHTML = ''; rowhdr.innerHTML = ''; wrap.innerHTML = '';
    for (var c = 0; c < cols; c++) { var ch = document.createElement('div'); ch.textContent = String(c + 1); colhdr.appendChild(ch); }
    for (var r = 0; r < rows; r++) { var rh = document.createElement('div'); rh.textContent = String.fromCharCode(65 + r); rowhdr.appendChild(rh); }
    for (var rr = 0; rr < rows; rr++) {
      for (var cc = 0; cc < cols; cc++) {
        var w = byPos[rr + ':' + cc], node = document.createElement('div');
        if (!w) { node.className = 'well empty'; wrap.appendChild(node); continue; }
        node.className = 'well';
        node.setAttribute('data-id', w.id);
        node.title = (w.pos ? w.pos + ' \u00b7 ' : '') + w.id;
        node.innerHTML = (QC.size === 96 ? '<span class="lbl">' + w.well + '</span>' : '') + '<span class="mark"></span>';
        wrap.appendChild(node);
      }
    }
    wrap.addEventListener('click', function (e) {
      var t = e.target;
      while (t && t !== wrap && !t.getAttribute('data-id')) t = t.parentNode;
      if (t && t.getAttribute && t.getAttribute('data-id')) select(t.getAttribute('data-id'));
    });
    wrap.addEventListener('mouseover', function (e) {
      var t = e.target;
      while (t && t !== wrap && !t.getAttribute('data-id')) t = t.parentNode;
      if (!t || !t.getAttribute || !t.getAttribute('data-id')) return;
      var w2 = BY_ID[t.getAttribute('data-id')];
      /* The droplet call is only meaningful when this plate has images; without them
         it would read "droplet no image" on all 384 wells. */
      el('hoverline').textContent = (w2.pos ? w2.pos + ' \u00b7 ' : '') + w2.id + ' \u00b7 ' + short(w2.reads)
        + ' reads \u00b7 gate ' + w2.status + (HAS_IMAGES ? ' \u00b7 droplet ' + w2.call : '');
    });
  }

  function buildDrawer() {
    var host = el('rgroups');
    if (!host) return;
    host.innerHTML = '';
    for (var g = 0; g < REASON_GROUPS.length; g++) {
      var box = document.createElement('div');
      box.className = 'rgroup';
      var head = document.createElement('div');
      head.className = 'micro';
      head.textContent = REASON_GROUPS[g].title;
      var chips = document.createElement('div');
      chips.className = 'rchips';
      for (var r = 0; r < REASON_GROUPS[g].items.length; r++) {
        var b = document.createElement('button');
        b.type = 'button'; b.className = 'rchip';
        b.setAttribute('data-reason', REASON_GROUPS[g].items[r]);
        b.textContent = REASON_GROUPS[g].items[r];
        chips.appendChild(b);
      }
      box.appendChild(head); box.appendChild(chips); host.appendChild(box);
    }
    host.addEventListener('click', function (e) {
      var r = e.target.getAttribute && e.target.getAttribute('data-reason');
      if (r) toggleReason(r);
    });
  }

  function wire() {
    el('dec-grid').addEventListener('click', function (e) {
      var t = e.target;
      while (t && !t.getAttribute('data-dec')) t = t.parentNode;
      if (t) setDecision(t.getAttribute('data-dec'));
    });
    el('confirm').addEventListener('click', confirmWell);
    el('caret').addEventListener('click', function () { S.drawer = !S.drawer; renderAll(); });
    el('notes').addEventListener('input', function (e) {
      if (!S.sel) return;
      S.decisions[S.sel].notes = e.target.value; save();
    });
    el('seg').addEventListener('click', function (e) {
      var ch = e.target.getAttribute && e.target.getAttribute('data-ch');
      if (ch) { S.channel = ch; renderAll(); }
    });
    el('subtabs').addEventListener('click', function (e) {
      if (!e.target.getAttribute) return;
      if (e.target.className.indexOf('tab') === -1) return;
      S.sub = e.target.getAttribute('data-sub') || null; renderAll();
    });
    el('search').addEventListener('input', function (e) { S.query = e.target.value; renderAll(); });
    el('copy-cmd').addEventListener('click', function () {
      /* writes to the plate's decisions path — the browser cannot, and the reviewer's
         terminal is rarely sitting in the right directory. Quoted heredoc: no expansion. */
      var path = QC.decisions_path || 'qc_decisions.csv';
      var dir = QC.decisions_dir || '.';
      copyText("mkdir -p '" + dir + "'\n"
        + "[ -e '" + path + "' ] && echo 'WARNING: overwriting existing " + path + "'\n"
        + "cat > '" + path + "' <<'QC_DECISIONS_EOF'\n" + csv() + 'QC_DECISIONS_EOF\n'
        + "echo 'Wrote " + path + "'\n", 'save command copied');
    });
    el('copy-csv').addEventListener('click', function () { copyText(csv(), 'CSV copied'); });
    el('download-csv').addEventListener('click', download);

    document.addEventListener('keydown', function (e) {
      var tag = (e.target.tagName || '').toUpperCase();
      if (tag === 'INPUT' || tag === 'TEXTAREA') return;
      if (e.key === 'j' || e.key === 'J') { step(1); e.preventDefault(); }
      else if (e.key === 'k' || e.key === 'K') { step(-1); e.preventDefault(); }
      else if (IS_REVIEW && '1234'.indexOf(e.key) !== -1) { setDecision(DECISIONS[Number(e.key) - 1]); e.preventDefault(); }
      else if (IS_REVIEW && e.key === 'Enter') { confirmWell(); e.preventDefault(); }
      else if (e.key === 'ArrowRight') { imgMoreTouched = true; S.imgMore = true; syncImgMore(); e.preventDefault(); }
      else if (e.key === 'ArrowLeft') { imgMoreTouched = true; S.imgMore = false; syncImgMore(); e.preventDefault(); }
    });
    el('img-more').addEventListener('click', function () {
      imgMoreTouched = true; S.imgMore = !S.imgMore; syncImgMore();
    });
    window.addEventListener('resize', function () { imgMoreTouched = false; renderAxis(); syncImgMore(); });
  }

  /* cn viewer: read-only. Hide every control that writes a decision, show the
     second-pass count and the genome-wide heatmap. */
  /* A plate with no CellenONE images loses the image row outright, so the Evidence
     section is just the copy-number profile and the bin histogram. Hiding it (rather
     than leaving an empty frame) is the whole no-cell-image mode on the page side —
     the pipeline simply never builds the cellenone/ directory for such a plate. */
  function applyImageLayer() {
    if (HAS_IMAGES) return;
    var row = el('img-row'); if (row) row.hidden = true;
  }

  function applyMode() {
    if (IS_REVIEW) return;
    var off = ['search', 'dec', 'export'];
    for (var i = 0; i < off.length; i++) { var n = el(off[i]); if (n) n.hidden = true; }
    var kept = 0;
    for (var j = 0; j < WELLS.length; j++) if (WELLS[j].included) kept++;
    var head = el('hdr-cn'); if (head) head.hidden = false;
    var ic = el('included-count'); if (ic) ic.textContent = kept;
    var cc = el('cn-cells'); if (cc) cc.textContent = kept;
    if (QC.heatmap) {
      el('cn-heatmap').src = QC.heatmap;
      el('cnmap').hidden = false;
    }
  }

  function start() {
    if (!WELLS.length) return;
    load();
    buildGrid();
    buildDrawer();
    wire();
    applyImageLayer();
    applyMode();
    var first = WELLS[0];
    S.sel = first.id;
    renderAll();
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start);
  else start();
})();
