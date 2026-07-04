/* scReportDE — Interactive Marker DotPlot Panel
 * ==============================================
 * Renders an interactive DotPlot using Plotly.react().
 * Data is precomputed in R and embedded as JSON in <script id="dotplot-data">.
 *
 * Controls:
 *   - Identity layer dropdown (x-axis grouping)
 *   - Gene panel mode (auto / selected / custom)
 *   - Top N / direction / max genes
 *   - Bubble size slider
 *   - Gene search / custom input
 */

(function(){
'use strict';

// ── State ────────────────────────────────────────────────
let allData    = [];       // full dotplot_data from JSON
let markerMeta = {};       // per-layer marker info
let identityLayers = [];   // available layers
let currentLayer = null;
let currentMode  = 'auto'; // auto | selected | custom
let currentTopN  = 10;
let currentDir   = 'up';
let currentMaxGenes = 80;
let currentBubbleMax = 14;
let customGenes  = [];
let selectedGroups = null; // for 'selected' mode

let sizeMin = 3;
let sizeMax = 14;

// ── Init ─────────────────────────────────────────────────
function initDotPlot() {
  var el = document.getElementById('dotplot-json-data');
  if (!el) return;
  try {
    var raw = JSON.parse(el.textContent);
    allData = raw.data || [];
    identityLayers = raw.identity_layers || [];
    markerMeta = raw.marker_meta || {};
    if (identityLayers.length > 0) {
      currentLayer = identityLayers[0];
    }
  } catch(e) {
    showError('Failed to parse DotPlot data.');
    return;
  }
  buildUI();
  renderPlot();
}

// ── Build Controls ──────────────────────────────────────
function buildUI() {
  var ctr = document.getElementById('dotplot-controls');
  if (!ctr) return;

  var html = '<div class="dp-row">';

  // Identity layer dropdown
  html += '<div class="dp-control">';
  html += '<label>Group by</label>';
  html += '<select id="dp-identity-layer">';
  identityLayers.forEach(function(lyr) {
    html += '<option value="' + escAttr(lyr) + '">' + escHtml(lyr) + '</option>';
  });
  html += '</select></div>';

  // Gene panel mode
  html += '<div class="dp-control">';
  html += '<label>Gene panel</label>';
  html += '<select id="dp-gene-mode">';
  html += '<option value="auto">Auto markers</option>';
  html += '<option value="selected">Selected groups</option>';
  html += '<option value="custom">Custom genes</option>';
  html += '</select></div>';

  // Top N
  html += '<div class="dp-control">';
  html += '<label>Top N</label>';
  html += '<select id="dp-top-n">';
  [5, 10, 20].forEach(function(n) {
    html += '<option value="' + n + '"' + (n === 10 ? ' selected' : '') + '>' + n + '</option>';
  });
  html += '</select></div>';

  // Direction
  html += '<div class="dp-control">';
  html += '<label>Direction</label>';
  html += '<select id="dp-direction">';
  ['up','down','both'].forEach(function(d) {
    html += '<option value="' + d + '">' + d + '</option>';
  });
  html += '</select></div>';

  html += '</div><div class="dp-row">';

  // Max genes
  html += '<div class="dp-control">';
  html += '<label>Max genes</label>';
  html += '<input type="number" id="dp-max-genes" value="80" min="10" max="200" style="width:60px">';
  html += '</div>';

  // Bubble max size
  html += '<div class="dp-control">';
  html += '<label>Bubble size: <span id="dp-size-val">14</span></label>';
  html += '<input type="range" id="dp-bubble-size" min="6" max="14" value="14" step="1">';
  html += '</div>';

  // Gene search / custom input
  html += '<div class="dp-control" style="flex:2">';
  html += '<label>Custom genes <span id="dp-gene-pool-hint">(from gene pool)</span></label>';
  html += '<input type="text" id="dp-custom-genes" placeholder="e.g. COL1A1, COL1A2, DCN" style="width:100%">';
  html += '</div>';

  html += '</div>';

  // Selected groups picker (hidden by default)
  html += '<div class="dp-row" id="dp-selected-groups-row" style="display:none">';
  html += '<div class="dp-control" style="flex:1">';
  html += '<label>Select groups to show markers for</label>';
  html += '<div id="dp-group-checkboxes"></div>';
  html += '</div></div>';

  // Messages
  html += '<div id="dp-message" class="dp-message"></div>';

  ctr.innerHTML = html;

  // Bind events
  document.getElementById('dp-identity-layer').addEventListener('change', onChangeLayer);
  document.getElementById('dp-gene-mode').addEventListener('change', onChangeMode);
  document.getElementById('dp-top-n').addEventListener('change', onChangeParam);
  document.getElementById('dp-direction').addEventListener('change', onChangeParam);
  document.getElementById('dp-max-genes').addEventListener('change', onChangeParam);
  document.getElementById('dp-bubble-size').addEventListener('input', function() {
    document.getElementById('dp-size-val').textContent = this.value;
    currentBubbleMax = parseInt(this.value);
    renderPlot();
  });
  document.getElementById('dp-custom-genes').addEventListener('input', debounce(onCustomGenes, 400));

  updateSelectedGroupsUI();
}

// ── Event Handlers ──────────────────────────────────────
function onChangeLayer() {
  currentLayer = document.getElementById('dp-identity-layer').value;
  updateSelectedGroupsUI();
  renderPlot();
}

function onChangeMode() {
  currentMode = document.getElementById('dp-gene-mode').value;
  updateSelectedGroupsUI();
  renderPlot();
}

function onChangeParam() {
  currentTopN  = parseInt(document.getElementById('dp-top-n').value) || 10;
  currentDir   = document.getElementById('dp-direction').value;
  currentMaxGenes = parseInt(document.getElementById('dp-max-genes').value) || 80;
  renderPlot();
}

function onCustomGenes() {
  var raw = document.getElementById('dp-custom-genes').value;
  customGenes = raw.split(/[,;\n]+/).map(function(g) { return g.trim(); }).filter(Boolean);
  if (currentMode === 'custom') renderPlot();
}

function updateSelectedGroupsUI() {
  var row = document.getElementById('dp-selected-groups-row');
  if (row) row.style.display = (currentMode === 'selected') ? '' : 'none';
  if (currentMode !== 'selected') return;

  // Get identity values for current layer
  var groups = [];
  allData.forEach(function(d) {
    if (d.identity_layer === currentLayer && groups.indexOf(d.identity_value) === -1) {
      groups.push(d.identity_value);
    }
  });
  groups.sort(naturalCompare);

  var ctr = document.getElementById('dp-group-checkboxes');
  if (!ctr) return;
  ctr.innerHTML = groups.map(function(g) {
    return '<label class="dp-checkbox"><input type="checkbox" value="' + escAttr(g) +
           '" checked>' + escHtml(g) + '</label>';
  }).join(' ');
  ctr.querySelectorAll('input').forEach(function(cb) {
    cb.addEventListener('change', function() {
      selectedGroups = getSelectedGroups();
      renderPlot();
    });
  });
  selectedGroups = getSelectedGroups();
}

function getSelectedGroups() {
  var cbs = document.querySelectorAll('#dp-group-checkboxes input:checked');
  return Array.from(cbs).map(function(cb) { return cb.value; });
}

// ── Filter Data ─────────────────────────────────────────
function getFilteredData() {
  return allData.filter(function(d) {
    return d.identity_layer === currentLayer;
  });
}

function getGeneList(filteredData) {
  if (currentMode === 'custom') {
    return customGenes;
  }

  // Build from marker meta for current layer
  var meta = markerMeta[currentLayer];
  if (!meta) return [];

  // Filter by direction
  var candidates = meta.ranked_genes || [];
  if (currentDir === 'up') {
    candidates = meta.up_genes || [];
  } else if (currentDir === 'down') {
    candidates = meta.down_genes || [];
  }

  // If selected mode, only genes from selected groups
  if (currentMode === 'selected' && selectedGroups && selectedGroups.length > 0) {
    var selSet = {};
    selectedGroups.forEach(function(g) {
      (meta.genes_by_group[g] || []).forEach(function(gene) {
        selSet[gene] = true;
      });
    });
    candidates = Object.keys(selSet);
  }

  // Take top N per group, merge, deduplicate
  var genes = [];
  var seen = {};
  // Follow group order
  var groups = meta.group_order || [];
  groups.forEach(function(g) {
    var gGenes = (meta.genes_by_group[g] || []).filter(function(gene) {
      return candidates.indexOf(gene) !== -1;
    });
    gGenes.slice(0, currentTopN).forEach(function(gene) {
      if (!seen[gene]) {
        seen[gene] = true;
        genes.push(gene);
      }
    });
  });

  // Cap
  if (genes.length > currentMaxGenes) {
    genes = genes.slice(0, currentMaxGenes);
  }

  return genes;
}

// ── Render Plot ─────────────────────────────────────────
function renderPlot() {
  var filtered = getFilteredData();
  if (filtered.length === 0) {
    showMessage('No data for this identity layer.');
    clearPlot();
    return;
  }

  var genes = getGeneList(filtered);
  if (genes.length === 0) {
    showMessage('No genes are available for the current DotPlot selection. Try lowering top_n, changing identity layer, or using custom genes.');
    clearPlot();
    return;
  }

  // Build plot data
  var geneSet = {};
  genes.forEach(function(g) { geneSet[g] = true; });

  var plotRows = filtered.filter(function(d) {
    return geneSet[d.gene];
  });

  if (plotRows.length === 0) {
    showMessage('No data matches the current gene selection.');
    clearPlot();
    return;
  }

  // Check for skipped custom genes
  var notFound = [];
  if (currentMode === 'custom') {
    customGenes.forEach(function(g) {
      if (!plotRows.some(function(r) { return r.gene === g; })) {
        notFound.push(g);
      }
    });
  }

  var msg = '';
  if (notFound.length > 0) {
    msg = 'The following genes are not available in the precomputed DotPlot gene pool and were skipped: ' + notFound.join(', ');
  }
  showMessage(msg);

  // Identity values (x-axis, natural sort)
  var identityValues = [];
  plotRows.forEach(function(d) {
    if (identityValues.indexOf(d.identity_value) === -1) {
      identityValues.push(d.identity_value);
    }
  });
  identityValues.sort(naturalCompare);
  var n_x = identityValues.length;
  var n_y = genes.length;

  // Scale bubble size
  var pctExprs = plotRows.map(function(d) { return d.pct_expr; });
  var capVal = quantile(pctExprs, 0.95);
  var adaptiveMax = sizeMax;
  if (n_x > 35) adaptiveMax = Math.min(adaptiveMax, 8);
  else if (n_x > 20) adaptiveMax = Math.min(adaptiveMax, 10);
  if (n_y > 60) adaptiveMax = Math.min(adaptiveMax, 9);
  adaptiveMax = Math.max(adaptiveMax, currentBubbleMax - 6);

  var sizes = plotRows.map(function(d) {
    var capped = Math.min(d.pct_expr, capVal);
    return sizeMin + Math.sqrt(capped) * (adaptiveMax - sizeMin);
  });

  // Build hover text
  var hoverTexts = plotRows.map(function(d) {
    var parts = [
      '<b>Gene:</b> ' + escHtml(d.gene),
      '<b>Group by:</b> ' + escHtml(d.identity_layer),
      '<b>Group:</b> ' + escHtml(d.identity_value),
      '<b>Average expression:</b> ' + fmtNum(d.avg_expr, 3),
      '<b>Scaled avg expression:</b> ' + fmtNum(d.avg_expr_scaled || 0, 3),
      '<b>Percent expressed:</b> ' + fmtPct(d.pct_expr)
    ];
    if (d.marker_avg_log2FC != null) {
      parts.push('<b>avg_log2FC:</b> ' + fmtNum(d.marker_avg_log2FC, 3));
    }
    if (d.marker_p_val_adj != null && d.marker_p_val_adj > 0) {
      parts.push('<b>adjusted p-value:</b> ' + fmtSci(d.marker_p_val_adj));
    }
    return parts.join('<br>');
  });

  // Plotly trace
  var trace = {
    x: plotRows.map(function(d) { return d.identity_value; }),
    y: plotRows.map(function(d) { return d.gene; }),
    mode: 'markers',
    type: 'scatter',
    marker: {
      size: sizes,
      sizemode: 'diameter',
      color: plotRows.map(function(d) { return d.avg_expr_scaled || 0; }),
      colorscale: [[0, '#dcdde1'], [1, '#E6194B']],
      showscale: true,
      colorbar: { title: 'Scaled avg expr', len: 0.5 },
      line: { width: 0.3, color: '#ffffff' }
    },
    text: hoverTexts,
    hoverinfo: 'text',
    hovertemplate: '%{text}',
  };

  // Layout
  var layout = {
    title: { text: 'Marker DotPlot — ' + currentLayer, font: { size: 13 } },
    xaxis: {
      title: currentLayer,
      type: 'category',
      categoryorder: 'array',
      categoryarray: identityValues,
      tickangle: identityValues.length > 8 ? -45 : 0,
      gridcolor: '#ecf0f1',
      zerolinecolor: '#dfe6e9'
    },
    yaxis: {
      title: 'Gene',
      type: 'category',
      categoryorder: 'array',
      categoryarray: genes,
      gridcolor: '#ecf0f1',
      zerolinecolor: '#dfe6e9'
    },
    plot_bgcolor: '#fafbfc',
    paper_bgcolor: '#ffffff',
    margin: { t: 40, r: 60, b: 60, l: Math.max(100, Math.min(200, n_y * 1.2)) },
    height: Math.max(500, n_y * 18 + 200),
    font: { family: '-apple-system, BlinkMacSystemFont, sans-serif', color: '#2d3436' }
  };

  var config = {
    displayModeBar: true,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    responsive: true
  };

  Plotly.react('dotplot-plot', [trace], layout, config);
}

function clearPlot() {
  Plotly.purge('dotplot-plot');
}

// ── Helpers ─────────────────────────────────────────────
function showMessage(msg) {
  var el = document.getElementById('dp-message');
  if (!el) return;
  el.textContent = msg;
  el.style.display = msg ? '' : 'none';
}

function showError(msg) {
  showMessage(msg);
}

function escHtml(s) {
  if (!s) return '';
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

function escAttr(s) {
  return String(s).replace(/&/g,'&amp;').replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}

function fmtNum(x, d) {
  if (x == null || isNaN(x)) return 'N/A';
  return parseFloat(x).toFixed(d || 2);
}

function fmtPct(x) {
  if (x == null || isNaN(x)) return 'N/A';
  return (x * 100).toFixed(1) + '%';
}

function fmtSci(x) {
  if (x == null) return 'N/A';
  if (x < 1e-4) return x.toExponential(2);
  return fmtNum(x, 4);
}

function naturalCompare(a, b) {
  var na = parseFloat(a), nb = parseFloat(b);
  if (!isNaN(na) && !isNaN(nb)) return na - nb;
  var sa = String(a), sb = String(b);
  return sa.localeCompare(sb, undefined, { numeric: true, sensitivity: 'base' });
}

function quantile(arr, q) {
  var sorted = arr.slice().sort(function(a,b){return a-b;});
  var pos = (sorted.length - 1) * q;
  var base = Math.floor(pos), rest = pos - base;
  if (sorted[base + 1] !== undefined) {
    return sorted[base] + rest * (sorted[base + 1] - sorted[base]);
  }
  return sorted[base];
}

function debounce(fn, ms) {
  var timer;
  return function() {
    var ctx = this, args = arguments;
    clearTimeout(timer);
    timer = setTimeout(function() { fn.apply(ctx, args); }, ms);
  };
}

// ── Start ────────────────────────────────────────────────
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', initDotPlot);
} else {
  initDotPlot();
}

})();
