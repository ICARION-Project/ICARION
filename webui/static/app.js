const statusEl = document.getElementById('status');
const fileInput = document.getElementById('fileInput');
const downloadBtn = document.getElementById('downloadBtn');
const domainSelect = document.getElementById('domainSelect');
const addDomainBtn = document.getElementById('addDomain');
const removeDomainBtn = document.getElementById('removeDomain');
const domainInstrument = document.getElementById('domainInstrument');
const envGasSelect = document.getElementById('envGasSelect');
const envGasCustomRow = document.getElementById('envGasCustomRow');
const envGasCustom = document.getElementById('envGasCustom');
const dbSpecies = document.getElementById('dbSpecies');
const dbReaction = document.getElementById('dbReaction');
const ionCloud = document.getElementById('ionCloud');
const applyDatabasesBtn = document.getElementById('applyDatabases');
const ionList = document.getElementById('ionList');
const addIonBtn = document.getElementById('addIon');
const removeIonBtn = document.getElementById('removeIon');
const applyIonsBtn = document.getElementById('applyIons');
const speciesDbInput = document.getElementById('speciesDbInput');
const downloadSpeciesDb = document.getElementById('downloadSpeciesDb');
const speciesList = document.getElementById('speciesList');
const addSpeciesBtn = document.getElementById('addSpecies');
const removeSpeciesBtn = document.getElementById('removeSpecies');
const spId = document.getElementById('spId');
const spName = document.getElementById('spName');
const spMass = document.getElementById('spMass');
const spCharge = document.getElementById('spCharge');
const spMobility = document.getElementById('spMobility');
const spCcs = document.getElementById('spCcs');
const spPolar = document.getElementById('spPolar');
const spGeomFile = document.getElementById('spGeomFile');
const speciesGeomInput = document.getElementById('speciesGeomInput');
const clearGeomBtn = document.getElementById('clearGeom');
const geomViewer = document.getElementById('geomViewer');
const openSpeciesDbPath = document.getElementById('openSpeciesDbPath');
const openReactionDbPath = document.getElementById('openReactionDbPath');
const loadGeomPathBtn = document.getElementById('loadGeomPath');
const reactionDbInput = document.getElementById('reactionDbInput');
const downloadReactionDb = document.getElementById('downloadReactionDb');
const reactionList = document.getElementById('reactionList');
const addReactionBtn = document.getElementById('addReaction');
const removeReactionBtn = document.getElementById('removeReaction');
const rxnId = document.getElementById('rxnId');
const rxnReactant = document.getElementById('rxnReactant');
const rxnProduct = document.getElementById('rxnProduct');
const rxnRate = document.getElementById('rxnRate');
const rxnModel = document.getElementById('rxnModel');
const rxnEa = document.getElementById('rxnEa');
const rxnN = document.getElementById('rxnN');
const rxnTref = document.getElementById('rxnTref');
const rxnOrder = document.getElementById('rxnOrder');
const waveformList = document.getElementById('waveformList');
const addWaveformBtn = document.getElementById('addWaveform');
const removeWaveformBtn = document.getElementById('removeWaveform');
const saveWaveformBtn = document.getElementById('saveWaveform');
const wfId = document.getElementById('wfId');
const wfType = document.getElementById('wfType');
const wfValue = document.getElementById('wfValue');
const wfStart = document.getElementById('wfStart');
const wfEnd = document.getElementById('wfEnd');
const wfStartTime = document.getElementById('wfStartTime');
const wfEndTime = document.getElementById('wfEndTime');
const wfClamp = document.getElementById('wfClamp');
const wfOffset = document.getElementById('wfOffset');
const wfAmplitude = document.getElementById('wfAmplitude');
const wfFrequency = document.getElementById('wfFrequency');
const wfPhase = document.getElementById('wfPhase');
const wfTimes = document.getElementById('wfTimes');
const wfValues = document.getElementById('wfValues');
const wfInterp = document.getElementById('wfInterp');
const dcAxialRef = document.getElementById('dcAxialRef');
const dcRadialRef = document.getElementById('dcRadialRef');
const dcQuadRef = document.getElementById('dcQuadRef');
const rfVoltageRef = document.getElementById('rfVoltageRef');
const rfFrequencyRef = document.getElementById('rfFrequencyRef');
const acVoltageRef = document.getElementById('acVoltageRef');
const acFrequencyRef = document.getElementById('acFrequencyRef');
const plotIonCloudBtn = document.getElementById('plotIonCloud');
const plotWaveformBtn = document.getElementById('plotWaveform');
const ionPlotOrientation = document.getElementById('ionPlotOrientation');
const ionPlotCanvas = document.getElementById('ionPlot');
const waveformPlotCanvas = document.getElementById('waveformPlot');

const GAS_OPTIONS = ['He', 'N2', 'O2', 'Ar', 'Air', 'CO2', 'H2'];

if (reactionDbInput) {
  reactionDbInput.addEventListener('change', async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    try {
      const text = await file.text();
      const parsed = JSON.parse(text);
      if (!Array.isArray(parsed.reactions)) {
        setStatus('Invalid reaction DB (expected { reactions: [ ... ] }).');
        return;
      }
      reactionDb = parsed;
      selectedReactionId = null;
      renderReactionList();
      populateReactionEditor(null);
      setStatus(`Loaded reaction DB: ${file.name}`);
    } catch (err) {
      setStatus(`Failed to load reaction DB: ${err.message || err}`);
    }
  });
}

if (downloadReactionDb) {
  downloadReactionDb.addEventListener('click', () => {
    const blob = new Blob([JSON.stringify(reactionDb, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'reactions_database.json';
    a.click();
    URL.revokeObjectURL(url);
    setStatus('Downloaded reaction database.');
  });
}

if (openReactionDbPath) {
  openReactionDbPath.addEventListener('click', async () => {
    const path = dbReaction.value.trim();
    if (!path) {
      setStatus('reaction_database path is empty.');
      return;
    }
    try {
      const parsed = await loadJsonFromPath(path);
      if (!Array.isArray(parsed.reactions)) {
        setStatus('Invalid reaction DB (expected { reactions: [ ... ] }).');
        return;
      }
      reactionDb = parsed;
      selectedReactionId = null;
      renderReactionList();
      populateReactionEditor(null);
      setStatus(`Loaded reaction DB from ${path}`);
    } catch (err) {
      setStatus(`Failed to load reaction DB: ${err.message || err}`);
    }
  });
}

const INSTRUMENT_FIELD_RULES = {
  IMS: { dc: true, dcAxial: true, dcRadial: false, dcQuad: true, rf: true, ac: false, geomDefault: true, orbitrap: false },
  LQIT: { dc: true, dcAxial: true, dcRadial: true, dcQuad: true, rf: true, ac: true, geomDefault: true, orbitrap: false },
  Quadrupole: { dc: true, dcAxial: true, dcRadial: false, dcQuad: true, rf: true, ac: false, geomDefault: true, orbitrap: false },
  TOF: { dc: true, dcAxial: true, dcRadial: false, dcQuad: false, rf: false, ac: false, geomDefault: true, orbitrap: false },
  FTICR: { dc: true, dcAxial: true, dcRadial: false, dcQuad: true, rf: true, ac: false, geomDefault: true, orbitrap: false },
  Orbitrap: { dc: true, dcAxial: false, dcRadial: true, dcQuad: false, rf: false, ac: false, geomDefault: true, orbitrap: true },
  NoFixedInstrument: { dc: true, dcAxial: true, dcRadial: true, dcQuad: true, rf: false, ac: false, geomDefault: true, orbitrap: false }
};

const fields = {
  simTotalTime: 'simulation.total_time_s',
  simDt: 'simulation.dt_s',
  simIntegrator: 'simulation.integrator',
  simWrite: 'simulation.write_interval',
  simSeed: 'simulation.rng_seed',
  outFolder: 'output.folder',
  outTrajectory: 'output.trajectory_file',
  outProgress: 'output.print_progress'
};

let config = null;
let speciesDb = { species: {} };
let selectedSpeciesId = null;
let reactionDb = { reactions: [] };
let selectedReactionId = null;
let selectedWaveformId = null;
let ionChart = null;
let waveformChart = null;
let viewer = null;
let viewerAtoms = [];

async function loadJsonFromPath(path) {
  if (!path) throw new Error('Empty path');
  const url = `/api/load-json?path=${encodeURIComponent(path)}`;
  const res = await fetch(url);
  if (!res.ok) {
    const detail = await res.json().catch(() => ({}));
    throw new Error(detail.detail || res.statusText);
  }
  return res.json();
}

const DEFAULTS = {
  simulation: {
    total_time_s: 1e-3,
    dt_s: 1e-9,
    integrator: 'RK4',
    write_interval: 100,
    rng_seed: 42
  },
  output: {
    folder: './results',
    trajectory_file: 'trajectories.h5',
    print_progress: true
  },
  databases: {
    species_database: '',
    reaction_database: '',
    ion_cloud: ''
  },
  ions: [
    {
      id: 'H3O+',
      count: 100,
      position: { type: 'gaussian', center: [0, 0, 0], std: [0.001, 0.001, 0.001] },
      velocity: { type: 'thermal', temperature_K: 300.0 }
    }
  ],
  domain: {
    instrument: 'IMS',
    env: {
      pressure_Pa: 101325.0,
      temperature_K: 300.0,
      gas_species: 'He'
    },
    geometry: {
      origin_m: [0, 0, 0],
      length_m: 0.05,
      radius_m: 0.01
    },
    fields: {
      DC: { axial_V: 0.0, radial_V: 0.0, quad_V: 0.0 },
      RF: { voltage_V: 0.0, frequency_Hz: 0.0 },
      AC: { voltage_V: 0.0, frequency_Hz: 0.0 }
    }
  }
};

function buildDefaultConfig() {
  return {
    simulation: { ...DEFAULTS.simulation },
    output: { ...DEFAULTS.output },
    species_database: DEFAULTS.databases.species_database,
    reaction_database: DEFAULTS.databases.reaction_database,
    ion_cloud: DEFAULTS.databases.ion_cloud,
    ions: { species: DEFAULTS.ions.map((s) => JSON.parse(JSON.stringify(s))) },
    domains: [
      {
        name: 'domain_0',
        instrument: DEFAULTS.domain.instrument,
        geometry: { ...DEFAULTS.domain.geometry },
        env: { ...DEFAULTS.domain.env },
        fields: {
          DC: { ...DEFAULTS.domain.fields.DC },
          RF: { ...DEFAULTS.domain.fields.RF },
          AC: { ...DEFAULTS.domain.fields.AC }
        }
      }
    ]
  };
}

function setStatus(message) {
  statusEl.textContent = message;
}

function parseNumber(value) {
  if (value === '' || value === null || value === undefined) return null;
  const num = Number(value);
  return Number.isFinite(num) ? num : null;
}

function getPathTokens(path) {
  return path.split('.').filter(Boolean).map(part => (part.match(/^\d+$/) ? Number(part) : part));
}

function setByPath(obj, path, value) {
  const tokens = getPathTokens(path);
  if (!tokens.length) return;
  let cur = obj;
  for (let i = 0; i < tokens.length - 1; i++) {
    const tok = tokens[i];
    const next = tokens[i + 1];
    if (typeof tok === 'string') {
      if (!cur[tok]) cur[tok] = typeof next === 'number' ? [] : {};
      cur = cur[tok];
    } else {
      if (!Array.isArray(cur)) throw new Error('Expected array for index path');
      while (cur.length <= tok) cur.push(typeof next === 'number' ? [] : {});
      cur = cur[tok];
    }
  }
  const last = tokens[tokens.length - 1];
  if (typeof last === 'string') cur[last] = value;
  else {
    if (!Array.isArray(cur)) throw new Error('Expected array for final index');
    while (cur.length <= last) cur.push(null);
    cur[last] = value;
  }
}

function getByPath(obj, path, fallback = '') {
  const tokens = getPathTokens(path);
  let cur = obj;
  for (const tok of tokens) {
    if (cur === null || cur === undefined) return fallback;
    cur = cur[tok];
  }
  return cur ?? fallback;
}

function updateTopLevelSection() {
  if (!config) return;
  Object.entries(fields).forEach(([id, path]) => {
    const el = document.getElementById(id);
    const value = el.value;
    if (value === '') return;
    if (id === 'outProgress') {
      setByPath(config, path, value === 'true');
      return;
    }
    if (id === 'simIntegrator' || id === 'outFolder' || id === 'outTrajectory') {
      setByPath(config, path, value);
      return;
    }
    const num = parseNumber(value);
    if (num !== null) setByPath(config, path, num);
  });
  setStatus('Updated simulation/output settings.');
}

function getSelectedDomain() {
  if (!config || !Array.isArray(config.domains) || config.domains.length === 0) return null;
  const idx = Number(domainSelect.value);
  return config.domains[idx];
}

function populateDomainSelect() {
  domainSelect.innerHTML = '';
  if (!config || !Array.isArray(config.domains)) return;
  config.domains.forEach((domain, index) => {
    const option = document.createElement('option');
    const name = domain.name || `domain_${index}`;
    const inst = domain.instrument || '?';
    option.value = index;
    option.textContent = `${index}: ${name} [${inst}]`;
    domainSelect.appendChild(option);
  });
}

function loadDomainValues() {
  const domain = getSelectedDomain();
  if (!domain) return;
  domainInstrument.value = domain.instrument || 'IMS';
  applyFieldVisibility(domainInstrument.value);
  const env = domain.env || domain.environment || {};
  const geom = domain.geometry || {};
  const fields = domain.fields || {};
  const dc = fields.DC || {};
  const rf = fields.RF || {};
  const ac = fields.AC || {};

  document.getElementById('envPressure').value = env.pressure_Pa ?? 101325.0;
  document.getElementById('envTemperature').value = env.temperature_K ?? 300.0;
  const gas = env.gas_species ?? 'He';
  if (GAS_OPTIONS.includes(gas)) {
    envGasSelect.value = gas;
    envGasCustomRow.classList.add('hidden');
    envGasCustom.value = '';
  } else {
    envGasSelect.value = '__custom__';
    envGasCustomRow.classList.remove('hidden');
    envGasCustom.value = gas;
  }

  const origin = geom.origin_m || [null, null, null];
  document.getElementById('geomOriginX').value = origin[0] ?? 0;
  document.getElementById('geomOriginY').value = origin[1] ?? 0;
  document.getElementById('geomOriginZ').value = origin[2] ?? 0;
  document.getElementById('geomLength').value = geom.length_m ?? '';
  document.getElementById('geomRadius').value = geom.radius_m ?? '';
  document.getElementById('geomInnerRadius').value = geom.radius_in_m ?? '';
  document.getElementById('geomOuterRadius').value = geom.radius_out_m ?? '';
  document.getElementById('geomCharacteristic').value = geom.radius_char_m ?? '';

  const dcAxialVal = dc.axial_V ?? '';
  dcAxialRef.value = typeof dcAxialVal === 'string' ? dcAxialVal : '';
  document.getElementById('dcAxial').value = typeof dcAxialVal === 'string' ? '' : dcAxialVal;
  const dcRadialVal = dc.radial_V ?? '';
  dcRadialRef.value = typeof dcRadialVal === 'string' ? dcRadialVal : '';
  document.getElementById('dcRadial').value = typeof dcRadialVal === 'string' ? '' : dcRadialVal;
  const dcQuadVal = dc.quad_V ?? '';
  dcQuadRef.value = typeof dcQuadVal === 'string' ? dcQuadVal : '';
  document.getElementById('dcQuad').value = typeof dcQuadVal === 'string' ? '' : dcQuadVal;

  const rfVoltageVal = rf.voltage_V ?? '';
  rfVoltageRef.value = typeof rfVoltageVal === 'string' ? rfVoltageVal : '';
  document.getElementById('rfVoltage').value = typeof rfVoltageVal === 'string' ? '' : rfVoltageVal;
  const rfFrequencyVal = rf.frequency_Hz ?? '';
  rfFrequencyRef.value = typeof rfFrequencyVal === 'string' ? rfFrequencyVal : '';
  document.getElementById('rfFrequency').value = typeof rfFrequencyVal === 'string' ? '' : rfFrequencyVal;

  const acVoltageVal = ac.voltage_V ?? '';
  acVoltageRef.value = typeof acVoltageVal === 'string' ? acVoltageVal : '';
  document.getElementById('acVoltage').value = typeof acVoltageVal === 'string' ? '' : acVoltageVal;
  const acFrequencyVal = ac.frequency_Hz ?? '';
  acFrequencyRef.value = typeof acFrequencyVal === 'string' ? acFrequencyVal : '';
  document.getElementById('acFrequency').value = typeof acFrequencyVal === 'string' ? '' : acFrequencyVal;
}

function applyFieldVisibility(instrument) {
  const container = document.getElementById('domainColumns');
  if (!container) return;
  const rules = INSTRUMENT_FIELD_RULES[instrument] || {
    dc: true, dcAxial: true, dcRadial: true, dcQuad: true, rf: true, ac: false, geomDefault: true, orbitrap: false
  };
  if (rules.ac) {
    container.classList.remove('ac-hidden');
  } else {
    container.classList.add('ac-hidden');
  }
  if (rules.rf) {
    container.classList.remove('rf-hidden');
  } else {
    container.classList.add('rf-hidden');
  }
  if (rules.dc) {
    container.classList.remove('dc-hidden');
  } else {
    container.classList.add('dc-hidden');
  }

  document.querySelectorAll('#domainColumns .dc-axial-only').forEach((el) => {
    el.classList.toggle('hidden', !rules.dc || !rules.dcAxial);
  });
  document.querySelectorAll('#domainColumns .dc-radial-only').forEach((el) => {
    el.classList.toggle('hidden', !rules.dc || !rules.dcRadial);
  });
  document.querySelectorAll('#domainColumns .dc-quad-only').forEach((el) => {
    el.classList.toggle('hidden', !rules.dc || !rules.dcQuad);
  });

  if (rules.orbitrap) {
    container.classList.remove('orbitrap-hidden');
    container.classList.remove('geom-default-hidden');
    container.classList.add('radius-hidden');
  } else {
    container.classList.add('orbitrap-hidden');
    container.classList.remove('geom-default-hidden');
    container.classList.remove('radius-hidden');
  }
}

function applyEnv() {
  const domain = getSelectedDomain();
  if (!domain) return;
  const env = domain.env || domain.environment || {};
  const pressure = parseNumber(document.getElementById('envPressure').value);
  const temperature = parseNumber(document.getElementById('envTemperature').value);
  let gas = envGasSelect.value;
  if (gas === '__custom__') {
    gas = envGasCustom.value.trim();
  }
  if (pressure !== null) env.pressure_Pa = pressure;
  if (temperature !== null) env.temperature_K = temperature;
  if (gas !== '') env.gas_species = gas;
  if (domain.env) domain.env = env; else domain.environment = env;
  setStatus('Environment updated.');
}

function applyGeometry() {
  const domain = getSelectedDomain();
  if (!domain) return;
  const geom = domain.geometry || {};
  const rules = INSTRUMENT_FIELD_RULES[domain.instrument] || { geomDefault: true, orbitrap: false };
  const ox = parseNumber(document.getElementById('geomOriginX').value);
  const oy = parseNumber(document.getElementById('geomOriginY').value);
  const oz = parseNumber(document.getElementById('geomOriginZ').value);
  const length = parseNumber(document.getElementById('geomLength').value);
  geom.origin_m = [ox ?? 0, oy ?? 0, oz ?? 0];
  if (length !== null) geom.length_m = length;

  if (rules.orbitrap) {
    const inner = parseNumber(document.getElementById('geomInnerRadius').value);
    const outer = parseNumber(document.getElementById('geomOuterRadius').value);
    const characteristic = parseNumber(document.getElementById('geomCharacteristic').value);
    if (inner !== null) geom.radius_in_m = inner;
    if (outer !== null) geom.radius_out_m = outer;
    if (characteristic !== null) geom.radius_char_m = characteristic;
    if ('radius_m' in geom) delete geom.radius_m;
  } else {
    const radius = parseNumber(document.getElementById('geomRadius').value);
    if (radius !== null) geom.radius_m = radius;
    if ('radius_in_m' in geom) delete geom.radius_in_m;
    if ('radius_out_m' in geom) delete geom.radius_out_m;
    if ('radius_char_m' in geom) delete geom.radius_char_m;
  }
  domain.geometry = geom;
  setStatus('Geometry updated.');
}

function applyFields() {
  const domain = getSelectedDomain();
  if (!domain) return;
  domain.instrument = domainInstrument.value;
  applyFieldVisibility(domain.instrument);
  domain.fields = domain.fields || {};
  domain.fields.DC = domain.fields.DC || {};
  domain.fields.RF = domain.fields.RF || {};
  domain.fields.AC = domain.fields.AC || {};

  const dcAxial = parseNumber(document.getElementById('dcAxial').value);
  const dcRadial = parseNumber(document.getElementById('dcRadial').value);
  const dcQuad = parseNumber(document.getElementById('dcQuad').value);

  const rfVoltage = parseNumber(document.getElementById('rfVoltage').value);
  const rfFrequency = parseNumber(document.getElementById('rfFrequency').value);

  const acVoltage = parseNumber(document.getElementById('acVoltage').value);
  const acFrequency = parseNumber(document.getElementById('acFrequency').value);
  const refs = {
    dcAxial: dcAxialRef.value.trim(),
    dcRadial: dcRadialRef.value.trim(),
    dcQuad: dcQuadRef.value.trim(),
    rfVoltage: rfVoltageRef.value.trim(),
    rfFrequency: rfFrequencyRef.value.trim(),
    acVoltage: acVoltageRef.value.trim(),
    acFrequency: acFrequencyRef.value.trim()
  };
  const rules = INSTRUMENT_FIELD_RULES[domain.instrument] || {
    dc: true, dcAxial: true, dcRadial: true, dcQuad: true, rf: true, ac: false
  };
  if (rules.dc) {
    if (rules.dcAxial) {
      if (refs.dcAxial) domain.fields.DC.axial_V = refs.dcAxial;
      else if (dcAxial !== null) domain.fields.DC.axial_V = dcAxial;
    } else {
      delete domain.fields.DC.axial_V;
    }
    if (rules.dcRadial) {
      if (refs.dcRadial) domain.fields.DC.radial_V = refs.dcRadial;
      else if (dcRadial !== null) domain.fields.DC.radial_V = dcRadial;
    } else {
      delete domain.fields.DC.radial_V;
    }
    if (rules.dcQuad) {
      if (refs.dcQuad) domain.fields.DC.quad_V = refs.dcQuad;
      else if (dcQuad !== null) domain.fields.DC.quad_V = dcQuad;
    } else {
      delete domain.fields.DC.quad_V;
    }
    if (Object.keys(domain.fields.DC).length === 0) delete domain.fields.DC;
  } else if (domain.fields.DC) {
    delete domain.fields.DC;
  }
  if (rules.rf) {
    if (refs.rfVoltage) domain.fields.RF.voltage_V = refs.rfVoltage;
    else if (rfVoltage !== null) domain.fields.RF.voltage_V = rfVoltage;
    if (refs.rfFrequency) domain.fields.RF.frequency_Hz = refs.rfFrequency;
    else if (rfFrequency !== null) domain.fields.RF.frequency_Hz = rfFrequency;
  } else if (domain.fields.RF) {
    delete domain.fields.RF;
  }
  if (rules.ac) {
    if (refs.acVoltage) domain.fields.AC.voltage_V = refs.acVoltage;
    else if (acVoltage !== null) domain.fields.AC.voltage_V = acVoltage;
    if (refs.acFrequency) domain.fields.AC.frequency_Hz = refs.acFrequency;
    else if (acFrequency !== null) domain.fields.AC.frequency_Hz = acFrequency;
  } else if (domain.fields.AC) {
    delete domain.fields.AC;
  }

  setStatus('Fields updated.');
}

function addDomain() {
  if (!config) {
    config = buildDefaultConfig();
  }
  if (!Array.isArray(config.domains)) config.domains = [];
  const newDomain = {
    name: `domain_${config.domains.length}`,
    instrument: 'IMS',
    geometry: {
      origin_m: [0, 0, 0],
      length_m: 0.05,
      radius_m: 0.01
    },
    env: {
      pressure_Pa: 101325.0,
      temperature_K: 300.0,
      gas_species: 'He'
    },
    fields: {
      DC: { axial_V: 0.0, radial_V: 0.0, quad_V: 0.0 },
      RF: { voltage_V: 0.0, frequency_Hz: 0.0 },
      AC: { voltage_V: 0.0, frequency_Hz: 0.0 }
    }
  };
  config.domains.push(newDomain);
  populateDomainSelect();
  domainSelect.value = String(config.domains.length - 1);
  loadDomainValues();
  setStatus('Added new domain.');
}

function removeDomain() {
  if (!config || !Array.isArray(config.domains) || config.domains.length === 0) return;
  const idx = Number(domainSelect.value);
  config.domains.splice(idx, 1);
  populateDomainSelect();
  if (config.domains.length > 0) {
    domainSelect.value = '0';
    loadDomainValues();
  }
  setStatus('Removed domain.');
}


function populateTopLevel() {
  if (!config) return;
  document.getElementById('simTotalTime').value = getByPath(config, fields.simTotalTime, 1e-3);
  document.getElementById('simDt').value = getByPath(config, fields.simDt, 1e-9);
  document.getElementById('simIntegrator').value = getByPath(config, fields.simIntegrator, 'RK4');
  document.getElementById('simWrite').value = getByPath(config, fields.simWrite, 100);
  document.getElementById('simSeed').value = getByPath(config, fields.simSeed, 42);

  document.getElementById('outFolder').value = getByPath(config, fields.outFolder, './results');
  document.getElementById('outTrajectory').value = getByPath(config, fields.outTrajectory, 'trajectories.h5');
  const prog = getByPath(config, fields.outProgress, true);
  document.getElementById('outProgress').value = prog ? 'true' : 'false';

  dbSpecies.value = config.species_database || '';
  dbReaction.value = config.reaction_database || '';
  ionCloud.value = config.ion_cloud || '';
}

function applyDatabases(silent = false) {
  if (!config) return;
  const speciesPath = dbSpecies.value.trim();
  const reactionPath = dbReaction.value.trim();
  const ionCloudPath = ionCloud.value.trim();
  if (speciesPath !== '') config.species_database = speciesPath; else delete config.species_database;
  if (reactionPath !== '') config.reaction_database = reactionPath; else delete config.reaction_database;
  if (ionCloudPath !== '') config.ion_cloud = ionCloudPath; else delete config.ion_cloud;
  if (!silent) setStatus('Database paths updated.');
}

function renderIonList() {
  if (!config) return;
  if (!config.ions) config.ions = { species: [] };
  if (!Array.isArray(config.ions.species)) config.ions.species = [];
  ionList.innerHTML = '';

  config.ions.species.forEach((ion, index) => {
    const card = document.createElement('div');
    card.className = 'ion-card';
    card.dataset.index = String(index);
    card.innerHTML = `
      <h4>Ion ${index + 1}</h4>
      <div class="form-grid">
        <label>id <input data-key="id" type="text" value="${ion.id ?? ''}" /></label>
        <label>count <input data-key="count" type="number" step="1" value="${ion.count ?? 0}" /></label>
        <label>mass_Da <input data-key="mass_Da" type="number" step="any" value="${ion.mass_Da ?? ''}" /></label>
        <label>charge <input data-key="charge" type="number" step="1" value="${ion.charge ?? ''}" /></label>
        <label>position.type
          <select data-key="pos_type">
            <option value="gaussian" ${ion.position?.type === 'gaussian' ? 'selected' : ''}>gaussian</option>
            <option value="point" ${ion.position?.type === 'point' ? 'selected' : ''}>point</option>
            <option value="uniform_box" ${ion.position?.type === 'uniform_box' ? 'selected' : ''}>uniform_box</option>
            <option value="uniform_sphere" ${ion.position?.type === 'uniform_sphere' ? 'selected' : ''}>uniform_sphere</option>
            <option value="uniform_cylinder" ${ion.position?.type === 'uniform_cylinder' ? 'selected' : ''}>uniform_cylinder</option>
          </select>
        </label>
        <label data-center>pos.center x <input data-key="pos_cx" type="number" step="any" value="${ion.position?.center?.[0] ?? 0}" /></label>
        <label data-center>pos.center y <input data-key="pos_cy" type="number" step="any" value="${ion.position?.center?.[1] ?? 0}" /></label>
        <label data-center>pos.center z <input data-key="pos_cz" type="number" step="any" value="${ion.position?.center?.[2] ?? 0}" /></label>
        <label data-std>pos.std x <input data-key="pos_sx" type="number" step="any" value="${ion.position?.std?.[0] ?? 0}" /></label>
        <label data-std>pos.std y <input data-key="pos_sy" type="number" step="any" value="${ion.position?.std?.[1] ?? 0}" /></label>
        <label data-std>pos.std z <input data-key="pos_sz" type="number" step="any" value="${ion.position?.std?.[2] ?? 0}" /></label>
        <label data-box>pos.min x <input data-key="pos_minx" type="number" step="any" value="${ion.position?.min?.[0] ?? -0.001}" /></label>
        <label data-box>pos.min y <input data-key="pos_miny" type="number" step="any" value="${ion.position?.min?.[1] ?? -0.001}" /></label>
        <label data-box>pos.min z <input data-key="pos_minz" type="number" step="any" value="${ion.position?.min?.[2] ?? -0.001}" /></label>
        <label data-box>pos.max x <input data-key="pos_maxx" type="number" step="any" value="${ion.position?.max?.[0] ?? 0.001}" /></label>
        <label data-box>pos.max y <input data-key="pos_maxy" type="number" step="any" value="${ion.position?.max?.[1] ?? 0.001}" /></label>
        <label data-box>pos.max z <input data-key="pos_maxz" type="number" step="any" value="${ion.position?.max?.[2] ?? 0.001}" /></label>
        <label data-sphere>pos.radius <input data-key="pos_radius" type="number" step="any" value="${ion.position?.radius ?? 0.001}" /></label>
        <label data-cyl>pos.radius <input data-key="pos_cyl_radius" type="number" step="any" value="${ion.position?.radius ?? 0.001}" /></label>
        <label data-cyl>pos.length <input data-key="pos_cyl_length" type="number" step="any" value="${ion.position?.length ?? 0.01}" /></label>
        <label>velocity.type
          <select data-key="vel_type">
            <option value="thermal" ${ion.velocity?.type === 'thermal' ? 'selected' : ''}>thermal</option>
            <option value="kinetic" ${ion.velocity?.type === 'kinetic' ? 'selected' : ''}>kinetic</option>
          </select>
        </label>
        <label data-thermal>temperature_K <input data-key="vel_temp" type="number" step="any" value="${ion.velocity?.temperature_K ?? ''}" /></label>
        <label data-kinetic>energy_eV <input data-key="vel_energy" type="number" step="any" value="${ion.velocity?.energy_eV ?? ''}" /></label>
        <label data-kinetic>direction x <input data-key="vel_dx" type="number" step="any" value="${ion.velocity?.direction?.[0] ?? 0}" /></label>
        <label data-kinetic>direction y <input data-key="vel_dy" type="number" step="any" value="${ion.velocity?.direction?.[1] ?? 0}" /></label>
        <label data-kinetic>direction z <input data-key="vel_dz" type="number" step="any" value="${ion.velocity?.direction?.[2] ?? 0}" /></label>
        <label data-kinetic>spread_angle_deg <input data-key="vel_spread" type="number" step="any" value="${ion.velocity?.spread_angle_deg ?? ''}" /></label>
      </div>
    `;

    card.addEventListener('click', () => {
      document.querySelectorAll('.ion-card').forEach((c) => c.classList.remove('selected'));
      card.classList.add('selected');
    });

    const posType = card.querySelector('[data-key="pos_type"]');
    const velType = card.querySelector('[data-key="vel_type"]');
    const syncVisibility = () => {
      const showStd = posType.value === 'gaussian';
      const showCenter = posType.value !== 'uniform_box';
      const showBox = posType.value === 'uniform_box';
      const showSphere = posType.value === 'uniform_sphere';
      const showCyl = posType.value === 'uniform_cylinder';
      card.querySelectorAll('[data-center]').forEach((el) => (el.style.display = showCenter ? '' : 'none'));
      card.querySelectorAll('[data-std]').forEach((el) => (el.style.display = showStd ? '' : 'none'));
      card.querySelectorAll('[data-box]').forEach((el) => (el.style.display = showBox ? '' : 'none'));
      card.querySelectorAll('[data-sphere]').forEach((el) => (el.style.display = showSphere ? '' : 'none'));
      card.querySelectorAll('[data-cyl]').forEach((el) => (el.style.display = showCyl ? '' : 'none'));
      const showThermal = velType.value === 'thermal';
      card.querySelectorAll('[data-thermal]').forEach((el) => (el.style.display = showThermal ? '' : 'none'));
      card.querySelectorAll('[data-kinetic]').forEach((el) => (el.style.display = showThermal ? 'none' : ''));
    };
    posType.addEventListener('change', syncVisibility);
    velType.addEventListener('change', syncVisibility);
    syncVisibility();

    ionList.appendChild(card);
  });
}

function collectIonsFromUI() {
  const cards = ionList.querySelectorAll('.ion-card');
  const species = [];
  cards.forEach((card) => {
    const getVal = (key) => card.querySelector(`[data-key="${key}"]`)?.value ?? '';
    const num = (v) => (v === '' ? null : Number(v));
    const ion = {
      id: getVal('id'),
      count: Number(getVal('count')) || 0
    };
    const mass = num(getVal('mass_Da'));
    const charge = num(getVal('charge'));
    if (mass !== null) ion.mass_Da = mass;
    if (charge !== null) ion.charge = charge;

    const posType = getVal('pos_type');
    if (posType === 'uniform_box') {
      ion.position = {
        type: 'uniform_box',
        min: [num(getVal('pos_minx')) ?? -0.001, num(getVal('pos_miny')) ?? -0.001, num(getVal('pos_minz')) ?? -0.001],
        max: [num(getVal('pos_maxx')) ?? 0.001, num(getVal('pos_maxy')) ?? 0.001, num(getVal('pos_maxz')) ?? 0.001]
      };
    } else if (posType === 'uniform_sphere') {
      ion.position = {
        type: 'uniform_sphere',
        center: [num(getVal('pos_cx')) ?? 0, num(getVal('pos_cy')) ?? 0, num(getVal('pos_cz')) ?? 0],
        radius: num(getVal('pos_radius')) ?? 0.001
      };
    } else if (posType === 'uniform_cylinder') {
      ion.position = {
        type: 'uniform_cylinder',
        center: [num(getVal('pos_cx')) ?? 0, num(getVal('pos_cy')) ?? 0, num(getVal('pos_cz')) ?? 0],
        radius: num(getVal('pos_cyl_radius')) ?? 0.001,
        length: num(getVal('pos_cyl_length')) ?? 0.01
      };
    } else {
      ion.position = {
        type: posType,
        center: [num(getVal('pos_cx')) ?? 0, num(getVal('pos_cy')) ?? 0, num(getVal('pos_cz')) ?? 0]
      };
      if (posType === 'gaussian') {
        ion.position.std = [num(getVal('pos_sx')) ?? 0, num(getVal('pos_sy')) ?? 0, num(getVal('pos_sz')) ?? 0];
      }
    }

    const velType = getVal('vel_type');
    ion.velocity = { type: velType };
    if (velType === 'thermal') {
      const temp = num(getVal('vel_temp'));
      if (temp !== null) ion.velocity.temperature_K = temp;
    } else {
      const energy = num(getVal('vel_energy'));
      if (energy !== null) ion.velocity.energy_eV = energy;
      ion.velocity.direction = [num(getVal('vel_dx')) ?? 0, num(getVal('vel_dy')) ?? 0, num(getVal('vel_dz')) ?? 0];
      const spread = num(getVal('vel_spread'));
      if (spread !== null) ion.velocity.spread_angle_deg = spread;
    }

    species.push(ion);
  });
  return species;
}

function applyIons() {
  if (!config) return;
  config.ions = { species: collectIonsFromUI() };
  setStatus('Ion species updated.');
}

function renderSpeciesList() {
  speciesList.innerHTML = '';
  const entries = Object.keys(speciesDb.species || {}).sort();
  entries.forEach((id) => {
    const item = document.createElement('div');
    item.className = 'species-item' + (id === selectedSpeciesId ? ' selected' : '');
    item.textContent = id;
    item.addEventListener('click', () => {
      selectedSpeciesId = id;
      renderSpeciesList();
      populateSpeciesEditor(speciesDb.species[id]);
    });
    speciesList.appendChild(item);
  });
}

function populateSpeciesEditor(spec) {
  spId.value = spec?.id ?? selectedSpeciesId ?? '';
  spName.value = spec?.name ?? '';
  spMass.value = spec?.mass_amu ?? '';
  spCharge.value = spec?.charge ?? '';
  spMobility.value = spec?.mobility_cm2Vs ?? '';
  spCcs.value = spec?.CCS_A2 ?? '';
  spPolar.value = spec?.polarizability_A3 ?? '';
  spGeomFile.value = spec?.geometry_file ?? '';
}

function addOrUpdateSpecies() {
  const id = spId.value.trim();
  if (!id) {
    setStatus('Species id is required.');
    return;
  }
  const spec = {
    id,
    name: spName.value.trim() || undefined,
    mass_amu: spMass.value === '' ? undefined : Number(spMass.value),
    charge: spCharge.value === '' ? undefined : Number(spCharge.value),
    mobility_cm2Vs: spMobility.value === '' ? undefined : Number(spMobility.value),
    CCS_A2: spCcs.value === '' ? undefined : Number(spCcs.value),
    polarizability_A3: spPolar.value === '' ? undefined : Number(spPolar.value),
    geometry_file: spGeomFile.value.trim() || undefined
  };
  Object.keys(spec).forEach((k) => spec[k] === undefined && delete spec[k]);
  if (!speciesDb.species) speciesDb.species = {};
  speciesDb.species[id] = spec;
  selectedSpeciesId = id;
  renderSpeciesList();
  setStatus('Species entry saved.');
}

function renderReactionList() {
  reactionList.innerHTML = '';
  (reactionDb.reactions || []).forEach((rxn) => {
    const item = document.createElement('div');
    const label = `${rxn.id || '(no id)'}: ${rxn.reactant || '?'} → ${rxn.product || '?'}`;
    item.className = 'species-item' + (rxn.id === selectedReactionId ? ' selected' : '');
    item.textContent = label;
    item.addEventListener('click', () => {
      selectedReactionId = rxn.id || null;
      renderReactionList();
      populateReactionEditor(rxn);
    });
    reactionList.appendChild(item);
  });
}

function renderWaveformList() {
  waveformList.innerHTML = '';
  if (!config) return;
  if (!config.waveforms) config.waveforms = {};
  Object.keys(config.waveforms).sort().forEach((id) => {
    const item = document.createElement('div');
    item.className = 'species-item' + (id === selectedWaveformId ? ' selected' : '');
    item.textContent = id;
    item.addEventListener('click', () => {
      selectedWaveformId = id;
      renderWaveformList();
      populateWaveformEditor(config.waveforms[id], id);
    });
    waveformList.appendChild(item);
  });
}

function setWaveformVisibility(type) {
  const setHidden = (selector, hidden) => {
    document.querySelectorAll(selector).forEach((el) => {
      el.classList.toggle('wf-hidden', hidden);
    });
  };
  setHidden('.wf-const', type !== 'constant');
  setHidden('.wf-linear', type !== 'linear');
  setHidden('.wf-sine', type !== 'sinusoidal');
  setHidden('.wf-arb', type !== 'arbitrary');
}

function populateWaveformEditor(wf, idOverride) {
  const wfIdVal = idOverride || '';
  wfId.value = wfIdVal;
  const type = wf?.type || 'constant';
  wfType.value = type;
  setWaveformVisibility(type);
  wfValue.value = wf?.value ?? '';
  wfStart.value = wf?.start ?? '';
  wfEnd.value = wf?.end ?? '';
  wfStartTime.value = wf?.start_time_s ?? '';
  wfEndTime.value = wf?.end_time_s ?? '';
  wfClamp.value = wf?.clamp === false ? 'false' : 'true';
  wfOffset.value = wf?.offset ?? '';
  wfAmplitude.value = wf?.amplitude ?? '';
  wfFrequency.value = wf?.frequency_Hz ?? '';
  wfPhase.value = wf?.phase_rad ?? '';
  wfTimes.value = Array.isArray(wf?.times) ? wf.times.join(',') : '';
  wfValues.value = Array.isArray(wf?.values) ? wf.values.join(',') : '';
  wfInterp.value = wf?.interpolation || 'linear';
}

function addOrUpdateWaveform() {
  if (!config) return;
  const id = wfId.value.trim();
  if (!id) {
    setStatus('Waveform id is required.');
    return;
  }
  const type = wfType.value;
  const wf = { type };
  if (type === 'constant') {
    wf.value = wfValue.value === '' ? 0 : Number(wfValue.value);
  } else if (type === 'linear') {
    wf.start = Number(wfStart.value);
    wf.end = Number(wfEnd.value);
    wf.start_time_s = wfStartTime.value === '' ? 0 : Number(wfStartTime.value);
    wf.end_time_s = Number(wfEndTime.value);
    wf.clamp = wfClamp.value !== 'false';
  } else if (type === 'sinusoidal') {
    wf.offset = wfOffset.value === '' ? 0 : Number(wfOffset.value);
    wf.amplitude = Number(wfAmplitude.value);
    wf.frequency_Hz = Number(wfFrequency.value);
    wf.phase_rad = wfPhase.value === '' ? 0 : Number(wfPhase.value);
  } else if (type === 'arbitrary') {
    const times = wfTimes.value.split(',').map((v) => Number(v.trim())).filter((v) => Number.isFinite(v));
    const values = wfValues.value.split(',').map((v) => Number(v.trim())).filter((v) => Number.isFinite(v));
    wf.times = times;
    wf.values = values;
    wf.interpolation = wfInterp.value;
  }

  if (!config.waveforms) config.waveforms = {};
  config.waveforms[id] = wf;
  selectedWaveformId = id;
  renderWaveformList();
  setStatus('Waveform saved.');
}

function removeWaveform() {
  if (!config || !config.waveforms) return;
  if (!selectedWaveformId) return;
  delete config.waveforms[selectedWaveformId];
  selectedWaveformId = null;
  renderWaveformList();
  populateWaveformEditor(null, '');
  setStatus('Waveform removed.');
}

function seededRng(seed) {
  let t = seed >>> 0;
  return () => {
    t += 0x6D2B79F5;
    let r = Math.imul(t ^ (t >>> 15), 1 | t);
    r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
    return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
  };
}

function randn(rng) {
  let u = 0;
  let v = 0;
  while (u === 0) u = rng();
  while (v === 0) v = rng();
  return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

function prepareCanvas(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  const attrWidth = Number(canvas.getAttribute('width')) || 300;
  const attrHeight = Number(canvas.getAttribute('height')) || 220;
  const cssWidth = rect.width > 0 ? rect.width : attrWidth;
  const cssHeight = rect.height > 0 ? rect.height : attrHeight;
  const width = Math.max(1, Math.floor(cssWidth * dpr));
  const height = Math.max(1, Math.floor(cssHeight * dpr));
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.scale(dpr, dpr);
  return { ctx, width: cssWidth, height: cssHeight };
}

function drawBaseAxes(ctx, width, height, xLabel, yLabel) {
  const m = { l: 56, r: 16, t: 16, b: 40 };
  const plotW = width - m.l - m.r;
  const plotH = height - m.t - m.b;
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = '#0b1220';
  ctx.fillRect(0, 0, width, height);
  ctx.strokeStyle = '#334155';
  ctx.lineWidth = 1;
  ctx.strokeRect(m.l, m.t, plotW, plotH);

  ctx.fillStyle = '#94a3b8';
  ctx.font = '12px Inter, sans-serif';
  ctx.fillText(xLabel, m.l + plotW / 2 - 20, height - 10);
  ctx.save();
  ctx.translate(14, m.t + plotH / 2 + 20);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();
  return { m, plotW, plotH };
}

function formatTick(v) {
  const av = Math.abs(v);
  if ((av > 0 && av < 1e-3) || av >= 1e4) return v.toExponential(2);
  return Number(v.toFixed(4)).toString();
}

function drawAxisTicks(ctx, m, plotW, plotH, xMin, xMax, yMin, yMax) {
  const ticks = 5;
  ctx.strokeStyle = '#334155';
  ctx.fillStyle = '#94a3b8';
  ctx.lineWidth = 1;
  ctx.font = '11px Inter, sans-serif';

  for (let i = 0; i <= ticks; i++) {
    const u = i / ticks;
    const x = m.l + u * plotW;
    const xv = xMin + u * (xMax - xMin);
    ctx.beginPath();
    ctx.moveTo(x, m.t + plotH);
    ctx.lineTo(x, m.t + plotH + 5);
    ctx.stroke();
    const txt = formatTick(xv);
    ctx.fillText(txt, x - 16, m.t + plotH + 18);
  }

  for (let i = 0; i <= ticks; i++) {
    const u = i / ticks;
    const y = m.t + (1 - u) * plotH;
    const yv = yMin + u * (yMax - yMin);
    ctx.beginPath();
    ctx.moveTo(m.l - 5, y);
    ctx.lineTo(m.l, y);
    ctx.stroke();
    const txt = formatTick(yv);
    ctx.fillText(txt, 6, y + 3);
  }
}

function drawDomainOverlay(ctx, mapX, mapY, orientation) {
  if (!config || !Array.isArray(config.domains) || config.domains.length === 0) return;
  const domain = config.domains[0];
  if (!domain || !domain.geometry) return;

  const geom = domain.geometry || {};
  const origin = geom.origin_m || [0, 0, 0];
  const ox = Number(origin[0]) || 0;
  const oy = Number(origin[1]) || 0;
  const oz = Number(origin[2]) || 0;
  const length = Math.max(0, Number(geom.length_m) || 0);
  const radius = Math.max(0, Number(geom.radius_m) || 0);
  const rIn = Math.max(0, Number(geom.radius_in_m) || 0);
  const rOut = Math.max(0, Number(geom.radius_out_m) || 0);
  const radialExtent = Math.max(radius, rOut);

  ctx.save();
  ctx.strokeStyle = 'rgba(16, 185, 129, 0.9)';
  ctx.fillStyle = 'rgba(16, 185, 129, 0.08)';
  ctx.setLineDash([6, 4]);
  ctx.lineWidth = 1.2;

  if (orientation === 'xy') {
    if (radius > 0) {
      const cpx = mapX(ox);
      const cpy = mapY(oy);
      const rx = Math.abs(mapX(ox + radius) - cpx);
      const ry = Math.abs(mapY(oy + radius) - cpy);
      ctx.beginPath();
      ctx.ellipse(cpx, cpy, rx, ry, 0, 0, 2 * Math.PI);
      ctx.fill();
      ctx.stroke();
    }
    if (rOut > 0) {
      const cpx = mapX(ox);
      const cpy = mapY(oy);
      const rxOut = Math.abs(mapX(ox + rOut) - cpx);
      const ryOut = Math.abs(mapY(oy + rOut) - cpy);
      ctx.beginPath();
      ctx.ellipse(cpx, cpy, rxOut, ryOut, 0, 0, 2 * Math.PI);
      ctx.stroke();
      if (rIn > 0) {
        const rxIn = Math.abs(mapX(ox + rIn) - cpx);
        const ryIn = Math.abs(mapY(oy + rIn) - cpy);
        ctx.beginPath();
        ctx.ellipse(cpx, cpy, rxIn, ryIn, 0, 0, 2 * Math.PI);
        ctx.stroke();
      }
    }
  } else if (orientation === 'xz') {
    if (radialExtent > 0 && length > 0) {
      const x0 = mapX(ox - radialExtent);
      const x1 = mapX(ox + radialExtent);
      const z0 = mapY(oz);
      const z1 = mapY(oz + length);
      const left = Math.min(x0, x1);
      const top = Math.min(z0, z1);
      const w = Math.abs(x1 - x0);
      const h = Math.abs(z1 - z0);
      ctx.fillRect(left, top, w, h);
      ctx.strokeRect(left, top, w, h);
    }
  } else if (orientation === 'yz') {
    if (radialExtent > 0 && length > 0) {
      const y0 = mapX(oy - radialExtent);
      const y1 = mapX(oy + radialExtent);
      const z0 = mapY(oz);
      const z1 = mapY(oz + length);
      const left = Math.min(y0, y1);
      const top = Math.min(z0, z1);
      const w = Math.abs(y1 - y0);
      const h = Math.abs(z1 - z0);
      ctx.fillRect(left, top, w, h);
      ctx.strokeRect(left, top, w, h);
    }
  }
  ctx.restore();
}

function plotIonCloud() {
  try {
  if (!config) {
    setStatus('No ions to plot.');
    return;
  }
  const species = ionList && ionList.children.length ? collectIonsFromUI() : (config.ions?.species || []);
  if (!Array.isArray(species) || species.length === 0) {
    setStatus('No ions to plot.');
    return;
  }
  const orientation = ionPlotOrientation?.value || 'xy';
  const axisMap = {
    xy: [0, 1, 'x (m)', 'y (m)'],
    xz: [0, 2, 'x (m)', 'z (m)'],
    yz: [1, 2, 'y (m)', 'z (m)']
  };
  const [axI, ayI, xLabel, yLabel] = axisMap[orientation] || axisMap.xy;
  const points = [];
  const maxPoints = 800;
  let xMin = Infinity;
  let xMax = -Infinity;
  let yMin = Infinity;
  let yMax = -Infinity;
  const sigmaScale = 3;

  species.forEach((ion, ionIndex) => {
    const count = Math.min(ion.count || 0, Math.floor(maxPoints / Math.max(1, species.length)));
    const pos = ion.position || {};
    const center = pos.center || [0, 0, 0];
    const std = pos.std || [0, 0, 0];
    const c3 = [Number(center[0]) || 0, Number(center[1]) || 0, Number(center[2]) || 0];
    const s3 = [Math.abs(Number(std[0]) || 0), Math.abs(Number(std[1]) || 0), Math.abs(Number(std[2]) || 0)];
    const cx = c3[axI];
    const cy = c3[ayI];
    const sx = s3[axI];
    const sy = s3[ayI];
    const seedOffset = (ionIndex + 1) * 7919;
    const localRng = seededRng(123456 + seedOffset);

    if (pos.type === 'uniform_box') {
      const pmin = pos.min || [0, 0, 0];
      const pmax = pos.max || [0, 0, 0];
      const min3 = [
        Math.min(Number(pmin[0]) || 0, Number(pmax[0]) || 0),
        Math.min(Number(pmin[1]) || 0, Number(pmax[1]) || 0),
        Math.min(Number(pmin[2]) || 0, Number(pmax[2]) || 0)
      ];
      const max3 = [
        Math.max(Number(pmin[0]) || 0, Number(pmax[0]) || 0),
        Math.max(Number(pmin[1]) || 0, Number(pmax[1]) || 0),
        Math.max(Number(pmin[2]) || 0, Number(pmax[2]) || 0)
      ];
      const minX = min3[axI];
      const maxX = max3[axI];
      const minY = min3[ayI];
      const maxY = max3[ayI];
      xMin = Math.min(xMin, minX);
      xMax = Math.max(xMax, maxX);
      yMin = Math.min(yMin, minY);
      yMax = Math.max(yMax, maxY);
      for (let i = 0; i < count; i++) {
        const xyz = [
          min3[0] + localRng() * (max3[0] - min3[0]),
          min3[1] + localRng() * (max3[1] - min3[1]),
          min3[2] + localRng() * (max3[2] - min3[2])
        ];
        points.push({ x: xyz[axI], y: xyz[ayI] });
      }
    } else if (pos.type === 'uniform_sphere') {
      const radius = Math.max(0, Number(pos.radius) || 0);
      xMin = Math.min(xMin, cx - radius);
      xMax = Math.max(xMax, cx + radius);
      yMin = Math.min(yMin, cy - radius);
      yMax = Math.max(yMax, cy + radius);
      for (let i = 0; i < count; i++) {
        const u = localRng();
        const v = localRng();
        const w = localRng();
        const rr = radius * Math.cbrt(u);
        const cosTheta = 2 * v - 1;
        const sinTheta = Math.sqrt(Math.max(0, 1 - cosTheta * cosTheta));
        const phi = 2 * Math.PI * w;
        const dx = rr * sinTheta * Math.cos(phi);
        const dy = rr * sinTheta * Math.sin(phi);
        const dz = rr * cosTheta;
        const xyz = [c3[0], c3[1], c3[2]];
        xyz[0] += dx;
        xyz[1] += dy;
        xyz[2] += dz;
        points.push({ x: xyz[axI], y: xyz[ayI] });
      }
    } else if (pos.type === 'uniform_cylinder') {
      const radius = Math.max(0, Number(pos.radius) || 0);
      const length = Math.max(0, Number(pos.length) || 0);
      const halfLen = length * 0.5;
      const axisExtent = (axisIndex) => (axisIndex === 2 ? halfLen : radius);
      xMin = Math.min(xMin, cx - axisExtent(axI));
      xMax = Math.max(xMax, cx + axisExtent(axI));
      yMin = Math.min(yMin, cy - axisExtent(ayI));
      yMax = Math.max(yMax, cy + axisExtent(ayI));
      for (let i = 0; i < count; i++) {
        const t = 2 * Math.PI * localRng();
        const r = radius * Math.sqrt(localRng());
        const dx = r * Math.cos(t);
        const dy = r * Math.sin(t);
        const dz = (localRng() - 0.5) * length;
        const xyz = [c3[0], c3[1], c3[2]];
        xyz[0] += dx;
        xyz[1] += dy;
        xyz[2] += dz;
        points.push({ x: xyz[axI], y: xyz[ayI] });
      }
    } else if (pos.type === 'point') {
      const boundX = 0;
      const boundY = 0;
      xMin = Math.min(xMin, cx - boundX);
      xMax = Math.max(xMax, cx + boundX);
      yMin = Math.min(yMin, cy - boundY);
      yMax = Math.max(yMax, cy + boundY);
      for (let i = 0; i < count; i++) {
        points.push({ x: cx, y: cy });
      }
    } else {
      const boundX = Math.abs(sx) * sigmaScale;
      const boundY = Math.abs(sy) * sigmaScale;
      xMin = Math.min(xMin, cx - boundX);
      xMax = Math.max(xMax, cx + boundX);
      yMin = Math.min(yMin, cy - boundY);
      yMax = Math.max(yMax, cy + boundY);
      for (let i = 0; i < count; i++) {
        const xyz = [
          c3[0] + randn(localRng) * s3[0],
          c3[1] + randn(localRng) * s3[1],
          c3[2] + randn(localRng) * s3[2]
        ];
        points.push({ x: xyz[axI], y: xyz[ayI] });
      }
    }
  });

  if (Array.isArray(config.domains) && config.domains.length > 0) {
    const d0 = config.domains[0] || {};
    const g = d0.geometry || {};
    const o = g.origin_m || [0, 0, 0];
    const ox = Number(o[0]) || 0;
    const oy = Number(o[1]) || 0;
    const oz = Number(o[2]) || 0;
    const len = Math.max(0, Number(g.length_m) || 0);
    const rad = Math.max(0, Number(g.radius_m) || 0);
    const rout = Math.max(0, Number(g.radius_out_m) || 0);
    const rEff = Math.max(rad, rout);
    const geomRanges = {
      xy: { min: [ox - rEff, oy - rEff], max: [ox + rEff, oy + rEff] },
      xz: { min: [ox - rEff, oz], max: [ox + rEff, oz + len] },
      yz: { min: [oy - rEff, oz], max: [oy + rEff, oz + len] }
    };
    const gr = geomRanges[orientation];
    if (gr) {
      xMin = Math.min(xMin, gr.min[0]);
      xMax = Math.max(xMax, gr.max[0]);
      yMin = Math.min(yMin, gr.min[1]);
      yMax = Math.max(yMax, gr.max[1]);
    }
  }

  if (!Number.isFinite(xMin) || !Number.isFinite(xMax)) {
    xMin = -1;
    xMax = 1;
  }
  if (!Number.isFinite(yMin) || !Number.isFinite(yMax)) {
    yMin = -1;
    yMax = 1;
  }
  if (xMin === xMax) {
    xMin -= 1e-6;
    xMax += 1e-6;
  }
  if (yMin === yMax) {
    yMin -= 1e-6;
    yMax += 1e-6;
  }
  const padX = (xMax - xMin) * 0.1;
  const padY = (yMax - yMin) * 0.1;

  const xLo = xMin - padX;
  const xHi = xMax + padX;
  const yLo = yMin - padY;
  const yHi = yMax + padY;
  const { ctx, width, height } = prepareCanvas(ionPlotCanvas);
  const { m, plotW, plotH } = drawBaseAxes(ctx, width, height, xLabel, yLabel);
  drawAxisTicks(ctx, m, plotW, plotH, xLo, xHi, yLo, yHi);
  const mapX = (x) => m.l + ((x - xLo) / (xHi - xLo)) * plotW;
  const mapY = (y) => m.t + (1 - (y - yLo) / (yHi - yLo)) * plotH;

  drawDomainOverlay(ctx, mapX, mapY, orientation);

  ctx.fillStyle = '#60a5fa';
  for (const p of points) {
    const px = mapX(p.x);
    const py = mapY(p.y);
    ctx.beginPath();
    ctx.arc(px, py, 1.8, 0, Math.PI * 2);
    ctx.fill();
  }
  setStatus(`Ion cloud plotted (${orientation.toUpperCase()}).`);
  } catch (err) {
    setStatus(`Ion cloud plot error: ${err.message || err}`);
  }
}

function sampleWaveform(wf) {
  const samples = 200;
  const ts = [];
  const ys = [];
  let tMin = 0;
  let tMax = 1e-3;
  const toNum = (value, fallback = 0) => {
    const num = Number(value);
    return Number.isFinite(num) ? num : fallback;
  };

  if (wf.type === 'linear') {
    tMin = toNum(wf.start_time_s, 0);
    tMax = toNum(wf.end_time_s, 1e-3);
    if (tMax <= tMin) tMax = tMin + 1e-9;
  } else if (wf.type === 'sinusoidal') {
    tMin = 0;
    const f = Math.max(1e-6, toNum(wf.frequency_Hz, 1));
    tMax = (1 / f) * 5;
  } else if (wf.type === 'arbitrary' && Array.isArray(wf.times) && wf.times.length > 1) {
    tMin = toNum(wf.times[0], 0);
    tMax = toNum(wf.times[wf.times.length - 1], 1e-3);
    if (tMax <= tMin) tMax = tMin + 1e-9;
  }

  for (let i = 0; i < samples; i++) {
    const t = tMin + (tMax - tMin) * (i / (samples - 1));
    let y = 0;
    if (wf.type === 'constant') {
      y = toNum(wf.value, 0);
    } else if (wf.type === 'linear') {
      const start = toNum(wf.start, 0);
      const end = toNum(wf.end, 0);
      const t0 = toNum(wf.start_time_s, 0);
      let t1 = toNum(wf.end_time_s, 1e-3);
      if (t1 <= t0) t1 = t0 + 1e-9;
      const u = Math.min(1, Math.max(0, (t - t0) / (t1 - t0)));
      y = start + (end - start) * u;
      if (wf.clamp === false && (t < t0 || t > t1)) y = start;
    } else if (wf.type === 'sinusoidal') {
      const offset = toNum(wf.offset, 0);
      const amp = toNum(wf.amplitude, 0);
      const f = Math.max(1e-6, toNum(wf.frequency_Hz, 1));
      const phase = toNum(wf.phase_rad, 0);
      y = offset + amp * Math.sin(2 * Math.PI * f * t + phase);
    } else if (wf.type === 'arbitrary' && Array.isArray(wf.times) && Array.isArray(wf.values)) {
      const times = wf.times.map((v) => toNum(v, 0));
      const values = wf.values.map((v) => toNum(v, 0));
      let idx = times.findIndex((tt) => tt >= t);
      if (idx <= 0) y = values[0];
      else if (idx === -1) y = values[values.length - 1];
      else {
        const t0 = times[idx - 1];
        const t1 = times[idx];
        const v0 = values[idx - 1];
        const v1 = values[idx];
        const u = (t - t0) / (t1 - t0);
        y = v0 + (v1 - v0) * u;
      }
    }
    if (Number.isFinite(t) && Number.isFinite(y)) {
      ts.push(t);
      ys.push(y);
    }
  }

  return { ts, ys };
}

function sampleWaveformInRange(wf, tMin, tMax, samples = 300) {
  const toNum = (value, fallback = 0) => {
    const num = Number(value);
    return Number.isFinite(num) ? num : fallback;
  };
  let a = toNum(tMin, 0);
  let b = toNum(tMax, 1e-3);
  if (b <= a) b = a + 1e-9;

  const ts = [];
  const ys = [];
  for (let i = 0; i < samples; i++) {
    const t = a + (b - a) * (i / (samples - 1));
    let y = 0;
    if (wf.type === 'constant') {
      y = toNum(wf.value, 0);
    } else if (wf.type === 'linear') {
      const start = toNum(wf.start, 0);
      const end = toNum(wf.end, 0);
      const t0 = toNum(wf.start_time_s, 0);
      let t1 = toNum(wf.end_time_s, 1e-3);
      if (t1 <= t0) t1 = t0 + 1e-9;
      const u = Math.min(1, Math.max(0, (t - t0) / (t1 - t0)));
      y = start + (end - start) * u;
      if (wf.clamp === false && (t < t0 || t > t1)) y = start;
    } else if (wf.type === 'sinusoidal') {
      const offset = toNum(wf.offset, 0);
      const amp = toNum(wf.amplitude, 0);
      const f = Math.max(1e-6, toNum(wf.frequency_Hz, 1));
      const phase = toNum(wf.phase_rad, 0);
      y = offset + amp * Math.sin(2 * Math.PI * f * t + phase);
    } else if (wf.type === 'arbitrary' && Array.isArray(wf.times) && Array.isArray(wf.values)) {
      const times = wf.times.map((v) => toNum(v, 0));
      const values = wf.values.map((v) => toNum(v, 0));
      let idx = times.findIndex((tt) => tt >= t);
      if (idx <= 0) y = values[0];
      else if (idx === -1) y = values[values.length - 1];
      else {
        const t0 = times[idx - 1];
        const t1 = times[idx];
        const v0 = values[idx - 1];
        const v1 = values[idx];
        const u = (t - t0) / (t1 - t0);
        y = v0 + (v1 - v0) * u;
      }
    }
    if (Number.isFinite(t) && Number.isFinite(y)) {
      ts.push(t);
      ys.push(y);
    }
  }
  return { ts, ys };
}

function plotWaveform() {
  try {
  if (!config || !config.waveforms) {
    config = config || buildDefaultConfig();
  }
  let wf = null;
  let wfLabel = selectedWaveformId || '';
  if (wfId && wfId.value.trim()) {
    wf = buildWaveformFromEditor();
    wfLabel = wfId.value.trim();
  }
  if (!wf && config.waveforms) {
    if (!selectedWaveformId) {
      const keys = Object.keys(config.waveforms);
      if (keys.length === 0) {
        setStatus('No waveforms available.');
        return;
      }
      selectedWaveformId = keys[0];
      renderWaveformList();
    }
    wf = config.waveforms[selectedWaveformId];
    wfLabel = selectedWaveformId;
  }
  if (!wf) {
    setStatus('Select a waveform to plot.');
    return;
  }
  const simTotalTime = Number(config?.simulation?.total_time_s);
  const hasSimWindow = Number.isFinite(simTotalTime) && simTotalTime > 0;
  const { ts, ys } = hasSimWindow
    ? sampleWaveformInRange(wf, 0, simTotalTime, 300)
    : sampleWaveform(wf);
  if (!ts.length || !ys.length) {
    setStatus('Waveform has no sample data.');
    return;
  }
  let xMin = ts[0];
  let xMax = ts[ts.length - 1];
  let yMin = Math.min(...ys);
  let yMax = Math.max(...ys);
  if (xMin === xMax) xMax = xMin + 1e-9;
  if (yMin === yMax) {
    yMin -= 1e-9;
    yMax += 1e-9;
  }
  const yPad = (yMax - yMin) * 0.1;
  yMin -= yPad;
  yMax += yPad;

  const { ctx, width, height } = prepareCanvas(waveformPlotCanvas);
  const { m, plotW, plotH } = drawBaseAxes(ctx, width, height, 'time (s)', 'value');
  drawAxisTicks(ctx, m, plotW, plotH, xMin, xMax, yMin, yMax);
  const mapX = (x) => m.l + ((x - xMin) / (xMax - xMin)) * plotW;
  const mapY = (y) => m.t + (1 - (y - yMin) / (yMax - yMin)) * plotH;

  ctx.strokeStyle = '#f59e0b';
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  ts.forEach((t, i) => {
    const px = mapX(t);
    const py = mapY(ys[i]);
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
  });
  ctx.stroke();

  ctx.fillStyle = '#e5e7eb';
  ctx.font = '12px Inter, sans-serif';
  ctx.fillText(wfLabel || 'waveform', m.l + 6, m.t + 14);
  setStatus('Waveform plotted.');
  } catch (err) {
    setStatus(`Waveform plot error: ${err.message || err}`);
  }
}

function buildWaveformFromEditor() {
  const type = wfType.value;
  const wf = { type };
  if (type === 'constant') {
    wf.value = wfValue.value === '' ? 0 : Number(wfValue.value);
  } else if (type === 'linear') {
    wf.start = wfStart.value === '' ? 0 : Number(wfStart.value);
    wf.end = wfEnd.value === '' ? 0 : Number(wfEnd.value);
    wf.start_time_s = wfStartTime.value === '' ? 0 : Number(wfStartTime.value);
    wf.end_time_s = wfEndTime.value === '' ? 1e-3 : Number(wfEndTime.value);
    wf.clamp = wfClamp.value !== 'false';
  } else if (type === 'sinusoidal') {
    wf.offset = wfOffset.value === '' ? 0 : Number(wfOffset.value);
    wf.amplitude = wfAmplitude.value === '' ? 0 : Number(wfAmplitude.value);
    wf.frequency_Hz = wfFrequency.value === '' ? 1 : Number(wfFrequency.value);
    wf.phase_rad = wfPhase.value === '' ? 0 : Number(wfPhase.value);
  } else if (type === 'arbitrary') {
    const times = wfTimes.value.split(',').map((v) => Number(v.trim())).filter((v) => Number.isFinite(v));
    const values = wfValues.value.split(',').map((v) => Number(v.trim())).filter((v) => Number.isFinite(v));
    wf.times = times;
    wf.values = values;
    wf.interpolation = wfInterp.value;
  }
  return wf;
}

function populateReactionEditor(rxn) {
  rxnId.value = rxn?.id ?? '';
  rxnReactant.value = rxn?.reactant ?? '';
  rxnProduct.value = rxn?.product ?? '';
  rxnRate.value = rxn?.rate_constant ?? '';
  rxnModel.value = rxn?.rate_model ?? 'Constant';
  rxnEa.value = rxn?.activation_energy_eV ?? '';
  rxnN.value = rxn?.temperature_exponent ?? '';
  rxnTref.value = rxn?.reference_temperature_K ?? '';
  rxnOrder.value = rxn?.order ? JSON.stringify(rxn.order, null, 2) : '';
}

function addOrUpdateReaction() {
  const id = rxnId.value.trim();
  if (!id) {
    setStatus('Reaction id is required.');
    return;
  }
  const rxn = {
    id,
    reactant: rxnReactant.value.trim(),
    product: rxnProduct.value.trim(),
    rate_constant: rxnRate.value === '' ? undefined : Number(rxnRate.value),
    rate_model: rxnModel.value || 'Constant',
    activation_energy_eV: rxnEa.value === '' ? undefined : Number(rxnEa.value),
    temperature_exponent: rxnN.value === '' ? undefined : Number(rxnN.value),
    reference_temperature_K: rxnTref.value === '' ? undefined : Number(rxnTref.value)
  };

  if (rxnOrder.value.trim() !== '') {
    try {
      rxn.order = JSON.parse(rxnOrder.value.trim());
    } catch (err) {
      setStatus('Order must be valid JSON array.');
      return;
    }
  }

  Object.keys(rxn).forEach((k) => rxn[k] === undefined && delete rxn[k]);
  if (!reactionDb.reactions) reactionDb.reactions = [];
  const existing = reactionDb.reactions.findIndex((r) => r.id === id);
  if (existing >= 0) reactionDb.reactions[existing] = rxn;
  else reactionDb.reactions.push(rxn);
  selectedReactionId = id;
  renderReactionList();
  setStatus('Reaction entry saved.');
}

function removeReaction() {
  if (!reactionDb.reactions) return;
  reactionDb.reactions = reactionDb.reactions.filter((r) => r.id !== selectedReactionId);
  selectedReactionId = null;
  renderReactionList();
  populateReactionEditor(null);
  setStatus('Reaction entry removed.');
}

function removeSpecies() {
  if (!selectedSpeciesId) return;
  delete speciesDb.species[selectedSpeciesId];
  selectedSpeciesId = null;
  renderSpeciesList();
  populateSpeciesEditor(null);
  setStatus('Species entry removed.');
}

function populateDefaultsUI() {
  document.getElementById('simTotalTime').value = DEFAULTS.simulation.total_time_s;
  document.getElementById('simDt').value = DEFAULTS.simulation.dt_s;
  document.getElementById('simIntegrator').value = DEFAULTS.simulation.integrator;
  document.getElementById('simWrite').value = DEFAULTS.simulation.write_interval;
  document.getElementById('simSeed').value = DEFAULTS.simulation.rng_seed;

  document.getElementById('outFolder').value = DEFAULTS.output.folder;
  document.getElementById('outTrajectory').value = DEFAULTS.output.trajectory_file;
  document.getElementById('outProgress').value = DEFAULTS.output.print_progress ? 'true' : 'false';

  domainInstrument.value = DEFAULTS.domain.instrument;
  applyFieldVisibility(domainInstrument.value);

  document.getElementById('envPressure').value = DEFAULTS.domain.env.pressure_Pa;
  document.getElementById('envTemperature').value = DEFAULTS.domain.env.temperature_K;
  envGasSelect.value = DEFAULTS.domain.env.gas_species;
  envGasCustomRow.classList.add('hidden');
  envGasCustom.value = '';

  const origin = DEFAULTS.domain.geometry.origin_m;
  document.getElementById('geomOriginX').value = origin[0];
  document.getElementById('geomOriginY').value = origin[1];
  document.getElementById('geomOriginZ').value = origin[2];
  document.getElementById('geomLength').value = DEFAULTS.domain.geometry.length_m;
  document.getElementById('geomRadius').value = DEFAULTS.domain.geometry.radius_m;

  document.getElementById('dcAxial').value = DEFAULTS.domain.fields.DC.axial_V;
  document.getElementById('dcRadial').value = DEFAULTS.domain.fields.DC.radial_V;
  document.getElementById('dcQuad').value = DEFAULTS.domain.fields.DC.quad_V;
  document.getElementById('rfVoltage').value = DEFAULTS.domain.fields.RF.voltage_V;
  document.getElementById('rfFrequency').value = DEFAULTS.domain.fields.RF.frequency_Hz;
  document.getElementById('acVoltage').value = DEFAULTS.domain.fields.AC.voltage_V;
  document.getElementById('acFrequency').value = DEFAULTS.domain.fields.AC.frequency_Hz;

  dbSpecies.value = DEFAULTS.databases.species_database;
  dbReaction.value = DEFAULTS.databases.reaction_database;
  ionCloud.value = DEFAULTS.databases.ion_cloud;

  config = buildDefaultConfig();
  populateDomainSelect();
  loadDomainValues();
  renderIonList();
}

fileInput.addEventListener('change', async (e) => {
  const file = e.target.files[0];
  if (!file) return;
  const text = await file.text();
  try {
    config = JSON.parse(text);
  } catch (err) {
    setStatus('Invalid JSON file.');
    return;
  }
  setStatus(`Loaded ${file.name}`);
  populateTopLevel();
  populateDomainSelect();
  loadDomainValues();
  renderIonList();
});

if (speciesDbInput) {
  speciesDbInput.addEventListener('change', async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    try {
      const text = await file.text();
      const parsed = JSON.parse(text);
      if (!parsed.species || typeof parsed.species !== 'object') {
        setStatus('Invalid species DB (expected { species: { ... } }).');
        return;
      }
      Object.entries(parsed.species).forEach(([id, spec]) => {
        if (spec && typeof spec === 'object' && !spec.id) spec.id = id;
      });
      speciesDb = parsed;
      selectedSpeciesId = null;
      renderSpeciesList();
      populateSpeciesEditor(null);
      setStatus(`Loaded species DB: ${file.name}`);
    } catch (err) {
      setStatus(`Failed to load species DB: ${err.message || err}`);
    }
  });
}

if (downloadSpeciesDb) {
  downloadSpeciesDb.addEventListener('click', () => {
    const blob = new Blob([JSON.stringify(speciesDb, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'species_database.json';
    a.click();
    URL.revokeObjectURL(url);
    setStatus('Downloaded species database.');
  });
}

if (openSpeciesDbPath) {
  openSpeciesDbPath.addEventListener('click', async () => {
    const path = dbSpecies.value.trim();
    if (!path) {
      setStatus('species_database path is empty.');
      return;
    }
    try {
      const parsed = await loadJsonFromPath(path);
      if (!parsed.species || typeof parsed.species !== 'object') {
        setStatus('Invalid species DB (expected { species: { ... } }).');
        return;
      }
      Object.entries(parsed.species).forEach(([id, spec]) => {
        if (spec && typeof spec === 'object' && !spec.id) spec.id = id;
      });
      speciesDb = parsed;
      selectedSpeciesId = null;
      renderSpeciesList();
      populateSpeciesEditor(null);
      setStatus(`Loaded species DB from ${path}`);
    } catch (err) {
      setStatus(`Failed to load species DB: ${err.message || err}`);
    }
  });
}

function initViewer() {
  if (!geomViewer || viewer) return;
  const width = geomViewer.clientWidth;
  const height = geomViewer.clientHeight;
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0b1220);
  const camera = new THREE.PerspectiveCamera(45, width / height, 0.01, 1000);
  camera.position.set(0, 0, 5);

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(width, height);
  geomViewer.innerHTML = '';
  geomViewer.appendChild(renderer.domElement);

  const ambient = new THREE.AmbientLight(0xffffff, 0.6);
  scene.add(ambient);
  const dirLight = new THREE.DirectionalLight(0xffffff, 0.8);
  dirLight.position.set(5, 5, 5);
  scene.add(dirLight);

  const controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;

  viewer = { scene, camera, renderer, controls, meshes: [] };

  const animate = () => {
    if (!viewer) return;
    requestAnimationFrame(animate);
    viewer.controls.update();
    viewer.renderer.render(viewer.scene, viewer.camera);
  };
  animate();
}

function clearViewer() {
  if (!viewer) return;
  viewer.meshes.forEach((m) => viewer.scene.remove(m));
  viewer.meshes = [];
  viewerAtoms = [];
}

function elementColor(el) {
  const map = {
    H: 0xffffff,
    C: 0x9ca3af,
    N: 0x3b82f6,
    O: 0xef4444,
    S: 0xf59e0b,
    P: 0x22c55e
  };
  return map[el] || 0xa855f7;
}

function elementRadius(el) {
  const map = {
    H: 0.2,
    C: 0.35,
    N: 0.35,
    O: 0.35,
    S: 0.45,
    P: 0.45
  };
  return map[el] || 0.35;
}

function loadGeometryJson(geom) {
  if (!geom?.molecule?.atoms || !Array.isArray(geom.molecule.atoms)) {
    setStatus('Invalid geometry JSON (expected molecule.atoms).');
    return;
  }
  initViewer();
  clearViewer();

  const atoms = geom.molecule.atoms;
  const positions = atoms.map((a) => a.pos || [0, 0, 0]);
  const center = positions.reduce(
    (acc, p) => [acc[0] + p[0], acc[1] + p[1], acc[2] + p[2]],
    [0, 0, 0]
  ).map((v) => v / positions.length);

  atoms.forEach((atom) => {
    const color = elementColor(atom.element);
    const radius = elementRadius(atom.element);
    const geometry = new THREE.SphereGeometry(radius, 24, 24);
    const material = new THREE.MeshStandardMaterial({ color });
    const mesh = new THREE.Mesh(geometry, material);
    const [x, y, z] = atom.pos || [0, 0, 0];
    mesh.position.set(x - center[0], y - center[1], z - center[2]);
    viewer.scene.add(mesh);
    viewer.meshes.push(mesh);
  });

  viewer.camera.position.set(0, 0, 6);
  viewer.controls.update();
  viewerAtoms = atoms;
  setStatus(`Loaded geometry with ${atoms.length} atoms.`);
}

if (domainSelect) {
  domainSelect.addEventListener('change', () => {
    loadDomainValues();
  });
}

if (domainInstrument) {
  domainInstrument.addEventListener('change', () => {
    applyFieldVisibility(domainInstrument.value);
  });
}

if (addDomainBtn) {
  addDomainBtn.addEventListener('click', addDomain);
}

if (removeDomainBtn) {
  removeDomainBtn.addEventListener('click', removeDomain);
}

if (applyDatabasesBtn) {
  applyDatabasesBtn.addEventListener('click', applyDatabases);
}

if (addIonBtn) {
  addIonBtn.addEventListener('click', () => {
    if (!config) return;
    if (!config.ions) config.ions = { species: [] };
    config.ions.species.push(JSON.parse(JSON.stringify(DEFAULTS.ions[0])));
    renderIonList();
  });
}

if (removeIonBtn) {
  removeIonBtn.addEventListener('click', () => {
    if (!config || !config.ions || !Array.isArray(config.ions.species)) return;
    const selected = ionList.querySelector('.ion-card.selected');
    if (!selected) return;
    const idx = Number(selected.dataset.index);
    config.ions.species.splice(idx, 1);
    renderIonList();
  });
}

if (applyIonsBtn) {
  applyIonsBtn.addEventListener('click', applyIons);
}

if (addSpeciesBtn) {
  addSpeciesBtn.addEventListener('click', addOrUpdateSpecies);
}

if (removeSpeciesBtn) {
  removeSpeciesBtn.addEventListener('click', removeSpecies);
}

if (addReactionBtn) {
  addReactionBtn.addEventListener('click', addOrUpdateReaction);
}

if (removeReactionBtn) {
  removeReactionBtn.addEventListener('click', removeReaction);
}

if (addWaveformBtn) {
  addWaveformBtn.addEventListener('click', () => {
    if (!config) return;
    if (!config.waveforms) config.waveforms = {};
    populateWaveformEditor({ type: 'constant', value: 0 }, '');
  });
}

if (removeWaveformBtn) {
  removeWaveformBtn.addEventListener('click', removeWaveform);
}

if (saveWaveformBtn) {
  saveWaveformBtn.addEventListener('click', addOrUpdateWaveform);
}

if (wfType) {
  wfType.addEventListener('change', () => setWaveformVisibility(wfType.value));
}

if (plotIonCloudBtn) {
  plotIonCloudBtn.addEventListener('click', plotIonCloud);
}

if (ionPlotOrientation) {
  ionPlotOrientation.addEventListener('change', () => {
    if (config) plotIonCloud();
  });
}

if (plotWaveformBtn) {
  plotWaveformBtn.addEventListener('click', plotWaveform);
}

if (speciesGeomInput) {
  speciesGeomInput.addEventListener('change', async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    try {
      const text = await file.text();
      const parsed = JSON.parse(text);
      loadGeometryJson(parsed);
    } catch (err) {
      setStatus('Failed to load geometry JSON.');
    }
  });
}

if (clearGeomBtn) {
  clearGeomBtn.addEventListener('click', () => {
    clearViewer();
    setStatus('Geometry viewer cleared.');
  });
}

if (loadGeomPathBtn) {
  loadGeomPathBtn.addEventListener('click', async () => {
    const path = spGeomFile.value.trim();
    if (!path) {
      setStatus('geometry_file path is empty.');
      return;
    }
    try {
      const parsed = await loadJsonFromPath(path);
      loadGeometryJson(parsed);
    } catch (err) {
      setStatus(`Failed to load geometry: ${err.message || err}`);
    }
  });
}

if (envGasSelect) {
  envGasSelect.addEventListener('change', () => {
    if (envGasSelect.value === '__custom__') {
      envGasCustomRow.classList.remove('hidden');
      envGasCustom.focus();
    } else {
      envGasCustomRow.classList.add('hidden');
      envGasCustom.value = '';
    }
  });
}

document.getElementById('applySimulation').addEventListener('click', updateTopLevelSection);
document.getElementById('applyOutput').addEventListener('click', updateTopLevelSection);
document.getElementById('applyEnv').addEventListener('click', applyEnv);
document.getElementById('applyGeometry').addEventListener('click', applyGeometry);
document.getElementById('applyFields').addEventListener('click', applyFields);

downloadBtn.addEventListener('click', () => {
  if (!config) {
    config = buildDefaultConfig();
    populateDomainSelect();
    loadDomainValues();
  }
  applyDatabases(true);
  if (typeof config.species_database === 'string' && config.species_database.trim() === '') delete config.species_database;
  if (typeof config.reaction_database === 'string' && config.reaction_database.trim() === '') delete config.reaction_database;
  if (typeof config.ion_cloud === 'string' && config.ion_cloud.trim() === '') delete config.ion_cloud;
  const blob = new Blob([JSON.stringify(config, null, 2)], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'icarion_config.json';
  a.click();
  URL.revokeObjectURL(url);
  setStatus('Downloaded JSON.');
});

populateDefaultsUI();
setStatus('Defaults loaded. Load a JSON to edit and download.');
renderSpeciesList();
renderReactionList();
renderWaveformList();
