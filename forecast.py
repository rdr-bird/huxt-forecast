"""
HUXt-ensemble-ennuste Maata kohti tuleville CME:ille.

Syötteet (kaikki julkisia):
  * CME:t:       NASA DONKI, WSA-Enlil-simulaatioiden kartio-CME-syötteet (cmeInputs)
  * Taustatuuli: NOAA SWPC, mitattu aurinkotuuli (plasma-7-day) + WSA-Enlil-ennuste Maan kohdalla
                 (enlil_time_series.json, CME-materiaali poistettu cloud-kentän avulla)
Tulos: docs/huxt_forecast.json

Ajo:  python forecast.py               (verkosta)
      python forecast.py --offline     (testidata kansiosta testdata/)
"""
import argparse, datetime as D, hashlib, json, os, sys, time, warnings
import numpy as np
import requests
import astropy.units as u
from astropy.time import Time

warnings.filterwarnings('ignore')
import huxt.huxt as H
import huxt.huxt_inputs as HIN

UTC = D.timezone.utc
SYN = 27.2753                       # synodinen kierto (vrk)
R_IN = 21.5 * u.solRad              # DONKI:n kartio-CME:t annetaan 21,5 Rs:n korkeudella
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'docs', 'huxt_forecast.json')

N_MEMBERS = int(os.environ.get('HUXT_MEMBERS', 500))
MAX_ABS_LON = 90.0                 # tätä kauempana Maan suunnasta (HEEQ) olevat CME:t eivät voi osua -> ei ajeta
FORECAST_DAYS = 5.0                 # kuinka pitkälle eteenpäin nykyhetkestä
CME_LOOKBACK_DAYS = 5.0             # mukaan CME:t, jotka ohittivat 21,5 Rs viimeisen 5 vrk aikana
RERUN_HOURS = 12                    # ajetaan uudelleen viimeistään näin usein, vaikka CME:t eivät muutu
# Ensemble-hajonta (tasajakauma ±). Valittu niin, että testitapaus (CME 11.–12.10.2026)
# toisti Readingin HUXt-ensemblen osumatodennäköisyydet (65 % / 25 %) mahdollisimman hyvin.
PERT = dict(lon=15.0, lat=15.0, width=0.20, v=0.15, amb_shift=15.0)

NASA_KEY = os.environ.get('NASA_API_KEY', 'DEMO_KEY')
PROXY = os.environ.get('DONKI_PROXY', 'https://revontulet-jade.vercel.app/api/proxy')
HDRS = {'User-Agent': 'huxt-forecast (github actions)'}


# ---------------------------------------------------------------- data
def get_json(url, **kw):
    r = requests.get(url, headers=HDRS, timeout=60, **kw)
    r.raise_for_status()
    return r.json()


def fetch_donki(now):
    start = (now - D.timedelta(days=8)).strftime('%Y-%m-%d')
    end = now.strftime('%Y-%m-%d')
    try:
        return get_json('https://api.nasa.gov/DONKI/WSAEnlilSimulations',
                        params=dict(startDate=start, endDate=end, api_key=NASA_KEY))
    except Exception as e:  # NASA:n rajapinta voi ohjautua tai rajoittaa -> sivuston oma proxy
        print('DONKI suoraan epäonnistui, käytetään proxya:', e)
        return get_json(PROXY, params=dict(source='donki-enlil', startDate=start, endDate=end))


def parse_t(s):
    s = s.replace('Z', '')
    return D.datetime.fromisoformat(s).replace(tzinfo=UTC)


def collect_cmes(sims, now):
    """Kaikki tuoreet kartio-CME:t simulaatioiden syötteistä; tuorein analyysi voittaa."""
    best = {}
    for s in sorted(sims or [], key=lambda s: s.get('modelCompletionTime') or ''):
        for c in s.get('cmeInputs') or []:
            try:
                t = parse_t(c['time21_5'])
                cid = c.get('CMEID') or c['cmeStartTime']
                if not c.get('isMostAccurate', True):
                    continue
                if not (now - D.timedelta(days=CME_LOOKBACK_DAYS) <= t <= now):
                    continue
                best[cid] = dict(id=cid, start=c.get('cmeStartTime'), t215=t.isoformat(),
                                 lon=float(c['longitude']), lat=float(c['latitude']),
                                 half=float(c['halfAngle']), v=float(c['speed']),
                                 sim=s.get('simulationID'), link=s.get('link'))
            except (KeyError, TypeError, ValueError):
                continue
    return sorted(best.values(), key=lambda c: c['t215'])


def fetch_ambient(now):
    """Maan kohdan nopeusaikasarja: mitattu (7 vrk taakse) + WSA-Enlil-ennuste (CME:t poistettu)."""
    pts, src = [], []
    try:
        rows = get_json('https://services.swpc.noaa.gov/products/solar-wind/plasma-7-day.json')
        hdr = rows[0]; it, iv = hdr.index('time_tag'), hdr.index('speed')
        obs = [(parse_t(r[it]), float(r[iv])) for r in rows[1:] if r[iv] not in (None, '')]
        # tuntikeskiarvot, jotta yksittäiset piikit eivät dominoi
        buckets = {}
        for t, v in obs:
            buckets.setdefault(t.replace(minute=0, second=0, microsecond=0), []).append(v)
        pts += [(t, float(np.median(v))) for t, v in buckets.items()]
        src.append('SWPC plasma-7-day (mitattu)')
    except Exception as e:
        print('plasma-7-day epäonnistui:', e)
    try:
        en = get_json('https://services.swpc.noaa.gov/json/enlil_time_series.json')
        last_obs = max((t for t, _ in pts), default=now - D.timedelta(days=1))
        fc = []
        for r in en:
            try:
                t = parse_t(r['time_tag']); v = float(r['v_r'])
                cloud = float(r.get('cloud') or 0)
            except (KeyError, TypeError, ValueError):
                continue
            if t > last_obs and cloud <= 0:   # cloud > 0 = CME-materiaalia Enlil-mallissa
                fc.append((t, v))
        pts += fc
        src.append(f'SWPC WSA-Enlil (enlil_time_series, {len(fc)} pistettä, CME:t poistettu)')
    except Exception as e:
        print('enlil_time_series epäonnistui:', e)
    if len(pts) < 24:
        raise RuntimeError('Taustatuulen dataa liian vähän')
    pts.sort()
    return pts, src


# ---------------------------------------------------------------- HUXt
class Model:
    def __init__(self, t0, ambient_pts, simtime_days):
        self.t0 = t0
        self.simtime = simtime_days * u.day
        tE = np.array([(t - t0).total_seconds() / 86400 for t, _ in ambient_pts])
        vE = np.array([v for _, v in ambient_pts])
        self.cr, self.cr_lon_init = HIN.datetime2huxtinputs(t0.replace(tzinfo=None))
        nlon = 128
        self.lon = np.arange(np.pi / nlon, 2 * np.pi, 2 * np.pi / nlon)
        # pituuspiiri L ohittaa Maan hetkellä t0 + (L_maa - L)/omega  (paikallaan pysyvä virtausrakenne)
        dt_lon = ((self.cr_lon_init.value - self.lon) % (2 * np.pi)) / (2 * np.pi) * SYN
        known_t = np.concatenate([tE - SYN, tE, tE + SYN])
        known_v = np.concatenate([vE, vE, vE])
        o = np.argsort(known_t)
        v1au = np.interp(dt_lon, known_t[o], known_v[o]) * u.km / u.s   # aukot: jaksollinen interpolointi
        self.coverage_days = float(tE.max() - tE.min())
        self.vb = HIN.map_v_boundary_inwards(v1au, 215 * u.solRad, R_IN)
        self.elat = HIN.get_earth_lat(t0.replace(tzinfo=None))

    def run(self, cmes, rng=None):
        vb = self.vb
        if rng is not None:
            sh = np.deg2rad(rng.uniform(-PERT['amb_shift'], PERT['amb_shift']))
            vb = np.interp(self.lon, (self.lon + sh) % (2 * np.pi), self.vb.value, period=2 * np.pi) * u.km / u.s
        m = H.HUXt(v_boundary=vb, cr_num=self.cr, cr_lon_init=self.cr_lon_init, latitude=self.elat.to(u.deg),
                   r_min=R_IN, lon_out=0 * u.deg, simtime=self.simtime, dt_scale=4, frame='synodic')
        cl = []
        for c in cmes:
            lo, la, w, v = c['lon'], c['lat'], 2 * c['half'], c['v']
            if rng is not None:
                lo += rng.uniform(-PERT['lon'], PERT['lon']); la += rng.uniform(-PERT['lat'], PERT['lat'])
                w *= rng.uniform(1 - PERT['width'], 1 + PERT['width']); v *= rng.uniform(1 - PERT['v'], 1 + PERT['v'])
            tl = (parse_t(c['t215']) - self.t0).total_seconds() * u.s
            cl.append(H.ConeCME(t_launch=tl, longitude=lo * u.deg, latitude=la * u.deg, v=v * u.km / u.s,
                                width=w * u.deg, thickness=0 * u.solRad, initial_height=R_IN, label=c['id']))
        m.solve(cl)
        out = []
        for cm in m.cmes:              # HUXt kopioi CME-listan; seurantatiedot ovat mallin kopioissa
            st = cm.compute_arrival_at_body('EARTH')
            out.append((bool(st['hit']), st['t_arrive'].unix if st['hit'] else None,
                        float(st['v'].value) if st['hit'] else None))
        return m, out


_MODEL = None; _CMES = None
def _member(seed):
    """Yksi ensemble-jäsen (ajetaan rinnakkain eri prosesseissa)."""
    return _MODEL.run(_CMES, np.random.default_rng(seed))[1]


def iso(unix):
    return D.datetime.fromtimestamp(unix, UTC).strftime('%Y-%m-%dT%H:%MZ')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--offline', action='store_true', help='testidata kansiosta testdata/')
    ap.add_argument('--now', help='nykyhetki (UTC, ISO) testausta varten')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    now = parse_t(a.now) if a.now else D.datetime.now(UTC).replace(microsecond=0)

    here = os.path.dirname(os.path.abspath(__file__))
    if a.offline:
        sims = json.load(open(os.path.join(here, 'testdata', 'donki_enlil_2026-10-12.json')))
        amb = [(parse_t(t + 'Z'), v) for t, v in json.load(open(os.path.join(here, 'testdata', 'ambient_digitized.json')))]
        amb_src = ['TESTI: Readingin taustatuuli digitoituna kuvasta']
    else:
        sims = fetch_donki(now)
        amb, amb_src = fetch_ambient(now)

    cmes = collect_cmes(sims, now)
    key = hashlib.sha1(json.dumps([[c['id'], c['t215'], c['lon'], c['lat'], c['half'], c['v']] for c in cmes]).encode()).hexdigest()[:12]

    prev = None
    if os.path.exists(OUT):
        try: prev = json.load(open(OUT))
        except Exception: prev = None
    if prev and not a.force and prev.get('input_key') == key:
        age_h = (now - parse_t(prev['generated'])).total_seconds() / 3600
        if age_h < RERUN_HOURS:
            print(f'Ei uusia CME:itä ja edellinen ajo {age_h:.1f} h vanha -> ei ajeta.')
            return

    result = dict(generated=now.strftime('%Y-%m-%dT%H:%M:%SZ'), model='HUXt ' + H.get_version(),
                  input_key=key, members=N_MEMBERS, perturbation=PERT, ambient_source=amb_src,
                  note='Kokeellinen ennuste. Tulos on yhtä hyvä kuin syötteet (DONKI-kartio-CME:t ja SWPC:n taustatuuli).',
                  cmes=[])
    if not cmes:
        print('Ei tuoreita CME:itä.')
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        json.dump(result, open(OUT, 'w'), indent=1, ensure_ascii=False)
        return

    # vain Maan puoleiset CME:t ajetaan; muut kirjataan osumatodennäköisyydellä 0
    # Lisäksi pois CME:t, jotka olisivat ehtineet Maahan jo ennen nykyhetkeä hitaallakin (0,8 x) nopeudella.
    AU_KM = (215 - 21.5) * 696000.0
    def still_relevant(c):
        transit_s = AU_KM / (0.8 * c['v'])
        return parse_t(c['t215']).timestamp() + transit_s > now.timestamp() - 12 * 3600
    run_cmes = [c for c in cmes if abs(((c['lon'] + 180) % 360) - 180) <= MAX_ABS_LON and still_relevant(c)]
    for c in cmes:
        if c in run_cmes: continue
        rec = dict(c, width=2 * c['half'], hit_pct=None, arrival_median=None,
                   status='ei kohti Maata' if abs(((c['lon'] + 180) % 360) - 180) > MAX_ABS_LON else 'jo ohi (ei ajettu)')
        rec.pop('half'); result['cmes'].append(rec)
    cmes = run_cmes
    if not cmes:
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        json.dump(result, open(OUT, 'w'), indent=1, ensure_ascii=False)
        print('Ei Maan puoleisia CME:itä.'); return
    t0 = (parse_t(cmes[0]['t215']) - D.timedelta(hours=12)).replace(minute=0, second=0)
    simdays = (now - t0).total_seconds() / 86400 + FORECAST_DAYS
    model = Model(t0, amb, simdays)
    result.update(model_start=t0.strftime('%Y-%m-%dT%H:%MZ'), ambient_coverage_days=round(model.coverage_days, 1))

    # deterministinen ajo (ilman häiriöitä) + taustatuuli Maan kohdalla kaaviota varten
    m0, det = model.run(cmes)
    e = m0.get_observer('earth'); ir = int(np.argmin(np.abs(m0.r - np.mean(e.r))))
    tt = m0.time_out.to(u.s).value; vv = m0.v_grid[:, ir, 0].value
    step = max(1, int(round(3 * 3600 / (tt[1] - tt[0]))))
    result['earth_speed_deterministic'] = [[iso(t0.timestamp() + t), round(float(v))] for t, v in zip(tt[::step], vv[::step])]

    global _MODEL, _CMES
    _MODEL, _CMES = model, cmes
    seeds = np.random.SeedSequence(int(now.timestamp()) // 3600).generate_state(N_MEMBERS)
    workers = max(1, min(int(os.environ.get('HUXT_WORKERS', os.cpu_count() or 1)), 8))
    tic = time.time()
    if workers > 1:
        import multiprocessing as mp
        with mp.get_context('fork').Pool(workers) as pool:
            members = pool.map(_member, seeds.tolist(), chunksize=4)
    else:
        members = [_member(sd) for sd in seeds.tolist()]
    result['workers'] = workers
    result['runtime_s'] = round(time.time() - tic, 1)

    for k, c in enumerate(cmes):
        hits = [mb[k] for mb in members if mb[k][0]]
        ts = np.array([h[1] for h in hits]); vs = np.array([h[2] for h in hits])
        pct = lambda p: iso(float(np.percentile(ts, p))) if len(ts) else None
        rec = dict(c, width=2 * c['half'], hit_pct=round(100 * len(hits) / N_MEMBERS, 1),
                   arrival_median=pct(50), arrival_p10=pct(10), arrival_p90=pct(90),
                   arrival_p5=pct(5), arrival_p95=pct(95),
                   speed_at_earth_median=round(float(np.median(vs))) if len(vs) else None,
                   deterministic_hit=det[k][0], deterministic_arrival=iso(det[k][1]) if det[k][0] else None)
        if not hits:
            rec['status'] = 'ei osumaa'
        elif float(np.percentile(ts, 95)) < now.timestamp():
            rec['status'] = 'saapunut / ohi'
        else:
            rec['status'] = 'tulossa'
        rec.pop('half')
        result['cmes'].append(rec)

    result['cmes'].sort(key=lambda c: c['t215'])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(result, open(OUT, 'w'), indent=1, ensure_ascii=False)
    print(json.dumps({k: v for k, v in result.items() if k != 'earth_speed_deterministic'}, indent=1, ensure_ascii=False))


if __name__ == '__main__':
    main()
